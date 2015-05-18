/*
 * InitialPartitionerBase.h
 *
 *  Created on: Apr 3, 2015
 *      Author: theuer
 */

#ifndef SRC_PARTITION_INITIAL_PARTITIONING_INITIALPARTITIONERBASE_H_
#define SRC_PARTITION_INITIAL_PARTITIONING_INITIALPARTITIONERBASE_H_

#include <vector>
#include <stack>
#include <map>
#include <algorithm>
#include <limits>

#include "lib/definitions.h"
#include "partition/Metrics.h"
#include "partition/initial_partitioning/InitialStatManager.h"
#include "partition/initial_partitioning/HypergraphPartitionBalancer.h"
#include "partition/Configuration.h"
#include "partition/refinement/KWayFMRefiner.h"
#include "partition/refinement/policies/FMImprovementPolicies.h"
#include "partition/refinement/policies/FMStopPolicies.h"
#include "partition/initial_partitioning/HypergraphPartitionBalancer.h"

using partition::Configuration;
using defs::Hypergraph;
using defs::HypernodeID;
using defs::HyperedgeID;
using defs::HypernodeWeight;
using defs::HyperedgeWeight;
using defs::HyperedgeIndexVector;
using defs::HyperedgeVector;
using defs::HyperedgeWeightVector;
using defs::HypernodeWeightVector;
using defs::HighResClockTimepoint;

using partition::HypergraphPartitionBalancer;
using partition::InitialStatManager;
using partition::StoppingPolicy;
using partition::NumberOfFruitlessMovesStopsSearch;
using partition::RandomWalkModelStopsSearch;
using partition::nGPRandomWalkStopsSearch;
using partition::KWayFMRefiner;

namespace partition {

struct node_assignment {
  HypernodeID hn;
  PartitionID from;
  PartitionID to;
};

class InitialPartitionerBase {

 public:

  InitialPartitionerBase(Hypergraph& hypergraph, Configuration& config) noexcept:
  _hg(hypergraph),
    _config(config),
    bisection_assignment_history(),
    refiner(_hg, _config),
    _balancer(_hg, _config) {

    for (const HypernodeID hn : _hg.nodes()) {
      total_hypergraph_weight += _hg.nodeWeight(hn);
      if(_hg.nodeWeight(hn) > heaviest_node) {
        heaviest_node = _hg.nodeWeight(hn);
      }
    }

    // TODO(heuer): This should be done in your application.
    // In general, the partitioner should only use the config parameters, not
    // set them.
    _config.partition.total_graph_weight = total_hypergraph_weight;

  }

  virtual ~InitialPartitionerBase() {}

  void recalculateBalanceConstraints() {
    // TODO(heuer): The calculation could be done once and the
    // values then only assigned.
    for (int i = 0; i < _config.initial_partitioning.k; i++) {
      _config.initial_partitioning.lower_allowed_partition_weight[i] =
          ceil(
              total_hypergraph_weight
              / static_cast<double>(_config.initial_partitioning.k))
          * (1.0 - _config.partition.epsilon);
      _config.initial_partitioning.upper_allowed_partition_weight[i] =
          ceil(
              total_hypergraph_weight
              / static_cast<double>(_config.initial_partitioning.k))
          * (1.0 + _config.partition.epsilon);
    }
  }

  void eraseConnectedComponents() {
    if(_config.initial_partitioning.erase_components) {
      HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
      _balancer.eraseSmallConnectedComponents();
      HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      InitialStatManager::getInstance().updateStat("Time Measurements", "Erasing Connected Components time",InitialStatManager::getInstance().getStat("Time Measurements", "Erasing Connected Components time") + static_cast<double>(elapsed_seconds.count()));
    }
  }

  void balancePartitions() {
    if(_config.initial_partitioning.balance) {
      _balancer.balancePartitions();
    }
  }

  void resetPartitioning(PartitionID unassigned_part) {
    _hg.resetPartitioning();
    // TODO(heuer): For efficiency: Do you want a reset-Method where you
    // can choose what unassigned means, i.e. reset everything and move
    // all nodes into a predefined part?
    if(unassigned_part != -1) {
      for(HypernodeID hn : _hg.nodes()) {
        _hg.setNodePart(hn, unassigned_part);
      }
    }
  }

  void performFMRefinement() {
    if(_config.initial_partitioning.refinement) {
      _config.partition.total_graph_weight = total_hypergraph_weight;
      refiner.initialize();

      std::vector<HypernodeID> refinement_nodes;
      for(HypernodeID hn : _hg.nodes()) {
        refinement_nodes.push_back(hn);
      }
      HyperedgeWeight cut_before = metrics::hyperedgeCut(_hg);
      HyperedgeWeight cut = cut_before;
      double imbalance = metrics::imbalance(_hg);

      // TODO(heuer): This is still an relevant issue! I think we should not test refinement as long as it is
      // not possible to give more than one upper bound to the refiner.
      // However, if I'm correct, the condition always evaluates to true if k=2^x right?
      //Only perform refinement if the weight of partition 0 and 1 is the same to avoid unexpected partition weights.
      if(_config.initial_partitioning.upper_allowed_partition_weight[0] == _config.initial_partitioning.upper_allowed_partition_weight[1]) {
        HypernodeWeight max_allowed_part_weight = _config.initial_partitioning.upper_allowed_partition_weight[0];
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        // TODO(heuer): If you look at the uncoarsening code that calls the refiner, you see, that
        // another idea is to restart the refiner as long as it finds an improvement on the current
        // level. This should also be evaluated. Actually, this is, what parameter --FM-reps is used
        //for.
        refiner.refine(refinement_nodes,_hg.numNodes(),max_allowed_part_weight,cut,imbalance);
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        InitialStatManager::getInstance().updateStat("Partitioning Results", "Cut increase during refinement",InitialStatManager::getInstance().getStat("Partitioning Results", "Cut increase during refinement") + (cut_before - metrics::hyperedgeCut(_hg)));
        InitialStatManager::getInstance().updateStat("Time Measurements", "Refinement time",InitialStatManager::getInstance().getStat("Time Measurements", "Refinement time") + static_cast<double>(elapsed_seconds.count()));
      }
    }
  }

  void rollbackToBestCut() {
    if (_config.initial_partitioning.rollback) {
      if(!bisection_assignment_history.empty()) {
        HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
        HyperedgeWeight cut_before = metrics::hyperedgeCut(_hg);
        HypernodeID current_node = std::numeric_limits<HypernodeID>::max();
        // TODO(heuer): Please use PartitionID from = ...;
        // PartitionID to = ...;
        PartitionID from, to;
        while (current_node != best_cut_node
               && !bisection_assignment_history.empty()) {
          if(current_node != std::numeric_limits<HypernodeID>::max()) {
            _hg.changeNodePart(current_node,to,from);
          }
          node_assignment assignment = bisection_assignment_history.top();
          bisection_assignment_history.pop();
          current_node = assignment.hn;
          from = assignment.from;
          to = assignment.to;
        }
        HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        // TODO(heuer): Try to restrict line length to 100.
        ASSERT(metrics::hyperedgeCut(_hg) == best_cut,
               "Calculation of the best seen cut failed. Should be " << best_cut << ", but is " << metrics::hyperedgeCut(_hg)<< ".");
        InitialStatManager::getInstance().updateStat("Partitioning Results", "Cut increase during rollback",InitialStatManager::getInstance().getStat("Partitioning Results", "Cut increase during rollback") + (cut_before-best_cut));
        InitialStatManager::getInstance().updateStat("Time Measurements", "Rollback time",InitialStatManager::getInstance().getStat("Time Measurements", "Rollback time") + static_cast<double>(elapsed_seconds.count()));
      }
    }
  }

  // TODO(heuer): Make parameters const, although they are passed by value. This way
  // you get a compiler error when modifying them unintentionally.
  bool assignHypernodeToPartition(HypernodeID hn, PartitionID target_part) {
    HypernodeWeight assign_partition_weight = _hg.partWeight(target_part)
                                              + _hg.nodeWeight(hn);
    PartitionID source_part = _hg.partID(hn);
    // TODO(heuer): Why do you need the second condition: && hn < _hg.numNodes()?
    if (assign_partition_weight
        <= _config.initial_partitioning.upper_allowed_partition_weight[target_part] && hn < _hg.numNodes()) {
      if (_hg.partID(hn) == -1) {
        _hg.setNodePart(hn, target_part);
      } else {
        if (_hg.partID(hn) != target_part) {
          _hg.changeNodePart(hn, _hg.partID(hn), target_part);
        }
        else {
          return false;
        }
      }
      if(_config.initial_partitioning.rollback) {
        calculateCutAfterAssignment(hn, source_part, target_part);
      }
      ASSERT(_hg.partID(hn) == target_part,
             "Assigned partition of Hypernode "<< hn<< " should be " << target_part << ", but currently is " << _hg.partID(hn));
      return true;
    } else {
      return false;
    }
  }

  // TODO(heuer): Are you sure that this code, in combination with GHG will never
  // run in an endless loop?
  HypernodeID getUnassignedNode(PartitionID unassigned_part = -1) {
    HypernodeID unassigned_node = Randomize::getRandomInt(0,
                                                          _hg.numNodes() - 1);
    while (_hg.partID(unassigned_node) != unassigned_part) {
      unassigned_node = Randomize::getRandomInt(0,
                                                _hg.numNodes() - 1);
    }
    return unassigned_node;
  }

  // TODO(heuer): Kill this and use my version :)
  void extractPartitionAsHypergraph(Hypergraph& hyper, PartitionID part,
                                    HypernodeID& num_hypernodes, HyperedgeID& num_hyperedges,
                                    HyperedgeIndexVector& index_vector, HyperedgeVector& edge_vector,
                                    HyperedgeWeightVector& hyperedge_weights,
                                    HypernodeWeightVector& hypernode_weights,
                                    std::vector<HypernodeID>& hgToExtractedPartitionMapping) {

    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();

    std::unordered_map<HypernodeID, HypernodeID> hypernodeMapper;
    for (HypernodeID hn : hyper.nodes()) {
      if (hyper.partID(hn) == part) {
        hypernodeMapper.insert(
            std::pair<HypernodeID, HypernodeID>(hn,
                                                hgToExtractedPartitionMapping.size()));
        hgToExtractedPartitionMapping.push_back(hn);
        hypernode_weights.push_back(hyper.nodeWeight(hn));
      }
    }

    ASSERT([&]() {
        for (const HypernodeID hn : hgToExtractedPartitionMapping) {
          if(hyper.partID(hn) != part)
            return false;
        }
        return true;
      }(),"There is a hypernode in the new hypergraph from a wrong partition.");

    ASSERT([&]() {
        for (unsigned int i = 0; i < hypernode_weights.size(); i++) {
          if(hyper.nodeWeight(hgToExtractedPartitionMapping[i]) != hypernode_weights[i])
            return false;
        }
        return true;
      }(),"Assign hypernode weights to the new hypergraph fail.");

    num_hypernodes = hgToExtractedPartitionMapping.size();
    index_vector.push_back(edge_vector.size());
    std::vector<HyperedgeID> hyperedgeMapper;
    for (HyperedgeID he : hyper.edges()) {
      if (hyper.connectivity(he) > 1) {
        continue;
      }
      bool is_part_edge = true;
      for (HypernodeID hn : hyper.pins(he)) {
        if (hyper.partID(hn) != part) {
          is_part_edge = false;
        }
        break;
      }
      if (is_part_edge) {
        for (HypernodeID hn : hyper.pins(he)) {
          edge_vector.push_back(hypernodeMapper[hn]);
        }
        index_vector.push_back(edge_vector.size());
        hyperedgeMapper.push_back(he);
        hyperedge_weights.push_back(hyper.edgeWeight(he));
      }
    }
    num_hyperedges = index_vector.size() - 1;

    ASSERT([&]() {
        for (unsigned int i = 0; i < hyperedgeMapper.size(); i++) {
          if((index_vector[i+1] - index_vector[i]) != hyper.edgeSize(hyperedgeMapper[i]))
            return false;
          if(hyperedge_weights[i] != hyper.edgeWeight(hyperedgeMapper[i]))
            return false;
        }
        return true;
      }(),"Size/Weight of a hyperedge in the extracted hypergraph differs from a hyperedge size/weight in original hypergraph.");

    ASSERT([&]() {
        for(int i = 0; i < hyperedgeMapper.size(); i++) {
          int startIndex = index_vector[i];
          for(HypernodeID hn : hyper.pins(hyperedgeMapper[i])) {
            if(hn != hgToExtractedPartitionMapping[edge_vector[startIndex]])
              return false;
            startIndex++;
          }
        }
        return true;
      }(),"Different hypernodes in hyperedge from the extracted hypergraph as in original hypergraph.");

    ASSERT([&]() {
        for(int i = 0; i < hyperedgeMapper.size(); i++) {
          for(HypernodeID hn : hyper.pins(hyperedgeMapper[i])) {
            if(hyper.partID(hn) != part)
              return false;
          }
        }
        return true;
      }(),"There are cut edges in the extracted hypergraph.");

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    InitialStatManager::getInstance().updateStat("Time Measurements", "Extracting Hypergraph time",InitialStatManager::getInstance().getStat("Time Measurements", "Extracting Hypergraph time") + static_cast<double>(elapsed_seconds.count()));
  }

 protected:
  Hypergraph& _hg;
  Configuration& _config;
  HypergraphPartitionBalancer _balancer;
  HypernodeWeight total_hypergraph_weight = 0;
  HypernodeWeight heaviest_node = 0;

 private:

  void calculateCutAfterAssignment(HypernodeID hn, PartitionID from, PartitionID to) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    for(HyperedgeID he : _hg.incidentEdges(hn)) {
      HyperedgeID pins_in_source_part_before = 0;
      // TODO(heuer): If from == -1, then pins_in_source_part should be sth like _hg.pinCountInPart(he,from)?
      // This might not be relevant right now for determining whether or not the hyperedge is a
      // cut hyperedge, however the current version only works for cut objective.
      if(from != -1) {
        pins_in_source_part_before = _hg.pinCountInPart(he,from) + 1;
      }
      // TODO(heuer): Make variables that are not intended to be modified const.
      HyperedgeID pins_in_target_part_after = _hg.pinCountInPart(he,to);
      HyperedgeID connectivity_before = _hg.connectivity(he);
      if(pins_in_source_part_before == 1) {
        connectivity_before++;
      }
      if(pins_in_target_part_after == 1) {
        connectivity_before--;
      }
      if(_hg.connectivity(he) == 1 && connectivity_before == 2) {
        current_cut -= _hg.edgeWeight(he);
      }
      if(_hg.connectivity(he) == 2 && connectivity_before == 1) {
        current_cut += _hg.edgeWeight(he);
      }
    }

    ASSERT(current_cut == metrics::hyperedgeCut(_hg),"Calculated hyperedge cut doesn't match with the real hyperedge cut!");

    bool isFeasibleSolution = true;
    // TODO(heuer): ??? Only to-part can become overloaded. It should not be necessary to check all parts
    // What about your lower bound? Shouldn't it be checked here as well?
    for(PartitionID k = 0; k < _config.initial_partitioning.k; k++) {
      if(_hg.partWeight(k) > _config.initial_partitioning.upper_allowed_partition_weight[k]) {
        isFeasibleSolution = false;
        break;
      }
    }
    if(isFeasibleSolution) {
      node_assignment assign;
      assign.hn = hn;
      assign.from = from;
      assign.to = to;
      if(current_cut < best_cut) {
        best_cut = current_cut;
        best_cut_node = hn;
      }
      // TODO(heuer): Potentially it is slightly more efficient to use a vector.
      // Additionally, using emplace, the node-assignment could be constructed inplace
      bisection_assignment_history.push(assign);
    }

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    InitialStatManager::getInstance().updateStat("Time Measurements", "Precalculation Rollback time",InitialStatManager::getInstance().getStat("Time Measurements", "Precalculation Rollback time") + static_cast<double>(elapsed_seconds.count()));
  }

  // TODO(heuer): naming conventions. Member variables are named _name!
  HypernodeID best_cut_node = std::numeric_limits<HypernodeID>::max();
  HyperedgeWeight best_cut = std::numeric_limits<HyperedgeWeight>::max();
  HyperedgeWeight current_cut = 0;
  std::stack<node_assignment> bisection_assignment_history;

  KWayFMRefiner<NumberOfFruitlessMovesStopsSearch,
                CutDecreasedOrInfeasibleImbalanceDecreased> refiner;

}
  ;

}

#endif /* SRC_PARTITION_INITIAL_PARTITIONING_INITIALPARTITIONERBASE_H_ */
