/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Nikolai Maas <nikolai.maas@student.kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include "kahypar/partition/bin_packing.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {
class BinPackingInitialPartitioner : public IInitialPartitioner,
                                 private InitialPartitionerBase<BinPackingInitialPartitioner>{
  using Base = InitialPartitionerBase<BinPackingInitialPartitioner>;
  friend Base;

 public:
  BinPackingInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context),
    _descending_nodes() { }

  ~BinPackingInitialPartitioner() override = default;

  BinPackingInitialPartitioner(const BinPackingInitialPartitioner&) = delete;
  BinPackingInitialPartitioner& operator= (const BinPackingInitialPartitioner&) = delete;

  BinPackingInitialPartitioner(BinPackingInitialPartitioner&&) = delete;
  BinPackingInitialPartitioner& operator= (BinPackingInitialPartitioner&&) = delete;

 private:
  void initializeNodes() {
    // exctract hypernodes in descending order
    if (_descending_nodes.empty()) {
      _descending_nodes = bin_packing::extract_nodes_with_descending_weight(_hg);
    }

    // shuffle nodes with equal weight
    size_t i = 0;
    for (size_t j = 1; j < _descending_nodes.size(); ++j) {
      if (_hg.nodeWeight(_descending_nodes[i]) != _hg.nodeWeight(_descending_nodes[j])) {
        if (j > i + 1) {
          Randomize::instance().shuffleVector(_descending_nodes, i, j);
        }
        i = j;
      }
    }

    ASSERT(_descending_nodes.size() == _hg.currentNumNodes(),
      "Number of extracted nodes does not match nodes in hypergraph.");

    ASSERT([&]() {
      for (const HypernodeID& hn : _hg.nodes()) {
        if (std::find(_descending_nodes.cbegin(), _descending_nodes.cend(), hn) == _descending_nodes.cend()) {
          return false;
        }
      }
      return true;
    } (),"The extracted nodes do not match the nodes in the hypergraph.");

    ASSERT([&]() {
      for (size_t i = 1; i < _descending_nodes.size(); ++i) {
        if (_hg.nodeWeight(_descending_nodes[i-1]) < _hg.nodeWeight(_descending_nodes[i])) {
          return false;
        }
      }
      return true;
    } (), "The hypernodes must be sorted in descending order of weight.");
  }

  void partitionImpl() override final {
    Base::multipleRunsInitialPartitioning();

    if (_context.initial_partitioning.balancing == WeightBalancingStrategy::restart_ip_algorithm) {
      Base::restartAtImbalancedBins();
    }
  }

  void initialPartition() {
    const PartitionID unassigned_part = _context.initial_partitioning.unassigned_part;
    _context.initial_partitioning.unassigned_part = -1;
    Base::resetPartitioning();
    initializeNodes();

    std::vector<PartitionID> partitions = bin_packing::apply_bin_packing_to_nodes(_hg, _context, _descending_nodes);

    for (size_t i = 0; i < _descending_nodes.size(); ++i) {
      HypernodeID hn = _descending_nodes[i];

      bool assigned = false;
      PartitionID p = partitions[i];
      do {
        if (Base::assignHypernodeToPartition(hn, p)) {
          assigned = true;
          break;
        }

        ++p;
        p %= _context.initial_partitioning.k;
      } while (p != partitions[i]);

      // If the current hypernode fits in no part of the partition, we have
      // to assign it to a part which violates the imbalance definition.
      if (!assigned) {
        _hg.setNodePart(hn, p);
      }

      ASSERT(_hg.partID(hn) == p, "Hypernode" << hn << "should be in part" << p
                                              << ", but is actually in" << _hg.partID(hn) << ".");
    }

    _context.initial_partitioning.unassigned_part = unassigned_part;

    ASSERT([&]() {
        for (const HypernodeID& hn : _hg.nodes()) {
          if (_hg.partID(hn) == -1) {
            return false;
          }
        }
        return true;
      } (), "There are unassigned hypernodes!");

    Base::performFMRefinement();
  }

  using Base::_hg;
  using Base::_context;

  std::vector<HypernodeID> _descending_nodes;
};
}  // namespace kahypar
