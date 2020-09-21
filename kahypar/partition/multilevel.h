/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partition/coarsening/i_coarsener.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/initial_partition.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/utils/timer.h"

namespace kahypar {
namespace multilevel {
using bin_packing::BalancingLevel;
using bin_packing::IBinPacker;

static constexpr bool debug = false;

static inline void partition(Hypergraph& hypergraph,
                             ICoarsener& coarsener,
                             IRefiner& refiner,
                             const Context& context,
                             const std::vector<HypernodeWeight>& adjusted_weight = std::vector<HypernodeWeight>()) {
  io::printCoarseningBanner(context);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  coarsener.coarsen(context.coarsening.contraction_limit);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::coarsening,
                        std::chrono::duration<double>(end - start).count());

  if (!context.partition.quiet_mode && context.partition.verbose_output && context.type == ContextType::main) {
    io::printHypergraphInfo(hypergraph, "Coarsened Hypergraph");
  }

  if (!context.partition_evolutionary || context.evolutionary.action.requires().initial_partitioning) {
    if (context.partition_evolutionary && context.evolutionary.action.requires().initial_partitioning) {
      hypergraph.reset();
    }
    io::printInitialPartitioningBanner(context);

    start = std::chrono::high_resolution_clock::now();
    initial::partition(hypergraph, context, adjusted_weight);
    end = std::chrono::high_resolution_clock::now();
    Timer::instance().add(context, Timepoint::initial_partitioning,
                          std::chrono::duration<double>(end - start).count());

    hypergraph.initializeNumCutHyperedges();
    if (!context.partition.quiet_mode && context.partition.verbose_output && context.type == ContextType::main) {
      LOG << "Initial Partitioning Result:";
      LOG << "Initial" << context.partition.objective << "      ="
          << (context.partition.objective == Objective::cut ? metrics::hyperedgeCut(hypergraph) :
          metrics::km1(hypergraph));
      LOG << "Initial imbalance =" << metrics::imbalance(hypergraph, context);
      LOG << "Initial part sizes and weights:";
      io::printPartSizesAndWeights(hypergraph);
      LLOG << "Target weights:";
      if (context.partition.mode == Mode::direct_kway) {
        LLOG << "w(*) =" << context.partition.max_part_weights[0] << "\n";
      } else {
        LLOG << "(RB): w(0)=" << context.partition.max_part_weights[0]
             << "w(1)=" << context.partition.max_part_weights[1] << "\n";
      }
    }
  }

  if (context.partition_evolutionary &&
      context.evolutionary.action.requires().evolutionary_parent_contraction) {
    hypergraph.reset();
    ASSERT(!context.evolutionary.action.requires().initial_partitioning);

    // There is currently no reason why an evolutionary contraction should be used
    // in conjunction with initial partitioning ... Yet


    hypergraph.setPartition(*context.evolutionary.parent1);


    const HyperedgeWeight parent_1_objective = metrics::correctMetric(hypergraph, context);

    hypergraph.setPartition(*context.evolutionary.parent2);
    const HyperedgeWeight parent_2_objective = metrics::correctMetric(hypergraph, context);

    if (parent_1_objective < parent_2_objective) {
      hypergraph.setPartition(*context.evolutionary.parent1);
    }
  }


  if (context.partition_evolutionary) {
    hypergraph.initializeNumCutHyperedges();
  }
  DBG << V(metrics::km1(hypergraph));
  DBG << V(metrics::imbalance(hypergraph, context));

  if (context.partition.verbose_output && context.type == ContextType::main) {
    io::printLocalSearchBanner(context);
  }

  start = std::chrono::high_resolution_clock::now();
  coarsener.uncoarsen(refiner);
  end = std::chrono::high_resolution_clock::now();

  Timer::instance().add(context, Timepoint::local_search,
                        std::chrono::duration<double>(end - start).count());

  io::printLocalSearchResults(context, hypergraph);
}

static inline void partitionRepeatedOnInfeasible(Hypergraph& hypergraph,
                                                 const Context& context,
                                                 Context::PartitioningStats& stats,
                                                 const BalancingLevel level,
                                                 const HypernodeWeight maxFeasibleBin,
                                                 bool repeat) {
  ASSERT((context.partition.rb_upper_k - context.partition.rb_lower_k + 1) > 2,
         "Prepacking is not allowed for k <= 2: " << V(context.partition.rb_upper_k) << " - " << context.partition.rb_lower_k);

  Context packing_context = initial::createContext(hypergraph, context);
  packing_context.setupInitialPartitioningPartWeights();
  BalancingLevel currLevel = level;

  do {
    if (currLevel != level) {
      hypergraph.reset();

      // TODO(maas) remove?
      std::string key("restarts_early_level_");
      key += std::to_string(static_cast<uint8_t>(currLevel));
      stats.add(StatTag::InitialPartitioning, key, 1.0);
    }

    // perform prepacking of heavy vertices
    std::unique_ptr<IBinPacker> bin_packer = bin_packing::createBinPacker(context.initial_partitioning.bp_algo);
    bin_packer->prepacking(hypergraph, packing_context, currLevel);

    std::unique_ptr<ICoarsener> coarsener(
      CoarsenerFactory::getInstance().createObject(context.coarsening.algorithm, hypergraph, context, hypergraph.weightOfHeaviestNode()));
    std::unique_ptr<IRefiner> refiner(RefinerFactory::getInstance().createObject(context.local_search.algorithm, hypergraph, context));
    ASSERT(coarsener.get() != nullptr, "coarsener not found");
    ASSERT(refiner.get() != nullptr, "refiner not found");

    partition(hypergraph, *coarsener, *refiner, context, packing_context.initial_partitioning.upper_allowed_partition_weight);
    hypergraph.resetFixedVertices();
    currLevel = bin_packing::increaseBalancingRestrictions(currLevel);
  } while (repeat && currLevel != BalancingLevel::STOP
           && bin_packing::resultingMaxBin(hypergraph, packing_context) > maxFeasibleBin);
}
}  // namespace multilevel
}  // namespace kahypar
