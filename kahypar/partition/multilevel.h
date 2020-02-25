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
static constexpr bool debug = false;

static inline void partition(Hypergraph& hypergraph,
                             ICoarsener& coarsener,
                             IRefiner& refiner,
                             const Context& context) {
  // TODO is this correct?!
  ASSERT(!hypergraph.containsFixedVertices(), "Fixed vertices not allowed here.");
  io::printCoarseningBanner(context);

  std::vector<HypernodeWeight> adjusted_weight;
  PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;
  // TODO this is fragile, as it depends on rb_range_k == 1 before top-level coarsening
    // perform prepacking of heavy vertices
  if ((context.initial_partitioning.balancing == WeightBalancingStrategy::prepacking_pessimistic
      || context.initial_partitioning.balancing == WeightBalancingStrategy::prepacking_optimistic) && (rb_range_k > 2)) {
    Context packing_context = initial::createContext(hypergraph, context);
    bool optimistic = context.initial_partitioning.balancing == WeightBalancingStrategy::prepacking_optimistic;
    bin_packing::prepack_heavy_vertices(hypergraph, packing_context, rb_range_k, optimistic);
  } else if (context.initial_partitioning.balancing == WeightBalancingStrategy::prepacking_dynamic && (rb_range_k > 2)) {
    Context packing_context = initial::createContext(hypergraph, context);
    HypernodeWeight max_bin_weight = floor(packing_context.initial_partitioning.current_max_bin * (1.0 + packing_context.initial_partitioning.bin_epsilon));
    for (size_t i = 0; i < packing_context.initial_partitioning.k; ++i) {
      HypernodeWeight lower = packing_context.initial_partitioning.perfect_balance_partition_weight[i];
      HypernodeWeight& border = packing_context.initial_partitioning.upper_allowed_partition_weight[i];
      HypernodeWeight upper = packing_context.initial_partitioning.num_bins_per_partition[i] * max_bin_weight;

      // TODO ugly magic numbers
      if (upper - border < (border - lower) / 10) {
        border = (lower + upper) / 2;
        packing_context.partition.epsilon = static_cast<double>(border) / static_cast<double>(lower) - 1.0;
      }
    }

    if (context.initial_partitioning.bp_algo ==  BinPackingAlgorithm::worst_fit) {
      bin_packing::apply_prepacking_pessimistic<bin_packing::WorstFit>(hypergraph, packing_context);
    } else {
      bin_packing::apply_prepacking_pessimistic<bin_packing::FirstFit>(hypergraph, packing_context);
    }
    adjusted_weight = packing_context.initial_partitioning.upper_allowed_partition_weight;
  }

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  coarsener.coarsen(context.coarsening.contraction_limit);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
  Timer::instance().add(context, Timepoint::coarsening,
                        std::chrono::duration<double>(end - start).count());

  if (context.partition.verbose_output && context.type == ContextType::main) {
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
    if (context.partition.verbose_output && context.type == ContextType::main) {
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

  if (hypergraph.containsFixedVertices()) {
    for (const HypernodeID& hn : hypergraph.fixedVertices()) {
      hypergraph.setVertexNotFixed(hn);
    }
  }

  io::printLocalSearchResults(context, hypergraph);
}
}  // namespace multilevel
}  // namespace kahypar
