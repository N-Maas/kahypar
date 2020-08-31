/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/definitions.h"
#include "kahypar/partition/bin_packing/bin_packing.h"

namespace kahypar {

namespace bin_packing {

class IBinPacker {
 public:
  IBinPacker(const IBinPacker&) = delete;
  IBinPacker(IBinPacker&&) = delete;
  IBinPacker& operator= (const IBinPacker&) = delete;
  IBinPacker& operator= (IBinPacker&&) = delete;

  void prepacking(Hypergraph& hg, Context& context, const BalancingLevel level) {
    prepackingImpl(hg, context, level);
  }

  std::vector<PartitionID> twoLevelPacking(const Hypergraph& hg, const Context& context, const std::vector<HypernodeID>& nodes) {
    twoLevelPackingImpl(hg, context, nodes);
  }

  virtual ~IBinPacker() = default;

 protected:
  IBinPacker() = default;

 private:
  virtual void prepackingImpl(Hypergraph& hg, Context& context, const BalancingLevel level) = 0;
  virtual std::vector<PartitionID> twoLevelPackingImpl(const Hypergraph& hg, const Context& context, const std::vector<HypernodeID>& nodes) = 0;
};

template< class BPAlgorithm = WorstFit >
class BinPacker final : public IBinPacker {
 public:
  BinPacker() = default;

 private:
  void prepackingImpl(Hypergraph& hg, Context& context, const BalancingLevel level) override {
    ASSERT(!hg.containsFixedVertices(), "Fixed vertices not allowed here.");
    ASSERT(level != BalancingLevel::STOP, "Invalid balancing level: STOP");

    PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;

    if (level == BalancingLevel::optimistic) {
        bin_packing::prepack_heavy_vertices<BPAlgorithm>(hg, context, rb_range_k);
    } else if (level == BalancingLevel::guaranteed) {
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1.0 + context.initial_partitioning.bin_epsilon));
        for (size_t i = 0; i < static_cast<size_t>(context.initial_partitioning.k); ++i) {
            HypernodeWeight lower = context.initial_partitioning.perfect_balance_partition_weight[i];
            HypernodeWeight& border = context.initial_partitioning.upper_allowed_partition_weight[i];
            HypernodeWeight upper = context.initial_partitioning.num_bins_per_partition[i] * max_bin_weight;

            // TODO ugly magic numbers
            if (upper - border < (border - lower) / 10) {
                border = (lower + upper) / 2;
                context.partition.epsilon = static_cast<double>(border) / static_cast<double>(lower) - 1.0;
            }
        }

       bin_packing::apply_prepacking_pessimistic<BPAlgorithm>(hg, context);
    }
  }

  std::vector<PartitionID> twoLevelPackingImpl(const Hypergraph& hg, const Context& context, const std::vector<HypernodeID>& nodes) override {

  }
};

static std::unique_ptr<IBinPacker> createBinPacker(const BinPackingAlgorithm& bp_algo) {
switch (bp_algo) {
    case BinPackingAlgorithm::worst_fit: return std::make_unique<BinPacker<WorstFit>>();
    case BinPackingAlgorithm::first_fit: return std::make_unique<BinPacker<FirstFit>>();
    case BinPackingAlgorithm::UNDEFINED: break;
      // omit default case to trigger compiler warning for missing cases
  }
  return std::unique_ptr<IBinPacker>();
}
} // namespace bin_packing
} // namespace kahypar
