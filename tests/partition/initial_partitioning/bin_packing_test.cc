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

#include "gmock/gmock.h"

#include "kahypar/partition/bin_packing.h"
#include "kahypar/utils/randomize.h"

using ::testing::Eq;
using ::testing::Test;

namespace kahypar {
namespace bin_packing {

class BinPackingTest : public Test {
    public:
        BinPackingTest() :
        hypergraph(0, 0,
                   HyperedgeIndexVector { 0 },
                   HyperedgeVector {}) {
        }

        void initializeWeights(const HypernodeWeightVector& weights) {
            hypergraph = Hypergraph(weights.size(), 0, HyperedgeIndexVector(weights.size() + 1, 0), HyperedgeVector {});

            size_t i = 0;
            for(const HypernodeWeight& w : weights) {
                hypergraph.setNodeWeight(i++, w);
            }
        }

        void createTestContext(Context& c,
                                      const std::vector<HypernodeWeight>& upper_weights,
                                      const std::vector<HypernodeWeight>& perfect_weights,
                                      const std::vector<PartitionID>& num_bins_per_part,
                                      const PartitionID k,
                                      const PartitionID rb_range_k,
                                      const HypernodeWeight max_bin) {
            c.initial_partitioning.upper_allowed_partition_weight = upper_weights;
            c.initial_partitioning.perfect_balance_partition_weight = perfect_weights;
            c.initial_partitioning.num_bins_per_partition = num_bins_per_part;
            c.partition.k = k;
            c.initial_partitioning.k = k;
            c.partition.rb_lower_k = 0;
            c.partition.rb_upper_k = rb_range_k - 1;
            c.initial_partitioning.current_max_bin = max_bin;
            c.initial_partitioning.bin_epsilon = 0.0;
        }

        template< class BPAlg >
        static inline std::vector<PartitionID> two_level_packing(const Hypergraph& hg,
                                                              const std::vector<HypernodeID>& hypernodes,
                                                              const std::vector<HypernodeWeight>& max_allowed_partition_weights,
                                                              const std::vector<PartitionID>& num_bins_per_partition,
                                                              PartitionID rb_range_k,
                                                              HypernodeWeight max_bin_weight,
                                                              std::vector<PartitionID>&& partitions = {}) {
        PartitionID num_partitions = static_cast<PartitionID>(max_allowed_partition_weights.size());
        ALWAYS_ASSERT(num_bins_per_partition.size() == max_allowed_partition_weights.size(),
            "max_allowed_partition_weights and num_bins_per_partition have different sizes: "
            << V(max_allowed_partition_weights.size()) << "; " << V(num_bins_per_partition.size()));
        ALWAYS_ASSERT((num_partitions > 0) && (rb_range_k >= num_partitions),
            "num_partitions or rb_range_k invalid: " << V(num_partitions) << "; " << V(rb_range_k));
        ALWAYS_ASSERT(partitions.empty() || (partitions.size() == hypernodes.size()),
            "Size of fixed vertice partition IDs does not match the number of hypernodes: "
            << V(partitions.size()) << "; " << V(hypernodes.size()));
        ASSERT([&]() {
            for (size_t i = 1; i < hypernodes.size(); ++i) {
                if (hg.nodeWeight(hypernodes[i-1]) < hg.nodeWeight(hypernodes[i])) {
                    return false;
                }
            }
            return true;
        } (), "The hypernodes must be sorted in descending order of weight.");

        TwoLevelPacker<BPAlg> packer(rb_range_k, max_bin_weight);

        if (partitions.empty()) {
            partitions.resize(hypernodes.size(), -1);
        } else {
            // If fixed vertices are specified, calculate the total weight and extract all fixed vertices...
            HypernodeWeight total_weight = 0;
            std::vector<size_t> fixed_vertices;

            for (size_t i = 0; i < hypernodes.size(); ++i) {
                ASSERT(partitions[i] >= -1, "Invalid partition ID.");
                total_weight += hg.nodeWeight(hypernodes[i]);

                if (partitions[i] >= 0) {
                    ALWAYS_ASSERT(partitions[i] < num_partitions, "Invalid partition ID for node: "
                                  << partitions[i] << "; " << V(num_partitions));

                    fixed_vertices.push_back(i);
                }
            }

            // ...to pre-pack the fixed vertices with a first fit packing
            HypernodeWeight avg_bin_weight = (total_weight + rb_range_k - 1) / rb_range_k;
            PartitionID kbins_per_partition = rb_range_k / num_partitions;

            for (const size_t& index : fixed_vertices) {
                HypernodeWeight weight = hg.nodeWeight(hypernodes[index]);
                PartitionID part_id = partitions[index];

                PartitionID start_index = part_id * kbins_per_partition;
                PartitionID assigned_bin = start_index;

                if (kbins_per_partition > 1) {
                    for (PartitionID i = start_index; i < start_index + kbins_per_partition; ++i) {
                        HypernodeWeight current_bin_weight = packer.binWeight(i);

                        // The node is assigned to the first fitting bin or, if none fits, the smallest bin.
                        if (current_bin_weight + weight <= avg_bin_weight) {
                            assigned_bin = i;
                            break;
                        } else if (current_bin_weight < packer.binWeight(assigned_bin)) {
                            assigned_bin = i;
                        }
                    }
                }

                packer.addFixedVertex(assigned_bin, part_id, weight);
                partitions[index] = assigned_bin;
            }
        }

        // At the first level, a packing with k bins is calculated ....
        for (size_t i = 0; i < hypernodes.size(); ++i) {
            if (partitions[i] == -1) {
                HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);
                partitions[i] = packer.insertElement(weight);
            }
        }

        PartitionMapping mapping = packer.applySecondLevel(max_allowed_partition_weights, num_bins_per_partition).first;
        mapping.applyMapping(partitions);
        return partitions;
    }

    Hypergraph hypergraph;
};

TEST_F(BinPackingTest, BaseCases) {
  initializeWeights({});

  ASSERT_TRUE(two_level_packing<WorstFit>(hypergraph, {}, {0, 1}, {1, 1}, 2, 1).empty());
  ASSERT_TRUE(two_level_packing<WorstFit>(hypergraph, {}, {1, 0}, {2, 2}, 3, 0).empty());
  ASSERT_TRUE(two_level_packing<WorstFit>(hypergraph, {}, {1, 1}, {2, 2}, 4, 1).empty());

  ASSERT_TRUE(two_level_packing<FirstFit>(hypergraph, {}, {0, 1}, {1, 1}, 2, 1).empty());
  ASSERT_TRUE(two_level_packing<FirstFit>(hypergraph, {}, {1, 0}, {2, 2}, 3, 0).empty());
  ASSERT_TRUE(two_level_packing<FirstFit>(hypergraph, {}, {1, 1}, {2, 2}, 4, 1).empty());

  initializeWeights({1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0}, {1}, {1}, 1, 1);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0}, {1}, {1}, 1, 0);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = two_level_packing<FirstFit>(hypergraph, {0}, {1}, {1}, 1, 1);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = two_level_packing<FirstFit>(hypergraph, {0}, {1}, {1}, 1, 0);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  initializeWeights({1, 1});

  result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {1, 1}, {1, 1}, 2, 0);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1) || (result.at(0) == 1 && result.at(1) == 0), true);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1}, {1, 1}, {1, 1}, 2, 0);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1) || (result.at(0) == 1 && result.at(1) == 0), true);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {1, 1}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1) || (result.at(0) == 1 && result.at(1) == 0), true);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1}, {1, 1}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {2, 2}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1) || (result.at(0) == 1 && result.at(1) == 0), true);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1}, {2, 2}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
}

TEST_F(BinPackingTest, ReverseIndizes) {
  initializeWeights({1, 3, 2});

  auto result = two_level_packing<WorstFit>(hypergraph, {1, 2, 0}, {3, 3}, {1, 1}, 2, 3);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1 && result.at(2) == 1)
    || (result.at(0) == 1 && result.at(1) == 0 && result.at(2) == 0), true);

  result = two_level_packing<FirstFit>(hypergraph, {1, 2, 0}, {3, 3}, {1, 1}, 2, 3);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ((result.at(0) == 0 && result.at(1) == 1 && result.at(2) == 1)
    || (result.at(0) == 1 && result.at(1) == 0 && result.at(2) == 0), true);
}

TEST_F(BinPackingTest, WFTwoBinPacking) {
  initializeWeights({5, 4, 3, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {8, 8}, {1, 1}, 2, 0);
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
}

TEST_F(BinPackingTest, WFMultiBinPacking) {
  initializeWeights({10, 7, 3, 3, 3, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, {1, 1, 1, 1, 1, 1}, 6, 0);
  ASSERT_EQ(result.size(), 6);
  bool contained[6] = {false, false, false, false, false, false};
  for(size_t i = 0; i < result.size(); ++i) {
    contained[result.at(i)] = true;
  }
  for(size_t i = 0; i < 6; ++i) {
    ASSERT_TRUE(contained[i]);
  }

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {10, 10, 10}, {1, 1, 1}, 3, 0);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 2);
  ASSERT_EQ(result.at(3), 2);
  ASSERT_EQ(result.at(4), 2);
  ASSERT_EQ(result.at(5), 1);
}

TEST_F(BinPackingTest, WFTwoLevelPackingBase) {
  initializeWeights({4, 3, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {5, 5}, {3, 3}, 4, 0);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
}

TEST_F(BinPackingTest, WFTwoLevelPackingComplex) {
  initializeWeights({9, 7, 6, 4, 4, 4, 4, 3, 3, 3, 3});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {25, 25}, {2, 2}, 4, 0);
  // The packing in 4 bins:
  // (0.)9 (8.)3           12
  // (1.)7 (6.)4 (10.)3    14
  // (2.)6 (5.)4 (9.)3     13
  // (3.)4 (4.)4 (7.)3     11

  ASSERT_EQ(result.size(), 11);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
  ASSERT_EQ(result.at(5), 1);
  ASSERT_EQ(result.at(6), 0);
  ASSERT_EQ(result.at(7), 0);
  ASSERT_EQ(result.at(8), 1);
  ASSERT_EQ(result.at(9), 1);
  ASSERT_EQ(result.at(10), 0);
}

TEST_F(BinPackingTest, WFFixedVerticesBase) {
  initializeWeights({1, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, {1}, 1, 0, {-1, -1});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, {1}, 1, 0, {0, -1});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, {1}, 1, 0, {0, 0});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
}

TEST_F(BinPackingTest, WFFixedVerticesOneLevel) {
  initializeWeights({4, 3, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {6, 6}, {1, 1}, 2, 0, {0, 1, 0, 1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7, 7}, {2, 2, 2}, 3, 0, {0, 0, 2, 2});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 2);
  ASSERT_EQ(result.at(3), 2);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, {2, 2}, 2, 0, {1, -1, -1, -1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, {2, 2}, 2, 0, {-1, -1, 0, 0});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, {2, 2}, 2, 0, {-1, 1, 0, -1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
}

TEST_F(BinPackingTest, WFFixedVerticesTwoLevel) {
  initializeWeights({7, 5, 4, 3, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, {2, 2}, 4, 0, {0, 1, 0, 1, 1, 0});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);
  ASSERT_EQ(result.at(5), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {12, 12, 12}, {3, 3, 3}, 9, 0, {0, 2, 0, 2, 1, 0});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 2);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 2);
  ASSERT_EQ(result.at(4), 1);
  ASSERT_EQ(result.at(5), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, {2, 2}, 4, 0, {-1, 1, 1, -1, -1, -1});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
  ASSERT_EQ(result.at(5), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, {2, 2}, 4, 0, {1, -1, -1, -1, 0, 0});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
  ASSERT_EQ(result.at(5), 0);
}

TEST_F(BinPackingTest, UnevenBase) {
  initializeWeights({4, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {4, 4}, {2, 2}, 3, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {4, 4}, {2, 2}, 3, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {4, 4, 4}, {2, 2, 2}, 5, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 2);
  ASSERT_EQ(result.at(2), 1);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {4, 4, 4}, {2, 2, 2}, 5, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {1, 6}, {2, 2}, 3, 4, {});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {1, 6}, {2, 2}, 3, 2, {});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {1, 6}, {2, 2}, 3, 4, {0, -1, -1});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {1, 6}, {2, 2}, 3, 4, {0, -1, -1});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
}

TEST_F(BinPackingTest, WFUnevenAndFixed) {
  initializeWeights({5, 4, 3, 3, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 10}, {2, 2}, 3, 0, {-1, 0, -1, 1, -1});
  // The packing in 3 bins:
  // (F0)4 (3.)1           5
  // (F1)3 (2.)3           6
  // (1.)5                 5

  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 0);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 10}, {2, 2}, 3, 0, {1, -1, -1, 1, 1});
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);

  result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 6}, {2, 2}, 3, 0, {1, -1, -1, 1, 1});
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);
}

TEST_F(BinPackingTest, BinLimit) {
  initializeWeights({3, 2, 2, 1});

  auto result = two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, {1, 3}, 4, 2);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 1);

  result = two_level_packing<FirstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, {1, 3}, 4, 2);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 1);
}

TEST_F(BinPackingTest, ExtractNodes) {
  initializeWeights({});

  ASSERT_EQ(extract_nodes_with_descending_weight(hypergraph).size(), 0);

  initializeWeights({2, 1, 3, 6, 4, 5});

  auto result = extract_nodes_with_descending_weight(hypergraph);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 3);
  ASSERT_EQ(result.at(1), 5);
  ASSERT_EQ(result.at(2), 4);
  ASSERT_EQ(result.at(3), 2);
  ASSERT_EQ(result.at(4), 0);
  ASSERT_EQ(result.at(5), 1);
}

TEST_F(BinPackingTest, TreshholdPessBase) {
  initializeWeights({6, 5, 4, 3, 2, 1});
  std::pair<size_t, HypernodeWeight> zero_pair(0, 0);

  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 1, 0), zero_pair);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 3, 0), zero_pair);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 4, 1), zero_pair);

  std::pair<size_t, HypernodeWeight> two_zero(2, 0);
  std::pair<size_t, HypernodeWeight> one_one(1, 1);
  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5}, 4, 0), two_zero);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5}, 4, 1), one_one);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5}, 4, 2), one_two);
}

TEST_F(BinPackingTest, TreshholdOptBase) {
  initializeWeights({6, 5, 4, 3, 2, 1});
  std::pair<size_t, HypernodeWeight> zero_pair(0, 0);

  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 1, 0), zero_pair);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 3, 0), zero_pair);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 4, 1), zero_pair);

  std::pair<size_t, HypernodeWeight> two_zero(2, 0);
  std::pair<size_t, HypernodeWeight> one_one(1, 1);
  std::pair<size_t, HypernodeWeight> zero_two(0, 2);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5}, 4, 0), two_zero);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5}, 4, 1), one_one);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5}, 4, 2), zero_two);
}

TEST_F(BinPackingTest, TreshholdPessComplex) {
  initializeWeights({22, 18, 17, 8, 7, 3, 1, 1});

  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  std::pair<size_t, HypernodeWeight> five_two(5, 2);
  std::pair<size_t, HypernodeWeight> three_four(3, 4);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5, 6, 7}, 4, 2), one_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 2), five_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 4), three_four);
}

TEST_F(BinPackingTest, TreshholdOptComplex) {
  initializeWeights({22, 18, 17, 8, 7, 3, 1, 1});

  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  std::pair<size_t, HypernodeWeight> five_two(5, 2);
  std::pair<size_t, HypernodeWeight> three_four(3, 4);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5, 6, 7}, 4, 2), one_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 2), five_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 4), three_four);
}

TEST_F(BinPackingTest, TreshholdDiff) {
  initializeWeights({8, 6, 5, 5, 5, 2, 2, 2, 2});

  std::pair<size_t, HypernodeWeight> zero_two(0, 2);
  std::pair<size_t, HypernodeWeight> four_two(4, 2);
  std::pair<size_t, HypernodeWeight> zero_three(0, 3);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8}, 4, 2), zero_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8}, 4, 2), four_two);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 5, 6, 7, 8}, 4, 3), zero_three);
  ASSERT_EQ(calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 5, 6, 7, 8}, 4, 3), zero_three);
}

TEST_F(BinPackingTest, PrepackingPessimisticBase) {
  Context c;

  initializeWeights({4});
  createTestContext(c, {2, 2}, {2, 2}, {1, 1}, 2, 2, 2);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), false);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), false);

  initializeWeights({4, 4, 4, 4});
  createTestContext(c, {12, 12}, {8, 8}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 8);

  initializeWeights({4, 4, 4, 4});
  createTestContext(c, {12, 12}, {8, 8}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 8);
}

TEST_F(BinPackingTest, PrepackingPessimisticExtended) {
  Context c;

  initializeWeights({4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1});
  createTestContext(c, {13, 13}, {12, 12}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), false);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 8);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(1), 8);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 13);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 13);

  initializeWeights({4, 4, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1});
  createTestContext(c, {13, 13}, {12, 12}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), false);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 8);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(1), 8); 
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 13);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 13);

  initializeWeights({7, 7, 5, 5, 3, 2, 1, 1, 1});
  createTestContext(c, {17, 17}, {16, 16}, {2, 2}, 2, 4, 10);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), false);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 12);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(1), 12);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 19);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 19);

  initializeWeights({7, 7, 5, 5, 3, 2, 1, 1, 1});
  createTestContext(c, {17, 17}, {16, 16}, {2, 2}, 2, 4, 10);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(0), 14);
  ASSERT_EQ(hypergraph.fixedVertexPartWeight(1), 0);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 17);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 17);

  // tests for end of range failure
  initializeWeights({50, 50, 50, 1});
  createTestContext(c, {140, 140}, {76, 76}, {2, 2}, 2, 4, 75);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
  ASSERT_EQ(hypergraph.isFixedVertex(3), false);

  initializeWeights({50, 50, 26, 1});
  createTestContext(c, {140, 140}, {51, 51}, {2, 2}, 2, 4, 75);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), false);

  // tests handling of full partition edge case
  initializeWeights({4, 4, 4, 3, 3, 3, 3, 3, 3});
  createTestContext(c, {15, 15}, {15, 15}, {3, 3}, 2, 6, 6);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), true);
  ASSERT_EQ(hypergraph.isFixedVertex(5), true);
  ASSERT_EQ(hypergraph.isFixedVertex(6), true);
  ASSERT_EQ(hypergraph.isFixedVertex(7), true);
  ASSERT_EQ(hypergraph.isFixedVertex(8), true);

  // tests optimization in full partition edge case
  initializeWeights({5, 5, 3, 3, 3, 2, 1, 1});
  createTestContext(c, {12, 12}, {12, 12}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
  ASSERT_EQ(hypergraph.isFixedVertex(3), false);
  ASSERT_EQ(hypergraph.isFixedVertex(4), false);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 12);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 12);

  // invalid packing - test for resulting upper part weight
  initializeWeights({8, 8, 8, 7, 7});
  createTestContext(c, {19, 19}, {19, 19}, {2, 2}, 2, 4, 10);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), true);
  ASSERT_LE(c.initial_partitioning.upper_allowed_partition_weight[0], 20);
  ASSERT_LE(c.initial_partitioning.upper_allowed_partition_weight[1], 20);

  // impossible packing test
  initializeWeights({2, 2, 2, 2, 2});
  createTestContext(c, {5, 5}, {5, 5}, {2, 2}, 2, 4, 3);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), true);

  // optimization and last element
  initializeWeights({8, 8, 8, 8, 4, 4, 3});
  createTestContext(c, {26, 26}, {22, 22}, {2, 2}, 2, 4, 14);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), true);
  ASSERT_EQ(hypergraph.isFixedVertex(5), true);
  ASSERT_EQ(hypergraph.isFixedVertex(6), false);

  initializeWeights({8, 8, 8, 8, 4, 4, 3});
  createTestContext(c, {26, 26}, {22, 22}, {2, 2}, 2, 4, 14);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), true);
  ASSERT_EQ(hypergraph.isFixedVertex(5), true);
  ASSERT_EQ(hypergraph.isFixedVertex(6), false);

  // invalid optimization - test for regression
  initializeWeights({4, 4, 4, 4});
  createTestContext(c, {10, 10}, {8, 8}, {2, 2}, 2, 4, 6);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);

  initializeWeights({4, 4, 4, 4});
  createTestContext(c, {10, 10}, {8, 8}, {2, 2}, 2, 4, 6);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);

  initializeWeights({6, 4, 4, 4, 2, 2});
  createTestContext(c, {12, 12}, {11, 11}, {2, 2}, 2, 4, 7);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), true);
  ASSERT_EQ(hypergraph.isFixedVertex(4), false);
  ASSERT_EQ(hypergraph.isFixedVertex(5), false);

  // small final optimization
  initializeWeights({8, 8, 3, 2, 2});
  createTestContext(c, {12, 12}, {12, 12}, {2, 2}, 2, 4, 8);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 14);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 14);

  initializeWeights({10, 10, 3, 2, 1});
  createTestContext(c, {16, 16}, {13, 13}, {2, 2}, 2, 4, 10);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[0], 19);
  ASSERT_EQ(c.initial_partitioning.upper_allowed_partition_weight[1], 19);
}

TEST_F(BinPackingTest, PrepackingPessimisticUnequal) {
  Context c;

  // bigger partition edge case
  initializeWeights({5, 5, 5, 5, 4, 3, 2, 1});
  createTestContext(c, {30, 15}, {20, 10}, {4, 2}, 2, 6, 9);

  apply_prepacking_pessimistic<WorstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), true);
  ASSERT_EQ(hypergraph.isFixedVertex(3), false);

  // smaller partition edge case
  initializeWeights({6, 4, 4, 4, 4, 4});
  createTestContext(c, {10, 20}, {10, 20}, {2, 4}, 2, 6, 7);

  apply_prepacking_pessimistic<FirstFit>(hypergraph, c);
  ASSERT_EQ(hypergraph.isFixedVertex(0), true);
  ASSERT_EQ(hypergraph.isFixedVertex(1), true);
  ASSERT_EQ(hypergraph.isFixedVertex(2), false);
}

TEST_F(BinPackingTest, TwoLevelPackingFuzzingTest) {
  // Constraints:
  // 1) Nodes must be sorted by weight
  // 2) hypernodes.size() == partitions.size()
  // 3) max_allowed_partition_weights.size() == num_bins_per_partition.size()
  // 4) num_partitions > 0 && rb_range_k >= num_partitions
  // 5) sum(num_bins_per_partition >= rb_range_k)
  // 6) num_bins_per_partition[i] > 0 for all i

  const size_t NUM_ITERATIONS = 1000;
  const HypernodeID MAX_NUM_HNS = 500;
  const size_t MAX_SELECTION_STEP = 10;
  const HypernodeWeight MAX_HN_WEIGHT = 100;
  const HypernodeWeight WEIGHT_VARIANCE = 2;
  const PartitionID MAX_NUM_PARTITIONS = 8;
  const PartitionID MAX_RANGE_K = 32;
  const PartitionID NUM_BIN_VARIANCE = 4;

  std::mt19937 gen(std::random_device{}());
  int seed = std::uniform_int_distribution<>{ 0 }(gen);
  Randomize& random = Randomize::instance();
  random.setSeed(seed);
  std::cout << "Seed for TwoLevelPackingFuzzingTest: " << seed << std::endl;

  for (size_t i = 0; i < NUM_ITERATIONS; ++i) {
    // create random weights
    HypernodeID num_hns_total = random.getRandomInt(0, MAX_NUM_HNS);
    HypernodeWeight min_hn_weight = random.getRandomInt(0, MAX_HN_WEIGHT);
    PartitionID num_partitions = random.getRandomInt(1, MAX_NUM_PARTITIONS);
    PartitionID rb_range_k = random.getRandomInt(num_partitions, MAX_RANGE_K);
    HypernodeWeight max_partition_weight = random.getRandomInt(0, WEIGHT_VARIANCE * num_hns_total * MAX_HN_WEIGHT / num_partitions);
    HypernodeWeight max_bin_weight = random.getRandomInt(0, WEIGHT_VARIANCE * num_hns_total * MAX_HN_WEIGHT / rb_range_k);

    HypernodeWeightVector weights;
    weights.reserve(num_hns_total);
    for (HypernodeID j = 0; j < num_hns_total; ++j) {
      weights.push_back(random.getRandomInt(min_hn_weight, MAX_HN_WEIGHT));
    }
    initializeWeights(weights);
    std::vector<HypernodeID> all_nodes = extract_nodes_with_descending_weight(hypergraph);

    // choose random subset of the nodes
    std::vector<HypernodeID> nodes;
    size_t max_step = random.getRandomInt(1, MAX_SELECTION_STEP);
    for (size_t j = 0; j < all_nodes.size(); j += random.getRandomInt(1, max_step)) {
      nodes.push_back(all_nodes[j]);
    }

    // choose fixed vertices
    float fixed_vertex_probability = random.getRandomFloat(0, 1);
    std::vector<PartitionID> partitions;
    partitions.reserve(nodes.size());
    for (size_t j = 0; j < nodes.size(); ++j) {
      if (random.getRandomFloat(0, 1) < fixed_vertex_probability) {
        partitions.push_back(random.getRandomInt(0, num_partitions - 1));
      } else {
        partitions.push_back(-1);
      }
    }

    // choose partition weights
    HypernodeWeight min_partition_weight = random.getRandomInt(0, max_partition_weight);
    std::vector<HypernodeWeight> max_allowed_partition_weights;
    max_allowed_partition_weights.reserve(num_partitions);
    for (PartitionID j = 0; j < num_partitions; ++j) {
      max_allowed_partition_weights.push_back(random.getRandomInt(min_partition_weight, max_partition_weight));
    }

    // choose num bins
    PartitionID max_num_bins = random.getRandomInt(1, std::max(NUM_BIN_VARIANCE * rb_range_k / num_partitions, 1));
    PartitionID min_num_bins = random.getRandomInt(1, max_num_bins);
    std::vector<PartitionID> num_bins_per_partition;
    num_bins_per_partition.reserve(num_partitions);
    for (PartitionID j = 0; j < num_partitions; ++j) {
      num_bins_per_partition.push_back(random.getRandomInt(min_num_bins, max_num_bins));
    }

    while (std::accumulate(num_bins_per_partition.cbegin(), num_bins_per_partition.cend(), 0) < rb_range_k) {
      for(PartitionID& num : num_bins_per_partition) {
        ++num;
      }
    }

    std::vector<PartitionID> partitions_copy(partitions);
    auto result = two_level_packing<WorstFit>(hypergraph, nodes, max_allowed_partition_weights,
      num_bins_per_partition, rb_range_k, max_bin_weight, std::move(partitions));
    ASSERT_EQ(result.size(), nodes.size());

    result = two_level_packing<WorstFit>(hypergraph, nodes, max_allowed_partition_weights,
      num_bins_per_partition, rb_range_k, max_bin_weight, std::move(partitions_copy));
    ASSERT_EQ(result.size(), nodes.size());
  }
}
}  // namespace bin_packing
}  // namespace kahypar
