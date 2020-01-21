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

#include <memory>

#include "gmock/gmock.h"

#include "kahypar/partition/bin_packing.h"

using ::testing::Eq;
using ::testing::Test;

namespace kahypar {
using bin_packing::WorstFit;
using bin_packing::FirstFit;

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

    Hypergraph hypergraph;
};

TEST_F(BinPackingTest, BaseCases) {
  initializeWeights({});

  ASSERT_TRUE(bin_packing::two_level_packing<WorstFit>(hypergraph, {}, {0, 0}, 1, 0).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<WorstFit>(hypergraph, {}, {0, 1}, 2, 1).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<WorstFit>(hypergraph, {}, {1, 0}, 3, 0).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<WorstFit>(hypergraph, {}, {1, 1}, 4, 1).empty());

  ASSERT_TRUE(bin_packing::two_level_packing<FirstFit>(hypergraph, {}, {0, 0}, 1, 0).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<FirstFit>(hypergraph, {}, {0, 1}, 2, 1).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<FirstFit>(hypergraph, {}, {1, 0}, 3, 0).empty());
  ASSERT_TRUE(bin_packing::two_level_packing<FirstFit>(hypergraph, {}, {1, 1}, 4, 1).empty());

  initializeWeights({1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0}, {1}, 1, 1);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0}, {1}, 1, 0);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0}, {1}, 1, 1);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0}, {1}, 1, 0);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result.at(0), 0);

  initializeWeights({1, 1});

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {1, 1}, 2, 0);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1}, {1, 1}, 2, 0);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1}, {1, 1}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {1, 1}, 1, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1}, {1, 1}, 1, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {2, 2}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1}, {2, 2}, 2, 2);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
}

TEST_F(BinPackingTest, ReverseIndizes) {
  initializeWeights({1, 3, 2});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {1, 2, 0}, {3, 3}, 2, 3);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {1, 2, 0}, {3, 3}, 2, 3);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
}

TEST_F(BinPackingTest, WFTwoBinPacking) {
  initializeWeights({5, 4, 3, 2, 1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {8, 8}, 2, 0);
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
}

TEST_F(BinPackingTest, WFMultiBinPacking) {
  initializeWeights({10, 7, 3, 3, 3, 1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {0, 0, 0, 0, 0, 0}, 6, 0);
  ASSERT_EQ(result.size(), 6);
  bool contained[6] = {false, false, false, false, false, false};
  for(size_t i = 0; i < result.size(); ++i) {
    contained[result.at(i)] = true;
  }
  for(size_t i = 0; i < 6; ++i) {
    ASSERT_TRUE(contained[i]);
  }

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {10, 10, 10}, 3, 0);
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

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {5, 5}, 4, 0);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
}

TEST_F(BinPackingTest, WFTwoLevelPackingComplex) {
  initializeWeights({9, 7, 6, 4, 4, 4, 4, 3, 3, 3, 3});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {25, 25}, 4, 0);
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

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, 1, 0, {-1, -1});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, 1, 0, {0, -1});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1}, {2}, 1, 0, {0, 0});
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
}

TEST_F(BinPackingTest, WFFixedVerticesOneLevel) {
  initializeWeights({4, 3, 2, 1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {6, 6}, 2, 0, {0, 1, 0, 1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7, 7}, 3, 0, {0, 0, 2, 2});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 2);
  ASSERT_EQ(result.at(3), 2);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, 2, 0, {1, -1, -1, -1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, 2, 0, {-1, -1, 0, 0});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3}, {7, 7}, 2, 0, {-1, 1, 0, -1});
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
}

TEST_F(BinPackingTest, WFFixedVerticesTwoLevel) {
  initializeWeights({7, 5, 4, 3, 2, 1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, 4, 0, {0, 1, 0, 1, 1, 0});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);
  ASSERT_EQ(result.at(5), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {12, 12, 12}, 9, 0, {0, 2, 0, 2, 1, 0});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 2);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 2);
  ASSERT_EQ(result.at(4), 1);
  ASSERT_EQ(result.at(5), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, 4, 0, {-1, 1, 1, -1, -1, -1});
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
  ASSERT_EQ(result.at(3), 0);
  ASSERT_EQ(result.at(4), 0);
  ASSERT_EQ(result.at(5), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4, 5}, {15, 15}, 4, 0, {1, -1, -1, -1, 0, 0});
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

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {4, 4}, 3, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {4, 4}, 3, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {4, 4, 4}, 5, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 2);
  ASSERT_EQ(result.at(2), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {4, 4, 4}, 5, 4);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {1, 6}, 3, 4, {});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {1, 6}, 3, 2, {});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 0);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2}, {1, 6}, 3, 4, {0, -1, -1});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);

  result = bin_packing::two_level_packing<FirstFit>(hypergraph, {0, 1, 2}, {1, 6}, 3, 4, {0, -1, -1});
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result.at(0), 0);
  ASSERT_EQ(result.at(1), 1);
  ASSERT_EQ(result.at(2), 1);
}

TEST_F(BinPackingTest, WFUnevenAndFixed) {
  initializeWeights({5, 4, 3, 3, 1});

  auto result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 10}, 3, 0, {-1, 0, -1, 1, -1});
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

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 10}, 3, 0, {1, -1, -1, 1, 1});
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);

  result = bin_packing::two_level_packing<WorstFit>(hypergraph, {0, 1, 2, 3, 4}, {10, 6}, 3, 0, {1, -1, -1, 1, 1});
  ASSERT_EQ(result.size(), 5);
  ASSERT_EQ(result.at(0), 1);
  ASSERT_EQ(result.at(1), 0);
  ASSERT_EQ(result.at(2), 0);
  ASSERT_EQ(result.at(3), 1);
  ASSERT_EQ(result.at(4), 1);
}

TEST_F(BinPackingTest, ExtractNodes) {
  initializeWeights({});

  ASSERT_EQ(bin_packing::extract_nodes_with_descending_weight(hypergraph).size(), 0);

  initializeWeights({2, 1, 3, 6, 4, 5});

  auto result = bin_packing::extract_nodes_with_descending_weight(hypergraph);
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

  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 1, 0), zero_pair);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 3, 0), zero_pair);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {}, 4, 1), zero_pair);

  std::pair<size_t, HypernodeWeight> two_zero(2, 0);
  std::pair<size_t, HypernodeWeight> one_one(1, 1);
  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5}, 4, 0), two_zero);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5}, 4, 1), one_one);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5}, 4, 2), one_two);
}

TEST_F(BinPackingTest, TreshholdOptBase) {
  initializeWeights({6, 5, 4, 3, 2, 1});
  std::pair<size_t, HypernodeWeight> zero_pair(0, 0);

  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 1, 0), zero_pair);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 3, 0), zero_pair);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {}, 4, 1), zero_pair);

  std::pair<size_t, HypernodeWeight> two_zero(2, 0);
  std::pair<size_t, HypernodeWeight> one_one(1, 1);
  std::pair<size_t, HypernodeWeight> zero_two(0, 2);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5}, 4, 0), two_zero);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5}, 4, 1), one_one);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5}, 4, 2), zero_two);
}

TEST_F(BinPackingTest, TreshholdPessComplex) {
  initializeWeights({22, 18, 17, 8, 7, 3, 1, 1});

  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  std::pair<size_t, HypernodeWeight> five_two(5, 2);
  std::pair<size_t, HypernodeWeight> three_four(3, 4);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 5, 6, 7}, 4, 2), one_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 2), five_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 4), three_four);
}

TEST_F(BinPackingTest, TreshholdOptComplex) {
  initializeWeights({22, 18, 17, 8, 7, 3, 1, 1});

  std::pair<size_t, HypernodeWeight> one_two(1, 2);
  std::pair<size_t, HypernodeWeight> five_two(5, 2);
  std::pair<size_t, HypernodeWeight> three_four(3, 4);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 5, 6, 7}, 4, 2), one_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 2), five_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7}, 4, 4), three_four);
}

TEST_F(BinPackingTest, TreshholdDiff) {
  initializeWeights({8, 6, 5, 5, 5, 2, 2, 2, 2});

  std::pair<size_t, HypernodeWeight> zero_two(0, 2);
  std::pair<size_t, HypernodeWeight> four_two(4, 2);
  std::pair<size_t, HypernodeWeight> zero_three(0, 3);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8}, 4, 2), zero_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 2, 3, 4, 5, 6, 7, 8}, 4, 2), four_two);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_optimistic(hypergraph, {0, 1, 5, 6, 7, 8}, 4, 3), zero_three);
  ASSERT_EQ(bin_packing::calculate_heavy_nodes_treshhold_pessimistic(hypergraph, {0, 1, 5, 6, 7, 8}, 4, 3), zero_three);
}

}  // namespace kahypar
