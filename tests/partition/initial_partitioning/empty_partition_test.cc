/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Nikolai Maas <nikolai.maas@student.kit.edu>
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
#include <vector>

#include "gmock/gmock.h"

#include "kahypar/partition/initial_partitioning/bfs_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/label_propagation_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/greedy_hypergraph_growing_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/random_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/bin_packing_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/partition/initial_partitioning/policies/ip_start_node_selection_policy.h"
#include "kahypar/partition/initial_partitioning/policies/ip_greedy_queue_selection_policy.h"
#include "kahypar/partition/initial_partitioning/policies/ip_gain_computation_policy.h"
#include "kahypar/partition/metrics.h"

using ::testing::Eq;
using ::testing::Test;

namespace kahypar {
using BFSInitialPartitionerBFS = BFSInitialPartitioner<BFSStartNodeSelectionPolicy<> >;
using LPInitialPartitionerBFS_FM =
  LabelPropagationInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                    FMGainComputationPolicy>;
using GHGInitialPartitionerBFS_FM_SEQ =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            FMGainComputationPolicy,
                                            SequentialQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_FM_GLO =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            FMGainComputationPolicy,
                                            GlobalQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_FM_RND =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            FMGainComputationPolicy,
                                            RoundRobinQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXP_SEQ =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxPinGainComputationPolicy,
                                            SequentialQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXP_GLO =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxPinGainComputationPolicy,
                                            GlobalQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXP_RND =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxPinGainComputationPolicy,
                                            RoundRobinQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXN_SEQ =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxNetGainComputationPolicy,
                                            SequentialQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXN_GLO =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxNetGainComputationPolicy,
                                            GlobalQueueSelectionPolicy>;
using GHGInitialPartitionerBFS_MAXN_RND =
  GreedyHypergraphGrowingInitialPartitioner<BFSStartNodeSelectionPolicy<>,
                                            MaxNetGainComputationPolicy,
                                            RoundRobinQueueSelectionPolicy>;

void initializeContext(Context& context, PartitionID k,
                       double epsilon, HypernodeWeight hypergraph_weight) {
  context.initial_partitioning.k = k;
  context.partition.k = k;
  context.partition.epsilon = epsilon;
  context.initial_partitioning.unassigned_part = -1;
  context.initial_partitioning.refinement = false;
  context.initial_partitioning.nruns = 20;
  context.initial_partitioning.upper_allowed_partition_weight.resize(
    context.initial_partitioning.k);
  context.initial_partitioning.perfect_balance_partition_weight.resize(
    context.initial_partitioning.k);
  context.initial_partitioning.num_bins_per_partition.resize(
    context.initial_partitioning.k, 1);
  context.partition.max_part_weights.resize(context.partition.k);
  context.partition.perfect_balance_part_weights.resize(context.partition.k);
  for (int i = 0; i < context.initial_partitioning.k; i++) {
    context.initial_partitioning.perfect_balance_partition_weight[i] = ceil(
      hypergraph_weight
      / static_cast<double>(context.initial_partitioning.k));
    context.initial_partitioning.upper_allowed_partition_weight[i] = ceil(
      hypergraph_weight
      / static_cast<double>(context.initial_partitioning.k)) * (1.0 + context.partition.epsilon);
  }
  for (int i = 0; i < context.initial_partitioning.k; ++i) {
    context.partition.perfect_balance_part_weights[i] =
      context.initial_partitioning.perfect_balance_partition_weight[i];
    context.partition.max_part_weights[i] = context.initial_partitioning.upper_allowed_partition_weight[i];
  }
}

class ATestPartitioner : public Test {
 public:
  ATestPartitioner() :
    partitioner(nullptr),
    hypergraph(2, 1, HyperedgeIndexVector { 0, 2 },
               HyperedgeVector { 0, 1 }),
    context() {
    PartitionID k = 2;
    double epsilon = 1.1;
    initializeContext(context, k, epsilon, 2);
  }

  void testPartitioningRepeated() {
    for (size_t i = 0; i < 10; ++i) {
      hypergraph.resetPartitioning();
      partitioner->partition();

      ASSERT_EQ(hypergraph.partWeight(0), 1);
      ASSERT_EQ(hypergraph.partWeight(1), 1);
    }
  }

  template <class PartitionerT>
  void initializePartitioner() {
    partitioner = std::make_shared<PartitionerT>(hypergraph, context);
  }

  std::shared_ptr<IInitialPartitioner> partitioner;
  Hypergraph hypergraph;
  Context context;
};

TEST_F(ATestPartitioner, Random) {
  initializePartitioner<RandomInitialPartitioner >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, BFS) {
  initializePartitioner<BFSInitialPartitionerBFS >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, LP) {
  initializePartitioner<LPInitialPartitionerBFS_FM >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, BinPacking) {
  initializePartitioner<BinPackingInitialPartitioner >();
  context.initial_partitioning.bp_algo = BinPackingAlgorithm::worst_fit;
  testPartitioningRepeated();
  context.initial_partitioning.bp_algo = BinPackingAlgorithm::first_fit;
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedySeq) {
  initializePartitioner<GHGInitialPartitionerBFS_FM_SEQ >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyGlo) {
  initializePartitioner<GHGInitialPartitionerBFS_FM_GLO >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyRnd) {
  initializePartitioner<GHGInitialPartitionerBFS_FM_RND >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxpSeq) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXP_SEQ >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxpGlo) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXP_GLO >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxpRnd) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXP_RND >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxnSeq) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXN_SEQ >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxnGlo) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXN_GLO >();
  testPartitioningRepeated();
}

TEST_F(ATestPartitioner, GreedyMaxnRnd) {
  initializePartitioner<GHGInitialPartitionerBFS_MAXN_RND >();
  testPartitioningRepeated();
}
}  // namespace kahypar
