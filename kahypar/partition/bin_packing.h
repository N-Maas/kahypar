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

#include <vector>

#include <boost/range/adaptor/reversed.hpp>

#include "kahypar/partition/context.h"
#include "kahypar/datastructure/binary_heap.h"

namespace kahypar {

namespace bin_packing {
    using kahypar::ds::BinaryMinHeap;
    using kahypar::ds::BinaryMaxHeap;

    // Bin packing algorithm classes
    //
    // The following operations must be available for a type BPAlg:
    // 1)
    //      PartitionID num_bins = ...;
    //      BPAlg alg(num_bins);
    //
    // 2)
    //      PartitionID bin = ...;
    //      HypernodeWeight weight = ...;
    //      alg.addWeight(bin, weight);
    //
    // 3)
    //      PartitionID bin = ...;
    //      PartitionID partition = ...;
    //      alg.setPartition(bin, partition);
    //
    // 4)
    //      HypernodeWeight weight = ...;
    //      PartitionID resulting_bin = alg.insertNode(weight);
    //
    // 5)
    //      std::vector<PartitionID> asc_bins = alg.extractBinsAscending();
    //
    //      * Where asc_bins contains the bin ids in ascending order of bin weight.
    //
    // 6)
    //      ParititionID bin = ...;
    //      bool isFixed = alg.isFixedBin(bin);
    //
    // 7)
    //      ParititionID bin = ...;
    //      PartitionId partition = alg.binPartition(bin);
    //
    // 8)
    //      ParititionID bin = ...;
    //      HypernodeWeight weight = alg.binWeight(bin);
    //
    // It can be assumed that 4 is only called once.

    class WorstFit {
        public:
            WorstFit(PartitionID num_bins, HypernodeWeight max) :
                bin_queue(num_bins),
                partitions(num_bins, -1),
                weights() {
                for (PartitionID i = 0; i < num_bins; ++i) {
                    bin_queue.push(i, 0);
                }
            }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                bin_queue.increaseKeyBy(bin, weight);
            }

            void setPartition(PartitionID bin, PartitionID partition) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");
                ASSERT(partitions[bin] == -1 || partitions[bin] == partition,
                       "Bin already assigned to other partition");

                partitions[bin] = partition;
            }

            PartitionID insertNode(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");

                // assign node to bin with lowest weight
                PartitionID bin = bin_queue.top();
                bin_queue.increaseKeyBy(bin, weight);
                return bin;
            }

            bool isFixedBin(PartitionID bin) {
                return binPartition(bin) != -1;
            }

            PartitionID binPartition(PartitionID bin) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                return partitions[bin];
            }

            HypernodeWeight binWeight(PartitionID bin) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                return weights.empty() ? bin_queue.getKey(bin) : weights[bin];
            }

            std::vector<PartitionID> extractBinsAscending() {
                ASSERT(bin_queue.size() == partitions.size());

                std::vector<PartitionID> asc_bins;
                asc_bins.reserve(partitions.size());
                weights.resize(partitions.size(), 0);

                while (bin_queue.size() > 0) {
                    PartitionID bin = bin_queue.top();
                    ASSERT(bin >= 0 && bin < weights.size() && weights[bin] == 0);

                    weights[bin] = bin_queue.topKey();
                    asc_bins.push_back(bin);
                    bin_queue.pop();
                }
                return std::move(asc_bins);
            }

        private:
            BinaryMinHeap<PartitionID, HypernodeWeight> bin_queue;
            std::vector<PartitionID> partitions;
            std::vector<HypernodeWeight> weights;
    };

    class FirstFit {
        public:
            FirstFit(PartitionID num_bins, HypernodeWeight max) :
                max_bin_weight(max),
                partitions(num_bins, -1),
                weights(num_bins, 0) { }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                weights[bin] += weight;
            }

            void setPartition(PartitionID bin, PartitionID partition) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");
                ASSERT(partitions[bin] == -1 || partitions[bin] == partition,
                       "Bin already assigned to other partition");

                partitions[bin] = partition;
            }

            PartitionID insertNode(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");

                size_t assigned_bin = 0;
                for (size_t i = 0; i < weights.size(); ++i) {
                    // The node is assigned to the first fitting bin or, if none fits, the smallest bin.
                    if (weights[i] + weight <= max_bin_weight) {
                        assigned_bin = i;
                        break;
                    } else if (weights[i] < weights[assigned_bin]) {
                        assigned_bin = i;
                    }
                }
                weights[assigned_bin] += weight;
                return assigned_bin;
            }

            bool isFixedBin(PartitionID bin) {
                return binPartition(bin) != -1;
            }

            PartitionID binPartition(PartitionID bin) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                return partitions[bin];
            }

            HypernodeWeight binWeight(PartitionID bin) {
                ASSERT(bin >= 0 && bin < partitions.size(), "Invalid bin id.");

                return weights[bin];
            }

            std::vector<PartitionID> extractBinsAscending() {
                ASSERT(partitions.size() == weights.size());

                std::vector<PartitionID> asc_bins;
                asc_bins.reserve(partitions.size());
                for (PartitionID i = 0; i < partitions.size(); ++i) {
                    asc_bins.push_back(i);
                }

                std::sort(asc_bins.begin(), asc_bins.end(), [&](PartitionID i, PartitionID j) {
                    return binWeight(i) <= binWeight(j);
                });
                return std::move(asc_bins);
            }

        private:
            HypernodeWeight max_bin_weight;
            std::vector<PartitionID> partitions;
            std::vector<HypernodeWeight> weights;
    };


    /**
     * The hypernodes must be sorted in descending order of weight. 
     * The partitions parameter can be used to specify fixed vertices.
     * With reserved_weights, weights for partitions can be added (e.g. to handle uneven k).
     *
     * Attention: the algorithm assumus that the partitions can be balanced
     * with rb_range_k bins of equal weight.
     */
    template< class BPAlg >
    static inline std::vector<PartitionID> two_level_packing(const Hypergraph& hg,
                                                              const std::vector<HypernodeID>& hypernodes,
                                                              PartitionID num_partitions,
                                                              PartitionID rb_range_k,
                                                              HypernodeWeight max_bin_weight,
                                                              std::vector<PartitionID>&& partitions = {},
                                                              const std::vector<HypernodeWeight>& reserved_weights = {}) {
        ALWAYS_ASSERT((num_partitions > 0) && (rb_range_k > 0),
            "num_partitions and rb_range_k must be positive.");
        ALWAYS_ASSERT(partitions.empty() || (partitions.size() == hypernodes.size()),
            "Size of fixed vertice partition IDs does not match the number of hypernodes.");
        ALWAYS_ASSERT(reserved_weights.empty() || (reserved_weights.size() == static_cast<size_t>(num_partitions)),
            "Size of reserved weights does not match the number of partitions.");
        ALWAYS_ASSERT(reserved_weights.empty() || (num_partitions != rb_range_k),
            "Reserved weights are currently not used for one level.");
        ASSERT([&]() {
            for (size_t i = 1; i < hypernodes.size(); ++i) {
                if (hg.nodeWeight(hypernodes[i-1]) < hg.nodeWeight(hypernodes[i])) {
                    return false;
                }
            }
            return true;
        } (), "The hypernodes must be sorted in descending order of weight.");

        BPAlg bin_packer(rb_range_k, max_bin_weight);

        if (partitions.empty()) {
            partitions.resize(hypernodes.size(), -1);
        } else {
            // If fixed vertices are specified, calculate the total weight and extract all fixed vertices...
            HypernodeWeight total_weight = 0;
            std::vector<size_t> fixed_vertices;

            for (size_t i = 0; i < hypernodes.size(); ++i) {
                total_weight += hg.nodeWeight(hypernodes[i]);

                if (partitions[i] >= 0) {
                    ALWAYS_ASSERT(partitions[i] < num_partitions, "Invalid partition ID for node.");

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
                        HypernodeWeight current_bin_weight = bin_packer.binWeight(i);

                        // The node is assigned to the first fitting bin or, if none fits, the smallest bin.
                        if (current_bin_weight + weight <= avg_bin_weight) {
                            assigned_bin = i;
                            break;
                        } else if (current_bin_weight < bin_packer.binWeight(assigned_bin)) {
                            assigned_bin = i;
                        }
                    }
                }

                bin_packer.addWeight(assigned_bin, weight);
                bin_packer.setPartition(assigned_bin, part_id);
                partitions[index] = assigned_bin;
            }
        }

        // At the first level, a packing with k bins is calculated ....
        for (size_t i = 0; i < hypernodes.size(); ++i) {
            if (partitions[i] == -1) {
                HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);
                partitions[i] = bin_packer.insertNode(weight);
            }
        }

        // ... and at the second level, the resulting k bins are packed into the final bins
        if (rb_range_k != num_partitions) {
            // initialize result queue
            BinaryMinHeap<PartitionID, HypernodeWeight> result_queue(num_partitions);
            for (PartitionID i = 0; i < num_partitions; ++i) {
                result_queue.push(i, 0);
            }

            if (!reserved_weights.empty()) {
                ASSERT(reserved_weights.size() == result_queue.size());

                for (size_t i = 0; i < reserved_weights.size(); ++i) {
                    result_queue.increaseKeyBy(i, reserved_weights[i]);
                }
            }

            // read bins from kbin queue
            std::vector<PartitionID> kbin_ascending_ordering = bin_packer.extractBinsAscending();

            ASSERT(kbin_ascending_ordering.size() == rb_range_k);

            for (const PartitionID& bin : kbin_ascending_ordering) {
                if (bin_packer.isFixedBin(bin)) {
                    result_queue.increaseKeyBy(bin_packer.binPartition(bin), bin_packer.binWeight(bin));
                }
            }

            // iteration must be in reverse order, because the vector is ascending
            for (const PartitionID& bin : boost::adaptors::reverse(kbin_ascending_ordering)) {
                // if not fixed, assign bin to partition with lowest weight
                if (!bin_packer.isFixedBin(bin)) {
                    PartitionID partition = result_queue.top();
                    ASSERT((partition >= 0) && (partition < num_partitions));

                    result_queue.increaseKeyBy(partition, bin_packer.binWeight(bin));
                    bin_packer.setPartition(bin, partition);
                }
            }

            // remap the partition IDs of the nodes
            for (PartitionID& id : partitions) {
                id = bin_packer.binPartition(id);
                ASSERT((id >= 0) && (id < rb_range_k));
            }
        }

        return partitions;
    }

    /**
     * Returns the current hypernodes sorted in descending order of weight.
     */
    static inline std::vector<HypernodeID> extract_nodes_with_descending_weight(const Hypergraph& hg) {
        std::vector<HypernodeID> nodes;
        nodes.reserve(hg.currentNumNodes());

        for (const HypernodeID& hn : hg.nodes()) {
            nodes.push_back(hn);
        }
        ASSERT(hg.currentNumNodes() == nodes.size());

        std::sort(nodes.begin(), nodes.end(), [&hg](HypernodeID a, HypernodeID b) {
            return hg.nodeWeight(a) > hg.nodeWeight(b);
        });

        return nodes;
    }

    /**
     * Calculates an index i for a set of n hypernodes with descending order of weight
     * s.t. for any (balanced) bipartition of i and the remaining lighter hypernodes, the
     * biparition does not worsen the imbalance of lower bisection levels by more
     * then allowed_imbalance. Returns i and the resulting imbalance.
     *
     * The hypernodes must be sorted in descending order of weight.
     *
     * i is calculated in the following way: Let i the index of a hypernode and j the
     * minimum index s.t. sum_{l=i}^j (a_l) >= sum_{l=j+1}^n (a_l). Choose the smallest i
     * s.t. max{a_l - 1/k * sum_{m=l}^j (a_l) | l in {i, ..., j}} is within allowed_imbalance.
     *
     * The algorithm is implemented in O(n*log(n)) by reusing the values of the previous
     * step and tracking the curent maximum with a priority queue.
     */
    static inline std::pair<size_t, HypernodeWeight> calculate_heavy_nodes_treshhold_pessimistic(const Hypergraph& hg,
                                                         const std::vector<HypernodeID>& hypernodes,
                                                         const PartitionID& rb_range_k,
                                                         const HypernodeWeight& allowed_imbalance) {
        ALWAYS_ASSERT(rb_range_k > 0, "rb_range_k must be positive.");
        ASSERT([&]() {
            for (size_t i = 1; i < hypernodes.size(); ++i) {
                if (hg.nodeWeight(hypernodes[i-1]) < hg.nodeWeight(hypernodes[i])) {
                    return false;
                }
            }
            return true;
        } (), "The hypernodes must be sorted in descending order of weight.");

        if (rb_range_k <= 2) {
            return {0, 0};
        }

        PartitionID k = (rb_range_k + 1) / 2;
        HypernodeWeight total_remaining_weight = 0;
        for (const HypernodeID& node : hypernodes) {
            total_remaining_weight += hg.nodeWeight(node);
        }

        BinaryMaxHeap<size_t, HypernodeWeight> ratings(hypernodes.size());
        size_t i = 0;
        size_t j = 0;
        HypernodeWeight sum_of_range = 0;
        HypernodeWeight sum_of_previous = 0;
        HypernodeWeight imbalance = 0;

        while (i < hypernodes.size()) {
            // increase j
            while (sum_of_range < (total_remaining_weight + 1) / 2) {
                ASSERT(j < hypernodes.size());

                HypernodeWeight weight = hg.nodeWeight(hypernodes[j]);
                ratings.push(j, sum_of_previous + k * weight);
                sum_of_range += weight;
                sum_of_previous += weight;
                ++j;
            }

            // test whether current i is sufficient
            imbalance = (ratings.topKey() - sum_of_previous + k - 1) / k;
            HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);
            if (imbalance <= allowed_imbalance || weight <= 0) {
                break;
            }

            ASSERT(j >= i);

            // increase i
            ratings.remove(i);
            total_remaining_weight -= weight;
            sum_of_range -= weight;
            ++i;
        }

        ASSERT(imbalance >= 0, "error calculating threshhold - negative imbalance: " << imbalance);

        return {i, i == hypernodes.size() ? 0 : imbalance};
    }

    /**
     * Calculates an index i for a set of n hypernodes with descending order of weight
     * s.t. for any (balanced) bipartition of i and the remaining lighter hypernodes, the
     * biparition does not worsen the imbalance of lower bisection levels by more
     * then allowed_imbalance. Returns i and the resulting imbalance.
     *
     * The hypernodes must be sorted in descending order of weight.
     *
     * TODO
     */
    static inline std::pair<size_t, HypernodeWeight> calculate_heavy_nodes_treshhold_optimistic(const Hypergraph& hg,
                                                         const std::vector<HypernodeID>& hypernodes,
                                                         const PartitionID& rb_range_k,
                                                         const HypernodeWeight& allowed_imbalance) {
        ALWAYS_ASSERT(rb_range_k > 0, "rb_range_k must be positive.");
        ASSERT([&]() {
            for (size_t i = 1; i < hypernodes.size(); ++i) {
                if (hg.nodeWeight(hypernodes[i-1]) < hg.nodeWeight(hypernodes[i])) {
                    return false;
                }
            }
            return true;
        } (), "The hypernodes must be sorted in descending order of weight.");

        if (rb_range_k <= 2) {
            return {0, 0};
        }

        HypernodeWeight current_lower_sum = 0;
        HypernodeWeight max_imbalance = 0;
        HypernodeWeight imbalance = 0;

        for (int i = hypernodes.size() - 1; i >= 0; --i) {
            HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);
            current_lower_sum += weight;
            imbalance = weight - current_lower_sum / rb_range_k;

            if (imbalance > allowed_imbalance) {
                return {i + 1, max_imbalance};
            }

            max_imbalance = std::max(imbalance, max_imbalance);
        }

        return {0, max_imbalance};
    }

    static inline std::vector<PartitionID> apply_bin_packing_to_nodes(const Hypergraph& hg,
                                                                      const Context& context,
                                                                      const std::vector<HypernodeID>& nodes) {
        std::vector<PartitionID> partitions(nodes.size(), -1);
        if (hg.containsFixedVertices()) {
            for (size_t i = 0; i < nodes.size(); ++i) {
                HypernodeID hn = nodes[i];

                if (hg.isFixedVertex(hn)) {
                    partitions[i] = hg.fixedVertexPartID(hn);
                }
            }
        }

        PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;
        std::vector<HypernodeWeight> reserved_weights;

        // handle uneven k
        if ((rb_range_k % context.initial_partitioning.k) != 0) {
            reserved_weights.reserve(context.initial_partitioning.k);

            for (PartitionID p = 0; p < context.initial_partitioning.k; ++p) {
                reserved_weights.push_back(context.initial_partitioning.upper_allowed_partition_weight[p]);
            }
            ASSERT(!reserved_weights.empty());
            HypernodeWeight max = *std::max_element(reserved_weights.cbegin(), reserved_weights.cend());

            // the bin packing algorithm requires the complement of the allowed weights
            for (HypernodeWeight& w : reserved_weights) {
                w = max - w;
            }
        }

        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * context.initial_partitioning.bin_epsilon);
        // TODO
        partitions = two_level_packing<WorstFit>(hg, nodes, context.initial_partitioning.k, rb_range_k,
                                       max_bin_weight, std::move(partitions), reserved_weights);

        ASSERT(nodes.size() == partitions.size());
        ASSERT([&]() {
            for (size_t i = 0; i < nodes.size(); ++i) {
                HypernodeID hn = nodes[i];

                if (hg.isFixedVertex(hn) && (hg.fixedVertexPartID(hn) != partitions[i])) {
                    return false;
                }
            }
            return true;
        } (), "Partition of fixed vertex changed.");

        return partitions;
    }


    static inline void prepack_heavy_vertices(Hypergraph& hg, const Context& context, const PartitionID& rb_range_k, bool optimistic) {
        std::vector<HypernodeID> nodes = extract_nodes_with_descending_weight(hg);

        const std::vector<HypernodeWeight>& allowed_weights = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<HypernodeWeight>& perfect_weights = context.initial_partitioning.perfect_balance_partition_weight;

        ALWAYS_ASSERT((allowed_weights.size() == 2) && (perfect_weights.size() == 2));
        ASSERT(context.initial_partitioning.current_max_bin >= hg.totalWeight() / rb_range_k);

        // TODO optimal value?
        const double FACTOR = 1.0;
        HypernodeWeight allowed_imbalance = FACTOR * floor(context.initial_partitioning.current_max_bin * context.initial_partitioning.bin_epsilon);

        std::pair<size_t, HypernodeWeight> treshhold = optimistic ? 
            calculate_heavy_nodes_treshhold_optimistic(hg, nodes, rb_range_k, allowed_imbalance)
            : calculate_heavy_nodes_treshhold_pessimistic(hg, nodes, rb_range_k, allowed_imbalance);

        nodes.resize(treshhold.first);
        std::vector<PartitionID> partitions = apply_bin_packing_to_nodes(hg, context, nodes);

        for (size_t i = 0; i < nodes.size(); ++i) {
            hg.setFixedVertex(nodes[i], partitions[i]);
        }
    }

    // TODO refactor with metrics::finalLevelBinImbalance
    static inline HypernodeWeight maxBinWeight(const Hypergraph& hypergraph, const PartitionID& rb_range_k) {
        // initialize queue
        BinaryMinHeap<PartitionID, HypernodeWeight> queue(rb_range_k);
        for (PartitionID j = 0; j < rb_range_k; ++j) {
                queue.push(j, 0);
            }

        // assign nodes
        std::vector<HypernodeID> hypernodes = kahypar::bin_packing::extract_nodes_with_descending_weight(hypergraph);
        for (const HypernodeID& hn : hypernodes) {
            PartitionID bin = queue.top();
            queue.increaseKeyBy(bin, hypergraph.nodeWeight(hn));
        }

        // the maximum is the biggest bin, i.e. the last element in the queue
        while (queue.size() > 1) {
            queue.pop();
        }
        return queue.getKey(queue.top());
    }
}
}