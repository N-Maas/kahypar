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
    //      HypernodeWeight max_bin_weight = ...;
    //      BPAlg alg(num_bins, max_bin_weight);
    //
    // 2)
    //      PartitionID bin = ...;
    //      HypernodeWeight weight = ...;
    //      alg.addWeight(bin, weight);
    //
    // 3)
    //      HypernodeWeight weight = ...;
    //      PartitionID resulting_bin = alg.insertElement(weight);
    //
    // 4)
    //      ParititionID bin = ...;
    //      alg.lockBin(bin);
    //
    // 5)
    //      ParititionID bin = ...;
    //      HypernodeWeight weight = alg.binWeight(bin);

    class PartitionMapping {
        public:
            PartitionMapping(PartitionID num_bins) : partitions() {
                ASSERT(num_bins > 0, "Number of bins must be positive.");

                partitions.resize(num_bins, -1);
            }

            void setPartition(PartitionID bin, PartitionID partition) {
                assertBinIsValid(bin);
                ASSERT(partitions[bin] == -1 || partitions[bin] == partition,
                       "Bin already assigned to other partition");

                partitions[bin] = partition;
            }

            bool isFixedBin(PartitionID bin) const {
                return binPartition(bin) != -1;
            }

            PartitionID binPartition(PartitionID bin) const {
                assertBinIsValid(bin);

                return partitions[bin];
            }

            void applyMapping(std::vector<PartitionID>& elements) {
                for (PartitionID& id : elements) {
                    ASSERT(binPartition(id) >= 0, "Bin not assigned: " << id);
                    id = binPartition(id);
                }
            }

        protected:
            void assertBinIsValid(PartitionID bin) const {
                ASSERT(bin >= 0 && static_cast<size_t>(bin) < partitions.size(), "Invalid bin id.");
            }

        private:
            std::vector<PartitionID> partitions;
    };

    class WorstFit {
        public:
            WorstFit(PartitionID num_bins, HypernodeWeight max) :
                bin_queue(num_bins),
                weights(),
                num_bins(num_bins) {
                for (PartitionID i = 0; i < num_bins; ++i) {
                    bin_queue.push(i, 0);
                }
            }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");

                bin_queue.increaseKeyBy(bin, weight);
            }

            PartitionID insertElement(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");
                ASSERT(!bin_queue.empty(), "All available bins are locked.");

                // assign node to bin with lowest weight
                PartitionID bin = bin_queue.top();
                bin_queue.increaseKeyBy(bin, weight);
                return bin;
            }

            void lockBin(PartitionID bin) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");
                ASSERT(bin_queue.contains(bin), "Bin already locked.");

                if (weights.empty()) {
                    weights.resize(num_bins, 0);
                }
                weights[bin] = bin_queue.getKey(bin);
                bin_queue.remove(bin);
            }

            HypernodeWeight binWeight(PartitionID bin) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");

                return bin_queue.contains(bin) ? bin_queue.getKey(bin) : weights[bin];
            }

            PartitionID numBins() {
                return _num_bins;
            }

        private:
            BinaryMinHeap<PartitionID, HypernodeWeight> bin_queue;
            std::vector<HypernodeWeight> weights;
            PartitionID num_bins;
    };

    class FirstFit {
        public:
            FirstFit(PartitionID num_bins, HypernodeWeight max) :
                max_bin_weight(max),
                bins(num_bins, {0, false}),
                num_bins(num_bins) { }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");

                bins[bin].first += weight;
            }

            PartitionID insertElement(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");

                size_t assigned_bin = 0;
                for (size_t i = 0; i < bins.size(); ++i) {
                    if (bins[i].second) {
                        continue;
                    }

                    // The node is assigned to the first fitting bin or, if none fits, the smallest bin.
                    if (bins[i].first + weight <= max_bin_weight) {
                        assigned_bin = i;
                        break;
                    } else if (bins[assigned_bin].second || bins[i].first < bins[assigned_bin].first) {
                        assigned_bin = i;
                    }
                }

                ASSERT(!bins[assigned_bin].second, "All available bins are locked.");
                bins[assigned_bin].first += weight;
                return assigned_bin;
            }

            void lockBin(PartitionID bin) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");
                ASSERT(!bins[bin].second, "Bin already locked.");

                bins[bin].second = true;
            }

            HypernodeWeight binWeight(PartitionID bin) {
                ASSERT(bin >= 0 && bin < num_bins, "Invalid bin id.");

                return bins[bin].first;
            }

            PartitionID numBins() {
                return _num_bins;
            }

        private:
            HypernodeWeight max_bin_weight;
            std::vector<std::pair<HypernodeWeight, bool>> bins;
            PartitionID num_bins;
    };

    template< class BPAlg >
    class BinCountWrapper {
        public:
            BinCountWrapper(PartitionID num_bins, HypernodeWeight max, const std::vector<PartitionID>& max_num_bins) :
                alg(num_bins, max),
                bin_counts() {
                    ASSERT(static_cast<size_t>(num_bins) == max_num_bins.size(), "Inconsistent number of bins.");

                    bin_counts.reserve(max_num_bins.size());
                    for (const PartitionID& num_bins : max_num_bins) {
                        ASSERT(num_bins > 0, "Every partition must allow at least one bin.");
                        bin_counts.push_back({0, num_bins});
                    }
                }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                alg.addWeight(bin, weight);
            }

            PartitionID insertElement(HypernodeWeight weight) {
                PartitionID bin = alg.insertElement(weight);

                ASSERT(bin >= 0 && static_cast<size_t>(bin) < bin_counts.size(), "Invalid bin id returned.");
                bin_counts[bin].first++;
                if (bin_counts[bin].first >= bin_counts[bin].second) {
                    lockBin(bin);
                }
                return bin;
            }

            void lockBin(PartitionID bin) {
                alg.lockBin(bin);
            }

            HypernodeWeight binWeight(PartitionID bin) {
                return alg.binWeight(bin);
            }

        private:
            BPAlg alg;
            std::vector<std::pair<PartitionID, PartitionID>> bin_counts;
    };

    template< class BPAlg >
    class TwoLevelPacker {
        public:
            TwoLevelPacker(PartitionID num_bins, HypernodeWeight max_bin) :
                _alg(num_bins, max_bin),
                _bins_to_parts(num_bins) {
                ASSERT(num_bins > 0, "Number of bins must be positive.");
            }

            // At the first level, a packing with k bins is calculated ....
            PartitionID insertElement(HypernodeWeight weight) {
                return _alg.insertElement(weight);
            }

            void addFixedVertex(PartitionID bin, PartitionID partition, HypernodeWeight weight) {
                _bins_to_parts.setPartition(bin, partition);
                _alg.addWeight(bin, weight);
            }

            HypernodeWeight binWeight(PartitionID bin) {
                return _alg.binWeight(bin);
            }

            // ... and at the second level, the resulting k bins are packed into the final partitions
            std::pair<PartitionMapping, HypernodeWeight> applySecondLevel(const std::vector<HypernodeWeight>& max_allowed_partition_weights,
                                                                          const std::vector<PartitionID>& num_bins_per_partition) {
                ALWAYS_ASSERT(num_bins_per_partition.size() == max_allowed_partition_weights.size(),
                    "max_allowed_partition_weights and num_bins_per_partition have different sizes: "
                    << V(max_allowed_partition_weights.size()) << "; " << V(num_bins_per_partition.size()));

                PartitionID num_partitions = static_cast<PartitionID>(max_allowed_partition_weights.size());
                std::vector<PartitionID> bin_counts(max_allowed_partition_weights.size(), 0);
                PartitionMapping mapping = _bins_to_parts;

                HypernodeWeight max_partition = *std::max_element(max_allowed_partition_weights.cbegin(),
                                                                  max_allowed_partition_weights.cend());
                BPAlg partition_packer(num_partitions, max_partition);
                for (PartitionID i = 0; i < num_partitions; ++i) {
                    partition_packer.addWeight(i, max_partition - max_allowed_partition_weights[i]);
                }

                for (PartitionID bin = 0; bin < _alg.numBins(); ++bin) {
                    if(mapping.isFixedBin(bin)) {
                        partition_packer.addWeight(mapping.binPartition(bin), _alg.binWeight(bin));
                    }
                }

                std::vector<PartitionID> kbins_descending(_alg.numBins());
                std::iota(kbins_descending.begin(), kbins_descending.end(), 0);
                std::sort(kbins_descending.begin(), kbins_descending.end(), [&](PartitionID i, PartitionID j) {
                    return _alg.binWeight(i) > _alg.binWeight(j);
                });

                for (const PartitionID& bin : kbins_descending) {
                    // if not fixed, assign bin to partition
                    if (!mapping.isFixedBin(bin)) {
                        PartitionID partition = partition_packer.insertElement(_alg.binWeight(bin));
                        ASSERT((partition >= 0) && (partition < num_partitions));

                        size_t part_idx = static_cast<size_t>(partition);
                        bin_counts[part_idx]++;
                        if (bin_counts[part_idx] >= num_bins_per_partition[part_idx]) {
                            partition_packer.lockBin(partition);
                        }

                        mapping.setPartition(bin, partition);
                    }
                }

                return {std::move(mapping), _alg.binWeight(kbins_descending[0])};
            }

        private:
            BPAlg _alg;
            PartitionMapping _bins_to_parts;
    };
    /**
     * The hypernodes must be sorted in descending order of weight. 
     * The partitions parameter can be used to specify fixed vertices.
     * With max_allowed_partition_weights, weights for partitions can
     * be added (e.g. to handle uneven k).
     *
     * Attention: the algorithm assumus that the partitions can be balanced
     * with rb_range_k bins of equal weight.
     */
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

        BPAlg bin_packer(rb_range_k, max_bin_weight);
        PartitionMapping bins_to_parts(rb_range_k);

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
                bins_to_parts.setPartition(assigned_bin, part_id);
                partitions[index] = assigned_bin;
            }
        }

        // At the first level, a packing with k bins is calculated ....
        for (size_t i = 0; i < hypernodes.size(); ++i) {
            if (partitions[i] == -1) {
                HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);
                partitions[i] = bin_packer.insertElement(weight);
            }
        }

        std::vector<PartitionID> kbins_ascending(rb_range_k);
        std::iota(kbins_ascending.begin(), kbins_ascending.end(), 0);
        std::sort(kbins_ascending.begin(), kbins_ascending.end(), [&](PartitionID i, PartitionID j) {
                return bin_packer.binWeight(i) < bin_packer.binWeight(j);
            });
        ASSERT(kbins_ascending.size() == static_cast<size_t>(rb_range_k));

        // ... and at the second level, the resulting k bins are packed into the final partitions
        HypernodeWeight max_partition = *std::max_element(max_allowed_partition_weights.cbegin(),
                                                    max_allowed_partition_weights.cend());
        BinCountWrapper<BPAlg> partition_packer(num_partitions, max_partition, num_bins_per_partition);

        for (PartitionID i = 0; i < num_partitions; ++i) {
            partition_packer.addWeight(i, max_partition - max_allowed_partition_weights[i]);
        }

        for (const PartitionID& bin : kbins_ascending) {
            if (bins_to_parts.isFixedBin(bin)) {
                partition_packer.addWeight(bins_to_parts.binPartition(bin), bin_packer.binWeight(bin));
            }
        }

        // iteration must be in reverse order, because the vector is ascending
        for (const PartitionID& bin : boost::adaptors::reverse(kbins_ascending)) {
            // if not fixed, assign bin to partition
            if (!bins_to_parts.isFixedBin(bin)) {
                PartitionID partition = partition_packer.insertElement(bin_packer.binWeight(bin));
                ASSERT((partition >= 0) && (partition < num_partitions));

                bins_to_parts.setPartition(bin, partition);
            }
        }

        // remap the partition IDs of the nodes
        for (PartitionID& id : partitions) {
            id = bins_to_parts.binPartition(id);
            ASSERT((id >= 0) && (id < rb_range_k));
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

        ASSERT(static_cast<size_t>(context.partition.k) == context.initial_partitioning.upper_allowed_partition_weight.size());

        PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1 + context.initial_partitioning.bin_epsilon));
        partitions = context.initial_partitioning.bp_algo == BinPackingAlgorithm::worst_fit ?
            two_level_packing<WorstFit>(hg, nodes, context.initial_partitioning.upper_allowed_partition_weight,
                                        context.initial_partitioning.num_bins_per_partition, rb_range_k, max_bin_weight, std::move(partitions)) :
            two_level_packing<FirstFit>(hg, nodes, context.initial_partitioning.upper_allowed_partition_weight,
                                        context.initial_partitioning.num_bins_per_partition, rb_range_k, max_bin_weight, std::move(partitions));

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
        ASSERT(allowed_imbalance > 0, "allowed_imbalance is zero!");

        std::pair<size_t, HypernodeWeight> treshhold = optimistic ? 
            calculate_heavy_nodes_treshhold_optimistic(hg, nodes, rb_range_k, allowed_imbalance) :
            calculate_heavy_nodes_treshhold_pessimistic(hg, nodes, rb_range_k, allowed_imbalance);

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