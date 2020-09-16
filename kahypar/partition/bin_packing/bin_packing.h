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
#include "kahypar/datastructure/segment_tree.h"

namespace kahypar {

namespace bin_packing {
    using kahypar::ds::BinaryMinHeap;
    using kahypar::ds::BinaryMaxHeap;

    // segment tree that used to ensure the balance of a prepacking
    template<typename S>
    S balance_max(const S& v1, const S& v2, const std::vector<std::pair<S, S>>& /*seq*/, const S& /*k*/) {
        return std::max(v1, v2);
    }

    template<typename S>
    S balance_base(const size_t& i, const std::vector<std::pair<S, S>>& seq, const S& k) {
        return k * seq[i].first - seq[i].second;
    }

    template<typename S>
    using BalanceSegTree = ds::ParametrizedSegmentTree<std::pair<S, S>, S, S, balance_max<S>, balance_base<S>>;

    // Bin packing algorithm classes
    //
    // The following operations must be available for a type BPAlgorithm:
    // 1)
    //      PartitionID num_bins = ...;
    //      HypernodeWeight max_bin_weight = ...;
    //      BPAlgorithm alg(num_bins, max_bin_weight);
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
    //
    // 6)
    //      PartitionID numBins = alg.numBins();

    class PartitionMapping {
        public:
            PartitionMapping(PartitionID num_bins) : _partitions() {
                ASSERT(num_bins > 0, "Number of bins must be positive.");

                _partitions.resize(static_cast<size_t>(num_bins), -1);
            }

            void setPartition(PartitionID bin, PartitionID partition) {
                size_t index = static_cast<size_t>(bin);
                ASSERT(bin >= 0 && index < _partitions.size(), "Invalid bin id: " << V(bin));
                ASSERT(_partitions[index] == -1 || _partitions[index] == partition,
                       "Bin already assigned to other partition");

                _partitions[index] = partition;
            }

            bool isFixedBin(PartitionID bin) const {
                return binPartition(bin) != -1;
            }

            PartitionID binPartition(PartitionID bin) const {
                ASSERT(bin >= 0 && static_cast<size_t>(bin) < _partitions.size(), "Invalid bin id: " << V(bin));

                return _partitions[static_cast<size_t>(bin)];
            }

            void applyMapping(std::vector<PartitionID>& elements) {
                for (PartitionID& id : elements) {
                    ASSERT(binPartition(id) >= 0, "Bin not assigned: " << id);
                    id = binPartition(id);
                }
            }

        private:
            std::vector<PartitionID> _partitions;
    };

    class WorstFit {
        public:
            WorstFit(PartitionID num_bins, HypernodeWeight /*max*/) :
                _bin_queue(num_bins),
                _weights(),
                _num_bins(num_bins) {
                for (PartitionID i = 0; i < num_bins; ++i) {
                    _bin_queue.push(i, 0);
                }
            }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));

                _bin_queue.increaseKeyBy(bin, weight);
            }

            PartitionID insertElement(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");
                ASSERT(!_bin_queue.empty(), "All available bins are locked.");

                // assign node to bin with lowest weight
                PartitionID bin = _bin_queue.top();
                _bin_queue.increaseKeyBy(bin, weight);
                return bin;
            }

            void lockBin(PartitionID bin) {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));
                ASSERT(_bin_queue.contains(bin), "Bin already locked.");

                if (_weights.empty()) {
                    _weights.resize(_num_bins, 0);
                }
                _weights[bin] = _bin_queue.getKey(bin);
                _bin_queue.remove(bin);
            }

            HypernodeWeight binWeight(PartitionID bin) const {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));

                return _bin_queue.contains(bin) ? _bin_queue.getKey(bin) : _weights[bin];
            }

            PartitionID numBins() const {
                return _num_bins;
            }

        private:
            BinaryMinHeap<PartitionID, HypernodeWeight> _bin_queue;
            std::vector<HypernodeWeight> _weights;
            PartitionID _num_bins;
    };

    class FirstFit {
        public:
            FirstFit(PartitionID num_bins, HypernodeWeight max) :
                _max_bin_weight(max),
                _bins(num_bins, {0, false}),
                _num_bins(num_bins) { }

            void addWeight(PartitionID bin, HypernodeWeight weight) {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));

                _bins[bin].first += weight;
            }

            PartitionID insertElement(HypernodeWeight weight) {
                ASSERT(weight >= 0, "Negative weight.");

                size_t assigned_bin = 0;
                for (size_t i = 0; i < _bins.size(); ++i) {
                    if (_bins[i].second) {
                        continue;
                    }

                    // The node is assigned to the first fitting bin or, if none fits, the smallest bin.
                    if (_bins[i].first + weight <= _max_bin_weight) {
                        assigned_bin = i;
                        break;
                    } else if (_bins[assigned_bin].second || _bins[i].first < _bins[assigned_bin].first) {
                        assigned_bin = i;
                    }
                }

                ASSERT(!_bins[assigned_bin].second, "All available bins are locked.");
                _bins[assigned_bin].first += weight;
                return assigned_bin;
            }

            void lockBin(PartitionID bin) {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));
                ASSERT(!_bins[bin].second, "Bin already locked.");

                _bins[bin].second = true;
            }

            HypernodeWeight binWeight(PartitionID bin) const {
                ASSERT(bin >= 0 && bin < _num_bins, "Invalid bin id: " << V(bin));

                return _bins[bin].first;
            }

            PartitionID numBins() const {
                return _num_bins;
            }

        private:
            HypernodeWeight _max_bin_weight;
            std::vector<std::pair<HypernodeWeight, bool>> _bins;
            PartitionID _num_bins;
    };

    template< class BPAlgorithm >
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

            HypernodeWeight binWeight(PartitionID bin) const {
                return _alg.binWeight(bin);
            }

            // ... and at the second level, the resulting k bins are packed into the final partitions.
            // Returns the partition mapping for the bins and a vector of the resulting partition weights
            std::pair<PartitionMapping, std::vector<HypernodeWeight>> applySecondLevel(const std::vector<HypernodeWeight>& max_allowed_partition_weights,
                                                                          const std::vector<PartitionID>& num_bins_per_partition) const {
                ALWAYS_ASSERT(num_bins_per_partition.size() == max_allowed_partition_weights.size(),
                    "max_allowed_partition_weights and num_bins_per_partition have different sizes: "
                    << V(max_allowed_partition_weights.size()) << "; " << V(num_bins_per_partition.size()));

                PartitionID num_partitions = static_cast<PartitionID>(max_allowed_partition_weights.size());
                std::vector<PartitionID> bin_counts(max_allowed_partition_weights.size(), 0);
                PartitionMapping mapping = _bins_to_parts;

                HypernodeWeight max_partition = *std::max_element(max_allowed_partition_weights.cbegin(),
                                                                  max_allowed_partition_weights.cend());
                BPAlgorithm partition_packer(num_partitions, max_partition);
                for (PartitionID i = 0; i < num_partitions; ++i) {
                    partition_packer.addWeight(i, max_partition - max_allowed_partition_weights[i]);
                }

                for (PartitionID bin = 0; bin < _alg.numBins(); ++bin) {
                    if (mapping.isFixedBin(bin)) {
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

                // collect partition weights
                std::vector<HypernodeWeight> part_weights;
                part_weights.reserve(num_partitions);
                for (PartitionID i = 0; i < num_partitions; ++i) {
                    part_weights.push_back(partition_packer.binWeight(i) + max_allowed_partition_weights[i] - max_partition);
                }

                return {std::move(mapping), std::move(part_weights)};
            }

        private:
            BPAlgorithm _alg;
            PartitionMapping _bins_to_parts;
    };

    /**
     * Returns the current hypernodes sorted in descending order of weight.
     */
    static inline std::vector<HypernodeID> extractNodesWithDescendingWeight(const Hypergraph& hg) {
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

    // preassigns the fixed vertices to the packer with a first fit packing
    template< class BPAlgorithm >
    static inline void preassignFixedVertices(const Hypergraph& hg, const std::vector<HypernodeID>& nodes, std::vector<PartitionID>& partitions,
                                                TwoLevelPacker<BPAlgorithm>& packer, PartitionID k, PartitionID rb_range_k) {
        HypernodeWeight avg_bin_weight = (hg.totalWeight() + rb_range_k - 1) / rb_range_k;
        PartitionID kbins_per_partition = (rb_range_k + k - 1) / k;

        for (size_t i = 0; i < nodes.size(); ++i) {
            HypernodeID hn = nodes[i];

            if (hg.isFixedVertex(hn)) {
                HypernodeWeight weight = hg.nodeWeight(hn);
                PartitionID part_id = hg.fixedVertexPartID(hn);

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
                partitions[i] = assigned_bin;
            }
        }

    }

    static inline size_t get_max_part_idx(const Context& context, const std::vector<HypernodeWeight>& part_weight,
                                          HypernodeWeight next_element, bool skip_full_parts) {
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1.0 + context.initial_partitioning.bin_epsilon));
        const std::vector<HypernodeWeight>& upper_weight = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<PartitionID>& num_bins_per_part = context.initial_partitioning.num_bins_per_partition;
        ASSERT(part_weight.size() == upper_weight.size() && part_weight.size() == num_bins_per_part.size());

        size_t max_imb_part_id = 0;
        HypernodeWeight min_remaining = std::numeric_limits<HypernodeWeight>::max();
        for (size_t j = 0; j < upper_weight.size(); ++j) {
            // skip valid full partitions - those can not provide an upper bound, because the weight sequence is empty
            HypernodeWeight remaining = upper_weight[j] - part_weight[j];
            HypernodeWeight allowed = num_bins_per_part[j] * max_bin_weight - upper_weight[j];
            if (skip_full_parts && (next_element > remaining) && ((num_bins_per_part[j] - 1) * remaining <= allowed)) {
                continue;
            }

            if (remaining < min_remaining) {
                max_imb_part_id = j;
                min_remaining = remaining;
            }
        }

        return max_imb_part_id;
    }


    template< class BPAlgorithm >
    static inline void calculateExactPrepacking(Hypergraph& hg, Context& context, PartitionID rb_range_k, HypernodeWeight max_bin_weight) {
        ALWAYS_ASSERT(rb_range_k > context.partition.k);
        ASSERT(context.initial_partitioning.current_max_bin >= hg.totalWeight() / rb_range_k);
        ASSERT(!hg.containsFixedVertices(), "No fixed vertices allowed before prepacking.");

        std::vector<HypernodeWeight>& upper_weight = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<PartitionID>& num_bins_per_part = context.initial_partitioning.num_bins_per_partition;

        ASSERT(upper_weight.size() == static_cast<size_t>(context.partition.k)
               && num_bins_per_part.size() == static_cast<size_t>(context.partition.k));
        ASSERT(std::accumulate(num_bins_per_part.cbegin(), num_bins_per_part.cend(), 0) >= rb_range_k);

        // initialization: exctract descending nodes, calculate weight sum table and initialize segment tree
        PartitionID max_k = *std::max_element(num_bins_per_part.cbegin(), num_bins_per_part.cend());
        TwoLevelPacker<BPAlgorithm> packer(rb_range_k, max_bin_weight);
        std::vector<HypernodeID> nodes = extractNodesWithDescendingWeight(hg);
        std::vector<std::pair<HypernodeWeight, HypernodeWeight>> weights(nodes.size() + 1);
        HypernodeWeight sum = 0;
        weights[nodes.size()] = {0, 0};
        for (size_t i = nodes.size(); i > 0; --i) {
            HypernodeWeight w = hg.nodeWeight(nodes[i - 1]);
            sum += w;
            weights[i - 1] = {w, sum};
        }
        BalanceSegTree<HypernodeWeight> seg_tree(weights, max_k);
        std::vector<PartitionID> partitions;

        ASSERT(nodes.size() > 0);
        ASSERT([&]() {
            for (size_t i = 0; i < static_cast<size_t>(context.partition.k); ++i) {
                if (upper_weight[i] < context.initial_partitioning.perfect_balance_partition_weight[i]
                    || upper_weight[i] >= num_bins_per_part[i] * max_bin_weight
                    ) {
                    return false;
                }
            }
            return true;
        } (), "The allowed partition or bin weights are to small to find a valid partition.");

        size_t i;
        std::pair<PartitionMapping, std::vector<HypernodeWeight>> packing_result = packer.applySecondLevel(upper_weight, num_bins_per_part);
        for (i = 0; i < nodes.size(); ++i) {
            ASSERT(upper_weight.size() == packing_result.second.size());

            size_t max_part_idx = get_max_part_idx(context, packing_result.second, weights[i].first, true);
            HypernodeWeight remaining = std::max(0, std::min(packing_result.second[max_part_idx] - upper_weight[max_part_idx] + weights[i].second,
                                        weights[i].second));

            // find subrange of specified weight
            size_t j = weights.crend() - std::lower_bound(weights.crbegin(), weights.crend() - i, std::make_pair(0, remaining),
                       [&](const auto& v1, const auto& v2) {
                           return v1.second < v2.second;
                       }) - 1;

            // calculate the bound for the subrange
            if (j > i) {
                HypernodeWeight imbalance = seg_tree.query(i, j - 1) + weights[j].second;
                HypernodeWeight partWeight = (j == nodes.size()) ? packing_result.second[max_part_idx] + weights[i].second : upper_weight[max_part_idx]
                                             * max_k / num_bins_per_part[max_part_idx];
                if (partWeight + imbalance <= max_k * max_bin_weight) {
                    break;
                }
            }

            partitions.push_back(packer.insertElement(weights[i].first));
            packing_result = packer.applySecondLevel(upper_weight, num_bins_per_part);
        }
        ASSERT(partitions.size() <= nodes.size());

        packing_result.first.applyMapping(partitions);
        for (size_t i = 0; i < partitions.size(); ++i) {
            hg.setFixedVertex(nodes[i], partitions[i]);
        }

        // calculate optimization for allowed weights
        size_t max_part_idx = get_max_part_idx(context, packing_result.second, 0, false);
        HypernodeWeight range_weight = packing_result.second[max_part_idx];
        HypernodeWeight imbalance = 0;
        HypernodeWeight optimized = 0;
        for (size_t j = i; j < nodes.size(); ++j) {
            HypernodeWeight weight = weights[j].first;
            imbalance = std::max(imbalance - weight, (max_k - 1) * weight);
            range_weight += weight;
            HypernodeWeight bin_weight = range_weight * max_k / num_bins_per_part[max_part_idx];

            if (bin_weight + imbalance > max_k * max_bin_weight) {
                break;
            }
            optimized = max_k * max_bin_weight - imbalance;
        }

        for (size_t i = 0; i < upper_weight.size(); ++i) {
            upper_weight[i] = std::max(upper_weight[i], num_bins_per_part[i] * optimized / max_k);
        }
    }

    template< class BPAlgorithm >
    static inline void calculateHeuristicPrepacking(Hypergraph& hg, const Context& context, const PartitionID& rb_range_k, HypernodeWeight max_bin_weight) {
        ALWAYS_ASSERT(rb_range_k > context.partition.k);
        ASSERT(context.initial_partitioning.current_max_bin >= hg.totalWeight() / rb_range_k);
        ASSERT(!hg.containsFixedVertices(), "No fixed vertices allowed before prepacking.");
        ASSERT((context.initial_partitioning.upper_allowed_partition_weight.size() == static_cast<size_t>(context.partition.k)) &&
               (context.initial_partitioning.perfect_balance_partition_weight.size() == static_cast<size_t>(context.partition.k)));

        std::vector<HypernodeID> nodes = extractNodesWithDescendingWeight(hg);
        HypernodeWeight allowed_imbalance = max_bin_weight - context.initial_partitioning.current_max_bin;
        ASSERT(allowed_imbalance > 0, "allowed_imbalance is zero!");

        // calculate heuristic threshold
        HypernodeID threshold = 0;
        HypernodeWeight current_lower_sum = 0;
        HypernodeWeight imbalance = 0;

        for (int i = nodes.size() - 1; i >= 0; --i) {
            HypernodeWeight weight = hg.nodeWeight(nodes[i]);
            current_lower_sum += weight;
            imbalance = weight - current_lower_sum / rb_range_k;

            if (imbalance > allowed_imbalance) {
                threshold = i + 1;
                break;
            }
        }

        // bin pack the heavy nodes
        nodes.resize(threshold);
        std::vector<PartitionID> partitions(nodes.size(), -1);
        TwoLevelPacker<BPAlgorithm> packer(rb_range_k, max_bin_weight);

        for (size_t i = 0; i < nodes.size(); ++i) {
            HypernodeWeight weight = hg.nodeWeight(nodes[i]);
            partitions[i] = packer.insertElement(weight);
        }

        PartitionMapping packing_result = packer.applySecondLevel(context.initial_partitioning.upper_allowed_partition_weight,
                                                                    context.initial_partitioning.num_bins_per_partition).first;
        packing_result.applyMapping(partitions);

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
        std::vector<HypernodeID> hypernodes = extractNodesWithDescendingWeight(hypergraph);
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