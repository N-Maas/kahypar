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

    enum class BalancingLevel : uint8_t {
        none,
        optimistic,
        guaranteed,
        STOP
    };

    BalancingLevel increaseBalancingRestrictions(BalancingLevel previous) {
        switch (previous) {
            case BalancingLevel::none:
                return BalancingLevel::optimistic;
            case BalancingLevel::optimistic:
                return BalancingLevel::guaranteed;
            case BalancingLevel::guaranteed:
                return BalancingLevel::STOP;
            case BalancingLevel::STOP:
                break;
                // omit default case to trigger compiler warning for missing cases
        }
        ASSERT(false, "Tried to increase invalid balancing level: " << static_cast<uint8_t>(previous));
        return previous;
    }

    // Segment tree definition for heuristic
    template<typename S>
    S heuristic_max(const S& v1, const S& v2, const std::vector<std::pair<S, S>>& seq, const S& k) {
        return std::max(v1, v2);
    }

    template<typename S>
    S heuristic_base(const size_t& i, const std::vector<std::pair<S, S>>& seq, const S& k) {
        return k * seq[i].first - seq[i].second;
    }

    template<typename S>
    using HeuristicSegTree = ds::ParametrizedSegmentTree<std::pair<S, S>, S, S, heuristic_max<S>, heuristic_base<S>>;

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
            WorstFit(PartitionID num_bins, HypernodeWeight max) :
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
                BPAlg partition_packer(num_partitions, max_partition);
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
            BPAlg _alg;
            PartitionMapping _bins_to_parts;
    };

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

    // this functions handles already existing fixed vertices correctly
    template< class BPAlg >
    static inline std::vector<PartitionID> apply_bin_packing_to_nodes(const Hypergraph& hg,
                                                                      const Context& context,
                                                                      const std::vector<HypernodeID>& nodes) {
        ASSERT(static_cast<size_t>(context.partition.k) == context.initial_partitioning.upper_allowed_partition_weight.size());

        PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1 + context.initial_partitioning.bin_epsilon));
        std::vector<PartitionID> partitions(nodes.size(), -1);
        TwoLevelPacker<BPAlg> packer(rb_range_k, max_bin_weight);

        if (hg.containsFixedVertices()) {
            // pre-pack the fixed vertices with a first fit packing
            HypernodeWeight avg_bin_weight = (hg.totalWeight() + rb_range_k - 1) / rb_range_k;
            PartitionID kbins_per_partition = rb_range_k / context.partition.k;

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

        for (size_t i = 0; i < nodes.size(); ++i) {
            HypernodeID hn = nodes[i];

            if(!hg.isFixedVertex(hn)) {
                HypernodeWeight weight = hg.nodeWeight(hn);
                partitions[i] = packer.insertElement(weight);
            }
        }

        PartitionMapping packing_result = packer.applySecondLevel(context.initial_partitioning.upper_allowed_partition_weight,
                                                                  context.initial_partitioning.num_bins_per_partition).first;
        packing_result.applyMapping(partitions);

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

    static inline size_t get_max_part_idx(const Context& context, const std::vector<HypernodeWeight>& part_weight,
                                          HypernodeWeight next_element, bool skip_full_parts) {
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1.0 + context.initial_partitioning.bin_epsilon));
        const std::vector<HypernodeWeight>& upper_weight = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<PartitionID>& num_bins_per_part = context.initial_partitioning.num_bins_per_partition;
        ASSERT(part_weight.size() == upper_weight.size() && part_weight.size() == num_bins_per_part.size());

        size_t max_imb_part_id = 0;
        HypernodeWeight min_remaining = std::numeric_limits<HypernodeWeight>::max();
        for (size_t j = 0; j < upper_weight.size(); ++j) {
            // skip valid full partitions - those can not provide an upper bound, as the heuristic is empty
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


    template< class BPAlg >
    static inline void apply_prepacking_pessimistic(Hypergraph& hg, Context& context) {
        ALWAYS_ASSERT(!hg.containsFixedVertices(), "No fixed vertices allowed before prepacking.");

        PartitionID rb_range_k = context.partition.rb_upper_k - context.partition.rb_lower_k + 1;
        if (rb_range_k <= context.partition.k) {
            return;
        }

        std::vector<HypernodeWeight>& upper_weight = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<PartitionID>& num_bins_per_part = context.initial_partitioning.num_bins_per_partition;

        ALWAYS_ASSERT(upper_weight.size() == static_cast<size_t>(context.partition.k)
                      && num_bins_per_part.size() == static_cast<size_t>(context.partition.k),
                      "Invalid context: " << V(context.partition.k));
        ASSERT(context.initial_partitioning.perfect_balance_partition_weight.size() == static_cast<size_t>(context.partition.k));
        ASSERT(context.initial_partitioning.current_max_bin >= hg.totalWeight() / rb_range_k);
        ASSERT(std::accumulate(num_bins_per_part.cbegin(), num_bins_per_part.cend(), 0) >= rb_range_k);

        // initialization: exctract descending nodes, calculate weight sum table and initialize segment tree
        HypernodeWeight max_bin_weight = floor(context.initial_partitioning.current_max_bin * (1.0 + context.initial_partitioning.bin_epsilon));
        PartitionID max_k = *std::max_element(num_bins_per_part.cbegin(), num_bins_per_part.cend());
        TwoLevelPacker<BPAlg> packer(rb_range_k, max_bin_weight);
        std::vector<HypernodeID> nodes = extract_nodes_with_descending_weight(hg);
        std::vector<std::pair<HypernodeWeight, HypernodeWeight>> weights(nodes.size() + 1);
        HypernodeWeight sum = 0;
        weights[nodes.size()] = {0, 0};
        for (size_t i = nodes.size(); i > 0; --i) {
            HypernodeWeight w = hg.nodeWeight(nodes[i - 1]);
            sum += w;
            weights[i - 1] = {w, sum};
        }
        HeuristicSegTree<HypernodeWeight> seg_tree(weights, max_k);
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

            // calculate the heuristic of the subrange
            if (j > i) {
                HypernodeWeight imbalance = seg_tree.query(i, j - 1) + weights[j].second;
                HypernodeWeight partWeight = (j == nodes.size()) ? packing_result.second[max_part_idx] + weights[i].second : upper_weight[max_part_idx]
                                             * max_k / num_bins_per_part[max_part_idx];
                if (partWeight + imbalance <= max_k * max_bin_weight) {
                    break;
                }
            }

            // insert element
            partitions.push_back(packer.insertElement(weights[i].first));
            packing_result = packer.applySecondLevel(upper_weight, num_bins_per_part);
        }
        ASSERT(partitions.size() <= nodes.size());

        // apply fixed vertices
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

        // apply allowed weight optimization
        for (size_t i = 0; i < upper_weight.size(); ++i) {
            upper_weight[i] = std::max(upper_weight[i], num_bins_per_part[i] * optimized / max_k);
        }
    }

    template< class BPAlgorithm >
    static inline void prepack_heavy_vertices(Hypergraph& hg, const Context& context, const PartitionID& rb_range_k) {
        std::vector<HypernodeID> nodes = extract_nodes_with_descending_weight(hg);

        const std::vector<HypernodeWeight>& allowed_weights = context.initial_partitioning.upper_allowed_partition_weight;
        const std::vector<HypernodeWeight>& perfect_weights = context.initial_partitioning.perfect_balance_partition_weight;

        ALWAYS_ASSERT((allowed_weights.size() == 2) && (perfect_weights.size() == 2));
        ASSERT(context.initial_partitioning.current_max_bin >= hg.totalWeight() / rb_range_k);

        HypernodeWeight allowed_imbalance = floor(context.initial_partitioning.current_max_bin * context.initial_partitioning.bin_epsilon);
        ASSERT(allowed_imbalance > 0, "allowed_imbalance is zero!");

        std::pair<size_t, HypernodeWeight> treshhold = calculate_heavy_nodes_treshhold_optimistic(hg, nodes, rb_range_k, allowed_imbalance);

        nodes.resize(treshhold.first);
        std::vector<PartitionID> partitions = apply_bin_packing_to_nodes<BPAlgorithm>(hg, context, nodes);

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