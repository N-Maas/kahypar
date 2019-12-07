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

#include "kahypar/definitions.h"
#include "kahypar/datastructure/binary_heap.h"

namespace kahypar {

namespace bin_packing {
    using kahypar::ds::BinaryMinHeap;

    /**
     * The hypernodes must be sorted in descending order of weight. 
     * The partitions parameter can be used to specify fixed vertices.
     */
    static inline std::vector<PartitionID> two_level_packing(const Hypergraph& hg,
                                                              const std::vector<HypernodeID>& hypernodes,
                                                              const PartitionID& num_partitions,
                                                              const PartitionID& rb_range_k,
                                                              std::vector<PartitionID> partitions = {}) {
        ALWAYS_ASSERT((num_partitions > 0) && (rb_range_k % num_partitions == 0),
            "num_partitions must be positive and k must be a multiple of the number of bins.");
        ALWAYS_ASSERT(partitions.empty() || (partitions.size() == hypernodes.size()),
            "Size of fixed vertice partition IDs does not match the number of hypernodes.");
        ASSERT([&]() {
            for (size_t i = 1; i < hypernodes.size(); ++i) {
                if (hg.nodeWeight(hypernodes[i-1]) < hg.nodeWeight(hypernodes[i])) {
                    return false;
                }
            }
            return true;
        } (), "The hypernodes must be sorted in descending order of weight.");

        BinaryMinHeap<PartitionID, HypernodeWeight> k_bin_queue(rb_range_k);
        for(PartitionID i = 0; i < rb_range_k; ++i) {
            k_bin_queue.push(i, 0);
        }
        std::vector<PartitionID> kbin_to_partition_mapping(rb_range_k, -1);

        if(partitions.empty()) {
            partitions.resize(hypernodes.size(), -1);
        } else {
            // If fixed vertices are specified, calculate the total weight and extract all fixed vertices...
            HypernodeWeight total_weight = 0;
            std::vector<size_t> fixed_vertices;

            for(size_t i = 0; i < hypernodes.size(); ++i) {
                total_weight += hg.nodeWeight(hypernodes[i]);

                if(partitions[i] >= 0) {
                    ALWAYS_ASSERT(partitions[i] < num_partitions, "Invalid partition ID for node.");

                    fixed_vertices.push_back(i);
                }
            }

            // ...to pre-pack the fixed vertices with a first fit packing
            HypernodeWeight avg_bin_weight = (total_weight + rb_range_k - 1) / rb_range_k;
            PartitionID kbins_per_partition = rb_range_k / num_partitions;

            for(const auto& index : fixed_vertices) {
                HypernodeWeight weight = hg.nodeWeight(hypernodes[index]);
                PartitionID part_id = partitions[index];

                PartitionID start_index = part_id * kbins_per_partition;
                PartitionID assigned_bin = start_index;

                if(kbins_per_partition > 1) {
                    for(PartitionID i = start_index; i < start_index + kbins_per_partition; ++i) {
                        HypernodeWeight current_bin_weight = k_bin_queue.getKey(i);

                        // The vertex is assigned to the first fitting bin or, if none fits, the smallest bin.
                        if(current_bin_weight + weight <= avg_bin_weight) {
                            assigned_bin = i;
                            break;
                        } else if (current_bin_weight < k_bin_queue.getKey(assigned_bin)) {
                            assigned_bin = i;
                        }
                    }
                }

                k_bin_queue.increaseKeyBy(assigned_bin, weight);
                kbin_to_partition_mapping[assigned_bin] = part_id;
                partitions[index] = assigned_bin;
            }
        }

        // At the first level, a packing with k bins is calculated ....
        for(size_t i = 0; i < hypernodes.size(); ++i) {
            if(partitions[i] == -1) {
                // assign node to bin with lowest weight
                PartitionID kbin = k_bin_queue.top();
                k_bin_queue.increaseKeyBy(kbin, hg.nodeWeight(hypernodes[i]));
                partitions[i] = kbin;
            }
        }

        BinaryMinHeap<PartitionID, HypernodeWeight> result_queue(num_partitions);
        for(PartitionID i = 0; i < num_partitions; ++i) {
            result_queue.push(i, 0);
        }

        // ... and at the second level, the resulting k bins are packed into the final bins
        if(rb_range_k != num_partitions) {
            std::vector<std::pair<PartitionID, HypernodeWeight>> kbin_ascending_ordering;

            // read bins from queue
            while(k_bin_queue.size() > 0) {
                PartitionID kbin = k_bin_queue.top();
                HypernodeWeight weight = k_bin_queue.topKey();
                k_bin_queue.pop();

                if(kbin_to_partition_mapping[kbin] == -1) {
                    kbin_ascending_ordering.push_back({kbin, weight});
                } else {
                    result_queue.increaseKeyBy(kbin_to_partition_mapping[kbin], weight);
                }
            }
            ASSERT(k_bin_queue.size() == 0);

            // iteration must be in reverse order, because the vector is ascending
            for(const auto& kbin : boost::adaptors::reverse(kbin_ascending_ordering)) {
                // if not fixed, assign bin to partition with lowest weight
                if(kbin_to_partition_mapping[kbin.first] == -1) {
                    PartitionID partition = result_queue.top();
                    ASSERT((partition >= 0) && (partition < num_partitions));

                    result_queue.increaseKeyBy(partition, kbin.second);
                    kbin_to_partition_mapping[kbin.first] = partition;
                }
            }

            // remap the partition IDs of the nodes
            for(PartitionID& id : partitions) {
                id = kbin_to_partition_mapping[id];
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

        for(const HypernodeID& hn : hg.nodes()) {
            nodes.push_back(hn);
        }
        ASSERT(hg.currentNumNodes() == nodes.size());

        std::sort(nodes.begin(), nodes.end(), [&hg](HypernodeID a, HypernodeID b) {
            return hg.nodeWeight(a) > hg.nodeWeight(b);
        });

        return nodes;
    }
}
}