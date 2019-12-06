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

#include "kahypar/definitions.h"
#include "kahypar/datastructure/binary_heap.h"

namespace kahypar {

namespace bin_packing {
    using kahypar::ds::BinaryMinHeap;

    /**
     * The hypernodes must be sorted in descending order of weight. 
     * The partitions parameter can be used to specify fixed vertices.
     *
     * Attention: the algorithm assumes that any fixed vertices are already partitioned in a balanced way.
     */
    static inline std::vector<PartitionID> two_level_packing(const Hypergraph& hg,
                                                              const std::vector<HypernodeID>& hypernodes,
                                                              const PartitionID& num_bins,
                                                              const PartitionID& rb_range_k,
                                                              std::vector<PartitionID> partitions = {}) {
        ALWAYS_ASSERT((num_bins > 0) && (rb_range_k % num_bins == 0), "num_bins must be positive and k must be a multiple of the number of bins.");
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

        if(partitions.empty()) {
            partitions.resize(hypernodes.size(), -1);
        }

        BinaryMinHeap<PartitionID, HypernodeWeight> result_queue(num_bins);
        for(PartitionID i = 0; i < num_bins; ++i) {
            result_queue.push(i, 0);
        }

        // At the first level, a packing with k bins is calculated ....
        BinaryMinHeap<PartitionID, HypernodeWeight> k_bin_queue(rb_range_k);
        for(PartitionID i = 0; i < rb_range_k; ++i) {
            k_bin_queue.push(i, 0);
        }
        bool fixed_vertex_contained = false;

        for(size_t i = 0; i < hypernodes.size(); ++i) {
            HypernodeWeight weight = hg.nodeWeight(hypernodes[i]);

            if(partitions[i] == -1) {
                // assign node to bin with lowest weight
                PartitionID kbin = k_bin_queue.top();
                // TODO remove
                ASSERT(kbin < rb_range_k);

                k_bin_queue.increaseKeyBy(kbin, weight);
                partitions[i] = kbin;
            } else {
                // fixed vertex: skip the k-bins level and instead assigned directly to the final bins
                ALWAYS_ASSERT(partitions[i] < num_bins, "Invalid partition ID for node.");

                fixed_vertex_contained = true;
                result_queue.increaseKeyBy(partitions[i], weight);

                // shift indizes for fixed vertices, so they are not remapped
                // TODO: this is rather ugly
                partitions[i] += rb_range_k;
            }
        }

        ASSERT(k_bin_queue.size() == rb_range_k);

        // ... and at the second level, the resulting k bins are packed into the final bins
        if((rb_range_k != num_bins) || fixed_vertex_contained) {
            std::vector<HypernodeWeight> kbin_weights(rb_range_k);
            std::vector<PartitionID> kbin_decreasing_ordering(rb_range_k);

            // read bins from queue
            // iteration must be in reverse order, because the smallest bins are on top
            for(ssize_t i = rb_range_k - 1; i >= 0; --i) {
                // TODO remove
                ASSERT(k_bin_queue.top() < rb_range_k);

                kbin_decreasing_ordering[i] = k_bin_queue.top();
                kbin_weights[i] = k_bin_queue.topKey();
                k_bin_queue.pop();
            }
            ASSERT(k_bin_queue.size() == 0);

            // pack bins to partitions
            std::vector<PartitionID> kbin_to_partition_mapping(rb_range_k, -1);

            for(size_t i = 0; i < rb_range_k; ++i) {
                // assign bin to partition with lowest weight
                PartitionID partition = result_queue.top();
                // TODO remove
                ASSERT(partition < num_bins);

                result_queue.increaseKeyBy(partition, kbin_weights[i]);

                PartitionID kbin = kbin_decreasing_ordering[i];
                kbin_to_partition_mapping[kbin] = partition;
            }

            // remap the partition IDs of the nodes
            for(PartitionID& id : partitions) {
                if(id < rb_range_k) {
                    id = kbin_to_partition_mapping[id];
                } else {
                    // reshift fixed vertex id
                    id -= rb_range_k;
                }

                ASSERT((id != -1) && (id < rb_range_k));
            }
        }

        return partitions;
    }
}
}