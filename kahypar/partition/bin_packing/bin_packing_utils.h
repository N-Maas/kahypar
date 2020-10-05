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

#include "kahypar/partition/context.h"
#include "kahypar/datastructure/binary_heap.h"

namespace kahypar {
namespace bin_packing {
  using kahypar::ds::BinaryMinHeap;

  class PartitionMapping {
    public:
      PartitionMapping(PartitionID num_bins) : _parts() {
        ASSERT(num_bins > 0, "Number of bins must be positive.");

        _parts.resize(static_cast<size_t>(num_bins), -1);
      }

      void setPart(PartitionID bin, PartitionID part) {
        size_t index = static_cast<size_t>(bin);
        ASSERT(bin >= 0 && index < _parts.size(), "Invalid bin id: " << V(bin));
        ASSERT(_parts[index] == -1 || _parts[index] == part,
             "Bin already assigned to different part");

        _parts[index] = part;
      }

      bool isFixedBin(PartitionID bin) const {
        return binPartition(bin) != -1;
      }

      PartitionID binPartition(PartitionID bin) const {
        ASSERT(bin >= 0 && static_cast<size_t>(bin) < _parts.size(), "Invalid bin id: " << V(bin));

        return _parts[static_cast<size_t>(bin)];
      }

      void applyMapping(std::vector<PartitionID>& elements) {
        for (PartitionID& id : elements) {
          ASSERT(binPartition(id) >= 0, "Bin not assigned: " << id);
          id = binPartition(id);
        }
      }

    private:
      std::vector<PartitionID> _parts;
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

      void addFixedVertex(PartitionID bin, PartitionID part, HypernodeWeight weight) {
        _bins_to_parts.setPart(bin, part);
        _alg.addWeight(bin, weight);
      }

      HypernodeWeight binWeight(PartitionID bin) const {
        return _alg.binWeight(bin);
      }

      // ... and at the second level, the resulting k bins are packed into the final parts.
      // Returns the partition mapping for the bins and a vector of the resulting part weights
      std::pair<PartitionMapping, std::vector<HypernodeWeight>> applySecondLevel(const std::vector<HypernodeWeight>& max_allowed_part_weights,
                                                                                 const std::vector<PartitionID>& num_bins_per_part) const {
        ASSERT(num_bins_per_part.size() == max_allowed_part_weights.size());

        const PartitionID num_parts = static_cast<PartitionID>(max_allowed_part_weights.size());
        std::vector<PartitionID> bin_counts(max_allowed_part_weights.size(), 0);
        PartitionMapping mapping = _bins_to_parts;

        const HypernodeWeight max_part_weight = *std::max_element(max_allowed_part_weights.cbegin(),
                                                max_allowed_part_weights.cend());
        BPAlgorithm partition_packer(num_parts, max_part_weight);
        for (PartitionID i = 0; i < num_parts; ++i) {
          partition_packer.addWeight(i, max_part_weight - max_allowed_part_weights[i]);
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
          if (!mapping.isFixedBin(bin)) {
            const PartitionID part = partition_packer.insertElement(_alg.binWeight(bin));
            ASSERT((part >= 0) && (part < num_parts));

            const size_t part_idx = static_cast<size_t>(part);
            bin_counts[part_idx]++;
            if (bin_counts[part_idx] >= num_bins_per_part[part_idx]) {
              partition_packer.lockBin(part);
            }

            mapping.setPart(bin, part);
          }
        }

        // collect part weights
        std::vector<HypernodeWeight> part_weights;
        part_weights.reserve(num_parts);
        for (PartitionID i = 0; i < num_parts; ++i) {
          part_weights.push_back(partition_packer.binWeight(i) + max_allowed_part_weights[i] - max_part_weight);
        }

        return {std::move(mapping), std::move(part_weights)};
      }

    private:
      BPAlgorithm _alg;
      PartitionMapping _bins_to_parts;
  };

  // Returns the current hypernodes sorted in descending order of weight.
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

  // Preassigns the fixed vertices to the packer with a first fit packing.
  template< class BPAlgorithm >
  static inline void preassignFixedVertices(const Hypergraph& hg, const std::vector<HypernodeID>& nodes, std::vector<PartitionID>& parts,
                                            TwoLevelPacker<BPAlgorithm>& packer, PartitionID k, PartitionID rb_range_k) {
    const HypernodeWeight avg_bin_weight = (hg.totalWeight() + rb_range_k - 1) / rb_range_k;
    const PartitionID bins_per_part = (rb_range_k + k - 1) / k;

    for (size_t i = 0; i < nodes.size(); ++i) {
      const HypernodeID hn = nodes[i];

      if (hg.isFixedVertex(hn)) {
        const HypernodeWeight weight = hg.nodeWeight(hn);
        const PartitionID part_id = hg.fixedVertexPartID(hn);
        const PartitionID start_index = part_id * bins_per_part;
        PartitionID assigned_bin = start_index;

        if (bins_per_part > 1) {
          for (PartitionID i = start_index; i < std::min(start_index + bins_per_part, rb_range_k); ++i) {
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
        parts[i] = assigned_bin;
      }
    }

  }

  static inline HypernodeWeight currentMaxBin(const Hypergraph& hypergraph, const PartitionID& rb_range_k) {
    BinaryMinHeap<PartitionID, HypernodeWeight> queue(rb_range_k);
    for (PartitionID j = 0; j < rb_range_k; ++j) {
        queue.push(j, 0);
      }

    const std::vector<HypernodeID> hypernodes = extractNodesWithDescendingWeight(hypergraph);
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

// Assumes that the final parts are of equal size (is this always true?).
static inline HypernodeWeight resultingMaxBin(const Hypergraph& hypergraph, const Context& context) {
  ASSERT(!context.partition.use_individual_part_weights);

  const PartitionID num_parts = context.initial_partitioning.k;

  // initialize queues
  std::vector<BinaryMinHeap<PartitionID, HypernodeWeight>> part_queues;
  part_queues.reserve(num_parts);
  for (PartitionID i = 0; i < num_parts; ++i) {
    const PartitionID current_k = context.initial_partitioning.num_bins_per_part[i];
    BinaryMinHeap<PartitionID, HypernodeWeight> queue(current_k);
    for (PartitionID j = 0; j < current_k; ++j) {
        queue.push(j, 0);
    }
    part_queues.push_back(std::move(queue));
  }

  // assign nodes
  const std::vector<HypernodeID> hypernodes = bin_packing::extractNodesWithDescendingWeight(hypergraph);
  for (const HypernodeID& hn : hypernodes) {
    const PartitionID part_id = hypergraph.partID(hn);

    ALWAYS_ASSERT(part_id >= 0 && part_id < num_parts,
                  "Node not assigned or part id " << part_id << " invalid: " << hn);

    const PartitionID bin = part_queues[part_id].top();
    part_queues[part_id].increaseKeyBy(bin, hypergraph.nodeWeight(hn));
  }

  HypernodeWeight max = 0;
  for(auto& queue : part_queues) {
    // the maximum is the biggest bin, i.e. the last element in the queue
    while (queue.size() > 1) {
      queue.pop();
    }
    max = std::max(queue.getKey(queue.top()), max);
  }

  return max;
}
} // namespace bin_packing
} // namespace kahypar
