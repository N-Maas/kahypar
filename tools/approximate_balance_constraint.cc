/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Nikolai Maas <nikolai.maas@kit.edu>
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

#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>

#include <boost/integer.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/macros.h"

using namespace kahypar;

typedef std::priority_queue<HypernodeWeight, std::vector<HypernodeWeight>,
                            std::greater<HypernodeWeight>>
    BinQueue;
typedef boost::int_t<sizeof(HypernodeWeight) * CHAR_BIT * 2>::least
    LongHypernodeWeight;

HypernodeWeight calculate_worst_fit_decreasing_bin_size(
    const std::map<HypernodeWeight, size_t, std::greater<HypernodeID>>
        &weight_distribution,
    HypernodeID k) {
  // init priority queue of bins with zeros
  BinQueue bins{std::greater<HypernodeWeight>(),
                std::vector<HypernodeWeight>(k)};

  for (const std::pair<HypernodeWeight, size_t> &weight_entry :
       weight_distribution) {
    HypernodeWeight node_weight = weight_entry.first;
    size_t num_nodes = weight_entry.second;
    ASSERT(bins.size() == k);
    ASSERT(num_nodes > 0);

    if (node_weight == 0) {
      continue;
    }

    for (size_t i = 0; i < num_nodes; ++i) {
      // for every node, insert the node to the smallest bin and update the
      // queue
      HypernodeWeight bin_weight = bins.top();
      bins.pop();
      bins.push(bin_weight + node_weight);
    }
  }

  // the result is the biggest bin, i.e. the last element in the queue
  while (bins.size() > 1) {
    bins.pop();
  }
  return bins.top();
}

HypernodeWeight calculate_recursive_wfd_bin_size(
    const std::map<HypernodeWeight, size_t, std::greater<HypernodeID>>
        &weight_distribution,
    HypernodeID k) {
  if(k == 1) {
    HypernodeWeight sum = 0;
    for (const std::pair<HypernodeWeight, size_t> &entry :
       weight_distribution) {
         sum += entry.first * entry.second;
       }
    return sum;
  }

  ASSERT(k % 2 == 0);
  std::map<HypernodeWeight, size_t, std::greater<HypernodeID>> left_part;
  HypernodeWeight left_sum = 0;
  std::map<HypernodeWeight, size_t, std::greater<HypernodeID>> right_part;
  HypernodeWeight right_sum = 0;

  for (const std::pair<HypernodeWeight, size_t> &weight_entry :
       weight_distribution) {
    HypernodeWeight node_weight = weight_entry.first;
    size_t num_nodes = weight_entry.second;

    if (node_weight == 0) {
      continue;
    }

    for (size_t i = 0; i < num_nodes; ++i) {
      // for every node, insert the node to the smaller partition and update the sum
      if(left_sum <= right_sum) {
        left_sum += node_weight;
        left_part[node_weight] += 1;
      } else {
        right_sum += node_weight;
        right_part[node_weight] += 1;
      }
    }
  }

  return std::max(calculate_recursive_wfd_bin_size(left_part, k / 2), calculate_recursive_wfd_bin_size(right_part, k / 2));
}

// max{a - 1 / (k - 1) * sum_{a_i < a}(a_i)}
template <class RatingFunction>
HypernodeWeight calculate_dominating_element_heuristic(
    RatingFunction rating,
    const std::map<HypernodeWeight, size_t, std::greater<HypernodeID>>
        &weight_distribution,
    HypernodeWeight total_weight, HypernodeID k) {
  HypernodeWeight lower_sum = 0;
  HypernodeWeight maximum = std::numeric_limits<HypernodeWeight>::min();

  for (const std::pair<HypernodeWeight, size_t> &weight_entry :
       boost::adaptors::reverse(weight_distribution)) {
    HypernodeWeight node_weight = weight_entry.first;
    size_t num_nodes = weight_entry.second;

    maximum =
        std::max(maximum, rating(node_weight, lower_sum, total_weight, k));
    lower_sum += node_weight * num_nodes;
  }

  return maximum;
}

void eval_file(std::string hgr_filename, HypernodeID num_blocks) {
  std::cout << "> " << hgr_filename << " <" << std::endl;

  Hypergraph hypergraph(io::createHypergraphFromFile(hgr_filename, num_blocks));

  std::map<HypernodeWeight, size_t, std::greater<HypernodeID>>
      weight_distribution;

  for (const auto hn : hypergraph.nodes()) {
    ++weight_distribution[hypergraph.nodeWeight(hn)];
  }

  HypernodeWeight total_weight = hypergraph.totalWeight();
  HypernodeWeight weight_per_block =
      (total_weight + num_blocks - 1) / num_blocks;
  HypernodeWeight max_weight = weight_distribution.cbegin()->first;
  ASSERT(max_weight == hypergraph.weightOfHeaviestNode());

  std::cout << "number of nodes:                    "
            << hypergraph.currentNumNodes() << std::endl;
  std::cout << "maximum node weight:                " << max_weight
            << std::endl;
  std::cout << "total weight:                       " << total_weight
            << std::endl;
  std::cout << "average block weight:               " << weight_per_block
            << std::endl;
  std::cout << "---------- ------------ ---------- ----------" << std::endl;

  HypernodeWeight calculated_border =
      calculate_worst_fit_decreasing_bin_size(weight_distribution, num_blocks);

  std::cout << "calculated border for block weight: " << calculated_border
            << std::endl;
  double e_calculated = static_cast<double>(calculated_border - weight_per_block) / static_cast<double>(weight_per_block);
  HypernodeWeight border_gamma_1 = std::max(weight_per_block, max_weight);
  double e_gamma_1 = static_cast<double>(border_gamma_1 - weight_per_block) / static_cast<double>(weight_per_block);
  HypernodeWeight border_gamma_2 = weight_per_block + max_weight;
  double e_gamma_2 = static_cast<double>(border_gamma_2 - weight_per_block) / static_cast<double>(weight_per_block);

  std::cout << "RESULT" << " graph=" << hgr_filename.substr(hgr_filename.find_last_of("/\\") + 1) << " k=" << num_blocks
            << " num_nodes=" << hypergraph.currentNumNodes() << " max=" << max_weight << " total=" << total_weight
            << " avg=" << weight_per_block << " border_calc=" << calculated_border << " e_calc=" << e_calculated
            << " border_gamma_1=" << border_gamma_1 << " e_gamma_1=" << e_gamma_1 << " border_gamma_2=" << border_gamma_2
            << " e_gamma_2=" << e_gamma_2 << std::endl;

  std::cout << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "Missing argument" << std::endl;
    std::cout << "Usage: BalanceConstraint -k <# blocks> <.hgr>" << std::endl;
    exit(0);
  }
  std::vector<std::string> hgr_filenames;
  HypernodeID num_blocks = 0;
  bool is_arg_num_blocks = false;
  for (int i = 1; i < argc; ++i) {
    if (is_arg_num_blocks) {
      num_blocks = std::stoul(argv[i]);
      is_arg_num_blocks = false;
    } else if (std::string(argv[i]) == std::string("-k")) {
      is_arg_num_blocks = true;
    } else {
      hgr_filenames.push_back(std::string(argv[i]));
    }
  }

  for(const auto& file : hgr_filenames) {
    eval_file(file, num_blocks);
  }

  return 0;
}
