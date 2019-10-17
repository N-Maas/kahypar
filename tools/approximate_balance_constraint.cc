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
#include <sstream>
#include <string>
#include <map>
#include <queue>

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/macros.h"

using namespace kahypar;

typedef std::priority_queue<HypernodeWeight, std::vector<HypernodeWeight>, std::greater<HypernodeWeight>> BinQueue;

HypernodeWeight calculate_worst_fit_decreasing_bin_size(const std::map<HypernodeWeight, size_t, std::greater<HypernodeID>>& weight_distribution, HypernodeID k) {
  // init priority queue of bins with zeros
  BinQueue bins{ std::greater<HypernodeWeight>(), std::vector<HypernodeWeight>(k) };

  for(const std::pair<HypernodeWeight, size_t>& weight_entry : weight_distribution) {
    HypernodeWeight node_weight = weight_entry.first;
    size_t num_nodes = weight_entry.second;
    ASSERT(bins.size() == k);
    ASSERT(num_nodes > 0);

    if(node_weight == 0) {
      continue;
    }

    for(size_t i = 0; i < num_nodes; ++i) {
      // for every node, insert the node to the smallest bin and update the queue
      HypernodeWeight bin_weight = bins.top();
      bins.pop();
      bins.push(bin_weight + node_weight);
    }
  }

  // the result is the biggest bin, e.g. the last element in the queue
  while(bins.size() > 1) {
    bins.pop();
  }
  return bins.top();
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cout << "Missing argument" << std::endl;
    std::cout << "Usage: BalanceConstraint -k <# blocks> <.hgr>" << std::endl;
    exit(0);
  }
  std::string hgr_filename;
  HypernodeID num_blocks = 0;
  bool is_arg_num_blocks = false;
  for(int i = 1; i < argc; ++i) {
    if(is_arg_num_blocks) {
      num_blocks = std::stoul(argv[i]);
      is_arg_num_blocks = false;
    } else if(std::string(argv[i]) == std::string("-k")) {
      is_arg_num_blocks = true;
    } else {
      hgr_filename = std::string(argv[i]);
    }
  }

  Hypergraph hypergraph(io::createHypergraphFromFile(hgr_filename, num_blocks));

  std::map<HypernodeWeight, size_t, std::greater<HypernodeID>> weight_distribution;

  for (const auto hn : hypergraph.nodes()){
    ++weight_distribution[hypergraph.nodeWeight(hn)];
  }

  HypernodeWeight total_weight = hypergraph.totalWeight();
  HypernodeWeight weight_per_block = (total_weight + num_blocks - 1) / num_blocks;
  HypernodeWeight max_weight = weight_distribution.cbegin()->first;
  ASSERT(max_weight == hypergraph.weightOfHeaviestNode());

  std::cout << "printing distribution:" << std::endl;
  for(const auto& p : weight_distribution) {
    std::cout << p.first << " - " << p.second << std::endl;
  }

  std::cout << "number of nodes: " << hypergraph.currentNumNodes() << std::endl;
  std::cout << "maximum node weight: " << max_weight << std::endl;

  HypernodeWeight calculated_border = calculate_worst_fit_decreasing_bin_size(weight_distribution, num_blocks);
  std::cout << "total weight: " << total_weight << std::endl;
  std::cout << "average block weight: " << weight_per_block << std::endl;
  std::cout << "calculated border for block weight: " << calculated_border << std::endl;

  std::cout << " ... done!" << std::endl;
  return 0;
}
