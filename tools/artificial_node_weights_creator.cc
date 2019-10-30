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

#include <random>
#include <iostream>

#include "tools/artificial_node_weights.h"

using namespace nodeweights;

int main(int argc, char *argv[]) {
  if (argc < 6) {
    std::cout << "Missing argument" << std::endl;
    std::cout << "Usage: ArtificialNodeWeights -k <# blocks> -f <max imbalance factor> <.hgr>" << std::endl;
    exit(0);
  }
  std::vector<std::string> hgr_filenames;
  uint32_t num_blocks = 0;
  int max_imbalance_factor = 0;

  bool is_arg_num_blocks = false;
  bool is_arg_factor = false;
  for (int i = 1; i < argc; ++i) {
    if (is_arg_num_blocks) {
      num_blocks = std::stoul(argv[i]);
      is_arg_num_blocks = false;
    } else if (is_arg_factor) {
      max_imbalance_factor = std::stoul(argv[i]);
      is_arg_factor = false;
    } else if (std::string(argv[i]) == std::string("-k")) {
      is_arg_num_blocks = true;
    } else if (std::string(argv[i]) == std::string("-f")) {
      is_arg_factor = true;
    } else {
      hgr_filenames.push_back(std::string(argv[i]));
    }
  }

  for(const auto& file : hgr_filenames) {
    std::string new_file = file.substr(0, file.size() - 4) + ".weighted.hgr";

    uint32_t num_nodes = parseNumNodes(file);
    std::vector<int32_t> weights = createWeights(num_nodes, num_blocks, (max_imbalance_factor + 2) / 3, max_imbalance_factor);

    appendWeightsToHgr(file, new_file, weights);
  }

  return 0;
}
