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
#include <sstream>

#include "kahypar/macros.h"
#include "tools/artificial_node_weights.h"

using namespace kahypar;

const int WEIGHT_VARIANCE = 3;

namespace nodeweights {

HypernodeID parseNumNodes(const std::string& hgr_filename) {
  std::ifstream file(hgr_filename);
  std::string line;
  std::getline(file, line);
  file.close();

  std::istringstream sstream(line);
  std::string num_edges, num_nodes;
  sstream >> num_edges >> num_nodes;
  return std::stoul(num_nodes);
}

std::vector<HypernodeWeight> createWeights(HypernodeID num_nodes, HypernodeID max_k, int imbalance_factor) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::vector<HypernodeWeight> weights;
  weights.reserve(num_nodes);

  HypernodeID num_heavy_nodes = std::uniform_int_distribution<>{ (3 * max_k + 1) / 2, 3 * max_k }(gen);
  HypernodeID num_light_nodes = num_nodes - num_heavy_nodes;
  // TODO: random weights?
  weights.resize(1, num_light_nodes);

  std::vector<HypernodeWeight> big_weights = createBigWeights(num_heavy_nodes, num_light_nodes, max_k, imbalance_factor);
  weights.insert(weights.end(), big_weights.cbegin(), big_weights.cend());

  return std::move(weights);
}

std::vector<HypernodeWeight> createBigWeights(HypernodeID num_heavy_nodes, HypernodeWeight lower_total_weight, HypernodeID k, int factor) {
  std::random_device rd;
  std::mt19937 gen(rd());

  // magic formula derived from heuristic
  HypernodeWeight min_weight = factor * ((WEIGHT_VARIANCE + 1) * num_heavy_nodes
    * lower_total_weight / k - lower_total_weight / (num_heavy_nodes - 1));

  std::uniform_int_distribution<> distribution{ min_weight, WEIGHT_VARIANCE * min_weight };
  std::vector<HypernodeWeight> weights;
  weights.reserve(num_heavy_nodes);
  for (HypernodeID i = 0; i < num_heavy_nodes; ++i) {
    weights.push_back(distribution(gen));
  }

  return std::move(weights);
}

void appendWeightsToHgr(const std::string& hgr_filename, const std::vector<kahypar::HypernodeWeight>& weights) {
  std::ofstream output(hgr_filename, std::ios::app);

  for (const HypernodeWeight &w : weights) {
    output << w << std::endl;
  }

  output.close();
}

}  // namespace nodeweights
