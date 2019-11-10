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
#include <algorithm>

#include <boost/integer.hpp>

#include "kahypar/macros.h"
#include "tools/artificial_node_weights.h"
#include "kahypar/definitions.h"

using namespace kahypar;

typedef boost::int_t<sizeof(HypernodeWeight) * CHAR_BIT * 2>::least
    LongHypernodeWeight;
typedef boost::uint_t<sizeof(HypernodeID) * CHAR_BIT * 2>::least
    LongHypernodeID;

const long WEIGHT_VARIANCE = 2;

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

std::vector<HypernodeWeight> createWeights(HypernodeID num_nodes, HypernodeID max_k, int min_imbalance_factor, int max_imbalance_factor) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::vector<HypernodeWeight> weights;
  weights.reserve(num_nodes);

  HypernodeID num_heavy_nodes = std::uniform_int_distribution<>{ (6 * max_k + 1) / 5, 2 * max_k }(gen);
  HypernodeID num_light_nodes = num_nodes - num_heavy_nodes;

  // Generate small weights with randomized distribution of 1, 3, 5, 17
  int one_probability = std::uniform_int_distribution<>{ 4, 12 }(gen);
  int three_probability = std::uniform_int_distribution<>{ 3, 6 }(gen);
  int five_probability = std::uniform_int_distribution<>{ 2, 5 }(gen);
  int seventeen_probability = std::uniform_int_distribution<>{ 1, 2 }(gen);
  std::discrete_distribution<HypernodeWeight> light_weights_distribution{double(one_probability), 0, double(three_probability),
    0, double(five_probability), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, double(seventeen_probability)};
  HypernodeWeight total_weight = 0;

  for(HypernodeID i = 0; i < num_light_nodes; ++i) {
    HypernodeWeight weight = light_weights_distribution(gen);
    total_weight += weight;
    weights.push_back(weight);
  }

  int imbalance_factor = std::uniform_int_distribution<>{ min_imbalance_factor, max_imbalance_factor }(gen);
  std::vector<HypernodeWeight> big_weights = createBigWeights(num_heavy_nodes, total_weight, max_k, imbalance_factor);
  std::for_each(big_weights.cbegin(), big_weights.cend(), [&](auto w) {weights.push_back(w);});

  std::random_shuffle(weights.begin(), weights.end());
  return std::move(weights);
}

std::vector<HypernodeWeight> createBigWeights(HypernodeID num_heavy_nodes, HypernodeWeight lower_total_weight, HypernodeID k, int factor) {
  std::random_device rd;
  std::mt19937 gen(rd());

  // magic formula derived from heuristic
  HypernodeWeight min_weight = long(factor) * ((WEIGHT_VARIANCE + 1) * LongHypernodeID(num_heavy_nodes)
    * LongHypernodeWeight(lower_total_weight) / LongHypernodeID(k) - LongHypernodeWeight(lower_total_weight)) / (LongHypernodeID(num_heavy_nodes) - 1);

  std::uniform_int_distribution<> distribution{ min_weight, WEIGHT_VARIANCE * min_weight };
  std::vector<HypernodeWeight> weights;
  weights.reserve(num_heavy_nodes);
  for (HypernodeID i = 0; i < num_heavy_nodes; ++i) {
    weights.push_back(distribution(gen));
  }

  return std::move(weights);
}

void appendWeightsToHgr(const std::string& hgr_filename, const std::string& out_filename, const std::vector<kahypar::HypernodeWeight>& weights) {
  std::ifstream input(hgr_filename);
  std::ofstream output(out_filename);

  // change 1. line to weighted format
  std::string line;
  std::getline(input, line);
  std::istringstream sstream(line);
  std::string num_edges, num_nodes;
  sstream >> num_edges >> num_nodes;
  output << num_edges << " " << num_nodes << " 10 " << std::endl;
  output << input.rdbuf();

  input.close();
  output.close();

  std::ofstream append_stream(out_filename, std::ios::app);

  for (const HypernodeWeight &w : weights) {
    append_stream << w << std::endl;
  }

  append_stream.close();
}

}  // namespace nodeweights
