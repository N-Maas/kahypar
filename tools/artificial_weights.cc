/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2020 Nikolai Maas <nikolai.maas@student.kit.edu>
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
#include <fstream>
#include <sstream>
#include <algorithm>

#include "kahypar/macros.h"
#include "kahypar/definitions.h"

using namespace kahypar;

static const double EXPECTED_NUM_HEAVY_NODES = 120;
static const double ECPECTED_RELATIVE_WEIGHT = 0.5;

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

HypernodeWeight uniformWeight(std::mt19937& rand, HypernodeID num_nodes) {
  HypernodeWeight max_heavy_weight = ceil(2.0 * ECPECTED_RELATIVE_WEIGHT / (1 - ECPECTED_RELATIVE_WEIGHT) * num_nodes / EXPECTED_NUM_HEAVY_NODES);
  bool is_heavy = std::uniform_int_distribution<>{ 1, static_cast<int>(num_nodes) }(rand) <= EXPECTED_NUM_HEAVY_NODES;
  if (is_heavy) {
    return std::uniform_int_distribution<>{ 1, max_heavy_weight }(rand);
  } else {
     return 1;
  }
}

std::vector<double> powerLawDistribution(HypernodeID num_nodes) {
    std::vector<double> probabilities;
    probabilities.push_back(0);
    HypernodeWeight max_heavy_weight = ceil(3.0 * ECPECTED_RELATIVE_WEIGHT / (1 - ECPECTED_RELATIVE_WEIGHT) * num_nodes / EXPECTED_NUM_HEAVY_NODES);
    for (int i = 1; i < max_heavy_weight; ++i) {
        probabilities.push_back(1.0 / (std::pow(static_cast<double>(i), 2.5)));
    }
    return std::move(probabilities);
}

void appendWeightsToHgr(const std::string& hgr_filename, const std::string& out_filename, std::vector<kahypar::HypernodeWeight>& weights) {
  std::random_shuffle(weights.begin(), weights.end());
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

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cout << "Usage: ArtificialNodeWeights <in.hgr> <out.hgr>" << std::endl;
    exit(0);
  }

  std::random_device rd;
  std::mt19937 rand(rd());
  std::string input_file(argv[1]);
  std::string output_file(argv[2]);

  HypernodeID num_nodes = parseNumNodes(input_file);
  std::vector<HyperedgeWeight> weights;
  weights.reserve(num_nodes);
  // auto prob_vec = powerLawDistribution(num_nodes);
  // std::discrete_distribution<HypernodeWeight> distribution(prob_vec.cbegin(), prob_vec.cend());
  for (HypernodeID i = 0; i < num_nodes; ++i) {
    // weights.push_back(distribution(rand));
    weights.push_back(uniformWeight(rand, num_nodes));
  }

  appendWeightsToHgr(input_file, output_file, weights);
}
