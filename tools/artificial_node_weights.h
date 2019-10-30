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

#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

// TODO types
// #include "kahypar/definitions.h"

namespace nodeweights {

uint32_t parseNumNodes(const std::string& hgr_filename);
std::vector<int32_t> createWeights(uint32_t num_nodes, uint32_t max_k, int min_imbalance_factor, int max_imbalance_factor);
std::vector<int32_t> createBigWeights(uint32_t num_heavy_nodes, int32_t lower_total_weight, uint32_t k, int factor);
void appendWeightsToHgr(const std::string& hgr_filename, const std::string& out_filename, const std::vector<int32_t>& weights);
}  // namespace nodeweights
