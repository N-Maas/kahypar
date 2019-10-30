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

#include "kahypar/definitions.h"

namespace nodeweights {

kahypar::HypernodeID parseNumNodes(const std::string& hgr_filename);
std::vector<kahypar::HypernodeWeight> createWeights(kahypar::HypernodeID num_nodes,
  kahypar::HypernodeID max_k, int imbalance_factor);
std::vector<kahypar::HypernodeWeight> createBigWeights(kahypar::HypernodeID num_heavy_nodes,
  kahypar::HypernodeWeight lower_total_weight, kahypar::HypernodeID k, int factor);
void appendWeightsToHgr(const std::string& hgr_filename, const std::vector<kahypar::HypernodeWeight>& weights);
}  // namespace nodeweights
