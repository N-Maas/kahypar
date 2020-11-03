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

#include "kahypar/definitions.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/macros.h"

using namespace kahypar;

HypernodeID removeHeavyNodes(Hypergraph& hypergraph, HypernodeID k, double epsilon) {
  HypernodeWeight allowed_node_weight = std::floor((1 + epsilon) * std::ceil(
                                        static_cast<double>(hypergraph.totalWeight()) / k));
  std::vector<HypernodeID> to_remove;

    for (const HypernodeID& hn : hypergraph.nodes()) {
      if (hypergraph.nodeWeight(hn) > allowed_node_weight) {
        to_remove.push_back(hn);
      }
    }
    for (const HypernodeID& hn : to_remove) {
      hypergraph.removeNode(hn);
    }

  return k - to_remove.size();
}

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cout << "Wrong number of arguments" << std::endl;
    std::cout << "Usage: RemoveHeavyNodes <.hgr> <output_file> <k> <e>" << std::endl;
    exit(0);
  }

  std::string hgr_filename(argv[1]);
  std::string output_filename(argv[2]);
  HypernodeID k = std::stoul(argv[3]);
  double epsilon = std::stod(argv[4]);

  Hypergraph hypergraph(io::createHypergraphFromFile(hgr_filename, k));

  HypernodeID new_k = removeHeavyNodes(hypergraph, k, epsilon);
  for (const HyperedgeID& he : hypergraph.edges()) {
    auto pins_start_end = hypergraph.pins(he);
    if (pins_start_end.first == pins_start_end.second) {
      hypergraph.removeEdge(he);
    }
  }

  auto modified_hg = ds::reindex(hypergraph);
  std::cout << "k=" << new_k << std::endl;
  std::cout << "removed=" << (k - new_k) << std::endl;

  // change hypergraph type to avoid writing edge weights unnecessarily
  bool hasEdgeWeights = false;
  for (const HyperedgeID& edge : hypergraph.edges()) {
    hasEdgeWeights |= hypergraph.edgeWeight(edge) != 1;
  }
  if (!hasEdgeWeights) {
    modified_hg.first->setType(HypergraphType::NodeWeights);
  }

  io::writeHypergraphFile(*modified_hg.first, output_filename);

  return 0;
}
