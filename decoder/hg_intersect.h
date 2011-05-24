#ifndef _HG_INTERSECT_H_
#define _HG_INTERSECT_H_

#include <vector>

#include "lattice.h"

class Hypergraph;
struct HG {
  static bool Intersect(const Lattice& target, Hypergraph* hg);
  static bool HighlightIntersection(const Lattice& target, const Hypergraph& hg, std::vector<bool>* correct_edges);
};

#endif
