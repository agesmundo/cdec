#ifndef _APPLY_MODELS_H_
#define _APPLY_MODELS_H_

#include <iostream>

//added for LG
#include "hg.h"
#include "ff.h"

struct ModelSet;
struct Hypergraph;
struct SentenceMetadata;

struct exhaustive_t {};

struct IntersectionConfiguration {
enum {
  FULL,
  CUBE,
  GREEDY_UNDIRECTED,
  GREEDY_UNDIRECTED_TRAINING,
  N_ALGORITHMS
};

  const int algorithm; // 0 = full intersection, 1 = cube pruning ,, 2 = greedy undirected , 3 greedy undirected training
  const int pop_limit; // max number of pops off the heap at each node
  const string weight_file_name_; //out train file name in case of greedy undirected training
  IntersectionConfiguration(int alg, int k, string fn="") : algorithm(alg), pop_limit(k), weight_file_name_(fn){}
  IntersectionConfiguration(exhaustive_t /* t */) : algorithm(0), pop_limit() {}
};

inline std::ostream& operator<<(std::ostream& os, const IntersectionConfiguration& c) {
  if (c.algorithm == 0) { os << "FULL"; }
  else if (c.algorithm == 1) { os << "CUBE:k=" << c.pop_limit; }
  else if (c.algorithm == 2) { os << "GREEDY_UNDIRECTED"; }
  else if (c.algorithm == 3) { os << "GREEDY_UNDIRECTED_TRAINING"; }
  else if (c.algorithm == 4) { os << "N_ALGORITHMS"; }
  else os << "OTHER";
  return os;
}

void ApplyModelSet(const Hypergraph& in,
                   const SentenceMetadata& smeta,
                   /*const*/ ModelSet& models,
                   const IntersectionConfiguration& config,
                   Hypergraph* out);

#endif
