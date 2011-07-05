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

  const int algorithm; // 0 = full intersection, 1 = cube pruning
  const int pop_limit; // max number of pops off the heap at each node
  IntersectionConfiguration(int alg, int k) : algorithm(alg), pop_limit(k) {}
  IntersectionConfiguration(exhaustive_t /* t */) : algorithm(0), pop_limit() {}
};

inline std::ostream& operator<<(std::ostream& os, const IntersectionConfiguration& c) {
  if (c.algorithm == 0) { os << "FULL"; }
  else if (c.algorithm == 1) { os << "CUBE:k=" << c.pop_limit; }
  else if (c.algorithm == 2) { os << "N_ALGORITHMS"; }
  else os << "OTHER";
  return os;
}

/////////////////////////////////////
//added here because model uses UCand as generalization of Edge
//therefore need to be imported in ff.cc
//using namespace std;
//
//struct UCandidate;
//typedef vector<UCandidate*> UCandidateHeap;
//typedef vector<UCandidate*> UCandidateList;

//struct UCandidate {
//  int node_index_;                     // -1 until incorporated
//                                       // into the +LM forest
//  const Hypergraph::Edge* in_edge_;    // in -LM forest
//  FFState state_;
//
//  FeatureVector est_vals_;
//  //TODO LG remove
//  Hypergraph::Edge out_edge_;
//
//  //TODO LG replace with pointer structure
//  const SmallVectorInt j_;
//
//  //vit_prob_ and est_prob_ are not updated in LG training
//  prob_t vit_prob_;            // these are fixed until the cand
//                               // is popped, then they may be updated
//  prob_t est_prob_;
//
//  prob_t action_prob_;
//
//  UCandidate(const Hypergraph::Edge& e,
//            const SmallVectorInt& j,
//            const vector<UCandidateList>& D,
//            const FFStates& node_states,
//            const SentenceMetadata& smeta,
//            const ModelSet& models,
//            bool is_goal) :
//      node_index_(-1),
//      in_edge_(&e),
//      j_(j) {
//    InitializeUCandidate(smeta, D, node_states, models, is_goal);
//  }
//
//  // used to query uniqueness
//  UCandidate(const Hypergraph::Edge& e,
//            const SmallVectorInt& j) : in_edge_(&e), j_(j) {}
//
//  bool IsIncorporatedIntoHypergraph() const {
//    return node_index_ >= 0;
//  }
//
//  void InitializeUCandidate(
//                           const SentenceMetadata& smeta,
//                           const vector<vector<UCandidate*> >& D,
//                           const FFStates& node_states,
//                           const ModelSet& models,
//                           const bool is_goal) {
//    const Hypergraph::Edge& in_edge = *in_edge_;
//    out_edge_.rule_ = in_edge.rule_;
//    out_edge_.feature_values_ = in_edge.feature_values_;
//    out_edge_.i_ = in_edge.i_;
//    out_edge_.j_ = in_edge.j_;
//    out_edge_.prev_i_ = in_edge.prev_i_;
//    out_edge_.prev_j_ = in_edge.prev_j_;
//    Hypergraph::TailNodeVector& tail = out_edge_.tail_nodes_;
//    tail.resize(j_.size());
//    prob_t p = prob_t::One();
//    // cerr << "\nEstimating application of " << in_edge.rule_->AsString() << endl;
//    for (int i = 0; i < tail.size(); ++i) {
//      const UCandidate& ant = *D[in_edge.tail_nodes_[i]][j_[i]];
//      assert(ant.IsIncorporatedIntoHypergraph());
//      tail[i] = ant.node_index_;
//      p *= ant.vit_prob_;
//    }
//    prob_t edge_estimate = prob_t::One();
//    if (is_goal) {
//      assert(tail.size() == 1);
//      const FFState& ant_state = node_states[tail.front()];
//      models.AddFinalFeatures(ant_state, &out_edge_, smeta);
//    } else {
//      models.AddFeaturesToUCandidate(smeta, node_states, /*this,*/ &out_edge_, &state_, &edge_estimate);
//    }
//    vit_prob_ = out_edge_.edge_prob_ * p;
//    est_prob_ = vit_prob_ * edge_estimate;
//    action_prob_ = out_edge_.edge_prob_ * edge_estimate;
//  }
//};
//
//std::ostream& operator<<(ostream& os, const UCandidate& cand);


////////////////////////////////////////////


void ApplyModelSet(const Hypergraph& in,
                   const SentenceMetadata& smeta,
                   /*const*/ ModelSet& models,
                   const IntersectionConfiguration& config,
                   Hypergraph* out);

#endif
