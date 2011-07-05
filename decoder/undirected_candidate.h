#ifndef _UNDIRECTED_CANDIDATE_H_
#define _UNDIRECTED_CANDIDATE_H_

#include <iostream>
#include "ff.h"

struct ModelSet;
struct Hypergraph;
struct SentenceMetadata;

/////////////////////////////////////
//added here because model uses UCand as generalization of Edge
//therefore need to be imported in ff.cc
using namespace std;

typedef ValueArray<uint8_t> FFState;
typedef vector<FFState> FFStates;

struct UCandidate;
typedef vector<UCandidate*> UCandidateHeap;
typedef vector<UCandidate*> UCandidateList;

struct UCandidate {
  int node_index_;                     // -1 until incorporated
                                       // into the +LM forest
  const Hypergraph::Edge* in_edge_;    // in -LM forest
  FFState state_;

  FeatureVector est_vals_;
  //TODO LG remove
  Hypergraph::Edge out_edge_;

  //TODO LG replace with pointer structure
  const SmallVectorInt j_;

  //vit_prob_ and est_prob_ are not updated in LG training
  prob_t vit_prob_;            // these are fixed until the cand
                               // is popped, then they may be updated
  prob_t est_prob_;

  prob_t action_prob_;

  UCandidate(const Hypergraph::Edge& e,
            const SmallVectorInt& j,
            const vector<UCandidateList>& D,
            const FFStates& node_states,
            const SentenceMetadata& smeta,
            const ModelSet& models,
            bool is_goal);

  // used to query uniqueness
  UCandidate(const Hypergraph::Edge& e,
            const SmallVectorInt& j);

  bool IsIncorporatedIntoHypergraph() const;

  void InitializeUCandidate(
                           const SentenceMetadata& smeta,
                           const vector<vector<UCandidate*> >& D,
                           const FFStates& node_states,
                           const ModelSet& models,
                           const bool is_goal);
};

ostream& operator<<(ostream& os, const UCandidate& cand);

////////////////////////////////////////////

#endif
