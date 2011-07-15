#ifndef _UNDIRECTED_CANDIDATE_H_
#define _UNDIRECTED_CANDIDATE_H_

#include <iostream>
#include "ff.h"

struct ModelSet;
struct Hypergraph;
struct SentenceMetadata;

// Define the following macro if you want to see lots of debugging output
// when running Guided Undirected Greedy decofing
#define DEBUG_GU
//#undef DEBUG_GU

/////////////////////////////////////
//added here because model uses UCand as generalization of Edge
//therefore need to be imported in ff.cc
using namespace std;

typedef ValueArray<uint8_t> FFState;
typedef vector<FFState> FFStates;

struct UCandidate;
typedef vector<UCandidate*> UCandidateHeap;//TODO GU create struct with methods
typedef vector<UCandidate*> UCandidateList;
typedef SmallVector<UCandidate*,3> LinksVector;
typedef pair<int,FFState*>  Node2State;

struct UCandidate {
  //int ucand_index_;                     // -1 until popped from queue

  const Hypergraph::Edge* in_edge_;    // in -LM forest
  Node2State** outgoing_states_;         //in_node_id 2 state seen from that node //array max 2 elements (simple map)
  int states_size_; //TODO make constant!
//  FFState state_;

  FeatureVector est_vals_;
  FeatureVector feature_values_;

  //links to context
  //NB!!:
  //0 HEAD,
  //1 FIRST CHILD (left)
  //2 SECOND CHILD
  LinksVector context_links_; //NULL if no link, *UCand if link, -1 for head link if Goal
  int source_link_; //context_links_(id) for the source UCand (-1 for starting leaf)

  //vit_prob_ and est_prob_ are not updated in LG training
  //prob_t vit_prob_;            // these are fixed until the cand
                               // is popped, then they may be updated
  //prob_t est_prob_;

  prob_t action_prob_;

  UCandidate(const Hypergraph::Edge& e,
		    const LinksVector& lv,
            //const vector<UCandidateList>& D,
            //const FFStates& ucands_states,
            const SentenceMetadata& smeta,
            const ModelSet& models,
            const int sl
            /*bool is_goal*/);

  void InitStates(size_t state_size);

  ~UCandidate();

//  bool IsSelected() const;

//  bool HasSingleMissingLink() const;

  bool HasMissingLink() const;

  UCandidate* GetSourceUCand();

  int GetSourceNodeId();

  bool IsHeadIncomingState();

  FFState* GetHeadIncomingState();

  FFState* GetOutgoingState(int node_id);

  bool HasSource();

  bool IsGoal();

  bool CreateLink(UCandidate* ucand);

//  void InitializeUCandidate(
//                           const SentenceMetadata& smeta,
//                           //const vector<vector<UCandidate*> >& D,
//                           //const FFStates& ucands_states,
//                           const ModelSet& models//,
//                           /*const bool is_goal*/);
};

ostream& operator<<(ostream& os, const UCandidate& cand);

////////////////////////////////////////////

#endif
