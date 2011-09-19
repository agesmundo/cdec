#ifndef _UNDIRECTED_CANDIDATE_H_
#define _UNDIRECTED_CANDIDATE_H_

#include <iostream>
#include <stack>
#include "ff.h"

struct ModelSet;
struct Hypergraph;
struct SentenceMetadata;

// Define the following macro if you want to see lots of debugging output
// when running Guided Undirected Greedy decoding
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
//typedef pair<int,FFState*>  Node2State;

struct UCandidate {
  //int ucand_index_;                     // -1 until popped from queue

  const Hypergraph::Edge* in_edge_;    // in -LM forest

  FFState** outgoing_states_;

  //  Node2State** outgoing_states_;         //in_node_id 2 state seen from that node //array max 2 elements (simple map)
//  int outgoing_states_size_; //no need, need a state for each link (also source to check if update)

//  FFState state_;

  //TODO? add pointer to in_edge feature for local features and avoid copy in costructor?
  FeatureVector est_vals_;
  FeatureVector feature_values_;

  //TODO? make this a pointer to avoid copy in costructor?
  //links to context
  //NB!!:
  //0 HEAD,
  //1 FIRST CHILD (left)
  //2 SECOND CHILD
  LinksVector context_links_; //NULL if no link, *UCand if link, -1 (goal_head_link_) for head link if Goal

  int source_link_; //context_links_(id) for the source UCand (-1 for starting leaf)

  //vit_prob_ and est_prob_ are not updated in LG training
  //prob_t vit_prob_;            // these are fixed until the cand
                               // is popped, then they may be updated
  //prob_t est_prob_;

  prob_t action_prob_;

  //needed to call update to features //TODO should be static
  const SentenceMetadata& smeta_;
  const ModelSet& models_;

  static UCandidate* goal_head_link_;

  UCandidate(const Hypergraph::Edge& e,
		    const LinksVector& lv,
            //const vector<UCandidateList>& D,
            //const FFStates& ucands_states,
            const SentenceMetadata& smeta,
            const ModelSet& models,
            const int sl
            /*bool is_goal*/);

  void InitStates(size_t state_size);

  void AllocStates();

  ~UCandidate();

  void DeleteStates(FFState** states);

//  bool IsSelected() const;

//  bool HasSingleMissingLink() const;

  int NLinks() const;

  void UpdateStates(stack<UCandidate*> &stck, vector < pair < UCandidate*, int > > &links_to_expand);

  bool HasMissingLink() const;

  UCandidate* GetSourceUCand();

  int GetSourceNodeId();

  int GetNodeIdFromLinkId(int link_id);

//  bool IsHeadIncomingState();
//
//  bool IsTailIncomingState(int tail_id);

  FFState* GetHeadIncomingState();

  FFState* GetTailIncomingState(int tail_id);

  FFState* GetHeadOutgoingState();

  FFState* GetTailOutgoingState(int node_id);

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
