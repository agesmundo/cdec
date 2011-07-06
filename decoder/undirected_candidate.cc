#include <iostream>

//added for LG
#include "undirected_candidate.h"

/////////////////////////////////////
//added here because model uses UCand as generalization of Edge
//therefore need to be imported in ff.cc
using namespace std;

  UCandidate::UCandidate(const Hypergraph::Edge& e,
            const SmallVectorInt& j,
            //const vector<UCandidateList>& D,
            const FFStates& node_states,
            const SentenceMetadata& smeta,
            const ModelSet& models,
            bool is_goal) :
      node_index_(-1),
      in_edge_(&e),
      j_(j) {
    InitializeUCandidate(smeta,/* D,*/ node_states, models, is_goal);
  }

  // used to query uniqueness
  UCandidate::UCandidate(const Hypergraph::Edge& e,
            const SmallVectorInt& j) : in_edge_(&e), j_(j) {}

  bool UCandidate::IsIncorporatedIntoHypergraph() const {
    return node_index_ >= 0;
  }

  void UCandidate::InitializeUCandidate(
                           const SentenceMetadata& smeta,
                           //const vector<vector<UCandidate*> >& D,
                           const FFStates& node_states,
                           const ModelSet& models,
                           const bool is_goal) {
    const Hypergraph::Edge& in_edge = *in_edge_;
    feature_values_ = in_edge.feature_values_;
    out_edge_.rule_ = in_edge.rule_;
    out_edge_.feature_values_ = in_edge.feature_values_;
    out_edge_.i_ = in_edge.i_;
    out_edge_.j_ = in_edge.j_;
    out_edge_.prev_i_ = in_edge.prev_i_;
    out_edge_.prev_j_ = in_edge.prev_j_;
    Hypergraph::TailNodeVector& tail = out_edge_.tail_nodes_;
    tail.resize(j_.size());

    //make a feature with score of links?
//    prob_t p = prob_t::One();
//    // cerr << "\nEstimating application of " << in_edge.rule_->AsString() << endl;
//    for (int i = 0; i < tail.size(); ++i) {
//      const UCandidate& ant = *D[in_edge.tail_nodes_[i]][j_[i]];
//      assert(ant.IsIncorporatedIntoHypergraph());
//      tail[i] = ant.node_index_;
//      p *= ant.vit_prob_;
//    }

    prob_t edge_estimate = prob_t::One();
    if (is_goal) {
      assert(tail.size() == 1);
      const FFState& ant_state = node_states[tail.front()];
      models.AddFinalFeatures(ant_state, &out_edge_, smeta);
    } else {
      models.AddFeaturesToUCandidate(smeta, node_states, this, &out_edge_, &state_, &edge_estimate);
    }

//    vit_prob_ = out_edge_.edge_prob_ * p;
//    est_prob_ = vit_prob_ * edge_estimate;

  }

ostream& operator<<(ostream& os, const UCandidate& cand) {
  os << "UCAND[";
  if (!cand.IsIncorporatedIntoHypergraph()) { os << "PENDING "; }
  else { os << "+LM_node=" << cand.node_index_; }
  //os << " edge=" << cand.in_edge_->id_; //printed by Edge<<
  os << " j=<";
  for (int i = 0; i < cand.j_.size(); ++i)
    os << (i==0 ? "" : " ") << cand.j_[i];
  os << "> ";
//  os << vit=" << log(cand.vit_prob_);
//  os << " est=" << log(cand.est_prob_);
  os << " act=" << log(cand.action_prob_);
  os << " in_edge_= " << *(cand.in_edge_)<< "; ";
  os << endl;
  os << "\tFEATS : " << cand.out_edge_.feature_values_;
  os << endl;
  os << "\tEST_F : " << cand.out_edge_.est_vals_;
  return os << ']';
}

////////////////////////////////////////////

