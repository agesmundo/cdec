#include <iostream>

//added for LG
#include "undirected_candidate.h"

/////////////////////////////////////
//added here because model uses UCand as generalization of Edge
//therefore need to be imported in ff.cc
using namespace std;

  UCandidate::UCandidate(const Hypergraph::Edge& e,
            const LinksVector& context,
            //const vector<UCandidateList>& D,
            //const FFStates& node_states,
            const SentenceMetadata& smeta,
            const ModelSet& models,
            const int sl
            /*bool is_goal*/) :
      ucand_index_(-1),
      in_edge_(&e),
      context_links_(context),
      source_link_(sl){
	    const Hypergraph::Edge& in_edge = *in_edge_;
	    feature_values_ = in_edge.feature_values_;
	    models.AddFeaturesToUCandidate(smeta, /*node_states,*/ this,/* &out_edge_,*/ &state_/*, &edge_estimate*/);

  }

  //GU TODO is this needed??
  // used to query uniqueness
//  UCandidate::UCandidate(const Hypergraph::Edge& e,
//            const LinksVector& context) : in_edge_(&e), context_links_(context) {}

  bool UCandidate::IsSelected() const {
    return ucand_index_ >= 0;
  }

//  void UCandidate::InitializeUCandidate(
//                           const SentenceMetadata& smeta,
//                           //const vector<vector<UCandidate*> >& D,
//                           //const FFStates& node_states,
//                           const ModelSet& models//,
//                           //const bool is_goal
//                           ) {
//    const Hypergraph::Edge& in_edge = *in_edge_;
//    feature_values_ = in_edge.feature_values_;
////    out_edge_.rule_ = in_edge.rule_;
////    out_edge_.feature_values_ = in_edge.feature_values_;
////    out_edge_.i_ = in_edge.i_;
////    out_edge_.j_ = in_edge.j_;
////    out_edge_.prev_i_ = in_edge.prev_i_;
////    out_edge_.prev_j_ = in_edge.prev_j_;
////    Hypergraph::TailNodeVector& tail = out_edge_.tail_nodes_;
////    tail.resize(in_edge.tail_nodes_.size());
//
//    //make a feature with score of links?
////    prob_t p = prob_t::One();
////    // cerr << "\nEstimating application of " << in_edge.rule_->AsString() << endl;
////    for (int i = 0; i < tail.size(); ++i) {
////      const UCandidate& ant = *D[in_edge.tail_nodes_[i]][j_[i]];
////      assert(ant.IsIncorporatedIntoHypergraph());
////      tail[i] = ant.node_index_;
////      p *= ant.vit_prob_;
////    }
//
//    //prob_t edge_estimate = prob_t::One();
////    if (is_goal) {//TODO GU use diff feats extraction per root?
////      assert(tail.size() == 1);//TODO GU adapt NB!! no head link!! context_links_!!!
////      const FFState& ant_state = node_states[tail.front()];
////      models.AddFinalFeatures(ant_state, &out_edge_, smeta);
////    } else {
//      models.AddFeaturesToUCandidate(smeta, /*node_states,*/ this,/* &out_edge_,*/ &state_/*, &edge_estimate*/);
////    }
//
////    vit_prob_ = out_edge_.edge_prob_ * p;
////    est_prob_ = vit_prob_ * edge_estimate;
//
//  }

ostream& operator<<(ostream& os, const UCandidate& cand) {
  os << "UCAND";
  os << "(" << &cand << ")";
  os <<  "[";
  if (!cand.IsSelected()) { os << "PENDING "; }
  else { os << "+LM_node=" << cand.ucand_index_; }
  //os << " edge=" << cand.in_edge_->id_; //printed by Edge<<
  os << " context_links_=<";
  for (int i = 0; i < cand.context_links_.size(); ++i)
    os << (i==0 ? "" : " ") << cand.context_links_[i];
  os << "> ";
//  os << vit=" << log(cand.vit_prob_);
//  os << " est=" << log(cand.est_prob_);
  os << " act=" << log(cand.action_prob_);
  os << " in_edge_= " << *(cand.in_edge_)<< "; ";
  os << endl;
  os << "\tFEATS : " << cand.feature_values_;
  os << endl;
  os << "\tEST_F : " << cand.est_vals_;
  return os << ']';
}

////////////////////////////////////////////

