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
            const int sl//,
            /*bool is_goal*/) :
      //ucand_index_(-1),
      in_edge_(&e),
      context_links_(context), //TODO why copy this?//TODO!!! use pointer to array like outgoing states//NLinks to variable!!!
      source_link_(sl),
	  smeta_(smeta),//TODO check is not making a copy!
	  models_(models){//TODO check is not making a copy!
	    feature_values_ = in_edge_->feature_values_;
	    AllocStates();
	    models_.AddFeaturesToUCandidate(smeta_, /*node_states,*/ this/*, &out_edge_, &state_, &edge_estimate*/);
  }

  UCandidate::~UCandidate(){
	  DeleteStates(outgoing_states_);
  }

  void UCandidate::DeleteStates(FFState** states){
	  for(int i=0;i<NLinks();i++){
		  delete states[i]; //pair
	  }
	  delete[] states; //array
  }


  void UCandidate::AllocStates(){
	  assert(NLinks()>0);//TODO on debug macro
	  outgoing_states_ = new FFState*[NLinks()];
	  for(int it=0;it<NLinks();it++){
		  outgoing_states_[it]=new FFState;
	  }
  }

  void UCandidate::InitStates(size_t state_size){

	  for(int i=0;i<NLinks();i++){
		  FFState* state = outgoing_states_[i];
#ifdef DEBUG_GU
		  assert(state!=NULL);
#endif
		  state->resize(state_size);
		  if (state_size > 0) {
		    memset(&(*state)[0], 0, state_size);
		  }
#ifdef DEBUG_GU
//		  cerr << " init outgoing_states_["<<i<<"] " << outgoing_states_[i]->first << " - " << outgoing_states_[i]->second << endl;
#endif
	  }
  }


  int UCandidate::NLinks() const {
//	  return in_edge_.Arity()+1; //children + head
	  return context_links_.size();
  }

  int UCandidate::ConnectedLinks() const{
	  int c=0;
	  for(int i=0; i<NLinks(); i++){
		  if(context_links_[i]!=NULL)c++;
	  }
	  return c;
  }

  void UCandidate::UpdateStates(stack<UCandidate*> &stck, vector < pair < UCandidate*, int > > &links_to_expand){
	  //store link to old states for comparison
	  FFState** old_outgoing_states=outgoing_states_;

	  //allocate new states
	  AllocStates();

	  //update states //TODO GU this also update scores that is kind of useless
	  models_.AddFeaturesToUCandidate(smeta_, this);

	  //compare old states to new
	  for(int i=0;i<NLinks();i++){
	  	if(i==0 && context_links_[i]==UCandidate::goal_head_link_)continue;
	  	if( !(*outgoing_states_[i] == *old_outgoing_states[i])){
	  		if(context_links_[i]){
	  			//TODO assert current source is not added //otherwise loop! anyway to ensure complexity must be one way propagation (that should never happen)
	  			stck.push(context_links_[i]); //TODO GU? avoid inserting twice items that already went to stack? avoid loops (that should never happen)
	  		}
	  		else{//in this case we need to signal an update in cands
	  			links_to_expand.push_back(pair < UCandidate*, int >(this, i));
#ifdef DEBUG_GU
	  			cerr <<"\tUCands queue update signaled: ("<< this<< " , "<< i << ")" <<endl;
#endif
	  		}
	  	}
	  }

	  //delete old states
	  DeleteStates(old_outgoing_states);

  }

//  bool UCandidate::HasSingleMissingLink() const{//TODO keep counter instead of computing?
//	  int count =0;
//	  for (int i=0; i<context_links_.size();i++){
//		  if(context_links_[i]==NULL){
//			  count++;
//			  if(count>1)return false;
//		  }
//	  }
//	  return count==1;
//  }

  bool UCandidate::HasMissingLink() const{
	  for (int i=0; i<context_links_.size();i++){
		  if(context_links_[i]==NULL){
			  return true;
		  }
	  }
	  return false;
  }

  UCandidate* UCandidate::GetSourceUCand(){
	  if(source_link_<0)return NULL;
	  return context_links_[source_link_];
  }

  int UCandidate::GetSourceNodeId(){
//	  if(source_link_<0)return -1;
	  if(source_link_==0)return in_edge_->head_node_;
	  if(source_link_==1)return in_edge_->tail_nodes_[0];
	  if(source_link_==2)return in_edge_->tail_nodes_[1];
//	  abort();
	  return -1;
  }

  int UCandidate::GetNodeIdFromLinkId(int link_id){
  		  if(link_id==0)return in_edge_->head_node_;
  		  if(link_id==1)return in_edge_->tail_nodes_[0];
  		  if(link_id==2)return in_edge_->tail_nodes_[1];
  		  return -1;
  }

  FFState* UCandidate::GetHeadOutgoingState(){
	  return outgoing_states_[0];
  }

  FFState* UCandidate::GetTailOutgoingState(int node_id){
	  for(int i=0;i<in_edge_->tail_nodes_.size();i++){
		  if (in_edge_->tail_nodes_[i]==node_id ) return outgoing_states_[i+1];
	  }
	  return NULL;
  }

  bool UCandidate::HasSource(){
	  return source_link_>=0;
  }

  //TODO GU store variable is_goal_ ?
  bool UCandidate::IsGoal(){
	  return context_links_[0]==goal_head_link_;
  }

//  bool UCandidate::IsHeadIncomingState(){
//	  if(IsGoal())return false;
//	  return context_links_[0]!=NULL;
//  }

//  bool UCandidate::IsTailIncomingState(int tail_id){
//#ifdef DEBUG_GU
//	  assert(tail_id+1 < context_links_.size());
//#endif
//	  return context_links_[tail_id+1]!=NULL;
//  }


  //returns NULL if head is not incoming state
  FFState* UCandidate::GetHeadIncomingState(){
#ifdef DEBUG_GU
	  //XXX	  if(!IsHeadIncomingState())return NULL; //commented for performance, should call this method after check
	  //XXX	  assert(IsHeadIncomingState());
	  assert(0 < context_links_.size());
#endif
	  if(IsGoal())return NULL;
	  if (context_links_[0]==NULL) return NULL;
	  return context_links_[0]->GetTailOutgoingState(in_edge_->head_node_);
  }

  //returns NULL if tail is not incoming state
  FFState* UCandidate::GetTailIncomingState(int tail_id){
#ifdef DEBUG_GU
	  //XXX	  if(!IsTailIncomingState(tail_id))return NULL; //commented for performance, should call this method after check
	  //XXX	  assert(IsTailIncomingState(tail_id));
	  assert(tail_id+1 < context_links_.size());
	  assert(tail_id < in_edge_->tail_nodes_.size());
	  if(context_links_[tail_id+1]!=NULL){
		  assert(in_edge_->tail_nodes_[tail_id]==context_links_[tail_id+1]->in_edge_->head_node_);
	  }
#endif
	  if (context_links_[tail_id+1]==NULL) return NULL;
	  return context_links_[tail_id+1]->GetHeadOutgoingState();
  }

  bool UCandidate::CreateLink(UCandidate* ucand){
	  int node_id=ucand->GetSourceNodeId();
#ifdef DEBUG_GU
	  assert(node_id>=0);
	  assert(context_links_.size()>0);
#endif
//	  if (ucand->source_link_==0)node_id =ucand->in_edge_->head_node_;
//	  else node_id = ucand->in_edge_->tail_nodes_[ucand->source_link_-1];

	  if(context_links_[0]==NULL && in_edge_->head_node_==node_id){
		  context_links_[0]=ucand;
		  return true;
	  }
	  if(context_links_.size()>1 && context_links_[1]==NULL && in_edge_->tail_nodes_[0]==node_id){
		  context_links_[1]=ucand;
		  return true;
	  }
	  if(context_links_.size()>2 && context_links_[2]==NULL && in_edge_->tail_nodes_[1]==node_id){
		  context_links_[2]=ucand;
		  return true;
	  }
	  return false;
  }

  //GU TODO is this needed??
  // used to query uniqueness
//  UCandidate::UCandidate(const Hypergraph::Edge& e,
//            const LinksVector& context) : in_edge_(&e), context_links_(context) {}

//  bool UCandidate::IsSelected() const {
//    return ucand_index_ >= 0;
//  }

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
//  if (!cand.IsSelected()) { os << "PENDING "; }
//  else { os << "+LM_node=" << cand.ucand_index_; }
  //os << " edge=" << cand.in_edge_->id_; //printed by Edge<<
  os << " context_links_=<";
  for (int i = 0; i < cand.context_links_.size(); ++i)
    os << (i==0 ? "" : " ") << cand.context_links_[i];
  os << ">";
//  os << " outgoing_states_=<";
//  for (int i = 0; i < cand.NLinks(); ++i)
//    os << (i==0 ? "" : " ") << cand.outgoing_states_[i];
//  os << ">";
  os<< " source_link_="<< cand.source_link_;
//  os << vit=" << log(cand.vit_prob_);
//  os << " est=" << log(cand.est_prob_);
  os << " act=" << log(cand.action_prob_);
  os << " in_edge_= " << *(cand.in_edge_)<< "; ";
  os << endl;
  os << "\tFEATS : " << cand.feature_values_;
//  os << endl;
//  os << "\tEST_F : " << cand.est_vals_;
  return os << ']';
}

UCandidate* UCandidate::goal_head_link_ =(UCandidate*)-1;

////////////////////////////////////////////

