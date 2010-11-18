#include "apply_models.h"

#include <vector>
#include <algorithm>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

#include <boost/functional/hash.hpp>
#include <boost/shared_array.hpp>

#include "hg.h"
#include "ff.h"

#define NORMAL_CP 1
#define FAST_CP   2
#define FAST_CP_2 3

#define DUMMY 0

// Define the following macro if you want to see lots of debugging output
// when you run the GuidedPruning
//#define DEBUG_GP
//#undef DEBUG_GP

using namespace std;
using namespace std::tr1;

struct Candidate;
typedef SmallVector JVector;
typedef vector<Candidate*> CandidateHeap;
typedef vector<Candidate*> CandidateList;

// default vector size (* sizeof string is memory used)
static const size_t kRESERVE_NUM_NODES = 500000ul;

// life cycle: candidates are created, placed on the heap
// and retrieved by their estimated cost, when they're
// retrieved, they're incorporated into the +LM hypergraph
// where they also know the head node index they are
// attached to.  After they are added to the +LM hypergraph
// vit_prob_ and est_prob_ fields may be updated as better
// derivations are found (this happens since the successor's
// of derivation d may have a better score- they are
// explored lazily).  However, the updates don't happen
// when a candidate is in the heap so maintaining the heap
// property is not an issue.
struct Candidate {
	int node_index_;                     // -1 until incorporated into the +LM forest
	const Hypergraph::Edge* in_edge_;    // in -LM forest
	Hypergraph::Edge out_edge_;
	string state_;
	const JVector j_;//position in D[node][] for the correspondent node in the tail
	prob_t vit_prob_;            // these are fixed until the cand
	// is popped, then they may be updated
	prob_t est_prob_;

	Candidate(const Hypergraph::Edge& e,
			const JVector& j,
			const Hypergraph& out_hg,
			const vector<CandidateList>& D,
			const vector<string>& node_states,
			const SentenceMetadata& smeta,
			const ModelSet& models,
			bool is_goal) :
				node_index_(-1),
				in_edge_(&e),
				j_(j) {
				InitializeCandidate(out_hg, smeta, D, node_states, models, is_goal);
			}

			// used to query uniqueness
			Candidate(const Hypergraph::Edge& e,
					const JVector& j) : in_edge_(&e), j_(j) {}

			bool IsIncorporatedIntoHypergraph() const {
				return node_index_ >= 0;
			}

			void InitializeCandidate(const Hypergraph& out_hg,
					const SentenceMetadata& smeta,
					const vector<vector<Candidate*> >& D,
					const vector<string>& node_states,
					const ModelSet& models,
					const bool is_goal) {
				const Hypergraph::Edge& in_edge = *in_edge_;
				out_edge_.rule_ = in_edge.rule_;
				out_edge_.feature_values_ = in_edge.feature_values_;
				out_edge_.i_ = in_edge.i_;
				out_edge_.j_ = in_edge.j_;
				out_edge_.prev_i_ = in_edge.prev_i_;
				out_edge_.prev_j_ = in_edge.prev_j_;
				Hypergraph::TailNodeVector& tail = out_edge_.tail_nodes_;
				tail.resize(j_.size());
				prob_t p = prob_t::One();
				// cerr << "\nEstimating application of " << in_edge.rule_->AsString() << endl;
				for (int i = 0; i < tail.size(); ++i) {
					const Candidate& ant = *D[in_edge.tail_nodes_[i]][j_[i]];
					assert(ant.IsIncorporatedIntoHypergraph());
					tail[i] = ant.node_index_;
					p *= ant.vit_prob_;
				}
				prob_t edge_estimate = prob_t::One();
				if (is_goal) {
					assert(tail.size() == 1);
					const string& ant_state = node_states[tail.front()];
					models.AddFinalFeatures(ant_state, &out_edge_);
				} else {
					models.AddFeaturesToEdge(smeta, out_hg, node_states, &out_edge_, &state_, &edge_estimate);
				}
				vit_prob_ = out_edge_.edge_prob_ * p;
				est_prob_ = vit_prob_ * edge_estimate;
			}
};

ostream& operator<<(ostream& os, const Candidate& cand) {
	os << "CAND[";
	if (!cand.IsIncorporatedIntoHypergraph()) { os << "PENDING "; }
	else { os << "+LM_node=" << cand.node_index_; }
	os << " edge=" << cand.in_edge_->id_;
	os << " j=<";
	for (int i = 0; i < cand.j_.size(); ++i)
		os << (i==0 ? "" : " ") << cand.j_[i];
	os << "> vit=" << log(cand.vit_prob_);
	os << " est=" << log(cand.est_prob_);
	return os << ']';
}


struct GCandidate;
typedef vector<GCandidate*> GCandidateHeap;
typedef vector<GCandidate*> GCandidateList;

struct SharedArrayIterator{
	boost::shared_array<GCandidate*> sr_; //TODO ?? do we need to use this?
	int length_;
	int current_id_;

	//	SharedArrayIterator(const GCandidate* arr[], int length){
	//		sr_.reset(new GCandidate*[length]);
	//		length_=length;
	//		current_id_=0;
	//		memcpy(sr_.get(),arr,length);
	//	}

	//	SharedArrayIterator(const vector<GCandidate*>& toCopy){
	//		length_=toCopy.size();
	//		sr_.reset(new GCandidate*[length_]);
	//		current_id_=0;
	//		copy(toCopy.begin(),toCopy.end(),sr_.get());
	//	}

	SharedArrayIterator(const SharedArrayIterator& toCopy){
		sr_=toCopy.sr_; //point to same dynamic array
		length_=toCopy.length_;
		current_id_=toCopy.current_id_;
	}

	SharedArrayIterator(const	boost::shared_array<GCandidate*> sr, int length){
		sr_=sr; //point to same dynamic array
		length_=length;
		current_id_=0;
	}


	SharedArrayIterator(GCandidate* singleCand){
		sr_.reset(new GCandidate*[1]);
		sr_[0]=singleCand;
		length_=1;
		current_id_=0;
	}

	//avoid future expansion from this iterator
	void Freeze(){
		length_=current_id_+1;
	}

	inline bool HasMore(){
		return (current_id_+1 < length_);
	}

	GCandidate* GetCurrent(){
		return  sr_[current_id_];
	}

	bool Advance(){
		if(!HasMore())return false;
		current_id_++;
		return true;
	}
};


struct GCandidate {
	int node_index_;                     // -1 until incorporated into the +LM forest
	const Hypergraph::Edge* in_edge_;    // in -LM forest
	Hypergraph::Edge out_edge_; //TODO can avoid allocate this? allocate few fields?
	string state_;
	prob_t vit_prob_;          // these are fixed until the cand is popped, then they may be updated (not true in GP, cause not using dyn prog trick for cands)
	prob_t est_prob_;

	//always test befor access, if dummy =0 
	SharedArrayIterator* head_iterator_;
	SharedArrayIterator** tail_iterators_;

	GCandidate(const Hypergraph::Edge& e,
			const Hypergraph& out_hg,
			//const vector<string>& node_states,
			const SentenceMetadata& smeta,
			const ModelSet& models,
			const bool is_goal,
			SharedArrayIterator* headIterator,
			SharedArrayIterator** tailIterators) :
				node_index_(-1),
				in_edge_(&e),
				head_iterator_(headIterator),
				tail_iterators_(tailIterators)
				{
		InitializeGCandidate(out_hg, smeta, /*node_states,*/ models, is_goal);
				}

	// used to query uniqueness
	//GCandidate(const Hypergraph::Edge& e): in_edge_(&e) {} //TODO init iterators

	~GCandidate(){
		delete head_iterator_;
		if (tail_iterators_){
			for (int i=0; i<in_edge_->tail_nodes_.size();i++){
				delete tail_iterators_[i];
			}
			delete[] tail_iterators_;
		}
	}

		bool HasDummyTail()const{
		for(int i=0; i<TailSize();i++){
			if(tail_iterators_[i]==DUMMY){
				return true;	
			}
		}
		return false;
	}

	bool HasNotIncorporatedTail()const{
		for(int i=0; i<TailSize();i++){
			if(IsTailNotIncorporated(i)){
				return true;	
			}
		}
		return false;
	}

	bool IsTailNotIncorporated(int i)const{
		if(tail_iterators_[i]==DUMMY || tail_iterators_[i]->GetCurrent()->node_index_<0){
			return true;
		}
		return false;
	}

	GCandidate* GetCurrentTailIterator(int i) const{
		if(tail_iterators_[i]==DUMMY){
			return NULL;
		}
		return tail_iterators_[i]->GetCurrent();
	}

	GCandidate* GetCurrentHeadIterator() const {
		if(head_iterator_==DUMMY){
			return NULL;
		}
		return head_iterator_->GetCurrent();
	}

	bool IsIncorporatedIntoHypergraph() const {
		return node_index_ >= 0;
	}

	void InitializeGCandidate(const Hypergraph& out_hg,
			const SentenceMetadata& smeta,
			//const vector<string>& node_states,
			const ModelSet& models,
			const bool is_goal) {
		const Hypergraph::Edge& in_edge = *in_edge_;
		out_edge_.rule_ = in_edge.rule_;
		out_edge_.feature_values_ = in_edge.feature_values_;
		out_edge_.i_ = in_edge.i_;
		out_edge_.j_ = in_edge.j_;
		out_edge_.prev_i_ = in_edge.prev_i_;
		out_edge_.prev_j_ = in_edge.prev_j_;
		Hypergraph::TailNodeVector& tail = out_edge_.tail_nodes_;
		tail.resize(in_edge.tail_nodes_.size());
		prob_t p = prob_t::One();
		prob_t edge_estimate = prob_t::One();

		// compute vit prob from child
		for (int i = 0; i < tail.size(); ++i) {
			if (tail_iterators_[i]){//is not dummy
				const GCandidate* ant=tail_iterators_[i]->GetCurrent();
				tail[i] = ant->node_index_;
				p *= ant->vit_prob_;
			}
			else{//if no cand mark tail as absent
				tail[i]=-1;//NB if missing child, is it going to be incorp in +LM? it should not if we want same search space
				//TODO!! edge_estimate*= <the default prob for the missing link>
			}
		}

		//TODO!!
		//if father link is there add its vit_prob_ to edge_estimate
		//else edge_estimate *= <the default prob for the missing head link>

		if (is_goal) {
			assert(tail.size() == 1);
			assert(tail_iterators_[0]);
			const string& ant_state = tail_iterators_[0]->GetCurrent()->state_;
			models.AddFinalFeatures(ant_state, &out_edge_);
		} else {

			//collect antecedent states
			vector<string*> ant_states(tail.size());
			for (int i = 0; i < tail.size(); ++i) {
				if (tail_iterators_[i]){
					ant_states[i]=&tail_iterators_[i]->GetCurrent()->state_;
				}
			}

			models.AddFeaturesToEdgeGP(smeta, out_hg, ant_states, &out_edge_, &state_, &edge_estimate);
		}
		vit_prob_ = out_edge_.edge_prob_ * p;
		est_prob_ = vit_prob_ * edge_estimate;
	}

	int TailSize() const{
		return in_edge_->tail_nodes_.size();
	}

	bool HasMoreForTail(int i) const{
		return tail_iterators_[i] && tail_iterators_[i]->HasMore();
	}

	bool HasMoreHead() const{
		return head_iterator_ && head_iterator_->HasMore();
	}
};

ostream& operator<<(ostream& os, const SharedArrayIterator& sai) {
	os << "SharedArrayIterator=[";
	os << "id:"<<&sai << "|";
	//os <<sai.current_id_<< "|";
	//os <<sai.length_<< "|";
	for (int i=0; i<sai.length_ ; i++){
		if(i==sai.current_id_){
			os<< "(( ";
		}
		os << sai.sr_[i]->in_edge_->id_;
		if(i==sai.current_id_){
			os<< " ))";
		}
		os << "; ";
	}
	return os << ']';
}

ostream& operator<<(ostream& os, const GCandidate& cand) {
	os << "CAND " << &cand <<" =[";
	if (!cand.IsIncorporatedIntoHypergraph()) { os << "PENDING "<< "; "; }
	else { os << "+LM_node=" << cand.node_index_<< "; "; }
	os << " in_edge_= " << *(cand.in_edge_)<< "; ";
	os << " vit_prob_= " << log(cand.vit_prob_)<< "; ";
	os << " est_prob_= " << log(cand.est_prob_)<< "; ";

	os << "\n\t\thead_iterator_=";
	if(cand.head_iterator_){
		os<< *cand.head_iterator_<< ";";
	}
	else{
		os<< "DUMMY;";
	}

	//os << "\n\t\ttail_iterator_="<<cand.tail_iterators_;
	for(int i =0; i<cand.TailSize();i++){
		os << "\n\t\ttail_iterator_["<<i<<"]=";
		if(cand.tail_iterators_[i]){
			os<< *cand.tail_iterators_[i]<< ";";
		}
		else{
			os<< "DUMMY;";
		}
	}

	return os << "]";
}

struct HeapCandCompare {
	bool operator()(const Candidate* l, const Candidate* r) const {
		return l->est_prob_ < r->est_prob_;
	}
	bool operator()(const GCandidate* l, const GCandidate* r) const {
		return l->est_prob_ < r->est_prob_;
	}
};

struct EstProbSorter {
	bool operator()(const Candidate* l, const Candidate* r) const {
		return l->est_prob_ > r->est_prob_;
	}
	bool operator()(const GCandidate* l, const GCandidate* r) const {
		return l->est_prob_ > r->est_prob_;
	}
};

// the same candidate <edge, j> can be added multiple times if
// j is multidimensional (if you're going NW in Manhattan, you
// can first go north, then west, or you can go west then north)
// this is a hash function on the relevant variables from
// Candidate to enforce this.
struct CandidateUniquenessHash {
	size_t operator()(const Candidate* c) const {
		size_t x = 5381;
		x = ((x << 5) + x) ^ c->in_edge_->id_;
		for (int i = 0; i < c->j_.size(); ++i)
			x = ((x << 5) + x) ^ c->j_[i];
		return x;
	}

	size_t operator()(const GCandidate* c) const {
		size_t x = 5381;
		x = ((x << 5) + x) ^ c->in_edge_->id_;
		//    x = ((x << 5) + x) ^ c->in_edge_->head_node_;
		for (int i = 0; i < c->TailSize(); ++i)
			x = ((x << 5) + x) ^ (c->in_edge_->tail_nodes_[i]+1);//-1 dummy tail, what if result is 0? should put +2?
		return x;
	}
};

struct CandidateUniquenessEquals {
	bool operator()(const Candidate* a, const Candidate* b) const {
		return (a->in_edge_ == b->in_edge_) && (a->j_ == b->j_);
	}

	bool operator()(const GCandidate* a, const GCandidate* b) const {
		if (a->in_edge_ != b->in_edge_) return false;
		if (a->vit_prob_ != b->vit_prob_) return false;
		//if (a->GetCurrentHeadIterator() != b->GetCurrentHeadIterator())return false;
		if (a->TailSize()!=b->TailSize()) return false;
		for (int i =0; i<a->TailSize();i++){
			if (a->GetCurrentTailIterator(i)== NULL && b->GetCurrentTailIterator(i)!=NULL )return false;
			if (a->GetCurrentTailIterator(i)!= NULL && b->GetCurrentTailIterator(i)==NULL )return false;
			if (a->GetCurrentTailIterator(i)== NULL && b->GetCurrentTailIterator(i)==NULL )continue;
			if (a->GetCurrentTailIterator(i)->state_!=  b->GetCurrentTailIterator(i)->state_ )return false;
		}
		return true;
	}
};

typedef unordered_set<const Candidate*, CandidateUniquenessHash, CandidateUniquenessEquals> UniqueCandidateSet;
typedef unordered_map<string, Candidate*, boost::hash<string> > State2Node;

//GP version of State2Node
//in_hg Node_id & state-> GCand
typedef unordered_map<pair<int,string>, GCandidate*, boost::hash<pair<int,string> > >InNodeAndState2GCand;

typedef unordered_set<const GCandidate*, CandidateUniquenessHash, CandidateUniquenessEquals> UniqueGCandidateSet;

class CubePruningRescorer {

public:
	CubePruningRescorer(const ModelSet& m,
			const SentenceMetadata& sm,
			const Hypergraph& i,
			int pop_limit,
			Hypergraph* o,
			int s = NORMAL_CP ) :
				models(m),
				smeta(sm),
				in(i),
				out(*o),
				D(in.nodes_.size()),
				pop_limit_(pop_limit),
				strategy_(s)	{
		cerr << "  Applying feature functions (cube pruning, pop_limit = " << pop_limit_ << ')' << endl;
		node_states_.reserve(kRESERVE_NUM_NODES);
	}

	void Apply() {
		int num_nodes = in.nodes_.size();
		int goal_id = num_nodes - 1;
		int pregoal = goal_id - 1;
		int every = 1;
		if (num_nodes > 100) every = 10;
		assert(in.nodes_[pregoal].out_edges_.size() == 1);
		for (int i = 0; i < in.nodes_.size(); ++i) {
			if (i % every == 0) cerr << '.';
			if (strategy_==NORMAL_CP){
				KBest(i, i == goal_id);
			}
			if (strategy_==FAST_CP){
				KBestFast(i, i == goal_id);
			}
			if (strategy_==FAST_CP_2){
				KBestFast2(i, i == goal_id);
			}
		}
		cerr << endl;
		cerr << "  Best path: " << log(D[goal_id].front()->vit_prob_)
		<< "\t" << log(D[goal_id].front()->est_prob_) << endl;
		out.PruneUnreachable(D[goal_id].front()->node_index_);
		FreeAll();

		//		for (int i=0;i<out.edges_.size();i++){
		//			cerr << i << ": " << out.edges_[i] << endl;
		//		}
	}

private:
	void FreeAll() {
		for (int i = 0; i < D.size(); ++i) {
			CandidateList& D_i = D[i];
			for (int j = 0; j < D_i.size(); ++j)
				delete D_i[j];
		}
		D.clear();
	}

	void IncorporateIntoPlusLMForest(Candidate* item, State2Node* s2n, CandidateList* freelist) {
		Hypergraph::Edge* new_edge = out.AddEdge(item->out_edge_.rule_, item->out_edge_.tail_nodes_);
		new_edge->feature_values_ = item->out_edge_.feature_values_;
		new_edge->edge_prob_ = item->out_edge_.edge_prob_;
		new_edge->i_ = item->out_edge_.i_;
		new_edge->j_ = item->out_edge_.j_;
		new_edge->prev_i_ = item->out_edge_.prev_i_;
		new_edge->prev_j_ = item->out_edge_.prev_j_;
		Candidate*& o_item = (*s2n)[item->state_];
		if (!o_item){
			o_item = item;
		}

		int& node_id = o_item->node_index_;
		if (node_id < 0) {
			Hypergraph::Node* new_node = out.AddNode(in.nodes_[item->in_edge_->head_node_].cat_);
			node_states_.push_back(item->state_);
			node_id = new_node->id_;
		}
		Hypergraph::Node* node = &out.nodes_[node_id];
		out.ConnectEdgeToHeadNode(new_edge, node);

		// update candidate if we have a better derivation
		// note: the difference between the vit score and the estimated
		// score is the same for all items with a common residual DP
		// state
		if (item->vit_prob_ > o_item->vit_prob_) { 
			assert(o_item->state_ == item->state_);    // sanity check!
			o_item->est_prob_ = item->est_prob_;
			o_item->vit_prob_ = item->vit_prob_;
		}
		if (item != o_item) freelist->push_back(item);
	}

	void KBest(const int vert_index, const bool is_goal) {
		// cerr << "KBest(" << vert_index << ")\n";
		CandidateList& D_v = D[vert_index];
		assert(D_v.empty());
		const Hypergraph::Node& v = in.nodes_[vert_index];
		// cerr << "  has " << v.in_edges_.size() << " in-coming edges\n";
		const vector<int>& in_edges = v.in_edges_;
		CandidateHeap cand;
		CandidateList freelist;
		cand.reserve(in_edges.size());
		UniqueCandidateSet unique_cands;
		for (int i = 0; i < in_edges.size(); ++i) {
			const Hypergraph::Edge& edge = in.edges_[in_edges[i]];
			const JVector j(edge.tail_nodes_.size(), 0);
			cand.push_back(new Candidate(edge, j, out, D, node_states_, smeta, models, is_goal));
			assert(unique_cands.insert(cand.back()).second);  // these should all be unique!
		}
		//    cerr << "  making heap of " << cand.size() << " candidates\n";
		make_heap(cand.begin(), cand.end(), HeapCandCompare());
		State2Node state2node;   // "buf" in Figure 2
		int pops = 0;
		while(!cand.empty() && pops < pop_limit_) {
			pop_heap(cand.begin(), cand.end(), HeapCandCompare());
			Candidate* item = cand.back();
			cand.pop_back();
			// cerr << "POPPED: " << *item << endl;
			PushSucc(*item, is_goal, &cand, &unique_cands);
			IncorporateIntoPlusLMForest(item, &state2node, &freelist);
			++pops;
		}
		D_v.resize(state2node.size());
		int c = 0;
		for (State2Node::iterator i = state2node.begin(); i != state2node.end(); ++i)
			D_v[c++] = i->second;
		sort(D_v.begin(), D_v.end(), EstProbSorter());
		// cerr << "  expanded to " << D_v.size() << " nodes\n";

		for (int i = 0; i < cand.size(); ++i)
			delete cand[i];
		// freelist is necessary since even after an item merged, it still stays in
		// the unique set so it can't be deleted til now
		for (int i = 0; i < freelist.size(); ++i)
			delete freelist[i];
	}

	void KBestFast(const int vert_index, const bool is_goal) {
		// cerr << "KBest(" << vert_index << ")\n";
		CandidateList& D_v = D[vert_index];
		assert(D_v.empty());
		const Hypergraph::Node& v = in.nodes_[vert_index];
		// cerr << "  has " << v.in_edges_.size() << " in-coming edges\n";
		const vector<int>& in_edges = v.in_edges_;
		CandidateHeap cand;
		CandidateList freelist;
		cand.reserve(in_edges.size());
		//init with j<0,0> for all rules-edges that lead to node-(NT-span)
		for (int i = 0; i < in_edges.size(); ++i) {
			const Hypergraph::Edge& edge = in.edges_[in_edges[i]];
			const JVector j(edge.tail_nodes_.size(), 0);
			cand.push_back(new Candidate(edge, j, out, D, node_states_, smeta, models, is_goal));
		}
		//    cerr << "  making heap of " << cand.size() << " candidates\n";
		make_heap(cand.begin(), cand.end(), HeapCandCompare());
		State2Node state2node;   // "buf" in Figure 2
		int pops = 0;
		while(!cand.empty() && pops < pop_limit_) {
			pop_heap(cand.begin(), cand.end(), HeapCandCompare());
			Candidate* item = cand.back();
			cand.pop_back();
			//			cerr << "POPPED: " << *item << endl;

			PushSuccFast(*item, is_goal, &cand);
			IncorporateIntoPlusLMForest(item, &state2node, &freelist);
			++pops;
		}
		D_v.resize(state2node.size());
		int c = 0;
		for (State2Node::iterator i = state2node.begin(); i != state2node.end(); ++i){
			D_v[c++] = i->second;
			//			cerr << "MERGED: " << *i->second << endl;
		}
		//cerr <<"Node id: "<< vert_index<< endl;
		//#ifdef MEASURE_CA
		//		cerr << "countInProcess (pop/tot): node id: " << vert_index << " (" << count_in_process_pop << "/" << count_in_process_tot << ")"<<endl;
		//		cerr << "countAtEnd (pop/tot): node id: " << vert_index  << " (" << count_at_end_pop << "/" << count_at_end_tot << ")"<<endl;
		//#endif
		sort(D_v.begin(), D_v.end(), EstProbSorter());

		// cerr << "  expanded to " << D_v.size() << " nodes\n";

		for (int i = 0; i < cand.size(); ++i)
			delete cand[i];
		// freelist is necessary since even after an item merged, it still stays in
		// the unique set so it can't be deleted til now
		for (int i = 0; i < freelist.size(); ++i)
			delete freelist[i];
	}

	void KBestFast2(const int vert_index, const bool is_goal) {
		// cerr << "KBest(" << vert_index << ")\n";
		CandidateList& D_v = D[vert_index];
		assert(D_v.empty());
		const Hypergraph::Node& v = in.nodes_[vert_index];
		// cerr << "  has " << v.in_edges_.size() << " in-coming edges\n";
		const vector<int>& in_edges = v.in_edges_;
		CandidateHeap cand;
		CandidateList freelist;
		cand.reserve(in_edges.size());
		UniqueCandidateSet unique_accepted;
		//init with j<0,0> for all rules-edges that lead to node-(NT-span)
		for (int i = 0; i < in_edges.size(); ++i) {
			const Hypergraph::Edge& edge = in.edges_[in_edges[i]];
			const JVector j(edge.tail_nodes_.size(), 0);
			cand.push_back(new Candidate(edge, j, out, D, node_states_, smeta, models, is_goal));
		}
		//    cerr << "  making heap of " << cand.size() << " candidates\n";
		make_heap(cand.begin(), cand.end(), HeapCandCompare());
		State2Node state2node;   // "buf" in Figure 2
		int pops = 0;
		while(!cand.empty() && pops < pop_limit_) {
			pop_heap(cand.begin(), cand.end(), HeapCandCompare());
			Candidate* item = cand.back();
			cand.pop_back();
			assert(unique_accepted.insert(item).second);  // these should all be unique!
			//			cerr << "POPPED: " << *item << endl;

			PushSuccFast2(*item, is_goal, &cand, &unique_accepted);
			IncorporateIntoPlusLMForest(item, &state2node, &freelist);
			++pops;
		}
		D_v.resize(state2node.size());
		int c = 0;
		for (State2Node::iterator i = state2node.begin(); i != state2node.end(); ++i){
			D_v[c++] = i->second;
			//			cerr << "MERGED: " << *i->second << endl;
		}
		//cerr <<"Node id: "<< vert_index<< endl;
		//#ifdef MEASURE_CA
		//		cerr << "countInProcess (pop/tot): node id: " << vert_index << " (" << count_in_process_pop << "/" << count_in_process_tot << ")"<<endl;
		//		cerr << "countAtEnd (pop/tot): node id: " << vert_index  << " (" << count_at_end_pop << "/" << count_at_end_tot << ")"<<endl;
		//#endif
		sort(D_v.begin(), D_v.end(), EstProbSorter());

		// cerr << "  expanded to " << D_v.size() << " nodes\n";

		for (int i = 0; i < cand.size(); ++i)
			delete cand[i];
		// freelist is necessary since even after an item merged, it still stays in
		// the unique set so it can't be deleted til now
		for (int i = 0; i < freelist.size(); ++i)
			delete freelist[i];
	}

	void PushSucc(const Candidate& item, const bool is_goal, CandidateHeap* pcand, UniqueCandidateSet* cs) {
		CandidateHeap& cand = *pcand;
		for (int i = 0; i < item.j_.size(); ++i) {
			JVector j = item.j_;
			++j[i];
			if (j[i] < D[item.in_edge_->tail_nodes_[i]].size()) {
				Candidate query_unique(*item.in_edge_, j);
				if (cs->count(&query_unique) == 0) {
					Candidate* new_cand = new Candidate(*item.in_edge_, j, out, D, node_states_, smeta, models, is_goal);
					cand.push_back(new_cand);
					push_heap(cand.begin(), cand.end(), HeapCandCompare());
					assert(cs->insert(new_cand).second);  // insert into uniqueness set, sanity check
				}
			}
		}
	}

	//PushSucc following unique ancestor generation function
	void PushSuccFast(const Candidate& item, const bool is_goal, CandidateHeap* pcand){
		CandidateHeap& cand = *pcand;
		for (int i = 0; i < item.j_.size(); ++i) {
			JVector j = item.j_;
			++j[i];
			if (j[i] < D[item.in_edge_->tail_nodes_[i]].size()) {
				Candidate* new_cand = new Candidate(*item.in_edge_, j, out, D, node_states_, smeta, models, is_goal);
				cand.push_back(new_cand);
				push_heap(cand.begin(), cand.end(), HeapCandCompare());
			}
			if(item.j_[i]!=0){
				return;
			}
		}
	}

	//PushSucc only if all ancest Cand are added
	void PushSuccFast2(const Candidate& item, const bool is_goal, CandidateHeap* pcand, UniqueCandidateSet* ps){
		CandidateHeap& cand = *pcand;
		for (int i = 0; i < item.j_.size(); ++i) {
			JVector j = item.j_;
			++j[i];
			if (j[i] < D[item.in_edge_->tail_nodes_[i]].size()) {
				Candidate query_unique(*item.in_edge_, j);
				if (HasAllAncestors(&query_unique,ps)) {
					Candidate* new_cand = new Candidate(*item.in_edge_, j, out, D, node_states_, smeta, models, is_goal);
					cand.push_back(new_cand);
					push_heap(cand.begin(), cand.end(), HeapCandCompare());
				}
			}
		}
	}

	bool HasAllAncestors(const Candidate* item, UniqueCandidateSet* cs){
		for (int i = 0; i < item->j_.size(); ++i) {
			JVector j = item->j_;
			--j[i];
			if (j[i] >=0) {
				Candidate query_unique(*item->in_edge_, j);
				if (cs->count(&query_unique) == 0) {
					return false;
				}
			}
		}
		return true;
	}

	const ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	Hypergraph& out;

	vector<CandidateList> D;   // maps nodes in in-HG to the
	// equivalent nodes (many due to state
	// splits) in the out-HG.
	vector<string> node_states_;  // for each node in the out-HG what is
	// its q function value?
	const int pop_limit_;
	const int strategy_;       //switch Cube Pruning strategy: 1 normal, 2 fast (alg 2), 3 fast_2 (alg 3). (see: Gesmundo A., Henderson J,. Faster Cube Pruning, IWSLT 2010)
};

struct GCandidateSmartList{

private:
	GCandidateList list_;
	bool needNewSA_; //true if shared array need to be updated
	boost::shared_array<GCandidate*> sr_;

public:
	static int pop_limit_;

	GCandidateSmartList():needNewSA_(true){};

	void push_back(GCandidate* item){
		list_.push_back(item);
		needNewSA_=true;
	}

	unsigned int size(){
		return list_.size();
	}

	void SetTailSA(){
		if(needNewSA_){
			GCandidateList new_list;
			unordered_map<string, GCandidate*, boost::hash<string> > s2g;

#ifdef DEBUG_GP
			cerr << " SetTailSA(): "<<endl << "init list size : " << list_.size() <<endl ;
#endif

			sort(list_.begin(),list_.end(),EstProbSorter());
			for(int i=0;i<list_.size() && new_list.size()<GCandidateSmartList::pop_limit_;i++){
				GCandidate* gc=list_[i];
				GCandidate*& o_gc = s2g[gc->state_];
				if (!o_gc){
					o_gc = gc;
					new_list.push_back(gc);
#ifdef DEBUG_GP
					cerr<< " ADDED: \n" << *gc << endl ;
					assert(!gc->head_iterator_);
				}
				else{
					cerr<< " MERGED: \ndel:" << *gc << "\nkeep:" << *o_gc <<endl;
					assert(!gc->head_iterator_);
#endif
				}
			}
			list_.swap(new_list);
			sr_.reset(new GCandidate*[list_.size()]);
			copy(list_.begin(),list_.end(),sr_.get());
			needNewSA_=false;

#ifdef DEBUG_GP
			cerr << "end SetTailSA(), list size : " << list_.size() <<endl ;
#endif

		}
	}

	//NB this works for binary rules
	void SetHeadSA(){
		if(needNewSA_){
			GCandidateList new_list;
			unordered_map<pair<int,string>, GCandidate*, boost::hash<pair<int,string> > > es2g;

#ifdef DEBUG_GP
			cerr << " SetHeadSA(): "<<endl << "init list size : " << list_.size() <<endl ;
#endif

			sort(list_.begin(),list_.end(),EstProbSorter());
			for(int i=0; i<list_.size() && new_list.size()<GCandidateSmartList::pop_limit_;i++){
				GCandidate* gc=list_[i];

				string sibling_state;

				assert(gc->in_edge_->tail_nodes_.size()==2);
				assert((gc->tail_iterators_[0]==DUMMY && gc->tail_iterators_[1]!=DUMMY) ||(gc->tail_iterators_[0]!=DUMMY && gc->tail_iterators_[1]==DUMMY) );

				if(gc->tail_iterators_[0]!=DUMMY){
					sibling_state=gc->tail_iterators_[0]->GetCurrent()->state_;
				}
				else{
					sibling_state=gc->tail_iterators_[1]->GetCurrent()->state_;
				}

				pair<int, string> query(gc->in_edge_->id_,sibling_state);
				GCandidate*& o_gc = es2g[pair<int, string>(query)];
				if (!o_gc){
					o_gc = gc;
					new_list.push_back(gc);
#ifdef DEBUG_GP
					cerr<< " ADDED: \n" << *gc << endl ;
					assert(!gc->head_iterator_);
				}
				else{
					cerr<< " MERGED: \ndel:" << *gc << "\nkeep:" << *o_gc <<endl;
					assert(!gc->head_iterator_);
#endif
				}
			}
			list_.swap(new_list);
			sr_.reset(new GCandidate*[list_.size()]);
			copy(list_.begin(),list_.end(),sr_.get());
			needNewSA_=false;

#ifdef DEBUG_GP
			cerr << "end SetHeadSA(), list size : " << list_.size() <<endl ;
#endif

		}
	}

	//Warning, this allocate dyn memory!
	SharedArrayIterator* GetHeadIterator(){
		SetHeadSA();
		assert(list_.size()!=0);
		return new SharedArrayIterator(sr_,list_.size());
	}

	//Warning, this allocate dyn memory!
	SharedArrayIterator* GetTailIterator(){
		SetTailSA();
		assert(list_.size()!=0);
		return new SharedArrayIterator(sr_,list_.size());
	}

	//	size_t size() const {
	//		return list_.size();
	//	}

	//	const GCandidate* operator[](size_t id){
	//		if(IsEmpty())return NULL;
	//		SetSA();
	//		return list_[id];
	//	}

	inline bool IsEmpty(){
		return list_.size()==0; 
	}

	bool IsNewBest(GCandidate* cand){
		if(IsEmpty()) return true;
		sort(list_.begin(),list_.end(),EstProbSorter());
		//XXX cerr << "======" << list_[0];
		return list_[0]->est_prob_ < cand->est_prob_;
	}

};
int GCandidateSmartList::pop_limit_;


class GuidedPruningRescorer {

public:
	GuidedPruningRescorer(const ModelSet& m,
			const SentenceMetadata& sm,
			const Hypergraph& i,
			int pop_limit,
			Hypergraph* o) :
				models(m),
				smeta(sm),
				in(i),
				goal_id_(in.nodes_.size()-1),
				out_goal_id_(-1),
				out(*o),
				H(in.nodes_.size()),
				D(in.nodes_.size()),
				Hc(in.nodes_.size()),
				Dc(in.nodes_.size()),
				pop_limit_(pop_limit),
				free_(in.nodes_.size()*pop_limit_) { //TODO speed+-
		GCandidateSmartList::pop_limit_=pop_limit;
		cerr << "  Applying feature functions (guided pruning, pop_limit = " << pop_limit_ << ')' << endl;
		//node_states_.reserve(kRESERVE_NUM_NODES);
	}

	void Apply() {
		int num_nodes = in.nodes_.size();
		int goal_id = num_nodes - 1;
		int pregoal = goal_id - 1;
		assert(in.nodes_[pregoal].out_edges_.size() == 1);
		GCandidateHeap cands; //contains cands
		//UniqueGCandidateSet unique_cands; //to check that cadidate is unique at insertion in cands TODO shouldn't be needed!!!we use trick of alg2
		InNodeAndState2GCand state2node;//to apply dyn. prog. trik to +LM

		InitCands(cands/*, unique_cands*/);

		for (;!cands.empty();) {
#ifdef DEBUG_GP
			cerr<< "/////////////////////////////////////////////////////////////////////\n";
			cerr << " L" << " | cands.size(): "<< cands.size()<<"\n";
#endif
			GCandidate* aCand = PopBest(cands);

			if(!GCandPoppable(*aCand) && !IsNewBest(aCand)){//TODO comment out not poppable here (and in add GCand) and compare number of nodes that arrived at the beam size, make sure has same performances
#ifdef DEBUG_GP
				cerr << "NotPoppable and NotNewBest" << endl;
#endif
				continue;
			}

			IncorporateIntoPlusLMForest(aCand,&state2node);

			PushSucc(*aCand, cands/*, unique_cands*/);

			HeadPropagation(*aCand, cands/*, unique_cands*/);

			//TailPropagation(*aCand, cands/*, unique_cands*/);

			UpdateSLists(aCand);

		}

		//		cerr << "Best path: " <<log(D[goal_id][0]->vit_prob_)<<endl;

		//free memory used by cands
		FreeAll(cands);

#ifdef DEBUG_GP
		assert(out_goal_id_>=0);
		for (int i=0;i<out.edges_.size();i++){
			cerr << i << ": " << out.edges_[i] << endl;
		}
#endif

		out.TopologicallySortNodesAndEdges(out_goal_id_);
		
		//Not needed if Topo.Sort is called already 
		//out.PruneUnreachable(out_goal_id_);
	}

private:

	void UpdateSLists(GCandidate* aCand){
		if(!aCand->head_iterator_ && !aCand->HasDummyTail()){
			D[aCand->in_edge_->head_node_].push_back(aCand);
#ifdef DEBUG_GP
			cerr << "added in D[" <<aCand->in_edge_->head_node_ <<"] | size = " << D[aCand->in_edge_->head_node_].size() << endl;
#endif
		}
		if(!aCand->head_iterator_){
			for(int i =0; i <aCand->TailSize();i++){
				if(!aCand->tail_iterators_[i]){
					H[aCand->in_edge_->tail_nodes_[i]].push_back(aCand);
#ifdef DEBUG_GP
					cerr << "added in H[" << aCand->in_edge_->tail_nodes_[i] <<"] | size = " << H[aCand->in_edge_->tail_nodes_[i]].size() << endl;
#endif
				}
			}
		}
	}
	
	void UpdateSListsC(GCandidate* aCand){
		if(!aCand->head_iterator_ && !aCand->HasDummyTail()){
			Dc[aCand->in_edge_->head_node_].push_back(aCand);
#ifdef DEBUG_GP
			cerr << "added in Dc[" <<aCand->in_edge_->head_node_ <<"] | size = " << Dc[aCand->in_edge_->head_node_].size() << endl;
#endif
		}
		if(!aCand->head_iterator_){
			for(int i =0; i <aCand->TailSize();i++){
				if(!aCand->tail_iterators_[i]){
					Hc[aCand->in_edge_->tail_nodes_[i]].push_back(aCand);
#ifdef DEBUG_GP
					cerr << "added in Hc[" << aCand->in_edge_->tail_nodes_[i] <<"] | size = " << Hc[aCand->in_edge_->tail_nodes_[i]].size() << endl;
#endif
				}
			}
		}
	}

	inline bool IsGoal (const Hypergraph::Edge& edge){
		return (edge.head_node_==goal_id_);
	}

	inline bool IsGoal (const GCandidate& cand){
		return IsGoal(*cand.in_edge_);
	}

	size_t TailSize(const Hypergraph::Edge& edge){
		return edge.tail_nodes_.size();
	}

	size_t TailSize(const GCandidate& cand){
		return cand.in_edge_->tail_nodes_.size();
	}

	//initialize candidate heap with leafs
	void InitCands(GCandidateHeap& cands/*, UniqueGCandidateSet& unique_cands*/)
	{
#ifdef DEBUG_GP
		cerr << "InintCands(): " << "\n"; 
#endif
		for (int i = 0; i < in.edges_.size(); ++i) {//loop edges
			const Hypergraph::Edge& currentEdge= in.edges_.at(i);
			if(currentEdge.tail_nodes_.size()==0){//leafs
				GCandidate* new_cand = CreateCandidate(currentEdge, DUMMY, DUMMY);
				AddGCandidate(new_cand,cands/*,unique_cands*/);
			}
		}

#ifdef DEBUG_GP
		cerr << "cands.size(): "<< cands.size()<<"\n";
#endif

	}

	//order cands and pop best
	inline GCandidate* PopBest(GCandidateHeap& cands){
		make_heap(cands.begin(), cands.end(), HeapCandCompare());
		pop_heap(cands.begin(), cands.end(), HeapCandCompare());
		GCandidate* aCand = cands.back(); //accepted cand
		cands.pop_back();
#ifdef DEBUG_GP
		cerr << "free_[" << free_.size() << "] = " << aCand <<endl;
		cerr << "PopBest(): " << *aCand << "\n";
#endif
		free_.push_back(aCand);
		return aCand;
	}

	void FreeAll(GCandidateHeap& cands) {
		for (int i = 0; i < cands.size(); ++i){
			delete cands[i];
		}
		for (int i = 0; i < free_.size(); ++i){
			delete free_[i];
		}
	}

	//	inline bool NodePoppable(int node_id){
	//		return D[node_id].size()>=pop_limit_;
	//	}

	inline bool GCandPoppable(const GCandidate& cand){
		if(D[cand.in_edge_->head_node_].size()>=pop_limit_) return false;
		for(int i=0; i<cand.TailSize();i++){
			if(cand.tail_iterators_[i]==DUMMY && H[cand.in_edge_->tail_nodes_[i]].size()>=pop_limit_){
				return false;
			}
		}
		return true;
	}
	
	inline bool GCandPoppableC(const GCandidate& cand){
		if(Dc[cand.in_edge_->head_node_].size()>=pop_limit_) return false;
		for(int i=0; i<cand.TailSize();i++){
			if(cand.tail_iterators_[i]==DUMMY && Hc[cand.in_edge_->tail_nodes_[i]].size()>=pop_limit_){
				return false;
			}
		}
		return true;
	}

	bool AlreadyIncorporated(int out_node_id,Hypergraph::Edge& new_edge){

		Hypergraph::Node& n = out.nodes_[out_node_id];
		const vector<int>& in_edges = n.in_edges_;
		for (int i = 0; i < in_edges.size(); ++i) {
			const Hypergraph::Edge& edge = out.edges_[in_edges[i]];
			if(
					new_edge.rule_ == edge.rule_ &&
					new_edge.tail_nodes_ == edge.tail_nodes_ 
			){
#ifdef DEBUG_GP
				assert(new_edge.edge_prob_ == edge.edge_prob_);
				cerr << "Not Added Matched with: " << edge<< endl;
#endif
				return true;
			}
		}
		return false;
	}

	void IncorporateIntoPlusLMForest(GCandidate* item, InNodeAndState2GCand* state2node) {

#ifdef DEBUG_GP
		cerr << "Incorporate(): \n"; 
#endif
		//do not incorporate if missing any tail
		if(item->HasNotIncorporatedTail() /*|| !NodePoppable(*item)*/){ //Could check: block notPoppable if not new Best 
#ifdef DEBUG_GP
			cerr<< "Not Inc. hasDummyTail "<</*"or NodeNotPoppable" <<*/endl;
#endif
			return;
		}

		pair<int, string> query(item->in_edge_->head_node_,item->state_);
		GCandidate*& o_item = (*state2node)[query];
		if (!o_item){
#ifdef DEBUG_GP
			assert( item->node_index_ <0);
#endif
			Hypergraph::Node* new_node = out.AddNode(in.nodes_[item->in_edge_->head_node_].cat_);
			o_item = item;
			o_item->node_index_= new_node->id_;
		}
		else {
#ifdef DEBUG_GP
			assert( o_item->node_index_ >=0);
#endif
			item->node_index_ = o_item->node_index_;
			if(AlreadyIncorporated(o_item->node_index_,item->out_edge_)){
#ifdef DEBUG_GP
				cerr<< "AlreadyIncorporated()" <<endl;
#endif
				return;
			}
		}

		//create new edge
		Hypergraph::Edge* new_edge = out.AddEdge(item->out_edge_.rule_, item->out_edge_.tail_nodes_);
		new_edge->feature_values_ = item->out_edge_.feature_values_;
		new_edge->edge_prob_ = item->out_edge_.edge_prob_;
		new_edge->i_ = item->out_edge_.i_;
		new_edge->j_ = item->out_edge_.j_;
		new_edge->prev_i_ = item->out_edge_.prev_i_;
		new_edge->prev_j_ = item->out_edge_.prev_j_;

		Hypergraph::Node* node = &out.nodes_[o_item->node_index_];
		out.ConnectEdgeToHeadNode(new_edge, node);

#ifdef DEBUG_GP
		cerr<<"Added edge in +LM: " << *new_edge<<endl;
		//cerr<<"Updated cand     : " << *item;
#endif

		if(IsGoal(*item)){
#ifdef DEBUG_GP
			assert(out_goal_id_==-1 || out_goal_id_==node->id_);	
#endif
			out_goal_id_=node->id_;
		}

		//		if(node_pop_limit_[item->in_edge_->head_node_]>pop_limit_){
		//			cerr << "Extra for node" 
		//			return false;
		//		}

		//TODO NB the part below is missing in GP should be done when updating D and H
		//... it's the dyn prog trick at cands level

		// update candidate if we have a better derivation
		// note: the difference between the vit score and the estimated
		// score is the same for all items with a common residual DP
		// state
		/*		if (item->vit_prob_ > o_item->vit_prob_) { 
			assert(o_item->state_ == item->state_);    // sanity check!
			o_item->est_prob_ = item->est_prob_;
			o_item->vit_prob_ = item->vit_prob_;
		}
		if (item != o_item) freelist->push_back(item);*/

	}

	inline GCandidate* CreateCandidate(const Hypergraph::Edge& edge,SharedArrayIterator* headIterator,SharedArrayIterator** tailIterators){
		return new GCandidate(edge,out,/*node_states_,*/ smeta, models,IsGoal(edge),headIterator,tailIterators);
	}

	bool IsNewBest(GCandidate* cand){
		for(int i=0; i<cand->TailSize();i++){
			if(cand->tail_iterators_[i]==DUMMY){
				return H[cand->in_edge_->tail_nodes_[i]].IsNewBest(cand);
			}
		}
		return D[cand->in_edge_->head_node_].IsNewBest(cand);
	}
	
	bool IsNewBestC(GCandidate* cand){
		for(int i=0; i<cand->TailSize();i++){
			if(cand->tail_iterators_[i]==DUMMY){
				return Hc[cand->in_edge_->tail_nodes_[i]].IsNewBest(cand);
			}
		}
		return Dc[cand->in_edge_->head_node_].IsNewBest(cand);
	}

	inline void AddGCandidate( GCandidate* cand,GCandidateHeap& cands/*, UniqueGCandidateSet& unique_cands*/){
//		if(!unique_cands.insert(cand).second){
//#ifdef DEBUG_GP
//			cerr << "Not unique GCand: "<< *cand <<endl;
//#endif
//			free_.push_back(cand);
//			return;
//		}
		if(!GCandPoppableC(*cand) && !IsNewBestC(cand)){
#ifdef DEBUG_GP
			cerr << " NotPop and NotNewBest"<< endl;
#endif
			free_.push_back(cand); //TODO delete here and avoid copy of delete in caller 
			return;
		}
		//		if(!NodePoppable(*cand)){
		//			cerr << "Not Poppable but added since new best"<< endl;;
		//		}

		cands.push_back(cand);
		push_heap(cands.begin(), cands.end(), HeapCandCompare());
		
		UpdateSListsC(cand);
		
		//assert(unique_cands.insert(cand).second);  // insert into uniqueness set, sanity check
#ifdef DEBUG_GP
		cerr << "\tAddCandidate(): " << *cand << "\n"; 
		//assert(NodePoppable(*cand)); //TODO for complex cases assertion may be false, but better keep it untill tested 
#endif
	}

	inline SharedArrayIterator** CopyTailPTArray(SharedArrayIterator** toCopy, int size, int dummyId=-1){
		SharedArrayIterator** newTailIterators = new SharedArrayIterator*[size];//or (SharedArrayIterator**) malloc (tailSize * sizeof(SharedArrayIterator*))
		for(int i=0; i< size ; i++){
			if(toCopy[i] == NULL || i==dummyId){
				newTailIterators[i]=NULL;
			}
			else{
				newTailIterators[i]=new SharedArrayIterator(*toCopy[i]);
			}
		}
		return newTailIterators;
	}

	inline SharedArrayIterator* CopyHeadPT(SharedArrayIterator* toCopy){
		if (toCopy==NULL){
			return NULL;
		}
		return new SharedArrayIterator(*toCopy);
	}

	void PushSucc(const GCandidate& aCand, GCandidateHeap& cands/*, UniqueGCandidateSet& unique_cands*/) {
#ifdef DEBUG_GP
		cerr << "PushSucc(): \n"; 
#endif
		//		if(!NodePoppable(aCand)){
		//			return;
		//		}
		int tailSize =aCand.TailSize();
		SharedArrayIterator** newTailIterators;
		for (int i = 0; i < tailSize; ++i) {//tail nodes
			if(aCand.HasMoreForTail(i)){
#ifdef DEBUG_GP
				cerr << "HasMoreForTail " <<i<< endl;
#endif
				newTailIterators = CopyTailPTArray(aCand.tail_iterators_,tailSize);
				//???newTailIterators[i]=new SharedArrayIterator(*aCand.tail_iterators_[i]);//copy (NB inside there is SP)
				assert(newTailIterators[i]->Advance());//and advance
				//TODO check compatibility of child's head if has any... while not compatible advance... IF we add all in D and not just dummy
				GCandidate* new_cand = CreateCandidate(*aCand.in_edge_,CopyHeadPT(aCand.head_iterator_),newTailIterators);
				AddGCandidate(new_cand,cands/*,unique_cands*/);
			}
		}
		if(aCand.HasMoreHead()){//head
#ifdef DEBUG_GP
			cerr << "HasMoreHead" << endl;
#endif
			SharedArrayIterator* newHeadIterator=new SharedArrayIterator(*aCand.head_iterator_);
			assert(newHeadIterator->Advance());
			GCandidate* new_cand = CreateCandidate(*aCand.in_edge_,newHeadIterator,CopyTailPTArray(aCand.tail_iterators_,tailSize));
			AddGCandidate(new_cand,cands/*,unique_cands*/);
		}
	}

	void TailPropagation(GCandidate& aCand, GCandidateHeap& cands/*, UniqueGCandidateSet& unique_cands*/){
#ifdef DEBUG_GP
		cerr << "TailPropagation(): \n"; 
#endif		
		for(int i=0;i<aCand.TailSize();i++){
			if(!aCand.tail_iterators_[i]){//TODO double check, do not propagate when have tail, this is not simetric with head
				int currentTailNodeID = aCand.in_edge_->tail_nodes_[i];
				//				if(!NodePoppable(currentTailNodeID)){
				//					continue;
				//				}
				const Hypergraph::Node& currentTailNode = in.nodes_[currentTailNodeID];
				for(int j=0;j<currentTailNode.in_edges_.size();j++){
					const Hypergraph::Edge& currentTailEdge = in.edges_[currentTailNode.in_edges_[j]];
					int tailSize = TailSize(currentTailEdge);
					SharedArrayIterator** newTailIterators;
					if(tailSize==0){
						newTailIterators=DUMMY;
					}
					else{
						newTailIterators = new SharedArrayIterator*[tailSize];
						for(int k=0;k<tailSize;k++){
							if(D[currentTailNodeID].IsEmpty()){
								newTailIterators[k]=NULL;
							}
							else{
								newTailIterators[k]=D[currentTailNodeID].GetTailIterator();//check head compatibility if we add also not dummy head in D
#ifdef DEBUG_GP
								cerr<< "Get D[" << currentTailNodeID << "]" << *newTailIterators[k] << endl;
#endif
							}
						}
					}
					SharedArrayIterator* newHeadIterator = new SharedArrayIterator(&aCand);
					GCandidate* newTailCand = CreateCandidate(currentTailEdge ,newHeadIterator ,newTailIterators);
					AddGCandidate(newTailCand,cands/*,unique_cands*/);
					//add cands with one dummy tails??
				}
			}
		}
	}

	void HeadPropagation(GCandidate& aCand, GCandidateHeap& cands/*, UniqueGCandidateSet& unique_cands*/){

#ifdef DEBUG_GP
		cerr << "HeadPropagation(): \n"; 
#endif

		//do not head propagate if notInc tail
		//TODO IMP this is first version try variation for improvements
		if(aCand.HasNotIncorporatedTail()){
			return;
		}

		assert (!aCand.head_iterator_);
		//		if(!aCand.head_iterator_){
		const Hypergraph::Node& headNode=in.nodes_[aCand.in_edge_->head_node_];

		//iterate ancestor edges
		for (int i=0;i<headNode.out_edges_.size();i++){
			const Hypergraph::Edge& currentHeadEdge = in.edges_[headNode.out_edges_[i]];
			const int tailSize = currentHeadEdge.tail_nodes_.size();

			//				if(!NodePoppable(currentHeadEdge.head_node_)){
			//					continue;
			//				}

			bool hasDummyTail = false;
			//iterate tails of ancestor edge to compute new tail
			for(int j=0;j<tailSize;j++){
				int currentTailNodeID = currentHeadEdge.tail_nodes_[j];
				if(currentTailNodeID!=aCand.in_edge_->head_node_ && D[currentTailNodeID].IsEmpty()){
					hasDummyTail=true;
					break;
				}
			}
			//				if(hasDummyTail)continue;

			SharedArrayIterator** newTailIterators=new SharedArrayIterator*[tailSize];

			//iterate tails of ancestor edge to compute new tail
			for(int j=0;j<tailSize;j++){
				int currentTailNodeID = currentHeadEdge.tail_nodes_[j];
				if(currentTailNodeID==aCand.in_edge_->head_node_){
					newTailIterators[j]=new SharedArrayIterator(&aCand);
				}
				else if(D[currentTailNodeID].IsEmpty()){
					newTailIterators[j]=NULL;
				}
				else{
					newTailIterators[j]=D[currentTailNodeID].GetTailIterator();//check head compatibility if we add also not dummy head in D
#ifdef DEBUG_GP
					cerr<< "Get D[" << currentTailNodeID << "]" << *newTailIterators[j] << endl;
#endif
				}
			}

			//put dummy head cand
#ifdef DEBUG_GP
			cerr << "\tput dummy head cand:\n"; 
#endif
			GCandidate* newCandWithDummyHead=CreateCandidate(currentHeadEdge,NULL,newTailIterators);
			AddGCandidate(newCandWithDummyHead,cands/*,unique_cands*/);

//			if(!hasDummyTail){
//				//link new cand to fathers if are there
//				GCandidate* newCandWithHead;
//				if(currentHeadEdge.head_node_<H.size() && !H[currentHeadEdge.head_node_].IsEmpty()){
//#ifdef DEBUG_GP
//					cerr << "put cand with head:\n"; 
//					cerr<< "Get H[" << currentHeadEdge.head_node_ << "]\n";
//#endif
//					newCandWithHead=CreateCandidate(currentHeadEdge,H[currentHeadEdge.head_node_].GetHeadIterator(),CopyTailPTArray(newTailIterators, tailSize));
//					AddGCandidate(newCandWithHead,cands/*,unique_cands*/);
//				}
//			}

			//add cands with dummy tails
			//* NB tail that ref to propagating node cannot be dummy!
			for(int j=0;j<tailSize;j++){
				if(newTailIterators[j] &&
						currentHeadEdge.tail_nodes_[j]!= aCand.in_edge_->head_node_){//*
					SharedArrayIterator** newTailIteratorsWithDummy=CopyTailPTArray(newTailIterators,tailSize,j);
#ifdef DEBUG_GP
					cerr << "\tput cand with dummy tail "<< j<<" :\n"; 
					assert (newTailIteratorsWithDummy[j]==NULL);
#endif
					GCandidate* newCandWithHeadAndADummyTail=CreateCandidate(currentHeadEdge,DUMMY,newTailIteratorsWithDummy);
					AddGCandidate(newCandWithHeadAndADummyTail,cands/*,unique_cands*/);
				}
			}

		}

		//		}
		//		else{
		//#ifdef DEBUG_GP
		//			cerr << "\tput updated head:\n"; 
		//#endif
		//			GCandidate* oldHead= aCand.head_iterator_->GetCurrent();
		//			//			if(!NodePoppable(*oldHead)){
		//			//				return;
		//			//			}
		//			SharedArrayIterator** newTailIterators=CopyTailPTArray(oldHead->tail_iterators_,oldHead->TailSize());
		//			for(int i=0;i<oldHead->TailSize();i++){
		//				if(oldHead->in_edge_->tail_nodes_[i]==aCand.in_edge_->head_node_){
		//					assert(!newTailIterators[i]);
		//					newTailIterators[i]=new SharedArrayIterator(&aCand);
		//				}
		//				else{
		//					if(newTailIterators[i]){
		//						newTailIterators[i]->Freeze(); //stop iterator to children saw by aCand
		//					}
		//				}
		//			}
		//			SharedArrayIterator* newHeadIterator=DUMMY;
		//			if(oldHead->head_iterator_){
		//				newHeadIterator = new SharedArrayIterator(*oldHead->head_iterator_);
		//				newHeadIterator->Freeze();
		//			}
		//			GCandidate* newHeadCand=CreateCandidate(*oldHead->in_edge_,newHeadIterator,newTailIterators);
		//			AddGCandidate(newHeadCand,cands/*,unique_cands*/);
		//		}

	}

	const ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	int goal_id_;
	int out_goal_id_;	
	Hypergraph& out; //TODO should end with _ , refactor with eclipse

	vector<GCandidateSmartList> H;//maps (cand head) nodeID to cands with same node as dummy tail for match
	vector<GCandidateSmartList> D;//maps (cand tail) nodeID to cands with same node as (NB also not dummy!!! check match) head for match//TODO check match
	vector<GCandidateSmartList> Hc;
	vector<GCandidateSmartList> Dc;
	// maps nodes in in-HG to the
	// equivalent nodes (many due to state
	// splits) in the out-HG.

	//In GP we don't incorporate all cands, we need state of not-inc. cand. use state stored in cand
	//vector<string> node_states_;  // for each node in the out-HG what is

	int pop_limit_;
	GCandidateList free_; //store GCands pointer to free mem
};

struct NoPruningRescorer {
	NoPruningRescorer(const ModelSet& m, const SentenceMetadata &sm, const Hypergraph& i, Hypergraph* o) :
		models(m),
		smeta(sm),
		in(i),
		out(*o),
		nodemap(i.nodes_.size()) {
		cerr << "  Rescoring forest (full intersection)\n";
		node_states_.reserve(kRESERVE_NUM_NODES);
	}

	typedef unordered_map<string, int, boost::hash<string> > State2NodeIndex;

	void ExpandEdge(const Hypergraph::Edge& in_edge, bool is_goal, State2NodeIndex* state2node) {
		const int arity = in_edge.Arity();
		Hypergraph::TailNodeVector ends(arity);
		for (int i = 0; i < arity; ++i)
			ends[i] = nodemap[in_edge.tail_nodes_[i]].size();

		Hypergraph::TailNodeVector tail_iter(arity, 0);
		bool done = false;
		while (!done) {
			Hypergraph::TailNodeVector tail(arity);
			for (int i = 0; i < arity; ++i)
				tail[i] = nodemap[in_edge.tail_nodes_[i]][tail_iter[i]];
			Hypergraph::Edge* new_edge = out.AddEdge(in_edge.rule_, tail);
			new_edge->feature_values_ = in_edge.feature_values_;
			new_edge->i_ = in_edge.i_;
			new_edge->j_ = in_edge.j_;
			new_edge->prev_i_ = in_edge.prev_i_;
			new_edge->prev_j_ = in_edge.prev_j_;
			string head_state;
			if (is_goal) {
				assert(tail.size() == 1);
				const string& ant_state = node_states_[tail.front()];
				models.AddFinalFeatures(ant_state, new_edge);
			} else {
				prob_t edge_estimate; // this is a full intersection, so we disregard this
				models.AddFeaturesToEdge(smeta, out, node_states_, new_edge, &head_state, &edge_estimate);
			}
			int& head_plus1 = (*state2node)[head_state];
			if (!head_plus1) {
				head_plus1 = out.AddNode(in_edge.rule_->GetLHS())->id_ + 1;
				node_states_.push_back(head_state);
				nodemap[in_edge.head_node_].push_back(head_plus1 - 1);
			}
			const int head_index = head_plus1 - 1;
			out.ConnectEdgeToHeadNode(new_edge->id_, head_index);

			int ii = 0;
			for (; ii < arity; ++ii) {
				++tail_iter[ii];
				if (tail_iter[ii] < ends[ii]) break;
				tail_iter[ii] = 0;
			}
			done = (ii == arity);
		}
	}

	void ProcessOneNode(const int node_num, const bool is_goal) {
		State2NodeIndex state2node;
		const Hypergraph::Node& node = in.nodes_[node_num];
		for (int i = 0; i < node.in_edges_.size(); ++i) {
			const Hypergraph::Edge& edge = in.edges_[node.in_edges_[i]];
			ExpandEdge(edge, is_goal, &state2node);
		}
	}

	void Apply() {
		int num_nodes = in.nodes_.size();
		int goal_id = num_nodes - 1;
		int pregoal = goal_id - 1;
		int every = 1;
		if (num_nodes > 100) every = 10;
		assert(in.nodes_[pregoal].out_edges_.size() == 1);
		cerr << "    ";
		for (int i = 0; i < in.nodes_.size(); ++i) {
			if (i % every == 0) cerr << '.';
			ProcessOneNode(i, i == goal_id);
		}
		cerr << endl;
	}

private:
	const ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	Hypergraph& out;

	vector<vector<int> > nodemap;
	vector<string> node_states_;  // for each node in the out-HG what is
	// its q function value?
};

// each node in the graph has one of these, it keeps track of
void ApplyModelSet(const Hypergraph& in,
		const SentenceMetadata& smeta,
		const ModelSet& models,
		const IntersectionConfiguration& config,
		Hypergraph* out) {
	// TODO special handling when all models are stateless
	if (config.algorithm == 0) {
		NoPruningRescorer ma(models, smeta, in, out);
		ma.Apply();
	} else {
		int pl = config.pop_limit;
		if (pl > 100 && in.nodes_.size() > 80000) {
			cerr << "  Note: reducing pop_limit to " << pl << " for very large forest\n";
			pl = 30;
		}
		if (config.algorithm ==1){
			CubePruningRescorer ma(models, smeta, in, pl, out);
			ma.Apply();
		} else if (config.algorithm ==2){
			CubePruningRescorer ma(models, smeta, in, pl, out, FAST_CP);
			ma.Apply();
		} else if (config.algorithm ==3){
			CubePruningRescorer ma(models, smeta, in, pl, out, FAST_CP_2);
			ma.Apply();
		} else if (config.algorithm ==4){
			GuidedPruningRescorer ma(models, smeta, in, pl, out);
			ma.Apply();
		} else {
			cerr << "Don't understand intersection algorithm " << config.algorithm << endl;
			exit(1);
		}
	}
	out->is_linear_chain_ = in.is_linear_chain_;  // TODO remove when this is computed
	// automatically
}

#undef DEBUG_GP
