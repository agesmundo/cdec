////TODO: keep model state in forest?

//TODO: (for many nonterminals, or multi-rescoring pass) either global
//best-first, or group by (NT,span) - use prev forest outside as a (admissable,
//if models are a subset and weights are same) heuristic

#include "apply_models.h"

#include <vector>
#include <algorithm>
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <stack>

#include <boost/functional/hash.hpp>

#include "verbose.h"
//#include "hg.h"
#include "undirected_candidate.h"
//#include "ff.h"
#include "hg_intersect.h"
#include "sentence_metadata.h"
#include "inside_outside.h"
#include "semirings.h"

using namespace std;
using namespace std::tr1;

struct Candidate;
typedef SmallVectorInt JVector;
typedef vector<Candidate*> CandidateHeap;
typedef vector<Candidate*> CandidateList;

// default vector size (* sizeof string is memory used)
static const size_t kRESERVE_NUM_NODES = 500000ul;

//Undirected Candidate
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
	int node_index_;                     // -1 until incorporated
	// into the +LM forest
	const Hypergraph::Edge* in_edge_;    // in -LM forest
	Hypergraph::Edge out_edge_;
	FFState state_;
	const JVector j_;
	prob_t vit_prob_;            // these are fixed until the cand
	// is popped, then they may be updated
	prob_t est_prob_;

	Candidate(const Hypergraph::Edge& e,
			const JVector& j,
			const Hypergraph& out_hg,
			const vector<CandidateList>& D,
			const FFStates& node_states,
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
			const FFStates& node_states,
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
			const FFState& ant_state = node_states[tail.front()];
			models.AddFinalFeatures(ant_state, &out_edge_, smeta);
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

struct HeapCandCompare {
	bool operator()(const Candidate* l, const Candidate* r) const {
		return l->est_prob_ < r->est_prob_;
	}
	bool operator()(const UCandidate* l, const UCandidate* r) const {
		return l->action_prob_ < r->action_prob_;
	}
};

struct EstProbSorter {
	bool operator()(const Candidate* l, const Candidate* r) const {
		return l->est_prob_ > r->est_prob_;
	}
	bool operator()(const UCandidate* l, const UCandidate* r) const {
		return l->action_prob_ > r->action_prob_;
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
};
//struct UCandidateUniquenessHash {//TODO GU is this correct??
//  size_t operator()(const UCandidate* c) const {
//    size_t x = 5381;
//    x = ((x << 5) + x) ^ c->in_edge_->id_;
//    for (int i = 0; i < c->context_links_.size(); ++i){
//    	assert(c->context_links_[i]->ucand_index_+1>=0);
//      x = ((x << 5) + x) ^ (c->context_links_[i]->ucand_index_+1);
//    }
//    return x;
//  }
//};

struct CandidateUniquenessEquals {
	bool operator()(const Candidate* a, const Candidate* b) const {
		return (a->in_edge_ == b->in_edge_) && (a->j_ == b->j_);
	}
};
struct UCandidateUniquenessEquals {
	bool operator()(const UCandidate* a, const UCandidate* b) const {
		return (a->in_edge_ == b->in_edge_) && (a->context_links_ == b->context_links_);
	}
};

typedef unordered_set<const Candidate*, CandidateUniquenessHash, CandidateUniquenessEquals> UniqueCandidateSet;
typedef unordered_map<FFState, Candidate*, boost::hash<FFState> > State2Node;

//TODO remove later, not merging states in GU
//typedef unordered_set<const UCandidate*, UCandidateUniquenessHash, UCandidateUniquenessEquals> UniqueUCandidateSet;
typedef unordered_map<FFState, UCandidate*, boost::hash<FFState> > UState2Node;

class GreedyUndirectedRescorer {

public:
	GreedyUndirectedRescorer /*const*/
	(ModelSet & m, const SentenceMetadata & sm, const Hypergraph & i, bool is_training, Hypergraph *o, string wf="")
	:models(m), smeta(sm), in(i), out(*o), //D(in.nodes_.size()),
	 is_training_(is_training), weight_file_name_(wf)
	{

		Inside<MaxSum<prob_t>, EdgeProb>(in, &inside_);
		UCandidate::inside_=&inside_;

		Outside<MaxSum<prob_t>, EdgeProb>(in, inside_, &outside_);
		UCandidate::outside_=&outside_;

		if(!SILENT)
			cerr << "  Applying feature functions (training = " << is_training_ << ')' << endl;

		//ucands_states_.reserve(kRESERVE_NUM_NODES);
	}

	void Apply()
	{
#ifdef DEBUG_GU
				cerr << "\tInit weight vector: ";
				models.PrintWeights(cerr);
#endif
		int wrong_count=0;
		int num_nodes = in.nodes_.size();
		assert(num_nodes >= 2);
		int goal_id = num_nodes - 1;
		int pregoal = goal_id - 1;
		int every = 1;
		if(num_nodes > 100)
			every = 10;

		assert(in.nodes_[pregoal].out_edges_.size() == 1);
		UCandidateHeap cands; //unique queue/heap of candidates

		//find nodes that intersect with reference lattice
		if (is_training_) {
			assert(smeta.HasReference()==true);
			correct_edges_mask_ = new vector<bool> (in.edges_.size(), false);
			HG::HighlightIntersection(smeta.GetReference(), in, correct_edges_mask_);
		}

		InitCands(cands); //put leafs candidates in the queue
		UCandidate *topCand;
		for (;!cands.empty();) {//TODO? borders should not be empty, first pass ok to be empty

#ifdef DEBUG_GU
			cerr<< "\n\n/////////////////////////////////////////////////////////////////////\n";
			cerr << " L | cands.size(): "<< cands.size()<<"\n";
#endif

			//best action candidate
			swap(*max_element(cands.begin(), cands.end(),HeapCandCompare()), cands.back());//TODO use heap?? NO!
			topCand=cands.back();

#ifdef DEBUG_GU
			cerr << "BEST IS: " << *topCand << "\n";
//			cerr << "FIRST :" << *topCand->in_edge_ << endl;
//			if (IsCorrect(*topCand)) cerr << " IS CORR" <<endl;
//			else cerr << " WRONGGGG! " <<endl;
#endif

			if(!is_training_ || IsCorrect(*topCand) || wrong_count>=1){

//				if(is_training_ &&   wrong_count==0 /*&& IsCorrect(*topCand) */){ //wrong count = 0 it is correct
//					models.UpAvgCnt();//TODO check if is averaged models.is_avg
//				}

				wrong_count=0;

				//store info about links expansions and re-expansions (cause updated source)
				vector < pair < UCandidate*, int > > links_to_expand;

				//POP BEST FROM QUEUE
				cands.pop_back();

				//IF CORRECT || !TRAINING:
#ifdef DEBUG_GU
				cerr << "\nCORRECT || !TRAINING:" << endl;
#endif

#ifdef DEBUG_GU
				cerr << "\nREMOVE ALTERNATIVE CANDIDATES\n" << "\tInit cands.size(): "<< cands.size()<<endl;
#endif
				if(topCand->HasSource()){//this is not first loop leaf,
					//delete alternatives cands with same source node
					DelCandsFromNodeId(cands,topCand->GetSourceNodeId());

#ifdef DEBUG_GU
//					cerr << "\tTopCand (still exists): "<< *topCand<<endl;
					cerr << "\tFinal cands.size(): "<< cands.size()<<endl;
#endif

					//UPDATE SOURCE UCAND LINK
#ifdef DEBUG_GU
					cerr << "\nUPDATE SOURCE UCAND LINK\n\tSource: "<< *topCand->GetSourceUCand()<<endl;
#endif
					assert(topCand->GetSourceUCand()->CreateLink(topCand));
#ifdef DEBUG_GU
					cerr << "\tSource updated: "<< *topCand->GetSourceUCand()<<endl;
#endif

					//STATE UPDATE BACKWARD PROPAGATION IN TREE
#ifdef DEBUG_GU
					cerr << "\nSTATE UPDATE PROPAGATION"<< endl;
#endif
					stack<UCandidate*> candsToUpdate;
					candsToUpdate.push(topCand->GetSourceUCand());//TODO skip all this if state for source cand is all 0
#ifdef DEBUG_GU
					assert(links_to_expand.size()==0);
#endif
					while(!candsToUpdate.empty()){
						UCandidate* curr = candsToUpdate.top();
						candsToUpdate.pop();
						//TODO? speedup: avoid propagation for tails and root (no missing links, no further porpagations)

						curr->UpdateStates(candsToUpdate,links_to_expand);
#ifdef DEBUG_GU
						cerr << "\tCurrent cand (updated): "<< *curr<< endl;
#endif
					}
#ifdef DEBUG_GU
	  			cerr <<"\tlinks_to_expand.size : "<< links_to_expand.size()<<endl;
	  			for(int i=0;i<links_to_expand.size();i++){
	  				assert(links_to_expand[i].first!=topCand);
	  			}
#endif

					//REMOVE FROM QUEUE CANDIDATE ACTIONS TO BE UPDATED WITH RE-EXPANSION
#ifdef DEBUG_GU
					cerr << "\nREMOVE FROM QUEUE CANDIDATE ACTIONS TO BE UPDATED"<< endl;
					cerr << "\tInit cands.size(): "<< cands.size()<<endl;
#endif
					for (int i=0;i<links_to_expand.size();i++){
						int node_id = links_to_expand[i].first->GetNodeIdFromLinkId(links_to_expand[i].second);
#ifdef DEBUG_GU
						assert(node_id>=0);
						cerr << "\tRemoving cands from link: (" << links_to_expand[i].first << " , "<< links_to_expand[i].second <<")" <<endl;
#endif
						DelCandsFromNodeId(cands,node_id);
					}
#ifdef DEBUG_GU
					cerr << "\tFinal cands.size(): "<< cands.size()<<endl;
#endif

				}
				else{//this is first loop, delete all queue, no propagation

#ifdef DEBUG_GU
					cerr << " \tFirst loop: empty candidates queue" <<endl;
#endif
					for (int i = 0; i < cands.size(); ++i){
#ifdef DEBUG_GU
						cerr << "\tDelete : (id:"<<i<<")  " <<cands[i]<<endl;
#endif
						delete cands[i];
					}

					cands.erase(cands.begin(),cands.end());

#ifdef DEBUG_GU
					cerr << "\tTopCand (still exists): "<< *topCand<<endl;
					cerr << "\tFinal cands.size(): "<< cands.size()<<endl;
#endif
				}

				//ADD IN QUEUE UCANDS FROM TOP UCAND
#ifdef DEBUG_GU
				cerr << "\nEXPAND UCANDS LINKS"<<endl;
#endif

				//insert topCand link expansions in the list
				for (int k=0;k<topCand->context_links_.size();k++){
					if(topCand->context_links_[k]==NULL){//missing link, available for expansion
						links_to_expand.push_back(pair< UCandidate*, int >(topCand,k));
#ifdef DEBUG_GU
	  			cerr <<"\tlinks_to_expand add: ("<< topCand<< " , "<< k << ")" <<endl;
#endif
					}
				}
#ifdef DEBUG_GU
	  			cerr <<"\tlinks_to_expand.size : "<< links_to_expand.size()<<endl;
#endif

				//execute link expansions in the list
#ifdef DEBUG_GU
	  			cerr << "\tcands.size(): "<< cands.size()<<endl;
#endif
				for(int k=0; k<links_to_expand.size();k++ ){
#ifdef DEBUG_GU
	  			cerr <<"\tlinks expansion: ("<< links_to_expand[k].first<< " , "<< links_to_expand[k].second << ")" <<endl;
#endif
					GenerateCandsFromLinkExpansion(links_to_expand[k].first,links_to_expand[k].second,cands);
#ifdef DEBUG_GU
	  			cerr << "\tcands.size(): "<< cands.size()<<endl;
#endif
				}

			}else{
				wrong_count++;
				//IF WRONG
#ifdef DEBUG_GU
				cerr<< "\nLEARN FROM WRONG PREDICTION:" << endl;
#endif

				//find first correct
#ifdef DEBUG_GU
				cerr<< "\nFIND FIRST CORRECT" << endl;
#endif
				UCandidate* correctCand=NULL;
				sort(cands.begin(), cands.end(), EstProbSorter()); //TODO skip last element (is wrong best)//TODO try iteration of pop heap (faster?) {make_heap(),pop_heap() [O(n+logn)]}
				for(int i = 1;i<cands.size() ; i++){//start from 1 to skip last that is topCand (wrong)
#ifdef DEBUG_GU
					cerr << "\tIs correct?: " << *cands[i] << endl;
					if(cands[i]->action_prob_ == topCand->action_prob_){
						cerr << "!!!DISCARTED BECAUSE HAS SAME ACTION PROB!!!"	<<endl;
					}
#endif
					if(IsCorrect(*cands[i])
					&& (cands[i]->action_prob_ != topCand->action_prob_)){ //TODO need to avoid items with same features! brute if have not same action score then surely not same feats
						correctCand = cands[i];
						break;
					}
				}
//				assert(correctCand!=NULL);//TODO GU what do in case there are no correct edges?

				if(correctCand==NULL){
					(*correct_edges_mask_)[topCand->in_edge_->id_]==true;
				}else{
				//update weights vector
#ifdef DEBUG_GU
				cerr<< "\nUPDATE WEIGHT VECTOR" << endl;
#endif

				double loss=1;
//PASSIVE AGGRESSIVE PA
//				double margin = 1;
//				double loss = log(topCand->action_prob_) - log(correctCand->action_prob_) + margin;
//				assert(loss>=0);
//#ifdef DEBUG_GU
//				cerr << "\tLoss: " << loss << endl;
//#endif

				SparseVector<Featval> diff (correctCand->feature_values_);
				//    		cerr << diff << endl;
//				diff +=correctCand->est_vals_;
				//    		cerr << diff << endl;
				diff -=topCand->feature_values_;
				//    		cerr << diff << endl;
//				diff -=topCand->est_vals_;
#ifdef DEBUG_GU
				cerr << "\tDiff vector: " << diff << endl;
				cerr << "\tWeight vector: ";
				models.PrintWeights(cerr);
#endif
				models.UpdateWeight(diff,loss);
#ifdef DEBUG_GU
				cerr << "\tUpdated weight vector: ";
				models.PrintWeights(cerr);
#endif

				//update queue
#ifdef DEBUG_GU
				cerr<< "\nRESCORE QUEUE" << endl;
#endif
				for(UCandidateHeap::iterator it = cands.begin();it!=cands.end();it++){
					UCandidate& ucand = **it;

//					prob_t estimate = prob_t::One();
//					estimate.logeq(models.ScoreVector(ucand.est_vals_));

//					prob_t local = prob_t::One();
//					local.logeq(models.ScoreVector(ucand.feature_values_));

//					ucand.action_prob_ = local;// * estimate; //sum exps
					ucand.action_prob_.logeq(models.ScoreVector(ucand.feature_values_));
#ifdef DEBUG_GU
					cerr << "\t-\n\tRescored " << ucand <<endl;
					IsCorrect(ucand);
#endif
				}
			}
			}
		}
		if(is_training_){
			delete correct_edges_mask_;
			models.WriteToFile(weight_file_name_);
		}
		//LG transform the UCands structure in the out_hg
		BuildOutHG(topCand);
	}

private:

	/**
	 * given link_id and associated ucand, generate ucands from link expansion and insert in cands queue
	 */
	void GenerateCandsFromLinkExpansion(UCandidate * topCand, int link_id, UCandidateHeap & cands)
	{

#ifdef DEBUG_GU
			cerr<< "\t-\n\tExpanding UCandidate ("<< topCand<<") via link [" << link_id<<"]:"<<endl;
#endif

		//head
		if(link_id==0){
			const Hypergraph::Node& head_node = in.nodes_[topCand->in_edge_->head_node_];
#ifdef DEBUG_GU
			cerr<< "\tHEAD PROPAGATION:"<<endl;
			cerr << "\thead_node.out_edges_.size(): " <<head_node.out_edges_.size()<<endl;
#endif
			for (int j=0; j<head_node.out_edges_.size();j++){
				const Hypergraph::Edge& currEdge = in.edges_[head_node.out_edges_[j]];
#ifdef DEBUG_GU
				cerr << "\tcurrEdge("<<j<<"): "<< currEdge<<endl;
#endif
				LinksVector context(currEdge.Arity()+1, NULL);
				assert(currEdge.Arity()<=2);//constraint to binary rules only  //TODO? add in debug
				assert(currEdge.Arity()>=1);//must have one child since reached via head propagation
				int source_link;
				if(currEdge.tail_nodes_[0]==head_node.id_){//link with source cand
					context[1]=topCand;
					source_link=1;
				}else{
					assert(currEdge.tail_nodes_[1]==head_node.id_);
					context[2]=topCand;
					source_link=2;
				}
				if(IsGoal(currEdge))context[0]=(UCandidate*)-1;//TODO speedup, if goal there is nothing to choose just link it score it, do not put in the queue
				cands.push_back(new UCandidate(currEdge, context,/* D, ucands_states_,*/ smeta, models, source_link/*, false*/));
#ifdef DEBUG_GU
				cerr << "\tPush UCand (" << cands.size() << ") :" << *cands.back() << endl;
#endif
			}
		}
		//tails
		else /*if(k==1 ||k==2)*/{
			const Hypergraph::Node& tail_node = in.nodes_[topCand->in_edge_->tail_nodes_[link_id-1]];
#ifdef DEBUG_GU
			assert (link_id==1 || link_id==2);
			if(link_id==1) cerr<< "\tLEFT CHILD PROPAGATION"<<endl;
			if(link_id==2) cerr<< "\tRIGHT CHILD PROPAGATION"<<endl;
			cerr << "\ttail_node.in_edges_.size(): " <<tail_node.in_edges_.size()<<endl;
#endif
			for (int j=0; j<tail_node.in_edges_.size();j++){
				const Hypergraph::Edge& currEdge = in.edges_[tail_node.in_edges_[j]];
#ifdef DEBUG_GU
				cerr << "currEdge("<<j<<"): "<< currEdge<<endl;
#endif
				LinksVector context(currEdge.Arity()+1, NULL);
				//XXX								assert(currEdge.Arity()<=2);//constraint to binary rules only
				assert(currEdge.head_node_==tail_node.id_);//TODO? add in debug
				context[0]=topCand;
				int	source_link=0;

				cands.push_back(new UCandidate(currEdge, context,/* D, ucands_states_,*/ smeta, models, source_link/*, false*/));
#ifdef DEBUG_GU
				cerr << "Push UCand (" << cands.size() << ") :" << *cands.back() << endl;
#endif
			}
		}
	}
	/**
	 * given cands list and node_id delete all cands that source from that node
	 */
	void DelCandsFromNodeId(UCandidateHeap & cands, int node_id)
	{

		for(int i = 0;i < cands.size();){
			if(cands[i]->GetSourceNodeId() == node_id){
				//TODO? try out algorithm with forward iterator(see ref remove)
#ifdef DEBUG_GU
				cerr << "\tDelete : (id:" << i << ")  " << cands[i] << endl;
#endif
				delete cands[i];
				swap(cands[i], cands.back());
				cands.pop_back();
			}else{
				i++;
			}
		}

	}

	//given any enrty point of the UCands structure builds the out_hg
	//and free memory
	void BuildOutHG(UCandidate *first)
	{
#ifdef DEBUG_GU
		cerr << "BUILD OUT HG" << endl;
#endif
		map<int,int> inNid2outNid; //in Node id to out Node id
		UCandidateList ucands_stack;
		ucands_stack.push_back(first);
		int goal_node_id = -1;
		while(!ucands_stack.empty()){
			UCandidate* ucand =ucands_stack.back();
			ucands_stack.pop_back();

			//copy in_edge in out_hg
			out.edges_.push_back(*ucand->in_edge_);
			Hypergraph::Edge* new_edge = &out.edges_.back();
			new_edge->id_ = out.edges_.size()-1;

#ifdef DEBUG_GU
			cerr << "/////////////////////////////////////////////////////////////////////\n";
			cerr << "NEXT UCAND: " << *ucand << endl;
			cerr << "IN  EDGE : " << *ucand->in_edge_ << endl;
			cerr << "OUT EDGE : " << *new_edge << endl;
#endif

			new_edge->edge_prob_ = ucand->in_edge_->edge_prob_;

			//link with head
			{
#ifdef DEBUG_GU
				assert(ucand->context_links_[0]!=NULL);
#endif
				UCandidate* head_ucand=ucand->context_links_[0];
				bool is_goal=(head_ucand==(UCandidate*)-1);
#ifdef DEBUG_GU
				if (is_goal) cerr<< "HEAD UCAND : " << "Goal"<<endl;
				else cerr<< "HEAD UCAND : " << head_ucand<<endl;
#endif
				map<int,int>::iterator it = inNid2outNid.find(new_edge->head_node_);

				int head_node_id;
				if(it!=inNid2outNid.end()){
					head_node_id = it->second;
				}else{
					head_node_id = out.AddNode(in.nodes_[ucand->in_edge_->head_node_].cat_)->id_;
					inNid2outNid[new_edge->head_node_]=head_node_id;
					if(is_goal){
						goal_node_id=head_node_id;
					}else{
						ucands_stack.push_back(head_ucand);
					}
				}
				new_edge->head_node_ = head_node_id;
				out.nodes_[head_node_id].in_edges_.push_back(new_edge->id_);
			}

			//link with left child
			if(ucand->context_links_.size()>=2){
#ifdef DEBUG_GU
				assert(ucand->context_links_[1]!=NULL);
#endif
				UCandidate* left_ucand=ucand->context_links_[1];
#ifdef DEBUG_GU
				cerr<< "LEFT UCAND : " << left_ucand<<endl;
#endif
				map<int,int>::iterator it = inNid2outNid.find(new_edge->tail_nodes_[0]);
				int left_node_id;
				if(it!=inNid2outNid.end()){
					left_node_id = it->second;
				}else{
					left_node_id = out.AddNode(in.nodes_[ucand->in_edge_->tail_nodes_[0]].cat_)->id_;
					inNid2outNid[new_edge->tail_nodes_[0]]=left_node_id;
					ucands_stack.push_back(left_ucand);
				}
				new_edge->tail_nodes_[0] = left_node_id;
				out.nodes_[left_node_id].out_edges_.push_back(new_edge->id_);
			}

			//link with right child
			if(ucand->context_links_.size()>=3){
#ifdef DEBUG_GU
				assert(ucand->context_links_[2]!=NULL);
#endif
				UCandidate* right_ucand=ucand->context_links_[2];
#ifdef DEBUG_GU
				cerr<< "RIGHT UCAND : " << right_ucand<<endl;
#endif
				map<int,int>::iterator it = inNid2outNid.find(new_edge->tail_nodes_[1]);
				int right_node_id;
				if(it!=inNid2outNid.end()){
					right_node_id = it->second;
				}else{
					right_node_id = out.AddNode(in.nodes_[ucand->in_edge_->tail_nodes_[1]].cat_)->id_;
					inNid2outNid[new_edge->tail_nodes_[1]]=right_node_id;
					ucands_stack.push_back(right_ucand);
				}
				new_edge->tail_nodes_[1] = right_node_id;
				out.nodes_[right_node_id].out_edges_.push_back(new_edge->id_);
			}
#ifdef DEBUG_GU
			cerr << "OUT EDGE : " << *new_edge << endl;
			cerr << "DELETE UCAND : " << *ucand <<endl;
#endif
			delete ucand;
		}
		assert(goal_node_id>=0);
		out.PruneUnreachable(goal_node_id);
	}

	//  //given any enrty point of the UCands structure builds the out_hg
	//  void BuildOutHG(UCandidate* first){
	//#ifdef DEBUG_GU
	//	  cerr << "BUILD OUT HG" << endl;
	//#endif
	//	  map<UCandidate*,int> ucand2edge;//TODO GU hash map is faster?
	//	  UCandidateList ucands_stack;
	//	  ucands_stack.push_back(first);
	//	  int goal_node_id=-1;
	//
	//	  while(!ucands_stack.empty()){
	//		  UCandidate* ucand =ucands_stack.back();
	//		  ucands_stack.pop_back();
	//
	//		  //copy in_edge in out_hg
	//		  out.edges_.push_back(*ucand->in_edge_);
	//		  Hypergraph::Edge* new_edge = &out.edges_.back();
	//		  new_edge->id_ = out.edges_.size()-1;
	//		  //clear links to nodes TODO this should be useless
	//		  new_edge->head_node_=-1;
	//		  for(int i=0;i<new_edge->tail_nodes_.size();i++){
	//			  new_edge->tail_nodes_[i]=-1;
	//		  }
	//
	//#ifdef DEBUG_GU
	//	  cerr << "/////////////////////////////////////////////////////////////////////\n";
	//	  cerr << "NEXT UCAND: " << *ucand << endl;
	//	  cerr << "IN  EDGE : " << *ucand->in_edge_ << endl;
	//	  cerr << "OUT EDGE : " << *new_edge << endl;
	//#endif
	//
	//		  ucand2edge[ucand]=new_edge->id_;
	//		  new_edge->edge_prob_ = ucand->in_edge_->edge_prob_;
	//
	//		  //link with head
	//		  {
	//#ifdef DEBUG_GU
	//			  assert(ucand->context_links_[0]!=NULL);
	//#endif
	//			  map<UCandidate*,int>::iterator it = ucand2edge.end();
	//			  UCandidate* head_ucand=ucand->context_links_[0];
	//			  bool is_goal=(head_ucand==(UCandidate*)-1);
	//#ifdef DEBUG_GU
	//			  if (is_goal) cerr<< "HEAD UCAND : " << "Goal"<<endl;
	//			  else cerr<< "HEAD UCAND : " << *head_ucand<<endl;
	//#endif
	//			  if(!is_goal){
	//				  it = ucand2edge.find(head_ucand);
	//			  }
	//			  int head_node_id;
	//			  if(it!=ucand2edge.end()){
	//				  const Hypergraph::Edge& head_edge = out.edges_[it->second];
	//				  if(head_ucand->context_links_[1]==ucand){
	//					  head_node_id=head_edge.tail_nodes_[0];
	//				  }else {
	//					  assert(head_ucand->context_links_[2]==ucand);
	//					  head_node_id=head_edge.tail_nodes_[1];
	//				  }
	//			  }else{
	//				  head_node_id = out.AddNode(in.nodes_[ucand->in_edge_->head_node_].cat_)->id_;
	//				  if(is_goal){
	//					  goal_node_id=head_node_id;
	//				  }else{
	//					  ucands_stack.push_back(head_ucand);
	//				  }
	//			  }
	//			  out.ConnectEdgeToHeadNode(new_edge, head_node_id);
	//		  }
	//
	//		  //link with left child
	//		  if(ucand->context_links_.size()>=2){
	//#ifdef DEBUG_GU
	//			  assert(ucand->context_links_[1]!=NULL);
	//#endif
	//			  UCandidate* left_ucand=ucand->context_links_[1];
	//#ifdef DEBUG_GU
	//			  cerr<< "LEFT UCAND : " << *left_ucand<<endl;
	//#endif
	//			  map<UCandidate*,int>::iterator it = ucand2edge.find(left_ucand);
	//			  int left_node_id;
	//			  if(it!=ucand2edge.end()){
	//				  left_node_id = out.edges_[it->second].head_node_;
	//			  }else{
	//				  left_node_id = out.AddNode(in.nodes_[ucand->in_edge_->tail_nodes_[0]].cat_)->id_;
	//				  ucands_stack.push_back(left_ucand);
	//			  }
	//			  new_edge->tail_nodes_[0] = left_node_id;
	//			  out.nodes_[left_node_id].out_edges_.push_back(new_edge->id_);
	//		  }
	//
	//		  //link with right child
	//		  if(ucand->context_links_.size()>=3){
	//#ifdef DEBUG_GU
	//			  assert(ucand->context_links_[2]!=NULL);
	//#endif
	//			  UCandidate* right_ucand=ucand->context_links_[2];
	//#ifdef DEBUG_GU
	//			  cerr<< "RIGHT UCAND : " << *right_ucand<<endl;
	//#endif
	//			  map<UCandidate*,int>::iterator it = ucand2edge.find(right_ucand);
	//			  int right_node_id;
	//			  if(it!=ucand2edge.end()){
	//				  right_node_id = out.edges_[it->second].head_node_;
	//			  }else{
	//				  right_node_id = out.AddNode(in.nodes_[ucand->in_edge_->tail_nodes_[1]].cat_)->id_;
	//				  ucands_stack.push_back(right_ucand);
	//			  }
	//			  new_edge->tail_nodes_[1] = right_node_id;
	//			  out.nodes_[right_node_id].out_edges_.push_back(new_edge->id_);
	//		  }
	//#ifdef DEBUG_GU
	//		  cerr << "OUT EDGE : " << *new_edge << endl;
	//#endif
	//	  }
	//
	//	  assert(goal_node_id>=0);
	//	  out.PruneUnreachable(goal_node_id);
	//  }
	//check if cand is correct (for training)
	bool IsCorrect(const UCandidate & ucand)
	{
#ifdef DEBUG_GU
		cerr << "\tIs correct? : (" << &ucand << ")";
		if((*correct_edges_mask_)[ucand.in_edge_->id_]){
			cerr << " true" << endl;
		}else{
			cerr << " false" << endl;
		}
#endif
		return (*correct_edges_mask_)[ucand.in_edge_->id_];
	}

	//initialize candidate heap with leafs
	void InitCands(UCandidateHeap & cands)
	{
#ifdef DEBUG_GU
		cerr << "InitCands(): " << "\n";
#endif
		for (int i = 0; i < in.edges_.size(); ++i) {//loop edges
			const Hypergraph::Edge& edge= in.edges_.at(i);
			if(edge.tail_nodes_.size()==0){//leafs
				//          		const Hypergraph::Edge& edge = in.edges_[i];
				const LinksVector context(/*edge.Arity()+*/1, NULL);//head link only
#ifdef DEBUG_GU
				//          		cerr << "=========================\tLEAF FROM :" <<edge<<endl;
#endif
				cands.push_back(new UCandidate(edge, context,/* D, ucands_states_,*/ smeta, models, -1/*, false*/));
#ifdef DEBUG_GU
				cerr << "Push Init UCand (" << cands.size() << ") :" << *cands.back() << endl;
#endif
			}
		}
		//          make_heap(cands.begin(), cands.end(), HeapCandCompare()); //useless not using heap
#ifdef DEBUG_GU
		cerr << "=========================\n cands.size(): " << cands.size() << "\n=========================\n";
#endif
	}

	bool IsGoal(const Hypergraph::Edge& edge){
		return edge.head_node_==in.nodes_.size()-1;
	}
	//  //order cands and pop best
	//  inline UCandidate* PopBest(UCandidateHeap& cands){
	//	  make_heap(cands.begin(), cands.end(), HeapCandCompare());
	//	  pop_heap(cands.begin(), cands.end(), HeapCandCompare());
	//	  UCandidate* ucand = cands.back(); //accepted cand
	//	  cands.pop_back();
	//#ifdef DEBUG_GU
	////	  cerr << "free_[" << free_.size() << "] = " << ucand <<endl;
	//	  cerr << "PopBest(): " << *ucand << "\n";
	//#endif
	////	  free_.push_back(ucand);
	//	  return ucand;
	//  }
	//
	//  void FreeAll(UCandidateHeap& cands) {
	//	  for (int i = 0; i < cands.size(); ++i){
	//		  cerr << "FREE : " <<*cands[i]<<endl;
	//		  delete cands[i];
	//	  }
	//  }

	//  void IncorporateIntoPlusLMForest(UCandidate* item, UState2Node* s2n, UCandidateList* freelist) {
	//    Hypergraph::Edge* new_edge = out.AddEdge(item->out_edge_);
	//    new_edge->edge_prob_ = item->out_edge_.edge_prob_;
	//    UCandidate*& o_item = (*s2n)[item->state_];
	//    if (!o_item) o_item = item;
	//
	//    int& node_id = o_item->node_index_;
	//    if (node_id < 0) {
	//      Hypergraph::Node* new_node = out.AddNode(in.nodes_[item->in_edge_->head_node_].cat_);
	//      node_states_.push_back(item->state_);
	//      node_id = new_node->id_;
	//    }
	//#if 0
	//    Hypergraph::Node* node = &out.nodes_[node_id];
	//    out.ConnectEdgeToHeadNode(new_edge, node);
	//#else
	//    out.ConnectEdgeToHeadNode(new_edge, node_id);
	//#endif
	//    // update candidate if we have a better derivation
	//    // note: the difference between the vit score and the estimated
	//    // score is the same for all items with a common residual DP
	//    // state
	//    if (item->vit_prob_ > o_item->vit_prob_) {
	//      assert(o_item->state_ == item->state_);    // sanity check!
	//      o_item->est_prob_ = item->est_prob_;
	//      o_item->vit_prob_ = item->vit_prob_;
	//    }
	//    if (item != o_item) freelist->push_back(item);
	//  }

	//  void KBest(const int vert_index, const bool is_goal) {
	//    // cerr << "KBest(" << vert_index << ")\n";
	//    UCandidateList& D_v = D[vert_index];
	//    assert(D_v.empty());
	//    const Hypergraph::Node& v = in.nodes_[vert_index];
	//    // cerr << "  has " << v.in_edges_.size() << " in-coming edges\n";
	//    const vector<int>& in_edges = v.in_edges_;
	//    UCandidateHeap cand;
	//    UCandidateList freelist;
	//    cand.reserve(in_edges.size());
	//    UniqueUCandidateSet unique_cands;
	//    for (int i = 0; i < in_edges.size(); ++i) {
	//      const Hypergraph::Edge& edge = in.edges_[in_edges[i]];
	//      const JVector j(edge.tail_nodes_.size(), 0);
	//      cand.push_back(new UCandidate(edge, j, D, node_states_, smeta, models, is_goal));
	//      assert(unique_cands.insert(cand.back()).second);  // these should all be unique!
	//    }
	////    cerr << "  making heap of " << cand.size() << " candidates\n";
	//    make_heap(cand.begin(), cand.end(), HeapCandCompare());
	//    UState2Node state2node;   // "buf" in Figure 2
	//    int pops = 0;
	//    while(!cand.empty() && pops < 100) {
	//      pop_heap(cand.begin(), cand.end(), HeapCandCompare());
	//      UCandidate* item = cand.back();
	//      cand.pop_back();
	//      // cerr << "POPPED: " << *item << endl;
	//      PushSucc(*item, is_goal, &cand, &unique_cands);
	//      IncorporateIntoPlusLMForest(item, &state2node, &freelist);
	//      ++pops;
	//    }
	//    D_v.resize(state2node.size());
	//    int c = 0;
	//    for (UState2Node::iterator i = state2node.begin(); i != state2node.end(); ++i)
	//      D_v[c++] = i->second;
	//    sort(D_v.begin(), D_v.end(), EstProbSorter());
	//    // cerr << "  expanded to " << D_v.size() << " nodes\n";
	//
	//    for (int i = 0; i < cand.size(); ++i)
	//      delete cand[i];
	//    // freelist is necessary since even after an item merged, it still stays in
	//    // the unique set so it can't be deleted til now
	//    for (int i = 0; i < freelist.size(); ++i)
	//      delete freelist[i];
	//  }

	//  void PushSucc(const UCandidate& item, const bool is_goal, UCandidateHeap* pcand, UniqueUCandidateSet* cs) {
	//    UCandidateHeap& cand = *pcand;
	//    for (int i = 0; i < item.j_.size(); ++i) {
	//      JVector j = item.j_;
	//      ++j[i];
	//      if (j[i] < D[item.in_edge_->tail_nodes_[i]].size()) {
	//        UCandidate query_unique(*item.in_edge_, j);
	//        if (cs->count(&query_unique) == 0) {
	//          UCandidate* new_cand = new UCandidate(*item.in_edge_, j, D, node_states_, smeta, models, is_goal);
	//          cand.push_back(new_cand);
	//          push_heap(cand.begin(), cand.end(), HeapCandCompare());
	//          assert(cs->insert(new_cand).second);  // insert into uniqueness set, sanity check
	//        }
	//      }
	//    }
	//  }

	/*const*/ ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	Hypergraph& out;

	//in GU is useless to rm
	//  vector<UCandidateList> D;   // maps nodes in in-HG to the
	// equivalent nodes (many due to state
	// splits) in the out-HG.

	//TODO delete and put states in the UCandidates
	//FFStates ucands_states_;  // for each node in the out-HG what is
	// its q function value?

	const bool is_training_;
	vector<bool>* correct_edges_mask_;//used only for training
	string& weight_file_name_;

	//UCandidateList free_; //store UCands pointer to free mem

	//NB! indexed with node ids ,size is hg.nodes_.size()
	vector<MaxSum<prob_t> > inside_;//inside scores of in_hg
	vector<MaxSum<prob_t> > outside_;//inside scores of in_hg
};

class CubePruningRescorer {

public:
	CubePruningRescorer(const ModelSet& m,
			const SentenceMetadata& sm,
			const Hypergraph& i,
			int pop_limit,
			Hypergraph* o) :
				models(m),
				smeta(sm),
				in(i),
				out(*o),
				D(in.nodes_.size()),
				pop_limit_(pop_limit) {
		if (!SILENT) cerr << "  Applying feature functions (cube pruning, pop_limit = " << pop_limit_ << ')' << endl;
		node_states_.reserve(kRESERVE_NUM_NODES);
	}

	void Apply() {
		int num_nodes = in.nodes_.size();
		assert(num_nodes >= 2);
		int goal_id = num_nodes - 1;
		int pregoal = goal_id - 1;
		int every = 1;
		if (num_nodes > 100) every = 10;
		assert(in.nodes_[pregoal].out_edges_.size() == 1);
		if (!SILENT) cerr << "    ";
		for (int i = 0; i < in.nodes_.size(); ++i) {
			if (!SILENT && i % every == 0) cerr << '.';
			KBest(i, i == goal_id);
		}
		if (!SILENT) {
			cerr << endl;
			cerr << "  Best path: " << log(D[goal_id].front()->vit_prob_)
          						 << "\t" << log(D[goal_id].front()->est_prob_) << endl;
		}
		out.PruneUnreachable(D[goal_id].front()->node_index_);
		FreeAll();
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
		Hypergraph::Edge* new_edge = out.AddEdge(item->out_edge_);
		new_edge->edge_prob_ = item->out_edge_.edge_prob_;
		Candidate*& o_item = (*s2n)[item->state_];
		if (!o_item) o_item = item;

		int& node_id = o_item->node_index_;
		if (node_id < 0) {
			Hypergraph::Node* new_node = out.AddNode(in.nodes_[item->in_edge_->head_node_].cat_);
			node_states_.push_back(item->state_);
			node_id = new_node->id_;
		}
#if 0
		Hypergraph::Node* node = &out.nodes_[node_id];
		out.ConnectEdgeToHeadNode(new_edge, node);
#else
		out.ConnectEdgeToHeadNode(new_edge, node_id);
#endif
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
		int pop_limit_eff=max(1,int(v.promise*pop_limit_));
		while(!cand.empty() && pops < pop_limit_eff) {
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

	const ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	Hypergraph& out;

	vector<CandidateList> D;   // maps nodes in in-HG to the
	// equivalent nodes (many due to state
	// splits) in the out-HG.
	FFStates node_states_;  // for each node in the out-HG what is
	// its q function value?
	const int pop_limit_;
};

struct NoPruningRescorer {
	NoPruningRescorer(const ModelSet& m, const SentenceMetadata &sm, const Hypergraph& i, Hypergraph* o) :
		models(m),
		smeta(sm),
		in(i),
		out(*o),
		nodemap(i.nodes_.size()) {
		if (!SILENT) cerr << "  Rescoring forest (full intersection)\n";
		node_states_.reserve(kRESERVE_NUM_NODES);
	}

	typedef unordered_map<FFState, int, boost::hash<FFState> > State2NodeIndex;

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
			Hypergraph::Edge* new_edge = out.AddEdge(in_edge, tail);
			FFState head_state;
			if (is_goal) {
				assert(tail.size() == 1);
				const FFState& ant_state = node_states_[tail.front()];
				models.AddFinalFeatures(ant_state, new_edge,smeta);
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
		if (!SILENT) cerr << "    ";
		for (int i = 0; i < in.nodes_.size(); ++i) {
			if (!SILENT && i % every == 0) cerr << '.';
			ProcessOneNode(i, i == goal_id);
		}
		if (!SILENT) cerr << endl;
	}

private:
	const ModelSet& models;
	const SentenceMetadata& smeta;
	const Hypergraph& in;
	Hypergraph& out;

	vector<vector<int> > nodemap;
	FFStates node_states_;  // for each node in the out-HG what is
	// its q function value?
};

// each node in the graph has one of these, it keeps track of
void ApplyModelSet(const Hypergraph& in,
		const SentenceMetadata& smeta,
		/*const*/ ModelSet& models,
		const IntersectionConfiguration& config,
		Hypergraph* out) {
	//force exhaustive if there's no state req. for model
	if (/*{LG TMP} models.stateless() ||*/ config.algorithm == IntersectionConfiguration::FULL) {
		NoPruningRescorer ma(models, smeta, in, out); // avoid overhead of best-first when no state
		ma.Apply();
	} else if (config.algorithm == IntersectionConfiguration::CUBE) {
		int pl = config.pop_limit;
		const int max_pl_for_large=50;
		if (pl > max_pl_for_large && in.nodes_.size() > 80000) {
			pl = max_pl_for_large;
			cerr << "  Note: reducing pop_limit to " << pl << " for very large forest\n";
		}
		CubePruningRescorer ma(models, smeta, in, pl, out);
		ma.Apply();
	} else if (config.algorithm == IntersectionConfiguration::GREEDY_UNDIRECTED) {
		//TODO config.pop_limit;
		GreedyUndirectedRescorer ma(models, smeta, in, false, out);
		ma.Apply();
	} else if (config.algorithm == IntersectionConfiguration::GREEDY_UNDIRECTED_TRAINING) {
		//TODO config.pop_limit;
		GreedyUndirectedRescorer ma(models, smeta, in, true, out, config.weight_file_name_);
		ma.Apply();
	} else {
		cerr << "Don't understand intersection algorithm " << config.algorithm << endl;
		exit(1);
	}
	out->is_linear_chain_ = in.is_linear_chain_;  // TODO remove when this is computed
	// automatically
}

#undef DEBUG_GU
