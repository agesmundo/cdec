#ifndef _MBR_SEMIRING_H_
#define _MBR_SEMIRING_H_

#include <iostream>
#include <vector>

#include "hg.h"
#include "tdict.h"

using namespace std;

// Define the following macros if you want to see lots of debugging output 
// ... for first pass
#define DEBUG_MBR_1
// ... for second pass
#define DEBUG_MBR_2

struct NGramScoresWeightType;
ostream& operator<<(ostream& os, NGramScoresWeightType& ng);

//////////////////////////////////////////
//////////////////////////////////////////
// First pass of MBR,
// Note that it is not possible to directly use Inside() because the 'x' operator in not binary

struct NGram{

	NGram(){}

	NGram(const vector<WordID>& buffer, const int start, const int end){
		for(int i=start; i<end; i++){
			word_Ids_.push_back(buffer[i]);
		}
	}

	int size(){
		return word_Ids_.size();
	}

	void append(WordID wid){
		word_Ids_.push_back(wid);
	}

	WordID operator[](unsigned int id){
		return  word_Ids_[id];
	}

	bool operator==(NGram other){
		if(size() != other.size()){
			return false;
		}
		for (int i=0; i<size();i++){
			if (word_Ids_[i]!=other[i]){
				return false;
			}
		}
		return true;
	}

	vector <WordID> word_Ids_; //ids of words that compose the ngram
};

#if defined DEBUG_MBR_1 || defined DEBUG_MBR_2
ostream& operator<<(ostream& os, const NGram& ng) {
	os << "NGram[";
	for(int i=0; i < ng.word_Ids_.size();i++){
		if (i!=0){
			os << ", ";
		}
		os << "\"" << TD::Convert(ng.word_Ids_[i])<< "\"";
	}
	return os << ']';
}
#endif 

struct compNGrams{
	bool operator()(const NGram& frst, const NGram& scnd) {
		if(frst.word_Ids_.size() < scnd.word_Ids_.size()){
			return true;
		}
		if(frst.word_Ids_.size() > scnd.word_Ids_.size()){
			return false;
		}
		for (int i=0; i<frst.word_Ids_.size();i++){
			if (frst.word_Ids_[i]<scnd.word_Ids_[i]){
				return true;
			}

			//TODO why??
			//removing ".word_Ids" generate error
			//passing «const NGram» as «this» argument of «WordID NGram::operator[](unsigned int)» discards qualifiers
			if(frst.word_Ids_[i]>scnd.word_Ids_[i]){
				return false;
			}
		}
		return false;
	}
};

const WordID kSTAR(TD::Convert("<{STAR}>"));

typedef set<NGram,compNGrams> NGramSet;

//TODO check Ngram.h ??
void ComputeNgramSets(const Hypergraph& hg, vector<NGramSet >& edgeToGeneratedNgrams, const int order){
	const int num_edges = hg.edges_.size();
	const int num_nodes = hg.nodes_.size();
	const int goal_id = num_nodes - 1;

	//output vector, map id of edge -> set of ngrams generated there
	edgeToGeneratedNgrams.resize(num_edges);

	//map id of node -> !(set) of boundary ngrams //TODO check that is true that there is only one state per node
	vector<NGram> nodeToBoundariesNgram;
	nodeToBoundariesNgram.resize(num_nodes);

	//loop nodes
	for (int node_index = 0; node_index < num_nodes; node_index++ ) {

#ifdef DEBUG_MBR_1
		cerr << "AT NODE: " << node_index << endl;
#endif

		const Hypergraph::Node& curr_node = hg.nodes_[node_index];
		const int num_in_edges = curr_node.in_edges_.size();
		NGram state;
		bool compute_state = (node_index != goal_id);

		//loop back-star edges for current node
		for (int tail_edge_index = 0; tail_edge_index < num_in_edges; ++tail_edge_index) {

			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[tail_edge_index]];
			
			//get starred sequence of terminals
			vector<WordID> buffer;
			const vector<WordID>& e = curr_edge.rule_->e_;

			//set of index of rule-terminal elements in buffer
			//needed to distinguish them from the child-state-terminals
			//since a generated NGram must contain at least one rule-terminal
			set<int> rule_terminal_ids;

#ifdef DEBUG_MBR_1
			cerr << "AT EDGE: " << tail_edge_index << " , " << curr_edge.id_ <<endl;
			cerr << "\tRule applied: ";
#endif

			for (int j = 0; j < e.size(); ++j) {

				//if target side of rule is <1, it is the index of the |e|-th non-terminal in the source side
				if (e[j] < 1) {

					//get antecedent state
					NGram antState = nodeToBoundariesNgram[curr_edge.tail_nodes_[-e[j]]];

#ifdef DEBUG_MBR_1
					cerr<< "[X," << (-e[j]) +1 <<": "; 
#endif

					//append antecendent State to buffer
					for (int k = 0; k < antState.size(); ++k){
						buffer.push_back(antState[k]);

#ifdef DEBUG_MBR_1		
						cerr<< "\"" << TD::Convert(antState[k])<< "\"" << " "; 
#endif

					}

#ifdef DEBUG_MBR_1		
					cerr<< "]"; 
#endif

				}

				//if target side of rule is >0, it is an index of word
				else {
					assert( (rule_terminal_ids.insert(buffer.size())).second );	
					buffer.push_back(e[j]);

#ifdef DEBUG_MBR_1		
					cerr<< "\"" << TD::Convert(e[j])<< "\"" << " ";
#endif

				}

			}

#ifdef DEBUG_MBR_1	
			cerr << endl;
#endif

			//extract state for node and extract set of generated ngrams
			//XXX see ff_lm.cc 266

			int start=0;
			for (int j=0; j< buffer.size();j++){
				if (buffer[j]==kSTAR){
					start=j+1;
				}
				else {
					if(j-start+1>=order){ //length span reached order
						start=j-order+1;
					}
					else if (start==0 && compute_state){
						state.append(buffer[j]);
					}

					//add all the ngrams that end at j
					NGramSet& currentNgramSet = edgeToGeneratedNgrams[curr_edge.id_];

#ifdef DEBUG_MBR_1			
					cerr<< "\tExtracting all n-grams generated at EDGE " << node_index << " , " << curr_edge.id_ <<" , ending in " << j << " , starting from "<< start <<endl ;
#endif

					bool pass =false;
					for (int k=j ;k>=start;k--){
						if(pass || ( rule_terminal_ids.find(k) != rule_terminal_ids.end()) ){
							NGram tmp = NGram(buffer,k,j+1);

#ifdef DEBUG_MBR_1	
							cerr<< "\t\t" << tmp << endl;
#endif

							currentNgramSet.insert(tmp);
						}
					}
				}

			}

			if(compute_state){

				if (buffer.size()>state.size()){
					state.append(kSTAR);

					if (buffer.size()- order +1 > start) start= buffer.size()- order+1;

					for(int w=start; w<buffer.size(); w++){
						state.append(buffer[w]);
					}

				}

				//compute state once
//#ifndef DEBUG_MBR_1 

				assert(nodeToBoundariesNgram[node_index].word_Ids_.size()==0); 
				nodeToBoundariesNgram[node_index] =state;
				compute_state=false;	
//#endif

//				//compute states for all incoming endges and check if matching
//#ifdef DEBUG_MBR_1 
//
//				cerr<< "\tState for NODE "<< node_index << " : " << state<< endl;
//
//				if(nodeToBoundariesNgram[node_index].word_Ids_.size()==0){
//					nodeToBoundariesNgram[node_index] =state;
//					cerr<< "INSERTED" << endl;
//				}
//				else{
//					assert(nodeToBoundariesNgram[node_index]==state);
//					cerr<< "MATCHING" << endl;
//				}
//				state=NGram();
//#endif



			}

		}

	}

#ifdef DEBUG_MBR_1

	//display states computed
	for (int node_index = 0; node_index < num_nodes; node_index++) {

		cerr << "NODE:  " <<node_index << " : " << nodeToBoundariesNgram[node_index] << endl ;

	}

	//display output, edge -> set of NGram
	for (int node_index = 0; node_index < num_nodes; node_index++) {
		const Hypergraph::Node& curr_node = hg.nodes_[node_index];
		for (int tail_edge_index =0 ; tail_edge_index <curr_node.in_edges_.size(); tail_edge_index++){
			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[tail_edge_index]];
			cerr << "EDGE:  " << curr_edge.id_ << " , Set of NGrams : { " << endl;
			const NGramSet& currentNgramSet = edgeToGeneratedNgrams[curr_edge.id_];
			for(NGramSet::iterator it = currentNgramSet.begin(); it != currentNgramSet.end(); it++){
				cerr << "\t " << *it << endl ;
			}
			cerr << "};" << endl;
		}
	}

#endif

}


//////////////////////////////////////////
//////////////////////////////////////////
// second pass
//


//WeightType for second pass
//map any possible ngram to score
typedef map<NGram, prob_t, compNGrams>  MapNGramScore;
struct NGramScoresWeightType{

	NGramScoresWeightType(){
		NGramScoresWeightType(prob_t());
	}

	NGramScoresWeightType(int x){
		NGramScoresWeightType(prob_t(x));
	}

	NGramScoresWeightType(prob_t x){
		tot_score_=x;
		default_score_=x;
	}

	static void AddFunction(const Hypergraph::Node& node,
			vector <NGramScoresWeightType>& child_edges_weights,
			NGramScoresWeightType* const cur_node_inside_score){

#ifdef DEBUG_MBR_2
		assert(node.in_edges_.size()==child_edges_weights.size());
#endif

		//*cur_node_inside_score =NGramScoresWeightType();
		for (int i =0; i<node.in_edges_.size(); i++){
			*cur_node_inside_score += child_edges_weights[i];
		}
	}

	NGramScoresWeightType& operator+=(const NGramScoresWeightType& o) {
		for(MapNGramScore::const_iterator it= o.ngram_score_.begin();it != o.ngram_score_.end(); it++){
			pair<MapNGramScore::iterator, bool> pointer = ngram_score_.insert(*it);
			if (!pointer.second){
				(*(pointer.first)).second += (*it).second;
			}
		}
		tot_score_ +=o.tot_score_;
		default_score_+=o.default_score_;
		return *this;
	}

	static void MultFunction(const Hypergraph::Edge& edge,
			vector <NGramScoresWeightType>& child_nodes_weights,
			NGramScoresWeightType* const cur_edge_inside_score){

#ifdef DEBUG_MBR_2
		assert(edge.tail_nodes_.size() == child_nodes_weights.size());
#endif

		NGramScoresWeightType& edge_scores = (*cur_edge_inside_score);
		NGramSet ngram_set = (*edge_to_ngram_set_)[edge.id_];
		int tail_size = child_nodes_weights.size() ; 

		///compute tot score
		//this could be done with regular inside
		edge_scores.tot_score_ = edge.edge_prob_;
		for(int i=0; i<tail_size; i++){
			edge_scores.tot_score_ *= child_nodes_weights[i].tot_score_;
		}

		///assign tot_score to ngram generated at current edge
		for (NGramSet::const_iterator it = ngram_set.begin(); it != ngram_set.end(); it++){
			assert( edge_scores.ngram_score_.insert(pair<NGram, prob_t>(*it, edge_scores.tot_score_)).second);
		}
		
#ifdef DEBUG_MBR_2
		cerr <<"-  base with gen ng :\n" << *cur_edge_inside_score << endl ;
#endif
		
		///assign score_of_paths_containing_ngram to ngrams generated in the yield
		//NB:
		//1) score_of_paths_containing_ngram = tot_score - score_of_paths_NOT_containing_ngram
		//2) score_of_paths_NOT_containing_ngram = rule_score * \MULT_{C \in child_nodes} score_of_paths_at_C_NOT_containing_ngram
		//3) score_of_paths_at_C_NOT_containing_ngram = tot_score_at_C - score_of_paths_at_C_containing_ngram

		NGramScoresWeightType tmp = NGramScoresWeightType(1);
		for(int i=0; i<tail_size; i++){
			NGramScoresWeightType child_update = NGramScoresWeightType(1);
			const NGramScoresWeightType child_weights= child_nodes_weights[i];

			for(MapNGramScore::const_iterator it = child_weights.ngram_score_.begin(); it != child_weights.ngram_score_.end(); it++){
				const NGram& key = (*it).first; //TODO try to test if ngram is already in set of gen ng
				prob_t child_notNg = child_nodes_weights[i].tot_score_ - child_nodes_weights[i].getScore(key); // <- 3)
				assert(child_update.insert(key, child_notNg));

			}
			
#ifdef DEBUG_MBR_2
			cerr <<"- update from child " <<  i << "  :\n" << child_update << endl ;
#endif 
			
			tmp *= child_update; // <- 2)a   ... MULT_{C \in child_nodes}
		}

		for(MapNGramScore::iterator it = tmp.ngram_score_.begin(); it != tmp.ngram_score_.end(); it++){
			(*it).second *= edge.edge_prob_; // <- 2)b rule_score * ...
		}
		
#ifdef DEBUG_MBR_2
		cerr <<"- mult all with rule score :\n" << tmp << endl ;
#endif
		
		for(MapNGramScore::const_iterator it = tmp.ngram_score_.begin(); it != tmp.ngram_score_.end(); it++){
			prob_t paths_with_ng = edge_scores.tot_score_ - (*it).second; // <- 3)
			edge_scores.insert((*it).first, paths_with_ng); //NB for ngrams generated at curr_edge there is already the value set and the insertion will fail
		}

	}

	bool insert (const NGram& ng, prob_t p){
		return ngram_score_.insert(pair <NGram, prob_t> (ng, p )).second;
	}

	prob_t getScore(NGram ng){
		MapNGramScore::const_iterator found = ngram_score_.find(ng);
		if (found == ngram_score_.end())
			return default_score_;
		else
			return found->second;
	}

	NGramScoresWeightType& operator*=(const NGramScoresWeightType& o) {
		for(MapNGramScore::const_iterator it= o.ngram_score_.begin();it != o.ngram_score_.end(); it++){
			pair<MapNGramScore::iterator, bool> pointer = ngram_score_.insert(*it);
			if (!pointer.second){
				(*(pointer.first)).second *= (*it).second;
			}
		}
		tot_score_ *=o.tot_score_;
		default_score_ *=o.default_score_;
		return *this;
	}

	prob_t tot_score_; //score all derivations
	MapNGramScore ngram_score_; //score of derivations generating NGram
	prob_t default_score_; //score for all NGram not explicitly cited in the map

	static vector< NGramSet >* edge_to_ngram_set_;	
};
vector< NGramSet >* NGramScoresWeightType::edge_to_ngram_set_;

ostream& operator<<(ostream& os, NGramScoresWeightType& ng) {
	os << "NGramScoresWeightType[\n";
	os << "\ttot_score_ : " << double(ng.tot_score_) << endl;
	os << "\tdefault_score_ : " << double(ng.default_score_) << endl ;
	os << "\tngram_score_ |" << ng.ngram_score_.size() << "| : {" ;
	for(MapNGramScore::const_iterator it= ng.ngram_score_.begin();it != ng.ngram_score_.end(); it++){		
		if (it!= ng.ngram_score_.begin()){
			os << ", ";
		}
		os << "(" << (*it).first << " | " << double((*it).second) << ")";
	}
	return os << "\t}\n]\n";
}


template<typename WeightType>
void GeneralizedInside(const Hypergraph& hg,
		std::vector<WeightType>* result = NULL){
	//const Hypergraph& hg, vector<NGramSet >& edgeToGeneratedNgrams, map <NGram,prob_t> ngramToPosterior)

	const int num_nodes = hg.nodes_.size();

	std::vector<WeightType> dummy;
	std::vector<WeightType>& inside_score =  result ? *result : dummy;
	inside_score.resize(num_nodes);
	std::fill(inside_score.begin(), inside_score.end(), WeightType(0));

	for (int i = 0; i < num_nodes; ++i) {

		const Hypergraph::Node& curr_node = hg.nodes_[i];
		WeightType* const cur_node_inside_score = &inside_score[i];
		const int num_in_edges = curr_node.in_edges_.size();

		if (num_in_edges == 0) {
			(*cur_node_inside_score)=WeightType(1);

			cerr << "XXXNODE : " << i << " , weights : \n" <<  *cur_node_inside_score << endl ;
			assert (false);
			continue;
		}

		vector <NGramScoresWeightType> child_edges_weights;
		child_edges_weights.resize(num_in_edges);
		std::fill(child_edges_weights.begin(), child_edges_weights.end(), WeightType());

		for (int j = 0; j < num_in_edges; ++j) {

			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[j]];
			int tail_size = curr_edge.tail_nodes_.size();

			vector <NGramScoresWeightType> child_nodes_weights;
			child_nodes_weights.resize(tail_size);

			for (int k = 0; k < tail_size; ++k) {
				child_nodes_weights[k] = inside_score[curr_edge.tail_nodes_[k]];
			}

			WeightType::MultFunction(curr_edge, child_nodes_weights, &child_edges_weights[j] );
#ifdef DEBUG_MBR_2
			cerr << "EDGE : " << curr_edge.id_ << " , rule score : " << curr_edge.edge_prob_ << " , tail:";
			for (int w=0 ; w<curr_edge.tail_nodes_.size();w++) cerr <<curr_edge.tail_nodes_[w]<< " ";
			cerr <<", weights : \n" <<  child_edges_weights[j] ;
			const NGramSet& currentNgramSet = (*NGramScoresWeightType::edge_to_ngram_set_)[curr_edge.id_];
			cerr <<"Set of NGrams |" <<currentNgramSet.size() << "| : { " ;
			for(NGramSet::iterator it = currentNgramSet.begin(); it != currentNgramSet.end(); it++){
				cerr << *it << " " ;
			}
			cerr << "};\n" << endl;
#endif
		}

		WeightType::AddFunction(curr_node, child_edges_weights, cur_node_inside_score);
#ifdef DEBUG_MBR_2
		cerr << "NODE : " << i << " , tail:";
		for (int w=0 ; w<curr_node.in_edges_.size();w++) cerr <<curr_node.in_edges_[w]<< " ";
		cerr << " , weights : \n" <<  *cur_node_inside_score << endl ;
#endif
	}

}





//WeightFunction for second pass
//struct MBR2WeightFunction {
//	MBR2WeightFunction(){
//		
//	}
//	
//	inline NGramScoresWeightType operator()(const Hypergraph::Edge& e) const {
//		//for each ngram
//		return e.edge_prob_;
//	}
//	
//	//shares the vector filled with results by Inside alg
//	//neede to know the set of ngrams at each node
//	vector<WeightType>& result_; 
//};

/*
// this file implements the first-order expectation semiring described
// in Li & Eisner (EMNLP 2009)

// requirements:
//   RType * RType ==> RType
//   PType * PType ==> PType
//   RType * PType ==> RType
// good examples:
//   PType scalar, RType vector
// BAD examples:
//   PType vector, RType scalar
template <typename PType, typename RType>
struct PRPair {
	PRPair() : p(), r() {}
	// Inside algorithm requires that T(0) and T(1)
	// return the 0 and 1 values of the semiring
	explicit PRPair(double x) : p(x), r() {}
	PRPair(const PType& p, const RType& r) : p(p), r(r) {}
	PRPair& operator+=(const PRPair& o) {
		p += o.p;
		r += o.r;
		return *this;
	}
	PRPair& operator*=(const PRPair& o) {
		r = (o.r * p) + (o.p * r);
		p *= o.p;
		return *this;
	}
	PType p;
	RType r;
};

template <typename P, typename R>
std::ostream& operator<<(std::ostream& o, const PRPair<P,R>& x) {
	return o << '<' << x.p << ", " << x.r << '>';
}

template <typename P, typename R>
const PRPair<P,R> operator+(const PRPair<P,R>& a, const PRPair<P,R>& b) {
	PRPair<P,R> result = a;
	result += b;
	return result;
}

template <typename P, typename R>
const PRPair<P,R> operator*(const PRPair<P,R>& a, const PRPair<P,R>& b) {
	PRPair<P,R> result = a;
	result *= b;
	return result;
}

template <typename P, typename PWeightFunction, typename R, typename RWeightFunction>
struct PRWeightFunction {
	explicit PRWeightFunction(const PWeightFunction& pwf = PWeightFunction(),
			const RWeightFunction& rwf = RWeightFunction()) :
		pweight(pwf), rweight(rwf) {}
	PRPair<P,R> operator()(const Hypergraph::Edge& e) const {
		const P p = pweight(e);
		const R r = rweight(e);
		return PRPair<P,R>(p, r * p);
	}
	const PWeightFunction pweight;
	const RWeightFunction rweight;
};
 */

#undef DEBUG_MBR_1
#undef DEBUG_MBR_2

#endif
