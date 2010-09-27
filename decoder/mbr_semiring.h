#ifndef _MBR_SEMIRING_H_
#define _MBR_SEMIRING_H_

#include <iostream>
#include <vector>

#include "hg.h"
#include "tdict.h"

using namespace std;


//////////////////////////////////////////
////////////////////////////////////////////
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

const int order = 3;//size of ngram is in language model feat.order


const WordID kSTAR(TD::Convert("<{STAR}>"));

typedef set<NGram,compNGrams> NGramSet;

//TODO check Ngram.h ??
void ComputeNgramSets(const Hypergraph& hg, vector<NGramSet >& edgeToGeneratedNgrams){
	const int num_nodes = hg.nodes_.size();
	const int num_edges = hg.edges_.size();

	//output vector, map id of edge -> set of ngrams generated there
	edgeToGeneratedNgrams.resize(num_edges);

	//map id of node -> !(set) of boundary ngrams //TODO check that is true that there is only one state per node
	vector<NGram> nodeToBoundariesNgram;

	//loop nodes
	for (int node_index = 0; node_index < num_nodes; ++node_index) {

		const Hypergraph::Node& curr_node = hg.nodes_[node_index];
		const int num_in_edges = curr_node.in_edges_.size();

		//DEBUG
		NGram DBGstate= NGram();
		//

		//loop back-star edges for current node
		for (int tail_edge_index = 0; tail_edge_index < num_in_edges; ++tail_edge_index) {

			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[tail_edge_index]];
			//			const int tail_size = curr_edge.tail_nodes_.size();

			//get starred sequence of terminals
			vector<WordID> buffer;
			const vector<WordID>& e = curr_edge.rule_->e_;
			for (int j = 0; j < e.size(); ++j) {

				//if target side of rule is <1, it is the index of the |e|-th non-terminal in the source side
				if (e[j] < 1) {

					//get antecedent state
					NGram antState = nodeToBoundariesNgram[curr_edge.tail_nodes_[-e[j]]];

					//append antecendent State to buffer
					for (int k = 0; k < antState.size(); ++k){
						buffer.push_back(antState[k]);
					}

				}

				//if target side of rule is >0, it is an index of word
				else {
					buffer.push_back(e[j]);
				}

			}


			//extract state for node and extract set of generated ngrams
			//XXX see ff_lm.cc 266
			//should do this just for first edge of head node
			//debug checking that all the edges create the same state


			NGram state;
			int start=0;
			for (int j=0; j< buffer.size();j++){
				if (buffer[j]==kSTAR){
					start=j+1;
				}
				else {
					if(j-start+1>=order){ //length span reached order
						start=j-order;
					}
					else if (start==0){
						state.append(buffer[j]);
					}
					//TODO add all the ngrams that end at j
					//					NGramSet& currentNgramSet = edgeToGeneratedNgrams[curr_edge.id_];
					NGramSet currentNgramSet;
					for (int k=start;k<=j;k++){
						NGram tmp = NGram(buffer,k,j+1);
						currentNgramSet.insert(tmp);
					}
				}

			}

			if (buffer.size()>state.size()){
				state.append(kSTAR);

				for(int w=start; w<buffer.size(); w++){
					state.append(buffer[w]);
				}

			}
			nodeToBoundariesNgram[curr_node.id_] =state;

			//DEBUG
			if(DBGstate.size()==0){
				DBGstate=state;
			}
			else{
				assert(DBGstate==state);
			}
			//


		}

	}

}

//void ComputeNgramPosteriors(const Hypergraph& hg, vector<NgramSet >& edgeToGeneratedNgrams, map <Ngram,prob_t> ngramToPosterior){
//	const int num_nodes = hg.nodes_.size();
//	for (int node_index = 0; node_index < num_nodes; ++node_index) {
//		const Hypergraph::Node& curr_node = hg.nodes_[node_index];
//		const vector<int>& in_edges = curr_node.in_edges_;
//		const int num_in_edges = in_edges.size();
//		for (int i = 0; i < num_in_edges; ++i) {
//			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[i]];
//			int tail_size = curr_edge->tail_nodes_.size();
//			prob_t w = curr_edge.edge_prob_;
//		}
//	}
//}


//////////////////////////////////////////
////////////////////////////////////////////
//TODO not true, reuse for second pass
//
//NB first pass cannot be computed as a semiring, computation is not local
//the 'x' is not a 
//WeightType for first MBR pass
/*
struct VectorSetIdWeightType{

	VectorSetIdWeightType()	{
		idSetsPtr_ = new std::set<int> ;
	}

	VectorSetIdWeightType(int x){

		//if id =0 init as identity for '+' operation
		//'+' operation is the union of the sets
		//therefore return empty list 

		//if x=1 init as identity for 'x' operation
		//'x' operation is the function that generates the ngrams from the boudary words and add them to the sets
		//therefore return empty list
		VectorSetIdWeightType();
	}

	VectorSetIdWeightType(const Hypergraph::Edge& e){
		VectorSetIdWeightType();
	}

	const VectorSetIdWeightType operator*=(const VectorSetIdWeightType& a) {
		//TODO implement here function that find
		return VectorSetIdWeightType();
	}

	const VectorSetIdWeightType operator+=(const VectorSetIdWeightType& a) {
		//TODO implement here union of two sets for each entry 
		//if one of two is empty vector return list
		return VectorSetIdWeightType();
	}

	std::set<int>* idSetsPtr_; //one set of ngram-Ids for each edge child
};*/

//WeightFunction for first MBR pass
//struct MBR1WeightFunction {
//	inline VectorSetIdWeightType operator()(const Hypergraph::Edge& e) const {
//		return VectorSetIdWeightType(e);
//	}
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
#endif
