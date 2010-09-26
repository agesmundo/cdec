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

struct Ngram{

	Ngram(){
		words_Id_ = new vector <WordID>();
	}

	Ngram(vector <WordID>* words_Id_to_Copy, int start_id){
		this();
		for(; start_id<word_Ids.size(); start_id++){
			words_Id_[start_id]=word_Ids_to_Copy[start_id];
		}
	}

	~Ngram(){
		delete word_Ids_;
	}

	int size(){
		return word_Ids_->size();
	}

	void append(WordID wid){
		word_Ids_->insert(wid);
	}

	vector <WordID>* word_Ids_; //ids of words that compose the ngram
};

const int order = 3;//size of ngram is in language model feat.order


const WordID kSTAR(TD::Convert("<{STAR}>"));

typedef set<Ngram> NgramSet;

//TODO check Ngram.h ??
void ComputeNgramSets(const Hypergraph& hg, vector<NgramSet >& edgeToGeneratedNgrams){
	const int num_nodes = hg.nodes_.size();
	edgeToGeneratedNgrams.resize(num_edges); //output vector, map id of edge -> set of ngrams generated there
	vector<Ngram> nodeToBoundariesNgram; //map id of node -> !(set) of boundary ngrams //TODO check that is true that there is only one state per node
	for (int node_index = 0; node_index < num_nodes; ++node_index) { //loop nodes
		const Hypergraph::Node& curr_node = hg.nodes_[node_index];
		const vector<int>& in_edges = curr_node.in_edges_;
		const int num_in_edges = in_edges.size();
		for (int i = 0; i < num_in_edges; ++i) { //loop back-star edges for current node
			const Hypergraph::Edge& curr_edge = hg.edges_[curr_node.in_edges_[j]];
			int tail_size = curr_edge->tail_nodes_.size();

			//get starred sequence of terminals
			vector<WordID> buffer;
			const vector<WordID>& e = curr_edge.rule_.e();
			for (int j = 0; j < e.size(); ++j) {
				if (e[j] < 1) {
					Ngram antState = nodeToBoundariesNgram[edge->tail_nodes_[-e[j]]];
					int slen = StateSize(astate);
					for (int k = 0; k < slen; ++k)
						buffer.insert(astate[k]);
				} else {
					buffer.insert(e[j]);
				}
			}

			//extract state for node and //extract set of generated ngrams
			//see ff_lm.cc 266
			//should do this just for first edge of head node
			//debug checking that all the edges create the same state
			Ngram state;
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
					for (int k=start;k<=j;k++){
						edgeToGeneratedNgrams.add(buffer.substring(k,j+1));
					}
				}

			}

			if (buffer.size()>state.size()){
				state.append(kSTAR);

				for(int w=start; w<buffer.size(); w++){
					state.append(buffer[w]);
				}

			}
			nodeToBoundariesNgram.insert(curr_node, state); //NB do not insert something that is in methond memory


		}
	}
}


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
