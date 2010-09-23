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

	Ngram(vector <WordID>* words_Id_to_Copy){
		this();
		for(int i=0; i<word_Ids.size(); i++){
			words_Id_[i]=word_Ids_to_Copy[i];
		}
	}
	
	~Ngram(){
		delete word_Ids_;
	}
	
	int size(){ //implemented depth first to reduce memory usage
		return word_Ids_->size();

	}

	void append(WordID wid){
		word_Ids_->insert(wid);
	}

	vector <WordID>* word_Ids_; //ids of words that compose the ngram
};

const int ngram_max_size = 3;
const WordID kSTAR(TD::Convert("<{STAR}>"));

typedef set<Ngram> NgramSet;

void ComputeNgramSets(const Hypergraph& hg, vector<NgramSet >& edgeToGeneratedNgrams){
	const int num_nodes = hg.nodes_.size();
	const int num_edges = hg.edges_.size();
	edgeToGeneratedNgrams.resize(num_edges); //output vector, map id of edge -> set of ngrams generated there
	vector<NgramSet> nodeToBoundariesNgram; //map id of node -> set of boundary ngrams
	for (int node_index = 0; node_index < num_nodes; ++node_index) { //loop nodes
		const Hypergraph::Node& cur_node = hg.nodes_[node_index];
		const vector<int>& in_edges = cur_node.in_edges_;
		const int num_in_edges = in_edges.size();
		for (int j = 0; j < num_in_edges; ++j) { //loop back-star edges for current node
			const Hypergraph::Edge& curr_edge = hg.edges_[cur_node.in_edges_[j]];
			int tail_size = curr_edge->tail_nodes_.size();
			vector<NgramSet> antsToBoundaryNgrams(tail_size);
			for (int k=0;k<tail_size;k++){
				antsToBoundaryNgrams[k]= nodeToBoundariesNgram[curr_edge->tail_nodes_[k]];
			}
			getSetOfGeneratedNgrams(curr_edge,antsToBoundaryNgrams, edgeToGeneratedNgrams[curr_edge.id_]);

			NgramSet toAdd=	getBoundariesNgrams(curr_edge,antsToBoundaryNgrams);
			nodeToBoundariesNgram[cur_node.id_].insert(toAdd.begin(),toAdd.end());
		}
	}
}

NgramSet getBoundariesNgrams(const Hypergraph::Edge& edge, const vector<NgramSet>& antsToBoundaryNgrams){
	return NULL;
}


struct NgramGenerator{

	NgramGenerator(){
		eId_(0);
		candidateBoundariesIterator_=NULL;
		ng_();
	}

	vector<NgramGenerator> step(NgramSet& ngSet){

	    for (int j = 0; j < e_.size(); ++j) {
	      if (e_[j] < 1) {
	        const int* astate = reinterpret_cast<const int*>(ant_states[-e_[j]]);
	        int slen = StateSize(astate);
	        for (int k = 0; k < slen; ++k)
	          buffer_[i--] = astate[k];
	      } else {
	    	  assert(ng_.size<ngram_max_size);
	    	  ng_.append(e_[j]);
	      }
	    }

	}

	static const vector<WordID>& e_;
	static const vector<NgramSet>& antsToBoundaryNgrams_;
	int eId_;
	Ngram ng_;
	set::iterator candidateBoundariesIterator_;//TODO pointer ref?? what is returned from set?

};


void getSetOfGeneratedNgrams(const Hypergraph::Edge& edge, const vector<NgramSet>& antsToBoundaryNgrams, NgramSet& ngSet){

	stack<NgramGenerator> stack;
	NgramGenerator::e_ = rule.e();
	NgramGenerator::antsToBoundaryNgrams_ = antsToBoundaryNgrams;
	stak.push(NgramGenerator());
	while(stack.size()!=0){
		NgramGenerator currGen = stack.pop();
		stack.pushAll(currGen.step(ngSet));
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
