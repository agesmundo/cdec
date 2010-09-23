#ifndef _MBR_SEMIRING_H_
#define _MBR_SEMIRING_H_

#include <iostream>
#include <vector>

#include "hg.h"

using namespace std;


// run the inside algorithm and return the inside score
// if result is non-NULL, result will contain the inside
// score for each node
// NOTE: WeightType()  must construct the semiring's additive identity
//       WeightType(1) must construct the semiring's multiplicative identity
//TODO test redoing one Inside call with this generalized version
template<typename WeightType, typename WeightFunctionMult, typename WeightFunctionAdd >
WeightType GenInside(const Hypergraph& hg,
		std::vector<WeightType>* result = NULL,
		const WeightFunctionMult& weight = WeightFunctionMult()) {
	const int num_nodes = hg.nodes_.size();
	std::vector<WeightType> dummy;
	std::vector<WeightType>& inside_score_nodes = result ? *result : dummy;
	inside_score_nodes.resize(num_nodes);
	std::fill(inside_score_nodes.begin(), inside_score_nodes.end(), WeightType());
	for (int i = 0; i < num_nodes; ++i) {
		const Hypergraph::Node& cur_node = hg.nodes_[i];
		WeightType* const cur_node_inside_score = &inside_score_nodes[i];
		const int num_in_edges = cur_node.in_edges_.size();
		if (num_in_edges == 0) {
			*cur_node_inside_score = WeightType(1);
			continue;
		}
		for (int j = 0; j < num_in_edges; ++j) {
			const Hypergraph::Edge& edge = hg.edges_[cur_node.in_edges_[j]];
			WeightType score = weight(edge);
			for (int k = 0; k < edge.tail_nodes_.size(); ++k) {
				const int tail_node_index = edge.tail_nodes_[k];
				score *= inside_score_nodes[tail_node_index];
			}
			*cur_node_inside_score += score;
		}
	}
	return inside_score_nodes.back();
}

//////////////////////////////////////////
////////////////////////////////////////////
// First pass of MBR,
// Note that it is not possible to directly use Inside() because the 'x' operator in not binary

struct Ngram{
	
	Ngram(vector <int>* words_Ids){
		words_Ids_=words_Ids;
	}
	
	~Ngram(){
		delete words_Ids_;
	}
	
	vector <int>* words_Ids_; //ids of words that compose the ngram
};

const int ngram_size = 3;

void ComputeNgramSets(const Hypergraph& hg, vector<set<Ngram> >& edgeToGeneratedNgrams){
	const int num_nodes = hg.nodes_.size();
	const int num_edges = hg.edges_.size();
	edgeToGeneratedNgrams.resize(num_edges); //output vector, map id of edge -> set of ngrams generated there
	vector<set<pair<Ngram,Ngram> > > nodeToBoundariesNgram; //map id of node -> set of boundary ngrams
	for (int node_index = 0; node_index < num_nodes; ++node_index) { //loop nodes
		const Hypergraph::Node& cur_node = hg.nodes_[node_index];
		const vector<int>& in_edges = cur_node.in_edges_;
		const int num_in_edges = in_edges.size();
		for (int j = 0; j < num_in_edges; ++j) { //loop back-star edges for current node
			
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
