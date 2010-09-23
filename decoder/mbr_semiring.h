#ifndef _MBR_SEMIRING_H_
#define _MBR_SEMIRING_H_

#include <iostream>
#include <vector>

#include "hg.h"

//////////////////////////////////////////
////////////////////////////////////////////
// First pass of MBR,
// Note that it is not possible to directly use Inside() because the 'x' operator in not binary

void ComputeNgramSets(const Hypergraph& in, std::vector<std::set<int> >& ngramTailSets , Dict& ngramDict){
    for (int i = 0; i < in.nodes_.size(); ++i) {
      ;
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
