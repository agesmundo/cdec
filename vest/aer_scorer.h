#ifndef _AER_SCORER_
#define _AER_SCORER_

#include <boost/shared_ptr.hpp>

#include "scorer.h"
#include "array2d.h"

class AERScorer : public SentenceScorer {
 public:
  AERScorer(const std::vector<std::vector<WordID> >& refs);
  Score* ScoreCandidate(const std::vector<WordID>& hyp) const;
  static Score* ScoreFromString(const std::string& in);
 private:
  boost::shared_ptr<Array2D<bool> > ref_;
};

#endif
