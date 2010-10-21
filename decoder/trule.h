#ifndef _RULE_H_
#define _RULE_H_

#include <algorithm>
#include <vector>
#include <cassert>
#include <boost/shared_ptr.hpp>

#include "sparse_vector.h"
#include "wordid.h"

using namespace std;

class TRule;
typedef boost::shared_ptr<TRule> TRulePtr;

struct NTSizeSummaryStatistics {
  NTSizeSummaryStatistics(int arity) : means(arity), vars(arity) {}
  std::vector<float> means;
  std::vector<float> vars;
};

// Translation rule
class TRule {
 public:
  TRule() : lhs_(0), prev_i(-1), prev_j(-1) { }
  TRule(WordID lhs, const WordID* src, int src_size, const WordID* trg, int trg_size, const int* feat_ids, const double* feat_vals, int feat_size, int arity) :
      e_(trg, trg + trg_size), f_(src, src + src_size), lhs_(lhs), arity_(arity), prev_i(-1), prev_j(-1) {
    for (int i = 0; i < feat_size; ++i)
      scores_.set_value(feat_ids[i], feat_vals[i]);
  }

  explicit TRule(const std::vector<WordID>& e) : e_(e), lhs_(0), prev_i(-1), prev_j(-1) {}
  TRule(const std::vector<WordID>& e, const std::vector<WordID>& f, const WordID& lhs) :
    e_(e), f_(f), lhs_(lhs), prev_i(-1), prev_j(-1) {}

  // deprecated - this will be private soon
  explicit TRule(const std::string& text, bool strict = false, bool mono = false) : prev_i(-1), prev_j(-1) {
    ReadFromString(text, strict, mono);
  }

  // deprecated, use lexer
  // make a rule from a hiero-like rule table, e.g.
  //    [X] ||| [X,1] DE [X,2] ||| [X,2] of the [X,1]
  // if misformatted, returns NULL
  static TRule* CreateRuleSynchronous(const std::string& rule);

  // deprecated, use lexer
  // make a rule from a phrasetable entry (i.e., one that has no LHS type), e.g:
  //    el gato ||| the cat ||| Feature_2=0.34
  static TRule* CreateRulePhrasetable(const std::string& rule);

  // deprecated, use lexer
  // make a rule from a non-synchrnous CFG representation, e.g.:
  //    [LHS] ||| term1 [NT] term2 [OTHER_NT] [YET_ANOTHER_NT]
  static TRule* CreateRuleMonolingual(const std::string& rule);

  static TRule* CreateLexicalRule(const WordID& src, const WordID& trg) {
    return new TRule(src, trg);
  }

  void ESubstitute(const std::vector<const std::vector<WordID>* >& var_values,
                   std::vector<WordID>* result) const {
    int vc = 0;
    result->clear();
    for (std::vector<WordID>::const_iterator i = e_.begin(); i != e_.end(); ++i) {
      const WordID& c = *i;
      if (c < 1) {
        ++vc;
        const std::vector<WordID>& var_value = *var_values[-c];
        std::copy(var_value.begin(),
                  var_value.end(),
                  std::back_inserter(*result));
      } else {
        result->push_back(c);
      }
    }
    assert(vc == var_values.size());
  }

  void FSubstitute(const std::vector<const std::vector<WordID>* >& var_values,
                   std::vector<WordID>* result) const {
    int vc = 0;
    result->clear();
    for (std::vector<WordID>::const_iterator i = f_.begin(); i != f_.end(); ++i) {
      const WordID& c = *i;
      if (c < 1) {
        const std::vector<WordID>& var_value = *var_values[vc++];
        std::copy(var_value.begin(),
                  var_value.end(),
                  std::back_inserter(*result));
      } else {
        result->push_back(c);
      }
    }
    assert(vc == var_values.size());
  }

  bool ReadFromString(const std::string& line, bool strict = false, bool monolingual = false);

  bool Initialized() const { return e_.size(); }

  std::string AsString(bool verbose = true) const;

  static TRule DummyRule() {
    TRule res;
    res.e_.resize(1, 0);
    return res;
  }

  const std::vector<WordID>& f() const { return f_; }
  const std::vector<WordID>& e() const { return e_; }

  int EWords() const { return ELength() - Arity(); }
  int FWords() const { return FLength() - Arity(); }
  int FLength() const { return f_.size(); }
  int ELength() const { return e_.size(); }
  int Arity() const { return arity_; }
  bool IsUnary() const { return (Arity() == 1) && (f_.size() == 1); }
  const SparseVector<double>& GetFeatureValues() const { return scores_; }
  double Score(int i) const { return scores_[i]; }
  WordID GetLHS() const { return lhs_; }
  void ComputeArity();

  // 0 = first variable, -1 = second variable, -2 = third ...
  std::vector<WordID> e_;
  // < 0: *-1 = encoding of category of variable
  std::vector<WordID> f_;
  WordID lhs_;
  SparseVector<double> scores_;
  
  char arity_;
  TRulePtr parent_rule_;  // usually NULL, except when doing constrained decoding

  // this is only used when doing synchronous parsing
  short int prev_i;
  short int prev_j;

  // may be null
  boost::shared_ptr<NTSizeSummaryStatistics> nt_size_summary_;

 private:
  TRule(const WordID& src, const WordID& trg) : e_(1, trg), f_(1, src), lhs_(), arity_(), prev_i(), prev_j() {}
  bool SanityCheck() const;
};

ostream& operator<<(ostream& os, const TRule& rule);

#endif
