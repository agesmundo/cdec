#ifndef _EXTRACT_H_
#define _EXTRACT_H_

#include <utility>
#include <vector>
#include "array2d.h"
#include "wordid.h"

struct AnnotatedParallelSentence;

// usually represents a consistent phrase, which may
// be annotated with a type (cat)
struct ParallelSpan {
  // i1 = i of f side
  // i2 = j of f side
  // j1 = i of e side
  // j2 = j of e side
  short i1,i2,j1,j2;
  // cat is set by AnnotatePhrasesWithCategoryTypes, otherwise it's 0
  WordID cat;  // category type of span (also overloaded by RuleItem class
               //                        to be a word ID)
  ParallelSpan() : i1(-1), i2(-1), j1(-1), j2(-1), cat() {}
  // used by Rule class to represent a terminal symbol:
  explicit ParallelSpan(WordID w) : i1(-1), i2(-1), j1(-1), j2(-1), cat(w) {}
  ParallelSpan(int pi1, int pi2, int pj1, int pj2) : i1(pi1), i2(pi2), j1(pj1), j2(pj2), cat() {}
  ParallelSpan(int pi1, int pi2, int pj1, int pj2, WordID c) : i1(pi1), i2(pi2), j1(pj1), j2(pj2), cat(c) {}

  // ParallelSpan is used in the Rule class where it is
  // overloaded to also represent terminal symbols
  inline bool IsVariable() const { return i1 != -1; }
};

struct Extract {
  struct RuleObserver {
    virtual void CountRule(WordID lhs,
                           const std::vector<WordID>& rhs_f,
                           const std::vector<WordID>& rhs_e,
                           const std::vector<std::pair<int, int> >& fe_terminal_alignments) = 0;
  };

  static void LoosenPhraseBounds(const AnnotatedParallelSentence& sentence,
                                 std::vector<ParallelSpan>* phrases);

  static void ExtractBasePhrases(const int max_base_phrase_size,
                        const AnnotatedParallelSentence& sentence,
                        std::vector<ParallelSpan>* phrases);

  // this uses the TARGET span (i,j) to annotate phrases, will copy
  // phrases if there is more than one annotation.
  // TODO: support source annotation
  static void AnnotatePhrasesWithCategoryTypes(const WordID default_cat,
                                      const Array2D<std::vector<WordID> >& types,
                                      std::vector<ParallelSpan>* phrases);

  static void ExtractConsistentRules(const AnnotatedParallelSentence& sentence,
                          const std::vector<ParallelSpan>& phrases,
                          const int max_vars,
                          const int max_syms,
                          const bool permit_adjacent_nonterminals,
                          RuleObserver* observer);
};

#endif
