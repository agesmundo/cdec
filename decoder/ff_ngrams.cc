#include "ff_ngrams.h"

#include <cstring>
#include <iostream>

#include <boost/scoped_ptr.hpp>

#include "filelib.h"
#include "stringlib.h"
#include "hg.h"
#include "tdict.h"

using namespace std;

static const unsigned char HAS_FULL_CONTEXT = 1;
static const unsigned char HAS_EOS_ON_RIGHT = 2;
static const unsigned char MASK             = 7;

namespace {
template <unsigned MAX_ORDER = 5>
struct State {
  explicit State() {
    memset(state, 0, sizeof(state));
  }
  explicit State(int order) {
    memset(state, 0, (order - 1) * sizeof(WordID));
  }
  State<MAX_ORDER>(char order, const WordID* mem) {
    memcpy(state, mem, (order - 1) * sizeof(WordID));
  }
  State(const State<MAX_ORDER>& other) {
    memcpy(state, other.state, sizeof(state));
  }
  const State& operator=(const State<MAX_ORDER>& other) {
    memcpy(state, other.state, sizeof(state));
  }
  explicit State(const State<MAX_ORDER>& other, unsigned order, WordID extend) {
    char om1 = order - 1;
    assert(om1 > 0);
    for (char i = 1; i < om1; ++i) state[i - 1]= other.state[i];
    state[om1 - 1] = extend;
  }
  const WordID& operator[](size_t i) const { return state[i]; }
  WordID& operator[](size_t i) { return state[i]; }
  WordID state[MAX_ORDER];
};
}

class NgramDetectorImpl {

  // returns the number of unscored words at the left edge of a span
  inline int UnscoredSize(const void* state) const {
    return *(static_cast<const char*>(state) + unscored_size_offset_);
  }

  inline void SetUnscoredSize(int size, void* state) const {
    *(static_cast<char*>(state) + unscored_size_offset_) = size;
  }

  inline State<5> RemnantLMState(const void* cstate) const {
    return State<5>(order_, static_cast<const WordID*>(cstate));
  }

  inline const State<5> BeginSentenceState() const {
    State<5> state(order_);
    state.state[0] = kSOS_;
    return state;
  }

  inline void SetRemnantLMState(const State<5>& lmstate, void* state) const {
    // if we were clever, we could use the memory pointed to by state to do all
    // the work, avoiding this copy
    memcpy(state, lmstate.state, (order_-1) * sizeof(WordID));
  }

  WordID IthUnscoredWord(int i, const void* state) const {
    const WordID* const mem = reinterpret_cast<const WordID*>(static_cast<const char*>(state) + unscored_words_offset_);
    return mem[i];
  }

  void SetIthUnscoredWord(int i, const WordID index, void *state) const {
    WordID* mem = reinterpret_cast<WordID*>(static_cast<char*>(state) + unscored_words_offset_);
    mem[i] = index;
  }

  inline bool GetFlag(const void *state, unsigned char flag) const {
    return (*(static_cast<const char*>(state) + is_complete_offset_) & flag);
  }

  inline void SetFlag(bool on, unsigned char flag, void *state) const {
    if (on) {
      *(static_cast<char*>(state) + is_complete_offset_) |= flag;
    } else {
      *(static_cast<char*>(state) + is_complete_offset_) &= (MASK ^ flag);
    }
  }

  inline bool HasFullContext(const void *state) const {
    return GetFlag(state, HAS_FULL_CONTEXT);
  }

  inline void SetHasFullContext(bool flag, void *state) const {
    SetFlag(flag, HAS_FULL_CONTEXT, state);
  }

  void FireFeatures(const State<5>& state, WordID cur, SparseVector<double>* feats) {
    FidTree* ft = &fidroot_;
    int n = 0;
    WordID buf[10];
    int ci = order_ - 1;
    WordID curword = cur;
    while(curword) {
      buf[n] = curword;
      int& fid = ft->fids[curword];
      ++n;
      if (!fid) {
        const char* code="_UBT456789";
        ostringstream os;
        os << code[n] << ':';
        for (int i = n-1; i >= 0; --i)
          os << (i != n-1 ? "_" : "") << TD::Convert(buf[i]);
        fid = FD::Convert(os.str());
      }
      feats->set_value(fid, 1);
      ft = &ft->levels[curword];
      --ci;
      if (ci < 0) break;
      curword = state[ci];
    }
  }

 public:
  void LookupWords(const TRule& rule, const vector<const void*>& ant_states, SparseVector<double>* feats, SparseVector<double>* est_feats, void* remnant) {
    double sum = 0.0;
    double est_sum = 0.0;
    int num_scored = 0;
    int num_estimated = 0;
    bool saw_eos = false;
    bool has_some_history = false;
    State<5> state;
    const vector<WordID>& e = rule.e();
    bool context_complete = false;
    for (int j = 0; j < e.size(); ++j) {
      if (e[j] < 1) {   // handle non-terminal substitution
        const void* astate = (ant_states[-e[j]]);
        int unscored_ant_len = UnscoredSize(astate);
        for (int k = 0; k < unscored_ant_len; ++k) {
          const WordID cur_word = IthUnscoredWord(k, astate);
          const bool is_oov = (cur_word == 0);
          SparseVector<double> p;
          if (cur_word == kSOS_) {
            state = BeginSentenceState();
            if (has_some_history) {  // this is immediately fully scored, and bad
              p.set_value(FD::Convert("Malformed"), 1.0);
              context_complete = true;
            } else {  // this might be a real <s>
              num_scored = max(0, order_ - 2);
            }
          } else {
            FireFeatures(state, cur_word, &p);
            const State<5> scopy = State<5>(state, order_, cur_word);
            state = scopy;
            if (saw_eos) { p.set_value(FD::Convert("Malformed"), 1.0); }
            saw_eos = (cur_word == kEOS_);
          }
          has_some_history = true;
          ++num_scored;
          if (!context_complete) {
            if (num_scored >= order_) context_complete = true;
          }
          if (context_complete) {
            (*feats) += p;
          } else {
            if (remnant)
              SetIthUnscoredWord(num_estimated, cur_word, remnant);
            ++num_estimated;
            (*est_feats) += p;
          }
        }
        saw_eos = GetFlag(astate, HAS_EOS_ON_RIGHT);
        if (HasFullContext(astate)) { // this is equivalent to the "star" in Chiang 2007
          state = RemnantLMState(astate);
          context_complete = true;
        }
      } else {   // handle terminal
        const WordID cur_word = e[j];
        SparseVector<double> p;
        if (cur_word == kSOS_) {
          state = BeginSentenceState();
          if (has_some_history) {  // this is immediately fully scored, and bad
            p.set_value(FD::Convert("Malformed"), -100);
            context_complete = true;
          } else {  // this might be a real <s>
            num_scored = max(0, order_ - 2);
          }
        } else {
          FireFeatures(state, cur_word, &p);
          const State<5> scopy = State<5>(state, order_, cur_word);
          state = scopy;
          if (saw_eos) { p.set_value(FD::Convert("Malformed"), 1.0); }
          saw_eos = (cur_word == kEOS_);
        }
        has_some_history = true;
        ++num_scored;
        if (!context_complete) {
          if (num_scored >= order_) context_complete = true;
        }
        if (context_complete) {
          (*feats) += p;
        } else {
          if (remnant)
            SetIthUnscoredWord(num_estimated, cur_word, remnant);
          ++num_estimated;
          (*est_feats) += p;
        }
      }
    }
    if (remnant) {
      SetFlag(saw_eos, HAS_EOS_ON_RIGHT, remnant);
      SetRemnantLMState(state, remnant);
      SetUnscoredSize(num_estimated, remnant);
      SetHasFullContext(context_complete || (num_scored >= order_), remnant);
    }
  }

  // this assumes no target words on final unary -> goal rule.  is that ok?
  // for <s> (n-1 left words) and (n-1 right words) </s>
  void FinalTraversal(const void* state, SparseVector<double>* feats) {
    if (add_sos_eos_) {  // rules do not produce <s> </s>, so do it here
      SetRemnantLMState(BeginSentenceState(), dummy_state_);
      SetHasFullContext(1, dummy_state_);
      SetUnscoredSize(0, dummy_state_);
      dummy_ants_[1] = state;
      LookupWords(*dummy_rule_, dummy_ants_, feats, NULL, NULL);
    } else {  // rules DO produce <s> ... </s>
#if 0
      double p = 0;
      if (!GetFlag(state, HAS_EOS_ON_RIGHT)) { p -= 100; }
      if (UnscoredSize(state) > 0) {  // are there unscored words
        if (kSOS_ != IthUnscoredWord(0, state)) {
          p -= 100 * UnscoredSize(state);
        }
      }
      return p;
#endif
    }
  }

 public:
  explicit NgramDetectorImpl(bool explicit_markers) :
      kCDEC_UNK(TD::Convert("<unk>")) ,
      add_sos_eos_(!explicit_markers) {
    order_ = 3;
    state_size_ = (order_ - 1) * sizeof(WordID) + 2 + (order_ - 1) * sizeof(WordID);
    unscored_size_offset_ = (order_ - 1) * sizeof(WordID);
    is_complete_offset_ = unscored_size_offset_ + 1;
    unscored_words_offset_ = is_complete_offset_ + 1;

    // special handling of beginning / ending sentence markers
    dummy_state_ = new char[state_size_];
    memset(dummy_state_, 0, state_size_);
    dummy_ants_.push_back(dummy_state_);
    dummy_ants_.push_back(NULL);
    dummy_rule_.reset(new TRule("[DUMMY] ||| [BOS] [DUMMY] ||| [1] [2] </s> ||| X=0"));
    kSOS_ = TD::Convert("<s>");
    kEOS_ = TD::Convert("</s>");
  }

  ~NgramDetectorImpl() {
    delete[] dummy_state_;
  }

  int ReserveStateSize() const { return state_size_; }

 private:
  const WordID kCDEC_UNK;
  WordID kSOS_;  // <s> - requires special handling.
  WordID kEOS_;  // </s>
  const bool add_sos_eos_; // flag indicating whether the hypergraph produces <s> and </s>
                     // if this is true, FinalTransitionFeatures will "add" <s> and </s>
                     // if false, FinalTransitionFeatures will score anything with the
                     // markers in the right place (i.e., the beginning and end of
                     // the sentence) with 0, and anything else with -100

  int order_;
  int state_size_;
  int unscored_size_offset_;
  int is_complete_offset_;
  int unscored_words_offset_;
  char* dummy_state_;
  vector<const void*> dummy_ants_;
  TRulePtr dummy_rule_;
  struct FidTree {
    map<WordID, int> fids;
    map<WordID, FidTree> levels;
  };
  mutable FidTree fidroot_;
};

NgramDetector::NgramDetector(const string& param) {
  string filename, mapfile, featname;
  bool explicit_markers = (param == "-x");
  pimpl_ = new NgramDetectorImpl(explicit_markers);
  SetStateSize(pimpl_->ReserveStateSize());
}

NgramDetector::~NgramDetector() {
  delete pimpl_;
}

void NgramDetector::TraversalFeaturesImpl(const SentenceMetadata& /* smeta */,
                                          const Hypergraph::Edge& edge,
                                          const vector<const void*>& ant_states,
                                          SparseVector<double>* features,
                                          SparseVector<double>* estimated_features,
                                          void* state) const {
  pimpl_->LookupWords(*edge.rule_, ant_states, features, estimated_features, state);
}

void NgramDetector::FinalTraversalFeatures(const void* ant_state,
                                           SparseVector<double>* features) const {
  pimpl_->FinalTraversal(ant_state, features);
}

