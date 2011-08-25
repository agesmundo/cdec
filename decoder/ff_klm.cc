#include "ff_klm.h"

#include <cstring>
#include <iostream>

#include <boost/scoped_ptr.hpp>

#include "filelib.h"
#include "stringlib.h"
#include "hg.h"
#include "tdict.h"
#include "lm/enumerate_vocab.hh"

using namespace std;

// Define the following macro if you want to see lots of debugging output
// when running Undirected KLM
//#define DEBUG_ULM
#undef DEBUG_ULM

static const unsigned char HAS_FULL_CONTEXT = 1;
static const unsigned char HAS_EOS_ON_RIGHT = 2;
static const unsigned char MASK             = 7;

// -x : rules include <s> and </s>
// -n NAME : feature id is NAME
bool ParseLMArgs(string const& in, string* filename, string* mapfile, bool* explicit_markers, string* featname) {
  vector<string> const& argv=SplitOnWhitespace(in);
  *explicit_markers = false;
  *featname="LanguageModel";
  *mapfile = "";
#define LMSPEC_NEXTARG if (i==argv.end()) {            \
    cerr << "Missing argument for "<<*last<<". "; goto usage; \
    } else { ++i; }

  for (vector<string>::const_iterator last,i=argv.begin(),e=argv.end();i!=e;++i) {
    string const& s=*i;
    if (s[0]=='-') {
      if (s.size()>2) goto fail;
      switch (s[1]) {
      case 'x':
        *explicit_markers = true;
        break;
      case 'm':
        LMSPEC_NEXTARG; *mapfile=*i;
        break;
      case 'n':
        LMSPEC_NEXTARG; *featname=*i;
        break;
#undef LMSPEC_NEXTARG
      default:
      fail:
        cerr<<"Unknown KLanguageModel option "<<s<<" ; ";
        goto usage;
      }
    } else {
      if (filename->empty())
        *filename=s;
      else {
        cerr<<"More than one filename provided. ";
        goto usage;
      }
    }
  }
  if (!filename->empty())
    return true;
usage:
  cerr << "KLanguageModel is incorrect!\n";
  return false;
}

template <class Model>
string KLanguageModel<Model>::usage(bool /*param*/,bool /*verbose*/) {
  return "KLanguageModel";
}

struct VMapper : public lm::ngram::EnumerateVocab {
  VMapper(vector<lm::WordIndex>* out) : out_(out), kLM_UNKNOWN_TOKEN(0) { out_->clear(); }
  void Add(lm::WordIndex index, const StringPiece &str) {
    const WordID cdec_id = TD::Convert(str.as_string());
    if (cdec_id >= out_->size())
      out_->resize(cdec_id + 1, kLM_UNKNOWN_TOKEN);
    (*out_)[cdec_id] = index;
  }
  vector<lm::WordIndex>* out_;
  const lm::WordIndex kLM_UNKNOWN_TOKEN;
};

template <class Model>
class KLanguageModelImpl {

  // returns the number of unscored words at the left edge of a span
  inline int UnscoredSize(const void* state) const {
    return *(static_cast<const char*>(state) + unscored_size_offset_);
  }

  inline void SetUnscoredSize(int size, void* state) const {
    *(static_cast<char*>(state) + unscored_size_offset_) = size;
  }

  static inline const lm::ngram::State& RemnantLMState(const void* state) {
    return *static_cast<const lm::ngram::State*>(state);
  }

  inline void SetRemnantLMState(const lm::ngram::State& lmstate, void* state) const {
    // if we were clever, we could use the memory pointed to by state to do all
    // the work, avoiding this copy
    memcpy(state, &lmstate, ngram_->StateSize());
  }

  lm::WordIndex IthUnscoredWord(int i, const void* state) const {
    const lm::WordIndex* const mem = reinterpret_cast<const lm::WordIndex*>(static_cast<const char*>(state) + unscored_words_offset_);
    return mem[i];
  }

  void SetIthUnscoredWord(int i, lm::WordIndex index, void *state) const {
    lm::WordIndex* mem = reinterpret_cast<lm::WordIndex*>(static_cast<char*>(state) + unscored_words_offset_);
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

 public:
  double LookupWords(const TRule& rule, const vector<const void*>& ant_states, double* pest_sum, double* oovs, double* est_oovs, void* remnant) {
    double sum = 0.0;
    double est_sum = 0.0;
    int num_scored = 0;
    int num_estimated = 0;
    if (oovs) *oovs = 0;
    if (est_oovs) *est_oovs = 0;
    bool saw_eos = false;
    bool has_some_history = false;
    lm::ngram::State state = ngram_->NullContextState();
    const vector<WordID>& e = rule.e();
    bool context_complete = false;
    for (int j = 0; j < e.size(); ++j) {
      if (e[j] < 1) {   // handle non-terminal substitution
        const void* astate = (ant_states[-e[j]]);
        int unscored_ant_len = UnscoredSize(astate);
        for (int k = 0; k < unscored_ant_len; ++k) {
          const lm::WordIndex cur_word = IthUnscoredWord(k, astate);
          const bool is_oov = (cur_word == 0);
          double p = 0;
          if (cur_word == kSOS_) {
            state = ngram_->BeginSentenceState();
            if (has_some_history) {  // this is immediately fully scored, and bad
              p = -100;
              context_complete = true;
            } else {  // this might be a real <s>
              num_scored = max(0, order_ - 2);
            }
          } else {
            const lm::ngram::State scopy(state);
            p = ngram_->Score(scopy, cur_word, state);
            if (saw_eos) { p = -100; }
            saw_eos = (cur_word == kEOS_);
          }
          has_some_history = true;
          ++num_scored;
          if (!context_complete) {
            if (num_scored >= order_) context_complete = true;
          }
          if (context_complete) {
            sum += p;
            if (oovs && is_oov) (*oovs)++;
          } else {
            if (remnant)
              SetIthUnscoredWord(num_estimated, cur_word, remnant);
            ++num_estimated;
            est_sum += p;
            if (est_oovs && is_oov) (*est_oovs)++;
          }
        }
        saw_eos = GetFlag(astate, HAS_EOS_ON_RIGHT);
        if (HasFullContext(astate)) { // this is equivalent to the "star" in Chiang 2007
          state = RemnantLMState(astate);
          context_complete = true;
        }
      } else {   // handle terminal
        const WordID cdec_word_or_class = ClassifyWordIfNecessary(e[j]);  // in future,
                                                                          // maybe handle emission
        const lm::WordIndex cur_word = MapWord(cdec_word_or_class); // map to LM's id
        double p = 0;
        const bool is_oov = (cur_word == 0);
        if (cur_word == kSOS_) {
          state = ngram_->BeginSentenceState();
          if (has_some_history) {  // this is immediately fully scored, and bad
            p = -100;
            context_complete = true;
          } else {  // this might be a real <s>
            num_scored = max(0, order_ - 2);
          }
        } else {
          const lm::ngram::State scopy(state);
          p = ngram_->Score(scopy, cur_word, state);
          if (saw_eos) { p = -100; }
          saw_eos = (cur_word == kEOS_);
        }
        has_some_history = true;
        ++num_scored;
        if (!context_complete) {
          if (num_scored >= order_) context_complete = true;
        }
        if (context_complete) {
          sum += p;
          if (oovs && is_oov) (*oovs)++;
        } else {
          if (remnant)
            SetIthUnscoredWord(num_estimated, cur_word, remnant);
          ++num_estimated;
          est_sum += p;
          if (est_oovs && is_oov) (*est_oovs)++;
        }
      }
    }
    if (pest_sum) *pest_sum = est_sum;
    if (remnant) {
      state.ZeroRemaining();
      SetFlag(saw_eos, HAS_EOS_ON_RIGHT, remnant);
      SetRemnantLMState(state, remnant);
      SetUnscoredSize(num_estimated, remnant);
      SetHasFullContext(context_complete || (num_scored >= order_), remnant);
    }
    return sum;
  }
//GU
  double UndirectedLookupWords(UCandidate& ucand, /*const vector<const void*>& ant_states, double* pest_sum,*/ double* oovs/*, double* est_oovs, void* remnant*/,int spos) {
    double sum = 0.0;
//    double est_sum = 0.0;
    int num_scored = 0;
    int num_estimated = 0;
    if (oovs) *oovs = 0;
//    if (est_oovs) *est_oovs = 0;
    bool saw_eos = false;
    bool has_some_history = false;
    bool context_complete = false;

    const Hypergraph::Edge& in_edge = *ucand.in_edge_;
    const TRule& rule=*in_edge.rule_.get();
    const int head_node_id = in_edge.head_node_;
    //const int source_node_id = ucand.GetSourceNodeId();
    const vector<WordID>& e = rule.e();
    lm::ngram::State state;
    //bool hole =false;
    void* current_outgoing_state;

    //head outgoing state
    current_outgoing_state=NULL;
    void* head_outgoing_state=NULL;
    FFState* ffs_head_out=ucand.GetOutgoingState(head_node_id);
    if(ffs_head_out!=NULL){
    	head_outgoing_state= FFS2LMS(ffs_head_out,spos);
    }
    current_outgoing_state=head_outgoing_state;

#ifdef DEBUG_ULM
    cerr << "-----------------------\nUNDIRECTED LOOKUP WORDS"<<endl;
#endif

    //TODO GU these alternatives are related to head_outgoing_state selection (above)
	FFState* ffs_head_in = ucand.GetHeadIncomingState();
	void* head_incoming_state=NULL;
    if(ffs_head_in!=NULL){ //head is incoming state
    	head_incoming_state = FFS2LMS(ffs_head_in,spos);
    	state = RemnantLMState(head_incoming_state);
    	num_scored = state.ValidLength();
    	context_complete =HasFullContext(head_incoming_state);
#ifdef DEBUG_ULM
    	cerr << " current_out : head_outgoing_state = ";
    	PringLMS(head_incoming_state);
#endif
    }
    else{
    	state= ngram_->NullContextState();
    	context_complete = false;
    }

    for (int j = 0; j < e.size(); ++j) {
      if (e[j] < 1) {   // handle non-terminal substitution
    	  int tail_id=-e[j];
		  FFState* ffs_tail_in = ucand.GetTailIncomingState(tail_id);
		  if(ffs_tail_in!=NULL){
			  void* tail_incoming_state = FFS2LMS(ffs_tail_in,spos);
#ifdef DEBUG_ULM
			  cerr << " tail_incoming_state = ";
			  PringLMS(tail_incoming_state);
#endif
			  int unscored_ant_len = UnscoredSize(tail_incoming_state);
			  for (int k = 0; k < unscored_ant_len; ++k) {
				  const lm::WordIndex cur_word = IthUnscoredWord(k, tail_incoming_state);
				  const bool is_oov = (cur_word == 0);
				  double p = 0;
				  if (cur_word == kSOS_) {
					  state = ngram_->BeginSentenceState();
					  if (has_some_history) {  // this is immediately fully scored, and bad
						  p = -100;
						  context_complete = true;
					  } else {  // this might be a real <s>
						  num_scored = max(0, order_ - 2);
					  }
				  } else {
					  const lm::ngram::State scopy(state);
					  p = ngram_->Score(scopy, cur_word, state);
#ifdef DEBUG_ULM
					  cerr << " scopy before Score() : " << scopy << endl;
					  cerr << " state after  Score() : " << state << endl;
#endif
					  if (saw_eos) { p = -100; }
					  saw_eos = (cur_word == kEOS_);
				  }
				  has_some_history = true;
				  ++num_scored;
				  if (!context_complete) {
					  if (num_scored >= order_) context_complete = true;
				  }
				  if (context_complete) {
					  sum += p;
					  if (oovs && is_oov) (*oovs)++;
				  } else {
					  if (current_outgoing_state /*&& !hole*/){
						  SetIthUnscoredWord(num_estimated, cur_word, current_outgoing_state);
#ifdef DEBUG_ULM
						  cerr << " curr_outgoing_state unscored word [" << num_estimated << "] to " << cur_word << endl;
#endif
					  }
					  ++num_estimated;
					  //					  est_sum += p;
					  //					  if (est_oovs && is_oov) (*est_oovs)++;
				  }
			  }
			  saw_eos = GetFlag(tail_incoming_state, HAS_EOS_ON_RIGHT);
			  if (HasFullContext(tail_incoming_state)) { // this is equivalent to the "star" in Chiang 2007
				  state = RemnantLMState(tail_incoming_state);
				  context_complete = true;
			  }
		  } else { //there is an hole

			  if (current_outgoing_state){
				  //set exiting sate details
				  state.ZeroRemaining();
				  SetFlag(saw_eos, HAS_EOS_ON_RIGHT, current_outgoing_state);//??correct? should be reset?
				  SetUnscoredSize(num_estimated, current_outgoing_state);
				  SetHasFullContext(context_complete || (num_scored >= order_), current_outgoing_state);
#ifdef DEBUG_ULM
				  cerr << " current_outgoing_state = ";
				  PringLMS(current_outgoing_state);
#endif
			  }

			  //set KLMState of tail out state
			  FFState* ffs_tail_out = ucand.GetOutgoingState(in_edge.tail_nodes_[tail_id]);
#ifdef DEBUG_ULM
			  assert(ffs_tail_out!=NULL);
#endif
			  current_outgoing_state = FFS2LMS(ffs_tail_out,spos);
			  SetRemnantLMState(state,current_outgoing_state);

			  state= ngram_->NullContextState();
			  SetUnscoredSize(num_estimated, current_outgoing_state);
			  num_estimated=0;
			  context_complete = false;
			  num_scored=0;
		  }
      } else {   // *1 handle terminal
    	  const WordID cdec_word_or_class = ClassifyWordIfNecessary(e[j]);  // in future,
    	  // maybe handle emission
    	  const lm::WordIndex cur_word = MapWord(cdec_word_or_class); // map to LM's id
    	  double p = 0;
    	  const bool is_oov = (cur_word == 0);
    	  if (cur_word == kSOS_) {
    		  state = ngram_->BeginSentenceState();
    		  if (has_some_history) {  // this is immediately fully scored, and bad
    			  p = -100;
    			  context_complete = true;
    		  } else {  // this might be a real <s>
    			  num_scored = max(0, order_ - 2);//actual order-1 since ++ later
    		  }
    	  } else {
    		  const lm::ngram::State scopy(state);
    		  p = ngram_->Score(scopy, cur_word, state);
#ifdef DEBUG_ULM
    		  cerr << " scopy before Score() : " << scopy << endl;
    		  cerr << " state after  Score() : " << state << endl;
#endif
    		  if (saw_eos) { p = -100; }
    		  saw_eos = (cur_word == kEOS_);
    	  }
    	  has_some_history = true;
    	  ++num_scored;
    	  if (!context_complete) {
    		  if (num_scored >= order_) context_complete = true;
    	  }
    	  if (context_complete) {
    		  sum += p;
    		  if (oovs && is_oov) (*oovs)++;
    	  } else {
    		  //set the head out state if needed
    		  if (current_outgoing_state/* && !hole*/){
    			  SetIthUnscoredWord(num_estimated, cur_word, current_outgoing_state);
#ifdef DEBUG_ULM
    			  cerr << " current_outgoing_state unscored word [" << num_estimated << "] to " << cur_word << endl;
#endif
    		  }
    		  ++num_estimated;
    		  //          est_sum += p;
    		  //          if (est_oovs && is_oov) (*est_oovs)++;
    	  }
      }
    }

    //read head_incoming_context unscored words (same as *1 handle terminal)
    if(head_incoming_state){
    	int unscored_size_his =UnscoredSize(head_incoming_state);
    	for(int i=0;i<unscored_size_his;i++){
    		const lm::WordIndex cur_word = IthUnscoredWord(i, head_incoming_state);
    		double p = 0;
    		const bool is_oov = (cur_word == 0);
    		if (cur_word == kSOS_) {
    			state = ngram_->BeginSentenceState();
    			if (has_some_history) {  // this is immediately fully scored, and bad
    				p = -100;
    				context_complete = true;
    			} else {  // this might be a real <s>
    				num_scored = max(0, order_ - 2);//actual order-1 since ++ later
    			}
    		} else {
    			const lm::ngram::State scopy(state);
    			p = ngram_->Score(scopy, cur_word, state);
#ifdef DEBUG_ULM
    			cerr << " scopy before Score() : " << scopy << endl;
    			cerr << " state after  Score() : " << state << endl;
#endif
    			if (saw_eos) { p = -100; }
    			saw_eos = (cur_word == kEOS_);
    		}
    		has_some_history = true;
    		++num_scored;
    		if (!context_complete) {
    			if (num_scored >= order_) context_complete = true;
    		}
    		if (context_complete) {
    			sum += p;
    			if (oovs && is_oov) (*oovs)++;
    		} else {
    			//set the head out state if needed
    			if (current_outgoing_state/* && !hole*/){
    				SetIthUnscoredWord(num_estimated, cur_word, current_outgoing_state);
#ifdef DEBUG_ULM
    				cerr << " current_outgoing_state unscored word [" << num_estimated << "] to " << cur_word << endl;
#endif
    			}
    			++num_estimated;
    			//          est_sum += p;
    			//          if (est_oovs && is_oov) (*est_oovs)++;
    		}
    	}
    }

    //    if (pest_sum) *pest_sum = est_sum;
    if (current_outgoing_state) {
    	state.ZeroRemaining(); //TODO try kMaxOrder to 3 for performances
    	SetFlag(saw_eos, HAS_EOS_ON_RIGHT, current_outgoing_state);//TODO? set this flag for head_outgoing_state?
    	SetUnscoredSize(num_estimated, current_outgoing_state);
    	if(head_outgoing_state) SetRemnantLMState(state, head_outgoing_state);
    	SetHasFullContext(context_complete || (num_scored/*+1*/ >= order_), current_outgoing_state);//?+1 since flag if next will have full context |order-1|
#ifdef DEBUG_ULM
    	if(head_outgoing_state){
    		cerr << " final head_outgoing_state = ";
    		PringLMS(head_outgoing_state);
    	}
#endif
    }

#ifdef DEBUG_ULM
    cerr << "-----------------------"<<endl;
#endif
    return sum;
  }
  //GU Utils
  void* FFS2LMS(FFState* ffs_state,int spos){ //FFState to LMState
	  return &(*ffs_state)[spos];
  }
  void PringLMS(void* lms){//Print LMS (both sides)
	  cerr << "LMSS (" << lms << ") = [ " ;
	  cerr << RemnantLMState(lms)<<endl;
	  int unscored_size=UnscoredSize(lms);
	  cerr << "\t\t unscored size    : " << unscored_size << endl;
	  cerr << "\t\t has eos on right : " << GetFlag(lms, HAS_EOS_ON_RIGHT) << endl;
	  cerr << "\t\t has full context : " << HasFullContext(lms) << endl;
	  cerr << "\t\t boundary words   : <";
	  for(int i=0; i<unscored_size; i++){
		  cerr << " W[" << i << "] : " << IthUnscoredWord(i, lms);
	  }
	  cerr << " >" << endl;
	  cerr << "]" << endl;
  }

  // this assumes no target words on final unary -> goal rule.  is that ok?
  // for <s> (n-1 left words) and (n-1 right words) </s>
  double FinalTraversalCost(const void* state, double* oovs) {
    if (add_sos_eos_) {  // rules do not produce <s> </s>, so do it here
      SetRemnantLMState(ngram_->BeginSentenceState(), dummy_state_);
      SetHasFullContext(1, dummy_state_);
      SetUnscoredSize(0, dummy_state_);
      dummy_ants_[1] = state;
      *oovs = 0;
      return LookupWords(*dummy_rule_, dummy_ants_, NULL, oovs, NULL, NULL);
    } else {  // rules DO produce <s> ... </s>
      double p = 0;
      if (!GetFlag(state, HAS_EOS_ON_RIGHT)) { p -= 100; }
      if (UnscoredSize(state) > 0) {  // are there unscored words
        if (kSOS_ != IthUnscoredWord(0, state)) {
          p -= 100 * UnscoredSize(state);
        }
      }
      return p;
    }
  }

  // if this is not a class-based LM, returns w untransformed,
  // otherwise returns a word class mapping of w,
  // returns TD::Convert("<unk>") if there is no mapping for w
  WordID ClassifyWordIfNecessary(WordID w) const {
    if (word2class_map_.empty()) return w;
    if (w >= word2class_map_.size())
      return kCDEC_UNK;
    else
      return word2class_map_[w];
  }

  // converts to cdec word id's to KenLM's id space, OOVs and <unk> end up at 0
  lm::WordIndex MapWord(WordID w) const {
    if (w >= cdec2klm_map_.size())
      return 0;
    else
      return cdec2klm_map_[w];
  }

 public:
  KLanguageModelImpl(const string& filename, const string& mapfile, bool explicit_markers) :
      kCDEC_UNK(TD::Convert("<unk>")) ,
      add_sos_eos_(!explicit_markers) {
    {
      VMapper vm(&cdec2klm_map_);
      lm::ngram::Config conf;
      conf.enumerate_vocab = &vm;
      ngram_ = new Model(filename.c_str(), conf);
    }
    order_ = ngram_->Order();
    cerr << "Loaded " << order_ << "-gram KLM from " << filename << " (MapSize=" << cdec2klm_map_.size() << ")\n";
    state_size_ = ngram_->StateSize() + 2 + (order_ - 1) * sizeof(lm::WordIndex);
    unscored_size_offset_ = ngram_->StateSize();
    is_complete_offset_ = unscored_size_offset_ + 1;
    unscored_words_offset_ = is_complete_offset_ + 1;

    // special handling of beginning / ending sentence markers
    dummy_state_ = new char[state_size_];
    memset(dummy_state_, 0, state_size_);
    dummy_ants_.push_back(dummy_state_);
    dummy_ants_.push_back(NULL);
    dummy_rule_.reset(new TRule("[DUMMY] ||| [BOS] [DUMMY] ||| [1] [2] </s> ||| X=0"));
    kSOS_ = MapWord(TD::Convert("<s>"));
    assert(kSOS_ > 0);
    kEOS_ = MapWord(TD::Convert("</s>"));
    assert(kEOS_ > 0);
    assert(MapWord(kCDEC_UNK) == 0); // KenLM invariant

    // handle class-based LMs (unambiguous word->class mapping reqd.)
    if (mapfile.size())
      LoadWordClasses(mapfile);
  }

  void LoadWordClasses(const string& file) {
    ReadFile rf(file);
    istream& in = *rf.stream();
    string line;
    vector<WordID> dummy;
    int lc = 0;
    cerr << "  Loading word classes from " << file << " ...\n";
    AddWordToClassMapping_(TD::Convert("<s>"), TD::Convert("<s>"));
    AddWordToClassMapping_(TD::Convert("</s>"), TD::Convert("</s>"));
    while(in) {
      getline(in, line);
      if (!in) continue;
      dummy.clear();
      TD::ConvertSentence(line, &dummy);
      ++lc;
      if (dummy.size() != 2) {
        cerr << "    Format error in " << file << ", line " << lc << ": " << line << endl;
        abort();
      }
      AddWordToClassMapping_(dummy[0], dummy[1]);
    }
  }

  void AddWordToClassMapping_(WordID word, WordID cls) {
    if (word2class_map_.size() <= word) {
      word2class_map_.resize((word + 10) * 1.1, kCDEC_UNK);
      assert(word2class_map_.size() > word);
    }
    if(word2class_map_[word] != kCDEC_UNK) {
      cerr << "Multiple classes for symbol " << TD::Convert(word) << endl;
      abort();
    }
    word2class_map_[word] = cls;
  }

  ~KLanguageModelImpl() {
    delete ngram_;
    delete[] dummy_state_;
  }

  int ReserveStateSize() const { return state_size_; }

 private:
  const WordID kCDEC_UNK;
  lm::WordIndex kSOS_;  // <s> - requires special handling.
  lm::WordIndex kEOS_;  // </s>
  Model* ngram_;
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
  vector<lm::WordIndex> cdec2klm_map_;
  vector<WordID> word2class_map_;        // if this is a class-based LM, this is the word->class mapping
  TRulePtr dummy_rule_;
};

template <class Model>
KLanguageModel<Model>::KLanguageModel(const string& param) {
  string filename, mapfile, featname;
  bool explicit_markers;
  if (!ParseLMArgs(param, &filename, &mapfile, &explicit_markers, &featname)) {
    abort();
  }
  try {
    pimpl_ = new KLanguageModelImpl<Model>(filename, mapfile, explicit_markers);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    abort();
  }
  fid_ = FD::Convert(featname);
  oov_fid_ = FD::Convert(featname+"_OOV");
  cerr << "FID: " << oov_fid_ << endl;
  SetStateSize(pimpl_->ReserveStateSize());
}

template <class Model>
Features KLanguageModel<Model>::features() const {
  return single_feature(fid_);
}

template <class Model>
KLanguageModel<Model>::~KLanguageModel() {
  delete pimpl_;
}

template <class Model>
void KLanguageModel<Model>::TraversalFeaturesImpl(const SentenceMetadata& /* smeta */,
                                          const Hypergraph::Edge& edge,
                                          const vector<const void*>& ant_states,
                                          SparseVector<double>* features,
                                          SparseVector<double>* estimated_features,
                                          void* state) const {
  double est = 0;
  double oovs = 0;
  double est_oovs = 0;
  features->set_value(fid_, pimpl_->LookupWords(*edge.rule_, ant_states, &est, &oovs, &est_oovs, state));
  estimated_features->set_value(fid_, est);
  if (oov_fid_) {
    if (oovs) features->set_value(oov_fid_, oovs);
    if (est_oovs) estimated_features->set_value(oov_fid_, est_oovs);
  }
}
//GU
template <class Model>
void KLanguageModel<Model>::TraversalUndirectedFeaturesImpl(const SentenceMetadata& /*smeta*/,
                                        UCandidate& ucand,
                                        int spos
                                        /*const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state*/) const {
	//add features like "worst&best scoring trigram" "number of trigrams"
	//fid_ add average score of trigrams generated by this action
	//oov_fid_ add number of out of vocabolary (LM) words
//	  double est = 0;//TODO GU meaning of estimated?
	  double oovs = 0;
//	  double est_oovs = 0;
	  ucand.feature_values_.set_value(fid_, pimpl_->UndirectedLookupWords(ucand, /*ant_states, &est,*/ &oovs/*, &est_oovs, state*/,spos));
//	  ucand.est_vals_.set_value(fid_, est);
	  if (oov_fid_) {
	    if (oovs) ucand.feature_values_.set_value(oov_fid_, oovs);
//	    if (est_oovs) ucand.est_vals_.set_value(oov_fid_, est_oovs);
	  }
}

template <class Model>
void KLanguageModel<Model>::FinalTraversalFeatures(const void* ant_state,
                                           SparseVector<double>* features) const {
  double oovs = 0;
  double lm = pimpl_->FinalTraversalCost(ant_state, &oovs);
  features->set_value(fid_, lm);
  if (oov_fid_ && oovs)
    features->set_value(oov_fid_, oovs);
}

// instantiate templates
template class KLanguageModel<lm::ngram::ProbingModel>;
template class KLanguageModel<lm::ngram::TrieModel>;
template class KLanguageModel<lm::ngram::QuantTrieModel>;

