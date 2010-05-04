#include <queue>
#include <iostream>
#include <vector>
#include <utility>
#include <tr1/unordered_map>
#include <set>

#include <boost/functional/hash.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include "tdict.h"
#include "wordid.h"
#include "array2d.h"
#include "filelib.h"

using namespace std;
using namespace std::tr1;
namespace po = boost::program_options;

static const size_t MAX_LINE_LENGTH = 100000;
WordID kBOS;
WordID kEOS;
WordID kDIVIDER;

struct ParallelSpan {
  // i1 = i of f side
  // i2 = j of f side
  // j1 = i of e side
  // j2 = j of e side
  short i1,i2,j1,j2;
  // cat is set by AnnotatePhrasesWithCategoryTypes, otherwise it's 0
  WordID cat;  // category type of span (also overloaded by Rule class)
  ParallelSpan() : i1(-1), i2(-1), j1(-1), j2(-1), cat() {}
  // used by Rule class to represent a terminal symbol:
  explicit ParallelSpan(WordID w) : i1(-1), i2(-1), j1(-1), j2(-1), cat(w) {}
  ParallelSpan(int pi1, int pi2, int pj1, int pj2) : i1(pi1), i2(pi2), j1(pj1), j2(pj2), cat() {}
  ParallelSpan(int pi1, int pi2, int pj1, int pj2, WordID c) : i1(pi1), i2(pi2), j1(pj1), j2(pj2), cat(c) {}

  // ParallelSpan is used in the Rule class where it is
  // overloaded to also represent terminal symbols
  bool IsVariable() const { return i1 != -1; }
};

// represents a parallel sentence with a word alignment and category
// annotations over subspans (currently in terms of f)
struct AnnotatedParallelSentence {
  vector<WordID> f, e;  // words in f and e
  vector<int> e_aligned, f_aligned; // counts the number of times column/row x is aligned
  Array2D<bool> aligned;

  Array2D<vector<WordID> > span_types;  // span_types(i,j) is the list of category
                               // types for a span (i,j) in the TARGET language.
  int f_len, e_len;
  vector<pair<int,int> > f_spans, e_spans;

  // read annotated parallel sentence from string
  void ParseInputLine(const char* buf);

 private:
  void HandleAlignment(const char* buf, int start, int end);
  void Reset() { f.clear(); e.clear(); e_aligned.clear(); f_aligned.clear(); aligned.clear(); span_types.clear(); f_spans.clear(); e_spans.clear(); }
  void HandleSpan(const char* buf, int start, int end);

  void AllocateForAlignment() {
    aligned.resize(f.size(), e.size(), false);
    f_aligned.resize(f.size(), 0);
    e_aligned.resize(e.size(), 0);
    span_types.resize(f.size(), f.size()+1);
    f_len = f.size();
    e_len = e.size();
  }
  static int ReadAlignmentPoint(const char* buf, int start, int end, bool permit_col, int* a, int* b);
};

void InitCommandLine(int argc, char** argv, po::variables_map* conf) {
  po::options_description opts("Configuration options");
  opts.add_options()
        ("input,i", po::value<string>()->default_value("-"), "Input file")
        ("default_category,d", po::value<string>(), "Default span type (use X for 'Hiero')")
        ("loose", "Use loose phrase extraction heuristic for base phrases")
        ("invert,I", "Invert the outputs")
        ("base_phrase,B", "Write base phrases")
        ("phrase_context,C", "Write base phrase contexts")
        ("phrase_context_size,S", po::value<int>()->default_value(2), "Use this many words of context on left and write when writing base phrase contexts")
        ("max_base_phrase_size,L", po::value<int>()->default_value(10), "Maximum starting phrase size")
        ("help,h", "Print this help message and exit");
  po::options_description clo("Command line options");
  po::options_description dcmdline_options;
  dcmdline_options.add(opts);

  po::store(parse_command_line(argc, argv, dcmdline_options), *conf);
  po::notify(*conf);

  if (conf->count("help") || conf->count("input") == 0) {
    cerr << "\nUsage: extractor [-options]\n";
    cerr << dcmdline_options << endl;
    exit(1);
  }
}

inline bool IsWhitespace(char c) { return c == ' ' || c == '\t'; }

inline void SkipWhitespace(const char* buf, int* ptr) {
  while (buf[*ptr] && IsWhitespace(buf[*ptr])) { ++(*ptr); }
}

int AnnotatedParallelSentence::ReadAlignmentPoint(const char* buf, int start, int end, bool permit_col, int* a, int* b) {
  if (end - start < 3) {
    cerr << "Alignment point badly formed: " << string(buf, start, end-start) << endl; abort();
  }
  int c = start;
  *a = 0;
  while(c < end && buf[c] != '-') {
    if (buf[c] < '0' || buf[c] > '9') {
      cerr << "Alignment point badly formed: " << string(buf, start, end-start) << endl;
      abort();
    }
    (*a) *= 10;
    (*a) += buf[c] - '0';
    ++c;
  }
  ++c;
  if (c >= end) {
    cerr << "Alignment point badly formed: " << string(buf, start, end-start) << endl; abort();
  }
  (*b) = 0;
  while(c < end && (!permit_col || (permit_col && buf[c] != ':'))) {
    if (buf[c] < '0' || buf[c] > '9') {
      cerr << "Alignment point badly formed: " << string(buf, start, end-start) << endl;
      abort();
    }
    (*b) *= 10;
    (*b) += buf[c] - '0';
    ++c;
  }
  return c;
}

void AnnotatedParallelSentence::HandleAlignment(const char* buf, int start, int end) {
  int a, b;
  ReadAlignmentPoint(buf, start, end, false, &a, &b);
  assert(a < f_len);
  assert(b < e_len);
  aligned(a,b) = true;
  ++f_aligned[a];
  ++e_aligned[b];
  // cerr << a << " " << b << endl;
}

void AnnotatedParallelSentence::HandleSpan(const char* buf, int start, int end) {
  int a,b;
  int c = ReadAlignmentPoint(buf, start, end, true, &a, &b) + 1;
  if (buf[c-1] != ':' || c >= end) {
    cerr << "Span badly formed: " << string(buf, start, end-start) << endl; abort();
  }
  // cerr << a << " " << b << " " << string(buf,c,end-c) << endl;
  span_types(a,b).push_back(-TD::Convert(string(buf, c, end-c)));
}

// INPUT FORMAT
// ein haus ||| a house ||| 0-0 1-1 ||| 0-0:DT 1-1:NN 0-1:NP
void AnnotatedParallelSentence::ParseInputLine(const char* buf) {
  Reset();
  int ptr = 0;
  SkipWhitespace(buf, &ptr);
  int start = ptr;
  int state = 0;  // 0 = French, 1 = English, 2 = Alignment, 3 = Spans
  while(char c = buf[ptr]) {
    if (!IsWhitespace(c)) { ++ptr; continue; } else {
      if (ptr - start == 3 && buf[start] == '|' && buf[start+1] == '|' && buf[start+2] == '|') {
        ++state;
        if (state == 4) { cerr << "Too many fields (ignoring):\n  " << buf << endl; return; }
        if (state == 2) {
          // cerr << "FLEN=" << f->size() << " ELEN=" << e->size() << endl;
          AllocateForAlignment();
        }
        SkipWhitespace(buf, &ptr);
        start = ptr;
        continue;
      }
      switch (state) {
        case 0:  f.push_back(TD::Convert(string(buf, start, ptr-start))); break;
        case 1:  e.push_back(TD::Convert(string(buf, start, ptr-start))); break;
        case 2:  HandleAlignment(buf, start, ptr); break;
        case 3:  HandleSpan(buf, start, ptr); break;
        default: cerr << "Can't happen\n"; abort();
      }
      SkipWhitespace(buf, &ptr);
      start = ptr;
    }
  }
  if (ptr > start) {
    switch (state) {
      case 0:  f.push_back(TD::Convert(string(buf, start, ptr-start))); break;
      case 1:  e.push_back(TD::Convert(string(buf, start, ptr-start))); break;
      case 2:  HandleAlignment(buf, start, ptr); break;
      case 3:  HandleSpan(buf, start, ptr); break;
      default: cerr << "Can't happen\n"; abort();
    }
  }
  if (state < 2) {
    cerr << "Not enough fields: " << buf << endl;
    abort();
  }
  if (e.empty() || f.empty()) {
    cerr << "Sentences must not be empty: " << buf << endl;
  }

  // for each alignment point in e, precompute the minimal consistent phrases in f
  // for each alignment point in f, precompute the minimal consistent phrases in e
  e_spans.resize(e_len, pair<int,int>(f_len, 0));
  f_spans.resize(f_len, pair<int,int>(e_len, 0));
  for (int i = 0; i < f_len; ++i) {
    for (int j = 0; j < e_len; ++j) {
      if (aligned(i,j)) {
        if (j < f_spans[i].first) f_spans[i].first = j;
        f_spans[i].second = j+1;
        if (i < e_spans[j].first) e_spans[j].first = i;
        e_spans[j].second = i+1;
      }
    }
  }
}

void LoosenPhraseBounds(const AnnotatedParallelSentence& sentence,
                        vector<ParallelSpan>* phrases) {
  assert(!"not implemented - TODO");
  (void) sentence;
  (void) phrases;
}

void ExtractBasePhrases(const int max_base_phrase_size,
                        const AnnotatedParallelSentence& sentence,
                        vector<ParallelSpan>* phrases) {
  phrases->clear();
  for (int i1 = 0; i1 < sentence.f_len; ++i1) {
    if (sentence.f_aligned[i1] == 0) continue;
    int j1 = sentence.e_len;
    int j2 = 0;
    const int i_limit = min(sentence.f_len, i1 + max_base_phrase_size);
    for (int i2 = i1 + 1; i2 <= i_limit; ++i2) {
      if (sentence.f_aligned[i2-1] == 0) continue;
      // cerr << "F has aligned span " << i1 << " to " << i2 << endl;
      j1 = min(j1, sentence.f_spans[i2-1].first);
      j2 = max(j2, sentence.f_spans[i2-1].second);
      if (j1 >= j2) continue;
      if (j2 - j1 > max_base_phrase_size) continue;
      int condition = 0;
      for (int j = j1; j < j2; ++j) {
        if (sentence.e_spans[j].first < i1) { condition = 1; break; }
        if (sentence.e_spans[j].second > i2) { condition = 2; break; }
      }
      if (condition == 1) break;
      if (condition == 2) continue;
      // category types added later!
      phrases->push_back(ParallelSpan(i1, i2, j1, j2));
      // cerr << i1 << " " << i2 << " : " << j1 << " " << j2 << endl;
    }
  }
}

// TODO how to handle alignment information?
void WriteBasePhrases(const AnnotatedParallelSentence& sentence,
                      const vector<ParallelSpan>& phrases) {
  vector<WordID> e,f;
  for (int it = 0; it < phrases.size(); ++it) {
    const ParallelSpan& phrase = phrases[it];
    e.clear();
    f.clear();
    for (int i = phrase.i1; i < phrase.i2; ++i)
      f.push_back(sentence.f[i]);
    for (int j = phrase.j1; j < phrase.j2; ++j)
      e.push_back(sentence.e[j]);
    cerr << TD::GetString(f) << " ||| " << TD::GetString(e) << endl;
  }
}

// TODO optional source context
void WritePhraseContexts(const AnnotatedParallelSentence& sentence,
                         const vector<ParallelSpan>& phrases,
                         const int ctx_size) {
  vector<WordID> context(ctx_size * 2 + 1);
  context[ctx_size] = kDIVIDER;
  for (int it = 0; it < phrases.size(); ++it) {
    const ParallelSpan& phrase = phrases[it];
    for (int i = 0; i < ctx_size; ++i) {
      int epos = phrase.j1 - 1 - i;
      const WordID left_ctx = (epos < 0) ? kBOS : sentence.e[epos];
      context[ctx_size - i - 1] = left_ctx;
      epos = phrase.j2 + i;
      const WordID right_ctx = (epos >= sentence.e_len) ? kEOS : sentence.e[epos];
      context[ctx_size + i + 1] = right_ctx;
    }
    cerr << TD::GetString(context) << endl;
  }
}

struct Rule {
  vector<ParallelSpan> f;
  int i,j,syms,vars;
  explicit Rule(int pi) : i(pi), j(pi), syms(), vars() {}
  void Extend(const WordID& fword) {
    f.push_back(ParallelSpan(fword));
    ++j;
    ++syms;
  }
  void Extend(const ParallelSpan& subphrase) {
    f.push_back(subphrase);
    j += subphrase.i2 - subphrase.i1;
    ++vars;
    ++syms;
  }
  bool RuleFEndsInVariable() const {
    if (f.size() > 0) {
      return f.back().IsVariable();
    } else { return false; }
  }
};

ostream& operator<<(ostream& os, const Rule& r) {
  os << "(" << r.i << "," << r.j << ") ";
  for (int i = 0; i < r.f.size(); ++i) {
    const ParallelSpan& x = r.f[i];
    if (x.IsVariable()) { os << "[" << TD::Convert(-x.cat) << "] "; } else
      os << TD::Convert(x.cat) << ' ';
  }
  return os;
}

void ExtractAndWriteRules(const AnnotatedParallelSentence& sentence,
                          const vector<ParallelSpan>& phrases,
                          const int max_vars,
                          const int max_syms) {
  unordered_map<pair<short, short>, vector<ParallelSpan>, boost::hash<pair<int, int> > > fspans;
  int max_len = -1;  // remove?
  set<int> starts;
  queue<Rule> q;
  vector<vector<ParallelSpan> > spans_by_start(sentence.f_len);
  for (int i = 0; i < phrases.size(); ++i) {
    fspans[make_pair(phrases[i].i1,phrases[i].i2)].push_back(phrases[i]);
    max_len = max(max_len, phrases[i].i2 - phrases[i].i1);
    if (starts.insert(phrases[i].i1).second)
      q.push(Rule(phrases[i].i1));
    spans_by_start[phrases[i].i1].push_back(phrases[i]);
  }
  // cerr << "MAX PHRASE: " << max_len << endl;
  bool adjacent_vars_permitted = false;
  vector<pair<int,int> > next_e(sentence.e_len);
  while(!q.empty()) {
    const Rule& rule = q.front();

    // extend the partial rule
    if (rule.j < sentence.f_len && (rule.j - rule.i) < max_len && rule.syms < max_syms) {
      Rule ew = rule;

      // extend with a word
      ew.Extend(sentence.f[ew.j]);
      q.push(ew);

      // with variables
      if (rule.vars < max_vars &&
          !spans_by_start[rule.j].empty() &&
          ((!rule.RuleFEndsInVariable()) || adjacent_vars_permitted)) {
        const vector<ParallelSpan>& sub_phrases = spans_by_start[rule.j];
        for (int it = 0; it < sub_phrases.size(); ++it) {
          if (sub_phrases[it].i2 - sub_phrases[it].i1 + rule.j - rule.i <= max_len) {
            Rule ev = rule;
            ev.Extend(sub_phrases[it]);
            q.push(ev);
            assert(ev.j <= sentence.f_len);
          }
        }
      }
    }
    // determine if rule is consistent
    if (rule.syms > 0 &&
        fspans.count(make_pair(rule.i,rule.j)) &&
        (!rule.RuleFEndsInVariable() || rule.syms > 1)) {
      const vector<ParallelSpan>& orig_spans = fspans[make_pair(rule.i,rule.j)];
      for (int s = 0; s < orig_spans.size(); ++s) {
        const ParallelSpan& orig_span = orig_spans[s];
        const WordID lhs = orig_span.cat;
        for (int j = orig_span.j1; j < orig_span.j2; ++j) next_e[j].first = -1;
        int nt_index_e = 0;
        for (int i = 0; i < rule.f.size(); ++i) {
          const ParallelSpan& cur = rule.f[i];
          if (cur.IsVariable())
            next_e[cur.j1] = pair<int,int>(cur.j2, ++nt_index_e);
        }
        cout << '[' << TD::Convert(-lhs) << "] |||";
        int nt_index_f = 0;
        for (int i = 0; i < rule.f.size(); ++i) {
          const ParallelSpan& cur = rule.f[i];
          if (cur.IsVariable()) {
            ++nt_index_f;
            cout << " [" << TD::Convert(-cur.cat) << ',' << nt_index_f << ']';
          } else {
            cout << ' ' << TD::Convert(cur.cat);
          }
        }
        cout << " |||";
        for (int j = orig_span.j1; j < orig_span.j2; ++j) {
          if (next_e[j].first < 0) {
            cout << ' ' << TD::Convert(sentence.e[j]);
          } else {
            cout << " [" << next_e[j].second << ']';
            j = next_e[j].first - 1;
          }
        }
        cout << endl;
      }
    }
    q.pop();
  }
}

// this uses the TARGET span (i,j) to annotate phrases, will copy
// phrases if there is more than one annotation.
// TODO: support source annotation
void AnnotatePhrasesWithCategoryTypes(const WordID default_cat,
                                      const Array2D<vector<WordID> >& types,
                                      vector<ParallelSpan>* phrases) {
  const int num_phrases = phrases->size();
  // have to use num_phrases since we may grow the size of phrases
  for (int i = 0; i < num_phrases; ++i) {
    ParallelSpan& phrase = (*phrases)[i];
    const vector<WordID>* pcats = &types(phrase.j1, phrase.j2);
    if (pcats->empty() && default_cat != 0) {
      static vector<WordID> s_default(1, default_cat);
      pcats = &s_default;
    }
    if (pcats->empty()) {
      cerr << "ERROR span " << phrase.i1 << "," << phrase.i2 << "-"
           << phrase.j1 << "," << phrase.j2 << " has no type. "
              "Did you forget --default_category?\n";
    }
    const vector<WordID>& cats = *pcats;
    phrase.cat = cats[0];
    for (int ci = 1; ci < cats.size(); ++ci) {
      ParallelSpan new_phrase = phrase;
      new_phrase.cat = cats[ci];
      phrases->push_back(new_phrase);
    }
  }
}

int main(int argc, char** argv) {
  po::variables_map conf;
  InitCommandLine(argc, argv, &conf);
  WordID default_cat = 0;  // 0 means no default- extraction will
                           // fail if a phrase is extracted without a
                           // category
  kBOS = TD::Convert("<s>");
  kEOS = TD::Convert("</s>");
  kDIVIDER = TD::Convert("|||");
  if (conf.count("default_category")) {
    string sdefault_cat = conf["default_category"].as<string>();
    default_cat = -TD::Convert(sdefault_cat);
    cerr << "Default category: " << sdefault_cat << endl;
  } else {
    cerr << "No default category (use --default_category if you want to set one)\n";
  }
  ReadFile rf(conf["input"].as<string>());
  istream& in = *rf.stream();

  char buf[MAX_LINE_LENGTH];
  AnnotatedParallelSentence sentence;
  vector<ParallelSpan> phrases;
  const int max_base_phrase_size = conf["max_base_phrase_size"].as<int>();
  const bool write_phrase_contexts = conf.count("phrase_context") > 0;
  const bool write_base_phrases = conf.count("base_phrase") > 0;
  const bool loose_phrases = conf.count("loose") > 0;
  const int ctx_size = conf["phrase_context_size"].as<int>();
  while(in) {
    in.getline(buf, MAX_LINE_LENGTH);
    if (buf[0] == 0) continue;
    sentence.ParseInputLine(buf);
    phrases.clear();
    ExtractBasePhrases(max_base_phrase_size, sentence, &phrases);
    if (loose_phrases)
      LoosenPhraseBounds(sentence, &phrases);
    if (phrases.empty()) {
      cerr << "WARNING no phrases extracted\n";
      continue;
    }
    if (write_phrase_contexts) {
      WritePhraseContexts(sentence, phrases, ctx_size);
      continue;
    }
    if (write_base_phrases) {
      WriteBasePhrases(sentence, phrases);
      continue;
    }
    AnnotatePhrasesWithCategoryTypes(default_cat, sentence.span_types, &phrases);
    int max_vars = 2;
    int max_syms = 5;
    ExtractAndWriteRules(sentence, phrases, max_vars, max_syms);
  }
  return 0;
}

