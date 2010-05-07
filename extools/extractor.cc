#include <iostream>
#include <vector>
#include <utility>

#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include "sentence_pair.h"
#include "extract.h"
#include "tdict.h"
#include "wordid.h"
#include "array2d.h"
#include "filelib.h"

using namespace std;
namespace po = boost::program_options;

static const size_t MAX_LINE_LENGTH = 100000;
WordID kBOS, kEOS, kDIVIDER;

void InitCommandLine(int argc, char** argv, po::variables_map* conf) {
  po::options_description opts("Configuration options");
  opts.add_options()
        ("input,i", po::value<string>()->default_value("-"), "Input file")
        ("default_category,d", po::value<string>(), "Default span type (use X for 'Hiero')")
        ("loose", "Use loose phrase extraction heuristic for base phrases")
        ("invert,I", "Invert the outputs")
        ("base_phrase,B", "Write base phrases")
        ("silent", "Write nothing to stderr except errors")
        ("phrase_context,C", "Write base phrase contexts")
        ("phrase_context_size,S", po::value<int>()->default_value(2), "Use this many words of context on left and write when writing base phrase contexts")
        ("max_base_phrase_size,L", po::value<int>()->default_value(10), "Maximum starting phrase size")
        ("max_syms,l", po::value<int>()->default_value(5), "Maximum number of symbols in final phrase size")
        ("max_vars,v", po::value<int>()->default_value(2), "Maximum number of nonterminal variables in final phrase size")
        ("permit_adjacent_nonterminals,A", "Permit adjacent nonterminals in source side of rules")
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

struct SimpleRuleWriter : public Extract::RuleObserver {
  virtual void CountRule(WordID lhs,
                         const vector<WordID>& rhs_f,
                         const vector<WordID>& rhs_e,
                         const vector<pair<int, int> >& fe_terminal_alignments) {
    cout << "[" << TD::Convert(-lhs) << "] |||";
    for (int i = 0; i < rhs_f.size(); ++i) {
      if (rhs_f[i] < 0) cout << " [" << TD::Convert(-rhs_f[i]) << ']';
      else cout << ' ' << TD::Convert(rhs_f[i]);
    }
    cout << " |||";
    for (int i = 0; i < rhs_e.size(); ++i) {
      if (rhs_e[i] <= 0) cout << " [" << (1-rhs_e[i]) << ']';
      else cout << ' ' << TD::Convert(rhs_e[i]);
    }
    cout << endl;
  }
};

int main(int argc, char** argv) {
  po::variables_map conf;
  InitCommandLine(argc, argv, &conf);
  kBOS = TD::Convert("<s>");
  kEOS = TD::Convert("</s>");
  kDIVIDER = TD::Convert("|||");

  WordID default_cat = 0;  // 0 means no default- extraction will
                           // fail if a phrase is extracted without a
                           // category
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
  const bool silent = conf.count("silent") > 0;
  const int max_syms = conf["max_syms"].as<int>();
  const int max_vars = conf["max_vars"].as<int>();
  const int ctx_size = conf["phrase_context_size"].as<int>();
  const bool permit_adjacent_nonterminals = conf.count("permit_adjacent_nonterminals") > 0;
  int line = 0;
  SimpleRuleWriter o;
  while(in) {
    ++line;
    in.getline(buf, MAX_LINE_LENGTH);
    if (buf[0] == 0) continue;
    if (!silent) {
      if (line % 200 == 0) cerr << '.';
      if (line % 8000 == 0) cerr << " [" << line << "]\n" << flush;
    }
    sentence.ParseInputLine(buf);
    phrases.clear();
    Extract::ExtractBasePhrases(max_base_phrase_size, sentence, &phrases);
    if (loose_phrases)
      Extract::LoosenPhraseBounds(sentence, &phrases);
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
    Extract::AnnotatePhrasesWithCategoryTypes(default_cat, sentence.span_types, &phrases);
    Extract::ExtractConsistentRules(sentence, phrases, max_vars, max_syms, permit_adjacent_nonterminals, &o);
  }
  if (!silent) cerr << endl;
  return 0;
}

