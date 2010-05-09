#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>

#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

#include "sentence_pair.h"
#include "fdict.h"

using namespace std;
namespace po = boost::program_options;

static const size_t MAX_LINE_LENGTH = 64000000;

void InitCommandLine(int argc, char** argv, po::variables_map* conf) {
  po::options_description opts("Configuration options");
  opts.add_options()
        ("normalize,n", "Normalize counts and write log relative frequencies")
        ("help,h", "Print this help message and exit");
  po::options_description clo("Command line options");
  po::options_description dcmdline_options;
  dcmdline_options.add(opts);

  po::store(parse_command_line(argc, argv, dcmdline_options), *conf);
  po::notify(*conf);

  if (conf->count("help")) {
    cerr << "\nUsage: mr_stripe_rule_reduce [-options]\n";
    cerr << dcmdline_options << endl;
    exit(1);
  }
}

namespace {
  inline bool IsWhitespace(char c) { return c == ' ' || c == '\t'; }

  inline void SkipWhitespace(const char* buf, int* ptr) {
    while (buf[*ptr] && IsWhitespace(buf[*ptr])) { ++(*ptr); }
  }
}

struct CountAndAlignment {
  CountAndAlignment() : count() {}
  double count;
  vector<pair<short,short> > aligns;
  CountAndAlignment& operator+=(const CountAndAlignment& rhs) {
    if (rhs.count > count) { aligns = rhs.aligns; }
    count += rhs.count;
    return *this;
  }
};

typedef map<int, CountAndAlignment> ID2CountAndAlignment;

void PlusEquals(const ID2CountAndAlignment& v, ID2CountAndAlignment* self) {
  for (ID2CountAndAlignment::const_iterator it = v.begin(); it != v.end(); ++it) {
    (*self)[it->first] += it->second;
  }
}

void ParseCountAndAlignment(const char* buf,
                            const int start,
                            const int end,
                            const bool has_alignment,
                            CountAndAlignment* canda) {
  canda->count = strtod(buf+start,NULL);
  canda->aligns.clear();
  if (has_alignment) {
    int ptr = start;
    while (buf[ptr] != '|') ++ptr;
    ++ptr;
    int cs;
    while(buf[ptr] != 0 && ptr < end) {
      SkipWhitespace(buf, &ptr);
      cs = ptr;
      while(ptr < end && !IsWhitespace(buf[ptr])) { ++ptr; }
      if (ptr > cs) {
        short a, b;
        AnnotatedParallelSentence::ReadAlignmentPoint(buf, cs, ptr, false, &a, &b);
        canda->aligns.push_back(make_pair(a,b));
      }
    }
  }
}

void ParseLine(const char* buf, int* cur_key, ID2CountAndAlignment* counts) {
  counts->clear();
  int ptr = 0;
  while(buf[ptr] != 0 && buf[ptr] != '\t') { ++ptr; }
  if (buf[ptr] != '\t') {
    cerr << "Missing tab separator between key and value!\n INPUT=" << buf << endl;
    exit(1);
  }
  *cur_key = FD::Convert(string(buf,0,ptr));
  ++ptr;
  int start = ptr;
  int end = ptr;
  int state = 0; // 0=reading label, 1=reading count
  int name = 0;
  bool seen_bar = false;
  while(buf[ptr] != 0) {
    while(buf[ptr] != 0 && buf[ptr] != '|') { ++ptr; }
    if (buf[ptr] == '|') {
      ++ptr;
      if (buf[ptr] != '|') {
        seen_bar = true;
      } else {
        ++ptr;
        if (buf[ptr] == '|') {
          ++ptr;
          end = ptr - 3;
          while (end > start && IsWhitespace(buf[end-1])) { --end; }
          if (start == end) {
            cerr << "Got empty token!\n  LINE=" << buf << endl;
            exit(1);
          }
          switch (state) {
            case 0: ++state; name=FD::Convert(string(buf,start,end-start)); break;
            case 1: --state; ParseCountAndAlignment(buf, start, end, seen_bar, &(*counts)[name]); break;
            default: cerr << "Can't happen\n"; abort();
          }
          SkipWhitespace(buf, &ptr);
          start = ptr;
          seen_bar = false;
        }
      }
    }
  }
  end=ptr;
  while (end > start && IsWhitespace(buf[end-1])) { --end; }
  if (end > start) {
    switch (state) {
      case 0: ++state; name=FD::Convert(string(buf,start,end-start)); break;
      case 1: --state; ParseCountAndAlignment(buf, start, end, seen_bar, &(*counts)[name]); break;
      default: cerr << "Can't happen\n"; abort();
    }
  }
}

void WriteValues(int key, const ID2CountAndAlignment& val) {
  const string& skey = FD::Convert(key);
  const bool emajor = (skey[0] == 'E');
  if (emajor) {
    const size_t divp = skey.find(" |||");
    const string lhs = skey.substr(2, divp - 2);
    const string rhs_e = skey.substr(divp + 5);
    for (ID2CountAndAlignment::const_iterator it = val.begin();
         it != val.end(); ++it) {
      cout << lhs << " ||| " << FD::Convert(it->first) << '\t' << rhs_e << " ||| FgivenE=" << it->second.count << endl;
    }
  } else {
    cout << skey.substr(2) << "\t";
    bool print_bar = false;
    for (ID2CountAndAlignment::const_iterator it = val.begin();
         it != val.end(); ++it) {
      if (print_bar) cout << " ||| "; else print_bar = true;
      cout << FD::Convert(it->first) << " ||| EgivenF=" << it->second.count;
      const vector<pair<short,short> >& aligns = it->second.aligns;
      if (aligns.size() > 0) {
        cout << " A=";
        for (int i = 0; i < aligns.size(); ++i) {
          if (i) cout << ',';
          cout << aligns[i].first << '-' << aligns[i].second;
        }
      }
    }
    cout << endl;
  }
}

void Normalize(ID2CountAndAlignment* v) {
  double sum = 0;
  for (ID2CountAndAlignment::iterator it = v->begin();
       it != v->end(); ++it)
    sum += it->second.count;
  if (!sum) return;

  sum = log(sum);
  for (ID2CountAndAlignment::iterator it = v->begin();
       it != v->end(); ++it)
    it->second.count = sum - log(it->second.count);
}

int main(int argc, char** argv) {
  po::variables_map conf;
  InitCommandLine(argc, argv, &conf);

  char* buf = new char[MAX_LINE_LENGTH];
  ID2CountAndAlignment acc, cur_counts;
  int key = -1, cur_key = -1;
  int line = 0;
  const bool normalize = conf.count("normalize") > 0;
  while(cin) {
    ++line;
    cin.getline(buf, MAX_LINE_LENGTH);
    if (buf[0] == 0) continue;
    ParseLine(buf, &cur_key, &cur_counts);
    if (cur_key != key) {
      if (key > 0) {
        if (normalize) Normalize(&acc);
        WriteValues(key, acc);
        acc.clear();
        if (FD::NumFeats() > 400000) {
          FD::dict_.clear();
          ParseLine(buf, &cur_key, &cur_counts);
        }
      }
      key = cur_key;
    }
    PlusEquals(cur_counts, &acc);
  }
  if (key > 0) {
    if (normalize) Normalize(&acc);
    WriteValues(key, acc);
  }
  return 0;
}

