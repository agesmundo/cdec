#include "trule.h"

#include <sstream>

#include "stringlib.h"
#include "tdict.h"

using namespace std;

static WordID ConvertTrgString(const string& w) {
  int len = w.size();
  WordID id = 0;
  // [X,0] or [0]
  // for target rules, we ignore the category, just keep the index
  if (len > 2 && w[0]=='[' && w[len-1]==']' && w[len-2] > '0' && w[len-2] <= '9' &&
      (len == 3 || (len > 4 && w[len-3] == ','))) {
    id = w[len-2] - '0';
    id = 1 - id;
  } else {
    id = TD::Convert(w);
  }
  return id;
}

static WordID ConvertSrcString(const string& w, bool mono = false) {
  int len = w.size();
  // [X,0]
  // for source rules, we keep the category and ignore the index (source rules are
  // always numbered 1, 2, 3...
  if (mono) {
    if (len > 2 && w[0]=='[' && w[len-1]==']') {
      if (len > 4 && w[len-3] == ',') {
        cerr << "[ERROR] Monolingual rules mut not have non-terminal indices:\n  "
             << w << endl;
        exit(1);
      }
      // TODO check that source indices go 1,2,3,etc.
      return TD::Convert(w.substr(1, len-2)) * -1;
    } else {
      return TD::Convert(w);
    }
  } else {
    if (len > 4 && w[0]=='[' && w[len-1]==']' && w[len-3] == ',' && w[len-2] > '0' && w[len-2] <= '9') {
      return TD::Convert(w.substr(1, len-4)) * -1;
    } else {
      return TD::Convert(w);
    }
  }
}

static WordID ConvertLHS(const string& w) {
  if (w[0] == '[') {
    int len = w.size();
    if (len < 3) { cerr << "Format error: " << w << endl; exit(1); }
    return TD::Convert(w.substr(1, len-2)) * -1;
  } else {
    return TD::Convert(w) * -1;
  }
}

TRule* TRule::CreateRuleSynchronous(const std::string& rule) {
  TRule* res = new TRule;
  if (res->ReadFromString(rule, true, false)) return res;
  cerr << "[ERROR] Failed to creating rule from: " << rule << endl;
  delete res;
  return NULL;
}

TRule* TRule::CreateRulePhrasetable(const string& rule) {
  // TODO make this faster
  // TODO add configuration for default NT type
  if (rule[0] == '[') {
    cerr << "Phrasetable rules shouldn't have a LHS / non-terminals:\n  " << rule << endl;
    return NULL;
  }
  TRule* res = new TRule("[X] ||| " + rule, true, false);
  if (res->Arity() != 0) {
    cerr << "Phrasetable rules should have arity 0:\n  " << rule << endl;
    delete res;
    return NULL;
  }
  return res;
}

TRule* TRule::CreateRuleMonolingual(const string& rule) {
  return new TRule(rule, false, true);
}

bool TRule::ReadFromString(const string& line, bool strict, bool mono) {
  e_.clear();
  f_.clear();
  scores_.clear();

  string w;
  istringstream is(line);
  int format = CountSubstrings(line, "|||");
  if (strict && format < 2) {
    cerr << "Bad rule format in strict mode:\n" << line << endl;
    return false;
  }
  if (format >= 2 || (mono && format == 1)) {
    while(is>>w && w!="|||") { lhs_ = ConvertLHS(w); }
    while(is>>w && w!="|||") { f_.push_back(ConvertSrcString(w, mono)); }
    if (!mono) {
      while(is>>w && w!="|||") { e_.push_back(ConvertTrgString(w)); }
    }
    int fv = 0;
    if (is) {
      string ss;
      getline(is, ss);
      //cerr << "L: " << ss << endl;
      int start = 0;
      const int len = ss.size();
      while (start < len) {
        while(start < len && (ss[start] == ' ' || ss[start] == ';'))
          ++start;
        if (start == len) break;
        int end = start + 1;
        while(end < len && (ss[end] != '=' && ss[end] != ' ' && ss[end] != ';'))
          ++end;
        if (end == len || ss[end] == ' ' || ss[end] == ';') {
          //cerr << "PROC: '" << ss.substr(start, end - start) << "'\n";
          // non-named features
          if (end != len) { ss[end] = 0; }
          string fname = "PhraseModel_X";
          if (fv > 9) { cerr << "Too many phrasetable scores - used named format\n"; abort(); }
          fname[12]='0' + fv;
          ++fv;
          // if the feature set is frozen, this may return zero, indicating an
          // undefined feature
          const int fid = FD::Convert(fname);
          if (fid)
            scores_.set_value(fid, atof(&ss[start]));
          //cerr << "F: " << fname << " VAL=" << scores_.value(FD::Convert(fname)) << endl;
        } else {
          const int fid = FD::Convert(ss.substr(start, end - start));
          start = end + 1;
          end = start + 1;
          while(end < len && (ss[end] != ' ' && ss[end] != ';'))
            ++end;
          if (end < len) { ss[end] = 0; }
	  assert(start < len);
          if (fid)
            scores_.set_value(fid, atof(&ss[start]));
          //cerr << "F: " << FD::Convert(fid) << " VAL=" << scores_.value(fid) << endl;
        }
        start = end + 1;
      }
    }
  } else if (format == 1) {
    while(is>>w && w!="|||") { lhs_ = ConvertLHS(w); }
    while(is>>w && w!="|||") { e_.push_back(ConvertTrgString(w)); }
    f_ = e_;
    int x = ConvertLHS("[X]");
    for (int i = 0; i < f_.size(); ++i)
      if (f_[i] <= 0) { f_[i] = x; }
  } else {
    cerr << "F: " << format << endl;
    cerr << "[ERROR] Don't know how to read:\n" << line << endl;
  }
  if (mono) {
    e_ = f_;
    int ci = 0;
    for (int i = 0; i < e_.size(); ++i)
      if (e_[i] < 0)
        e_[i] = ci--;
  }
  ComputeArity();
  return SanityCheck();
}

bool TRule::SanityCheck() const {
  vector<int> used(f_.size(), 0);
  int ac = 0;
  for (int i = 0; i < e_.size(); ++i) {
    int ind = e_[i];
    if (ind > 0) continue;
    ind = -ind;
    if ((++used[ind]) != 1) {
      cerr << "[ERROR] e-side variable index " << (ind+1) << " used more than once!\n";
      return false;
    }
    ac++;
  }
  if (ac != Arity()) {
    cerr << "[ERROR] e-side arity mismatches f-side\n";
    return false;
  }
  return true;
}

void TRule::ComputeArity() {
  int min = 1;
  for (vector<WordID>::const_iterator i = e_.begin(); i != e_.end(); ++i)
    if (*i < min) min = *i;
  arity_ = 1 - min;
}

static string AnonymousStrVar(int i) {
  string res("[v]");
  if(!(i <= 0 && i >= -8)) {
    cerr << "Can't handle more than 9 non-terminals: index=" << (-i) << endl;
    abort();
  }
  res[1] = '1' - i;
  return res;
}

string TRule::AsString(bool verbose) const {
  ostringstream os;
  int idx = 0;
  if (lhs_ && verbose) {
    os << '[' << TD::Convert(lhs_ * -1) << "] |||";
    for (int i = 0; i < f_.size(); ++i) {
      const WordID& w = f_[i];
      if (w < 0) {
        int wi = w * -1;
        ++idx;
        os << " [" << TD::Convert(wi) << ',' << idx << ']';
      } else {
        os << ' ' << TD::Convert(w);
      }
    }
    os << " ||| ";
  }
  if (idx > 9) {
    cerr << "Too many non-terminals!\n partial: " << os.str() << endl;
    exit(1);
  }
  for (int i =0; i<e_.size(); ++i) {
    if (i) os << ' ';
    const WordID& w = e_[i];
    if (w < 1)
      os << AnonymousStrVar(w);
    else
      os << TD::Convert(w);
  }
  if (!scores_.empty() && verbose) {
    os << " ||| " << scores_;
  }
  return os.str();
}

ostream& operator<<(ostream& os, const TRule& rule){
	return os << rule.AsString();
}
