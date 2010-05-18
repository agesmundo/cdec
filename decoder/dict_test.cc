#include "dict.h"

#include "fdict.h"

#include <iostream>
#include <gtest/gtest.h>
#include <cassert>
#include "filelib.h"

#include "tdict.h"
#include <google/sparse_hash_map>

using namespace google;
using namespace std;

class DTest : public testing::Test {
 public:
  DTest() {}
 protected:
  virtual void SetUp() { }
  virtual void TearDown() { }
};

struct LMNode {
  LMNode(size_t x) : ptr_(new sparse_hash_map<WordID, LMNode>(x)) {}
  LMNode() : ptr_(NULL) {
  //  ptr_.set_empty_key(0);
  }
  float prob_;
  float bo_prob_;
  sparse_hash_map<WordID, LMNode>* ptr_;
};

// -4.294035       ! teulade       -0.2394543
// -3.759696       ! thank

inline void ConvertARPA(const string& line, const int n, vector<WordID>* gram, float* prob, float* bo) {
  int p = 1;
  const size_t size = line.size();
  while(p < size && line[p] != '\t' && line[p] != ' ') { ++p; }
  *prob = strtof(&line[0], NULL);
  int wi = 0;
  int wstart = p+1;
  while(wi < n) {
    wstart = ++p;
    while(p < size && line[p] != ' ' && line[p] != '\t') { ++p; }
    (*gram)[wi] = TD::Convert(line.substr(wstart, p - wstart));
    ++wi;
  }
  if (p < size) {
    *bo = strtof(&line[p+1], NULL);
  }
}

inline void AddGram(const vector<WordID>& gram, LMNode* root, const float prob, float bo) {
  LMNode* cur = root;
  for (int i = 0; i < gram.size(); ++i) {
    sparse_hash_map<WordID, LMNode>*& ptr = cur->ptr_;
    if (!ptr) {
      ptr = new sparse_hash_map<WordID, LMNode>;
    }
    cur = &(*ptr)[gram[i]];
  }
  cur->prob_ = prob;
  cur->bo_prob_ = bo;
}

TEST_F(DTest, MyLM) {
  ReadFile rf("/Users/redpony/cdyer-svn-root/cdec/src/test_data/c2e.3gram.lm.gz");
  istream& in = *rf.stream();
  string x;
  int lc = 0;
  while (x.empty()) { getline(in, x); ++lc; }
  assert(x == "\\data\\");
  vector<unsigned long long int> counts(1,0);
  while(1){
    getline(in, x); ++lc;
    if (x.empty()) break;
    if (x[0]=='n' && x[1]=='g' && x[2]=='r' && x[3]=='a' && x[4]=='m' && x[5]==' ' && x[7]=='=') {
      int order = x[6] - '0';
      assert(order == counts.size());
      counts.push_back(strtoll(&x[8], NULL, 10));
      cerr << "NUM " << order << " : " << counts.back() << endl;
    }
  }
  const int order = counts.size() - 1;
  cerr << "ORDER = " << order << endl;
  LMNode root(counts[1]);
  const float fmark = -999999.0f;
  for (int n = 1; n <= order; ++n) {
    x.clear();
    while(x.empty() && in) { getline(in, x); ++lc; }
    string expected = "\\n-grams:";
    expected[1] = n + '0';
    if (x != expected) {
      cerr << "FORMAT ERROR: " << lc << ": expected " << expected << " got " << x << endl;
      abort();
    }
    const unsigned long long num_exp = counts[n];
    cerr << "Reading " << num_exp << " " << n << "-grams...\n";
    vector<WordID> gram(n);
    for (unsigned long long int c = 0; c < num_exp; ++c) {
      getline(in, x); ++lc;
      if (x.empty()) {
        if (!in) { cerr << "Unexpected end of file!\n"; abort(); }
        cerr << "Unexpected blank line!\n"; abort();
      }
      float prob;
      float bo = fmark;
      ConvertARPA(x, n, &gram, &prob, &bo);
      AddGram(gram, &root, prob, bo);
    }
  }
}

TEST_F(DTest, Convert) {
  Dict d;
  WordID a = d.Convert("foo");
  WordID b = d.Convert("bar");
  std::string x = "foo";
  WordID c = d.Convert(x);
  EXPECT_NE(a, b);
  EXPECT_EQ(a, c);
  EXPECT_EQ(d.Convert(a), "foo");
  EXPECT_EQ(d.Convert(b), "bar");
}

TEST_F(DTest, FDictTest) {
  int fid = FD::Convert("First");
  EXPECT_GT(fid, 0);
  EXPECT_EQ(FD::Convert(fid), "First");
  string x = FD::Escape("=");
  cerr << x << endl;
  EXPECT_NE(x, "=");
  x = FD::Escape(";");
  cerr << x << endl;
  EXPECT_NE(x, ";");
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

