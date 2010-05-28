/*
 * Build lexical translation table from alignment file to use for lexical translation probabilties when scoring a grammar
 *
 * Ported largely from the train-factored-phrase-model.perl script by Philipp Koehn
 */
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <cstdlib>
#include <fstream>
#include <tr1/unordered_map>

#include "sentence_pair.h"
#include "extract.h"
#include "fdict.h"
#include "tdict.h"

#include <boost/functional/hash.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>


using namespace std;
using namespace std::tr1;

static const size_t MAX_LINE_LENGTH = 64000000;

typedef unordered_map<vector<WordID>, RuleStatistics, boost::hash<vector<WordID> > > ID2RuleStatistics;


namespace {
  inline bool IsWhitespace(char c) { return c == ' ' || c == '\t'; }
  inline bool IsBracket(char c){return c == '[' || c == ']';}
  inline void SkipWhitespace(const char* buf, int* ptr) {
    while (buf[*ptr] && IsWhitespace(buf[*ptr])) { ++(*ptr); }
  }
}

int ReadPhraseUntilDividerOrEnd(const char* buf, const int sstart, const int end, vector<WordID>* p) {
  static const WordID kDIV = TD::Convert("|||");
  int ptr = sstart;
  while(ptr < end) {
    while(ptr < end && IsWhitespace(buf[ptr])) { ++ptr; }
    int start = ptr;
    while(ptr < end && !IsWhitespace(buf[ptr])) { ++ptr; }
    if (ptr == start) {cerr << "Warning! empty token.\n"; return ptr; }
    const WordID w = TD::Convert(string(buf, start, ptr - start));
    if (w == kDIV) return ptr;
    p->push_back(w);
  }
  return ptr;
}


void ParseLine(const char* buf, vector<WordID>* cur_key, ID2RuleStatistics* counts) {
  static const WordID kDIV = TD::Convert("|||");
  counts->clear();
  int ptr = 0;
  while(buf[ptr] != 0 && buf[ptr] != '\t') { ++ptr; }
  if (buf[ptr] != '\t') {
    cerr << "Missing tab separator between key and value!\n INPUT=" << buf << endl;
    exit(1);
  }
  cur_key->clear();
  // key is: "[X] ||| word word word"
  int tmpp = ReadPhraseUntilDividerOrEnd(buf, 0, ptr, cur_key);
  cur_key->push_back(kDIV);
  ReadPhraseUntilDividerOrEnd(buf, tmpp, ptr, cur_key);
  ++ptr;
  int start = ptr;
  int end = ptr;
  int state = 0; // 0=reading label, 1=reading count
  vector<WordID> name;
  while(buf[ptr] != 0) {
    while(buf[ptr] != 0 && buf[ptr] != '|') { ++ptr; }
    if (buf[ptr] == '|') {
      ++ptr;
      if (buf[ptr] == '|') {
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
	  case 0: ++state; name.clear(); ReadPhraseUntilDividerOrEnd(buf, start, end, &name); break;
	  case 1: --state; (*counts)[name].ParseRuleStatistics(buf, start, end); break;
	  default: cerr << "Can't happen\n"; abort();
          }
          SkipWhitespace(buf, &ptr);
          start = ptr;
        }
      }
    }
  }
  end=ptr;
  while (end > start && IsWhitespace(buf[end-1])) { --end; }
  if (end > start) {
    switch (state) {
    case 0: ++state; name.clear(); ReadPhraseUntilDividerOrEnd(buf, start, end, &name); break;
    case 1: --state; (*counts)[name].ParseRuleStatistics(buf, start, end); break;
    default: cerr << "Can't happen\n"; abort();
    }
  }
}

int main(int argc, char* argv[]){

  bool DEBUG = false;

  map < pair<WordID,WordID>,int > word_translation;
  map <WordID, int> total_foreign;
  map <WordID, int> total_english;

  AnnotatedParallelSentence sent;
  char* buf = new char[MAX_LINE_LENGTH];
  while(cin) 
    {
      cin.getline(buf, MAX_LINE_LENGTH);
      if (buf[0] == 0) continue;
      
      sent.ParseInputLine(buf);
      
      map <WordID, int> foreign_aligned;
      map <WordID, int> english_aligned;

      //iterate over the alignment to compute aligned words
            
      for(int i =0;i<sent.aligned.width();i++)
	{
	  for (int j=0;j<sent.aligned.height();j++)
	    {
	      if (DEBUG) cout << sent.aligned(i,j) << " ";
	      if( sent.aligned(i,j))
		{
		  if (DEBUG) cout << TD::Convert(sent.f[i])  << " aligned to " << TD::Convert(sent.e[j]);
		  //local counts
		  ++foreign_aligned[sent.f[i]];
		  ++english_aligned[sent.e[j]];

		  //global counts
		  ++word_translation[pair<WordID,WordID> (sent.f[i], sent.e[j])];
		  ++total_foreign[sent.f[i]];
		  ++total_english[sent.e[j]];
		}
	    }
	  if (DEBUG)  cout << endl;
	}
      if (DEBUG) cout << endl;
      
      static const WordID NULL_ = TD::Convert("NULL");
      //handle unaligned words - align them to null
      for (int j =0; j < sent.e_len; j++)
	{
	  if (english_aligned.count(sent.e[j])) continue;
	  ++word_translation[pair<WordID,WordID> (NULL_, sent.e[j])];
	  ++total_foreign[NULL_];
	  ++total_english[sent.e[j]];
	}

      for (int i =0; i < sent.f_len; i++)
	{
	  if (foreign_aligned.count(sent.f[i])) continue;
	  ++word_translation[pair<WordID,WordID> (sent.f[i], NULL_)];
	  ++total_english[NULL_];
	  ++total_foreign[sent.f[i]];
	}
      
    }

  for(map < pair<WordID,WordID>,int >::iterator it = word_translation.begin(); it != word_translation.end(); ++it)
    {
      cout <<  TD::Convert(it->first.first) <<  "," << TD::Convert(it->first.second) << "=" << it->second << "/" << total_foreign[it->first.first] << endl;
    }


  return 0;
}
