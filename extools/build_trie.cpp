/*
 * Build suffix trie representation of a data set for grammar filtering
 * ./build_trie <test set> <rule-table> <output filtered grammar>
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
#include "suffix_trie.h"
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
    //look in the buffer and see if its a nonterminal marker before integerizing it to wordID-anything with [...] or |||

    const WordID w = TD::Convert(string(buf, start, ptr - start));

    if((IsBracket(buf[start]) and IsBracket(buf[ptr-1])) or( w == kDIV))
      p->push_back(-1);
    else {
	if (w == kDIV) return ptr;
	p->push_back(w);
    }
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

  ifstream testSet (argv[1]);
  ofstream filter_grammar_;
  bool DEBUG = false;


  AnnotatedParallelSentence sent;
  char* buf = new char[MAX_LINE_LENGTH];
  cout << "Build suffix trie..." << endl;
  //root of the suffix trie
  Node<int> *root = new Node<int>;
  root->label_=-1;
  int line=0;

  /* process the data set to build trie
   */
  while(!testSet.eof())
    //  while(cin) 
    {
    ++line;
    testSet.getline(buf, MAX_LINE_LENGTH);
    if (buf[0] == 0) continue;

    //hack to read in the test set using the alignedparallelsentence methods
    strcat(buf," ||| fake ||| 0-0");   
    sent.ParseInputLine(buf);

    if (DEBUG)cout << line << "||| " << buf << " -- " << sent.f_len << endl;

    //add each successive suffix to the tree
    for(int i =0;i<sent.f_len;i++)
      {
	root->InsertEdge(sent.f, i, sent.f_len - 1);
      }
    if(DEBUG)cout<<endl;

  }  

  cout << "Filter grammar..." << endl;
  //process the unfiltered, unscored grammar    
  ifstream grammar (argv[2]);
  filter_grammar_.open(argv[3]);

  ID2RuleStatistics acc, cur_counts;
  vector<WordID> key, cur_key,temp_key;
  line = 0;

  
  while(!grammar.eof())
    //    while(cin)
      {
		++line;
	grammar.getline(buf, MAX_LINE_LENGTH);
	if (buf[0] == 0) continue;
	bool already_in = false, exists = false;
	ParseLine(buf, &cur_key, &cur_counts);	
	for(int i=0;i<cur_key.size()-1; i++)
	  {
	    if(!already_in)
	      {
		if (DEBUG)cout << line << " " << cur_key[i] << " ::: ";
		if(cur_key[i] == -1)
		  {
		    
		    if(!temp_key.empty())
		      {

			if (DEBUG){
			  cout << endl <<  "Going with : ";
			  for(int jj=0;jj<temp_key.size();jj++)
			    cout << temp_key[jj] << " ";
			}
			
			//go through trie here		
			exists = root->CheckEdge(temp_key);
			if(exists){ if (DEBUG) cout << "EXISTS "; already_in= true;  }
			if (DEBUG)cout<< endl;
			temp_key.clear();
		      }
		    
		  }
		else{
		  temp_key.push_back(cur_key[i]);      
		  
		}
	      }
	    
	  }

	if(!already_in)
	  if(!temp_key.empty())
	    { 
	      if (DEBUG){
		cout << endl <<  "Going with : ";
		for(int jj=0;jj<temp_key.size();jj++)
		  cout << temp_key[jj] << " ";
	      }

	      exists = root->CheckEdge(temp_key);
	      if(exists) if(DEBUG) cout << "EXISTS ";
	      if (DEBUG)cout << endl;
	      temp_key.clear();
	    }

	//associative array over keys 
	//	cout << "[" << cur_counts[] <<"]";
       
	if (DEBUG) cout << "|||" << endl;
	if(exists)
	  {
	    filter_grammar_ << buf << endl;
	  }
      }
    //get back word by calling td with string     
 
  filter_grammar_.close(); 

  return 0;
}
