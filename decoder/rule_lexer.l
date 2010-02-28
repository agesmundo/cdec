%{
#include <string>
#include <iostream>
#include <cstring>
#include "tdict.h"
#include "fdict.h"

#define YY_SKIP_YYWRAP 1
int num_rules = 0;
int yywrap() { return 1; }

#define MAX_TOKEN_SIZE 255
std::string scfglex_tmp_token(MAX_TOKEN_SIZE, '\0');

#define MAX_RULE_SIZE 48
int scfglex_src_rhs[MAX_RULE_SIZE];
int scfglex_trg_rhs[MAX_RULE_SIZE];
int scfglex_src_rhs_size;
int scfglex_trg_rhs_size;
int scfglex_lhs;
int scfglex_src_arity;
int scfglex_trg_arity;

#define MAX_FEATS 20
int scfglex_feat_ids[MAX_FEATS];
double scfglex_feat_vals[MAX_FEATS];
int scfglex_num_feats;

#define MAX_ARITY 20
int scfglex_nt_sanity[MAX_ARITY];
int scfglex_src_nts[MAX_ARITY];

void sanity_check_trg_symbol(int index) {
  // TODO once we populate scfglex_src_nts;
}

void sanity_check_trg_index(int index) {
  if (index > scfglex_src_arity) {
    std::cerr << "Target index " << index << " exceeds source arity " << scfglex_src_arity << std::endl;
    abort();
  }
  int& flag = scfglex_nt_sanity[index - 1];
  if (flag) {
    std::cerr << "Target index " << index << " used multiple times!" << std::endl;
    abort();
  }
  flag = 1;
}

void scfglex_reset() {
  scfglex_src_arity = 0;
  scfglex_trg_arity = 0;
  scfglex_num_feats = 0;
  scfglex_src_rhs_size = 0;
  scfglex_trg_rhs_size = 0;
}

%}

REAL -?[0-9]+(\.[0-9]*([eE][-+]*[0-9]+)?)?|inf|-inf
NT [A-Z_][A-Z+/=_0-9!:@]*

%x LHS_END SRC TRG FEATS FEATVAL
%%

<INITIAL>[ \t]	;

<INITIAL>\[{NT}\]   {
		scfglex_tmp_token.assign(yytext + 1, yyleng - 2);
		TD::Convert(scfglex_tmp_token);
		// std::cerr << scfglex_tmp_token << "\n";
  		BEGIN(LHS_END);
		}

<SRC>\[{NT}\]   {
		scfglex_tmp_token.assign(yytext + 1, yyleng - 2);
		TD::Convert(scfglex_tmp_token);
		++scfglex_src_arity;
		int index = scfglex_src_arity;
		}

<SRC>\[{NT},[1-9][0-9]?\]   {
		int index = yytext[yyleng - 2] - '0';
		if (yytext[yyleng - 3] == ',') {
		  scfglex_tmp_token.assign(yytext + 1, yyleng - 4);
		} else {
		  scfglex_tmp_token.assign(yytext + 1, yyleng - 5);
		  index += 10 * (yytext[yyleng - 3] - '0');
		}
		++scfglex_src_arity;
		if (scfglex_src_arity != index) {
			std::cerr << "Src indices must go in order: expected " << scfglex_src_arity << " but got " << index << std::endl;
			abort();
		}
		TD::Convert(scfglex_tmp_token);
}

<TRG>\[{NT},[1-9][0-9]?\]   {
		int index = yytext[yyleng - 2] - '0';
		if (yytext[yyleng - 3] == ',') {
		  scfglex_tmp_token.assign(yytext + 1, yyleng - 4);
		} else {
		  scfglex_tmp_token.assign(yytext + 1, yyleng - 5);
		  index += 10 * (yytext[yyleng - 3] - '0');
		}
		++scfglex_trg_arity;
		// std::cerr << "TRG INDEX: " << index << std::endl;
		sanity_check_trg_symbol(index);
		sanity_check_trg_index(index);
		TD::Convert(scfglex_tmp_token);
}

<TRG>\[[1-9][0-9]?\]   {
		int index = yytext[yyleng - 2] - '0';
		if (yyleng == 4) {
		  index += 10 * (yytext[yyleng - 3] - '0');
		}
		++scfglex_trg_arity;
		sanity_check_trg_index(index);
}

<LHS_END>[ \t] { ; }
<LHS_END>\|\|\|	{
		scfglex_reset();
		BEGIN(SRC);
		}
<INITIAL,LHS_END>.	{
		std::cerr << "Unexpected input in LHS: " << yytext << std::endl;
		abort();
		}

<SRC>\|\|\|	{
		memset(scfglex_nt_sanity, 0, scfglex_src_arity * sizeof(int));
		BEGIN(TRG);
		}
<SRC>[^ \t]+	{ 
		scfglex_tmp_token.assign(yytext, yyleng);
		TD::Convert(scfglex_tmp_token);
		}
<SRC>[ \t]+	{ ; }

<TRG>\|\|\|	{
		BEGIN(FEATS);
		}
<TRG>[^ \t]+	{
		scfglex_tmp_token.assign(yytext, yyleng);
		TD::Convert(scfglex_tmp_token);
		}
<TRG>[ \t]+	{ ; }

<TRG,FEATS>\n	{
		num_rules++;
		BEGIN(INITIAL);
		}

<FEATS>[ \t]	{ ; }
<FEATS>[^ \t=;]+=	{
		scfglex_tmp_token.assign(yytext, yyleng - 1);
		FD::Convert(scfglex_tmp_token);
		BEGIN(FEATVAL);
		}
<FEATVAL>{REAL}	{
		double fval = strtod(yytext, NULL);
		BEGIN(FEATS);
		}
<FEATVAL>.	{
		std::cerr << "Unexpected input in feature value: " << yytext << std::endl;
		abort();
		}
<FEATS>{REAL} 	{ 
		double fval = strtod(yytext, NULL);
		}
<FEATS>.	{
		std::cerr << "Unexpected input in features: " << yytext << std::endl;
		abort();
		}

%%

// int main() {
//        yylex();
//	printf("%d rules\n", num_rules);
//}


