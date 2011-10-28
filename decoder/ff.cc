//TODO: non-sparse vector for all feature functions?  modelset applymodels keeps track of who has what features?  it's nice having FF that could generate a handful out of 10000 possible feats, though.

//TODO: actually score rule_feature()==true features once only, hash keyed on rule or modify TRule directly?  need to keep clear in forest which features come from models vs. rules; then rescoring could drop all the old models features at once

#include "fast_lexical_cast.hpp"
#include <stdexcept>
#include "ff.h"

#include "tdict.h"
#include "hg.h"
#include "weights.h"

using namespace std;

FeatureFunction::~FeatureFunction() {}

void FeatureFunction::PrepareForInput(const SentenceMetadata&) {}

void FeatureFunction::FinalTraversalFeatures(const void* /* ant_state */,
                                             SparseVector<double>* /* features */) const {
}

string FeatureFunction::usage_helper(std::string const& name,std::string const& params,std::string const& details,bool sp,bool sd) {
  string r=name;
  if (sp) {
    r+=": ";
    r+=params;
  }
  if (sd) {
    r+="\n";
    r+=details;
  }
  return r;
}

Features FeatureFunction::single_feature(WordID feat) {
  return Features(1,feat);
}

Features ModelSet::all_features(std::ostream *warn,bool warn0) {
  return ::all_features(models_,weights_,warn,warn0);
}

void show_features(Features const& ffs,DenseWeightVector const& weights_,std::ostream &out,std::ostream &warn,bool warn_zero_wt) {
  out << "Weight  Feature\n";
  for (unsigned i=0;i<ffs.size();++i) {
    WordID fid=ffs[i];
    string const& fname=FD::Convert(fid);
    double wt=weights_[fid];
    if (warn_zero_wt && wt==0)
      warn<<"WARNING: "<<fname<<" has 0 weight."<<endl;
    out << wt << "  " << fname<<endl;
  }
}

void ModelSet::show_features(std::ostream &out,std::ostream &warn,bool warn_zero_wt)
{
//  ::show_features(all_features(),weights_,out,warn,warn_zero_wt);
  show_all_features(models_,weights_,out,warn,warn_zero_wt,warn_zero_wt);
}

// Hiero and Joshua use log_10(e) as the value, so I do to
WordPenalty::WordPenalty(const string& param) :
  fid_(FD::Convert("WordPenalty")),
    value_(-1.0 / log(10)) {

//	string featname = "WordPenalty";
//	bin_threshold_=7;
//	bin_fids_=new int[bin_threshold_];//TODO remember destroyer
//	  for(int i=0; i<bin_threshold_; i++){
//		  string currFeat;
//		  stringstream ss;
//		  string id;
//		  ss << i;
//		  ss >> id;
//
//		  string suff= "_WRP-BIN_" ;
//		  currFeat = featname+suff+id;
//		  bin_fids_[i] = FD::Convert(currFeat);
//		  cerr << currFeat << " ; FID: " << bin_fids_[i] << endl;
//	  }
  if (!param.empty()) {
    cerr << "Warning WordPenalty ignoring parameter: " << param << endl;
  }
}

void FeatureFunction::TraversalFeaturesImpl(const SentenceMetadata& smeta,
                                        const Hypergraph::Edge& edge,
                                        const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state) const {
  throw std::runtime_error("TraversalFeaturesImpl not implemented - override it or TraversalFeaturesLog.\n");
  abort();
}
//GU
void FeatureFunction::TraversalUndirectedFeaturesImpl(const SentenceMetadata& smeta,
                                        UCandidate& ucand,
                                        int spos
                                        /*const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state*/) const {
  throw std::runtime_error("TraversalUndirectedFeaturesImpl not implemented - override it.\n");
  abort();
}

void WordPenalty::TraversalFeaturesImpl(const SentenceMetadata& smeta,
                                        const Hypergraph::Edge& edge,
                                        const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state) const {
  (void) smeta;
  (void) ant_states;
  (void) state;
  (void) estimated_features;
  features->set_value(fid_, edge.rule_->EWords() * value_);
}
//GU
void WordPenalty::TraversalUndirectedFeaturesImpl(const SentenceMetadata& smeta,
                                        UCandidate& ucand,
                                        int spos/*
                                        const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state*/) const {
  (void) smeta;
//  (void) ant_states;
//  (void) state;
//  (void) estimated_features;

  //lin
  ucand.feature_values_.set_value(fid_, ucand.in_edge_->rule_->EWords() * value_);

//  //bin
//  int target_words=ucand.in_edge_->rule_->EWords();
//  if(target_words<bin_threshold_-1){
//	  ucand.feature_values_.set_value(bin_fids_[target_words], target_words* value_);
//  }else{
//	  ucand.feature_values_.set_value(bin_fids_[bin_threshold_-1], target_words* value_); //last is lin end
//  }

}

SourceWordPenalty::SourceWordPenalty(const string& param) :
    fid_(FD::Convert("SourceWordPenalty")),
    value_(-1.0 / log(10)) {
  if (!param.empty()) {
    cerr << "Warning SourceWordPenalty ignoring parameter: " << param << endl;
  }
}

Features SourceWordPenalty::features() const {
  return single_feature(fid_);
}

Features WordPenalty::features() const {
  return single_feature(fid_);
}


void SourceWordPenalty::TraversalFeaturesImpl(const SentenceMetadata& smeta,
                                        const Hypergraph::Edge& edge,
                                        const std::vector<const void*>& ant_states,
                                        SparseVector<double>* features,
                                        SparseVector<double>* estimated_features,
                                        void* state) const {
  (void) smeta;
  (void) ant_states;
  (void) state;
  (void) estimated_features;
  features->set_value(fid_, edge.rule_->FWords() * value_);
}

ArityPenalty::ArityPenalty(const std::string& param) :
    value_(-1.0 / log(10)) {
  string fname = "Arity_";
  unsigned MAX=DEFAULT_MAX_ARITY;
  using namespace boost;
  if (!param.empty())
    MAX=lexical_cast<unsigned>(param);
  for (unsigned i = 0; i <= MAX; ++i) {
    WordID fid=FD::Convert(fname+lexical_cast<string>(i));
    fids_.push_back(fid);
  }
  while (!fids_.empty() && fids_.back()==0) fids_.pop_back(); // pretty up features vector in case FD was frozen.  doesn't change anything
}

Features ArityPenalty::features() const {
  return Features(fids_.begin(),fids_.end());
}

void ArityPenalty::TraversalFeaturesImpl(const SentenceMetadata& smeta,
                                         const Hypergraph::Edge& edge,
                                         const std::vector<const void*>& ant_states,
                                         SparseVector<double>* features,
                                         SparseVector<double>* estimated_features,
                                         void* state) const {
  (void) smeta;
  (void) ant_states;
  (void) state;
  (void) estimated_features;
  unsigned a=edge.Arity();
  features->set_value(a<fids_.size()?fids_[a]:0, value_);
}

ModelSet::ModelSet(const vector<double>& w, const vector<const FeatureFunction*>& models) :
    models_(models),
    weights_(w),
    weights_avg_(w),
    is_avg_(true),
    count_avg_(0),
    state_size_(0),
    model_state_pos_(models.size()) {
  for (int i = 0; i < models_.size(); ++i) {
    model_state_pos_[i] = state_size_;
    state_size_ += models_[i]->NumBytesContext();
  }
  //TODO fix this in initialization find method
  for(int i=0; i<weights_avg_.size();i++){
	  weights_avg_[i]=0;
  }
}

void ModelSet::PrepareForInput(const SentenceMetadata& smeta) {
  for (int i = 0; i < models_.size(); ++i)
    const_cast<FeatureFunction*>(models_[i])->PrepareForInput(smeta);
}

void ModelSet::AddFeaturesToEdge(const SentenceMetadata& smeta,
                                 const Hypergraph& /* hg */,
                                 const FFStates& node_states,
                                 Hypergraph::Edge* edge,
                                 FFState* context,
                                 prob_t* combination_cost_estimate) const {
  edge->reset_info();
  context->resize(state_size_);
  if (state_size_ > 0) {
    memset(&(*context)[0], 0, state_size_);
  }
  SparseVector<double> est_vals;  // only computed if combination_cost_estimate is non-NULL
  if (combination_cost_estimate) *combination_cost_estimate = prob_t::One();
  for (int i = 0; i < models_.size(); ++i) {
    const FeatureFunction& ff = *models_[i];
    void* cur_ff_context = NULL;
    vector<const void*> ants(edge->tail_nodes_.size());
    bool has_context = ff.NumBytesContext() > 0;
    if (has_context) {
      int spos = model_state_pos_[i];
      cur_ff_context = &(*context)[spos];
      for (int i = 0; i < ants.size(); ++i) {
        ants[i] = &node_states[edge->tail_nodes_[i]][spos];
      }
    }
    ff.TraversalFeatures(smeta, *edge, ants, &edge->feature_values_, &est_vals, cur_ff_context);
  }
  if (combination_cost_estimate)
    combination_cost_estimate->logeq(est_vals.dot(weights_));
  edge->edge_prob_.logeq(edge->feature_values_.dot(weights_));
}
//GU
void ModelSet::AddFeaturesToUCandidate(const SentenceMetadata& smeta,
                                 //const FFStates& node_states,
                                 UCandidate* ucand//,
                                 //Hypergraph::Edge* edge,
                                 //FFState* context,
                                 /*prob_t* combination_cost_estimate*/) const {
	ucand->InitStates(state_size_);

  //SparseVector<double> est_vals;  // only computed if combination_cost_estimate is non-NULL
  //if (combination_cost_estimate) *combination_cost_estimate = prob_t::One();
  for (int i = 0; i < models_.size(); ++i) {
    const FeatureFunction& ff = *models_[i];
//    void* cur_ff_context = NULL;
//    /*rm*/assert(ucand->in_edge_->tail_nodes_.size()==ucand->in_edge_->rule_->Arity());//debug TODO RM
//    vector<const void*> ants(ucand->in_edge_->tail_nodes_.size());
    bool has_context = ff.NumBytesContext() > 0;
    int spos=0;
        if (has_context) {
          int spos = model_state_pos_[i];
//          cur_ff_context = &(*context)[spos];
        }
//    if (has_context) {
//      int spos = model_state_pos_[i];
//      cur_ff_context = &(*context)[spos];
//      for (int i = 0; i < ants.size(); ++i) {//TODO NB ants useles each UCand ha its own state, not merging per out node as CP
//        ants[i] = &node_states[edge->tail_nodes_[i]][spos];//TODO then remove ants when know how to adapt LM (context space issue)
//      }
//    }
  //TODO pass spos (and smeta) all the time seems like a waste of time
    ff.TraversalUndirectedFeatures(smeta, *ucand, spos/*ants, &ucand->feature_values_, &ucand->est_vals_, cur_ff_context*/);
  }
  prob_t estimate = prob_t::One();
  estimate.logeq(ucand->est_vals_.dot(weights_));

  prob_t local = prob_t::One();
  local.logeq(ucand->feature_values_.dot(weights_));

  ucand->action_prob_ = local * estimate; //sum exps
}

void ModelSet::AddFinalFeatures(const FFState& state, Hypergraph::Edge* edge,SentenceMetadata const& smeta) const {
  assert(1 == edge->rule_->Arity());
  edge->reset_info();
  for (int i = 0; i < models_.size(); ++i) {
    const FeatureFunction& ff = *models_[i];
    const void* ant_state = NULL;
    bool has_context = ff.NumBytesContext() > 0;
    if (has_context) {
      int spos = model_state_pos_[i];
      ant_state = &state[spos];
    }
    ff.FinalTraversalFeatures(smeta, *edge, ant_state, &edge->feature_values_);
  }
  edge->edge_prob_.logeq(edge->feature_values_.dot(weights_));
}

ostream& ModelSet::PrintWeights(ostream& os) {
	os<< "WEIGHTS: ";
	for(int i=0; i< weights_.size();i++){
		os << FD::Convert(i)<<"("<<i<<")"<<"="<<weights_[i]<< " ";
	}
	os<<endl;

	if(is_avg_){
		os<<"AVG_WEIGHTS: ";
		std::vector<double> avg_w=ComputeAvgWeight();
		for(int i=0; i< avg_w.size();i++){
			os << FD::Convert(i)<<"("<<i<<")"<<"="<<avg_w[i]<< " ";
		}
		os<<endl;
		os<<"COUNT_AVG: "<< count_avg_<<endl;
	}

}

void ModelSet::UpdateWeight(SparseVector<Featval> vector, double loss){
//PASSIVE AGGRESSIVE
//	double norm = vector.l2norm_sq();
//	assert(norm!=0);
//	double alpha =  loss / norm;

//	cerr << "\tALPHA= "<<alpha<<endl;

	for (SparseVector<Featval>::const_iterator i=vector.begin(),e=vector.end();i!=e;++i) {
		if (weights_.size() <= i->first) weights_.resize(i->first+1);
		weights_[i->first] += i->second /* * alpha*/;
		if(is_avg_){
			if (weights_avg_.size() <= i->first) weights_avg_.resize(i->first+1);
			weights_avg_[i->first] += i->second * (count_avg_+1.0)/* * alpha*/;
		}
	}
	count_avg_++;
}

double ModelSet::ScoreVector(SparseVector<Featval> vector){
	return vector.dot(weights_);
}

std::vector<double> ModelSet::ComputeAvgWeight(){
	std::vector<double> avg_w (weights_);
	assert(weights_.size()==weights_avg_.size());
	for(int i =0 ; i<weights_.size(); i++){
		avg_w[i] -= (weights_avg_[i] / (count_avg_+1));
	}
	return avg_w;
}

void ModelSet::WriteToFile(const std::string& fname){
	Weights w;

	if(is_avg_){
		w.InitFromVector(ComputeAvgWeight());
	}
	else{
		w.InitFromVector(weights_);
	}
	w.WriteToFile(fname);
}
