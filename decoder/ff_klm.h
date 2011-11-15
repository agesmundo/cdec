#ifndef _KLM_FF_H_
#define _KLM_FF_H_

#include <vector>
#include <string>

#include "ff.h"
#include "lm/model.hh"

template <class Model> struct KLanguageModelImpl;

// the supported template types are instantiated explicitly
// in ff_klm.cc.
template <class Model>
class KLanguageModel : public FeatureFunction {
 public:
  // param = "filename.lm [-o n]"
  KLanguageModel(const std::string& param);
  ~KLanguageModel();
  virtual void FinalTraversalFeatures(const void* context,
                                      SparseVector<double>* features) const;
  static std::string usage(bool param,bool verbose);
  Features features() const;
 protected:
  virtual void TraversalFeaturesImpl(const SentenceMetadata& smeta,
                                     const Hypergraph::Edge& edge,
                                     const std::vector<const void*>& ant_contexts,
                                     SparseVector<double>* features,
                                     SparseVector<double>* estimated_features,
                                     void* out_context) const;

  virtual void TraversalUndirectedFeaturesImpl(const SentenceMetadata& smeta,
                                          UCandidate& ucand,
                                          int spos
                                          /*const std::vector<const void*>& ant_states,
                                          SparseVector<double>* features,
                                          SparseVector<double>* estimated_features,
                                          void* state*/) const;

 private:
  int fid_; // conceptually const; mutable only to simplify constructor
  int oov_fid_; // will be zero if extra OOV feature is not configured by decoder
  //  int est_fid_;
  //  int lnk_fid_; //number of edge's links (1+arity) //binary version win, see lnk_bin_fids_
  int cln_fid_; //connected links: number of connected edge's links

  //feats # bounded by edge's links #
  static const int max_lnks_ =3; //max number of edge's links (head+2 children) for binary grammar
  int* lnk_bin_fids_; //binary version of lnk_fid_
//  int* cln_bin_fids_; //binary version of cln_fid_ //lin version win, only 2 values effective

  //feats # bounded by order
  int* ngram_avg_fids_;//for each 0<=i<order fid of the feat storing the average score of all ngrams of size i+1
  int* ngram_cnt_fids_;//for each 0<=i<order fid of the feat storing the count all ngrams of size i+1
//  int* ngram_oov_fids_;//

//  //cnt bin-lin
//  int** ngram_cnt_bin_fids_;
//  int ngram_cnt_bin_threshold_;

  //inside-outside
  int inout_fid_;

  KLanguageModelImpl<Model>* pimpl_;
};

#endif
