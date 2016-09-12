#ifndef EM_ANALYZER_H_
#define EM_ANALYZER_H_

#include <assert.h>
#include "math.h"

#include <iostream>
#include <string>
#include <vector>

#include "bamtools/include/api/BamMultiReader.h"
#include "mathops.h"
#include "snp.h"
#include "vcf_reader.h"

class AlignmentEntry{
 private:
  int num_samples_;
  double* log_sample_posteriors_;
  std::vector<SNP> snps_;
  std::vector<char> bases_;
  std::vector<char> qualities_;

 public:
  AlignmentEntry(int num_samples){
    log_sample_posteriors_ = new double[num_samples_];
  }
  
  ~AlignmentEntry(){
    delete [] log_sample_posteriors_;
  }

  void add_snp(SNP& snp, char base, char qual){
    snps_.push_back(snp);
    bases_.push_back(base);
    qualities_.push_back(qual);
  }

  double recalc_log_sample_posteriors(double* log_sample_priors, double log_correct, double log_error);

  double* get_log_sample_posteriors(){ return log_sample_posteriors_; }
};





class EMAnalyzer {
 private:
  VCF::VCFReader snp_vcf_;
  int num_samples_;
  double* log_sample_priors_;
  std::vector<AlignmentEntry> aln_entries_;
  bool init_randomly_;

  // Parameters for EM converge
  int MAX_EM_ITER;
  double ABS_LL_CONVERGE;   // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
  double FRAC_LL_CONVERGE;  // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE
  double LOG_READ_PRIOR_PC; // Pseudocount used when recomputing sample priors

  int DIST_BASE_FROM_INDEL; // Only consider a SNP if the base in the read is at least this many
                            // base pairs from the nearest indel in its alignment
  
  void init_log_sample_priors(){
    if (init_randomly_){
      printErrorAndDie("Random initialization not yet implemented");
    }
    else
      std::fill(log_sample_priors_, log_sample_priors_+num_samples_, -1);
    
    double total_LL = log_sum_exp(log_sample_priors_, log_sample_priors_+num_samples_);
    for (unsigned int i = 0; i < num_samples_; i++)
      log_sample_priors_[i] -= total_LL;
  }

  bool run_EM(double error_rate, double& LL);

 public:
 EMAnalyzer(std::string& snp_vcf_file, bool init_randomly): snp_vcf_(snp_vcf_file){
    num_samples_         = snp_vcf_.get_samples().size();
    log_sample_priors_   = new double[num_samples_];
    MAX_EM_ITER          = 100;
    ABS_LL_CONVERGE      = 0.01;
    FRAC_LL_CONVERGE     = 0.001;
    LOG_READ_PRIOR_PC    = log(0.1);
    DIST_BASE_FROM_INDEL = 7;
    init_log_sample_priors();
  }

  ~EMAnalyzer(){
    delete [] log_sample_priors_;
  }

  bool analyze(BamTools::BamMultiReader& reader, double error_rate);

  void print_sample_priors(std::ostream& out){
    const std::vector<std::string>& samples = snp_vcf_.get_samples();
    for (unsigned int i = 0; i < num_samples_; i++)
      out << samples[i] << exp(log_sample_priors_[i]) << "\n";
    out << "\n";
  }
};




#endif
