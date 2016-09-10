#include <assert.h>

#include <cfloat>

#include "em_analyzer.h"
#include "mathops.h"

double AlignmentEntry::recalc_log_sample_posteriors(double* log_sample_priors, double log_correct, double log_error){
  std::memcpy(log_sample_posteriors_, log_sample_priors, num_samples_*sizeof(double));
  for(unsigned int i = 0; i < snps_.size(); i++)
    for (unsigned int j = 0; j < num_samples_; j++)
      log_sample_posteriors_[j] += snps_[i].get_base_log_likelihood(bases_[i], qualities_[i], j, log_correct, log_error);

  double total_LL = log_sum_exp(log_sample_posteriors_, log_sample_posteriors_+num_samples_);
  for (int i = 0; i < num_samples_; i++)
    log_sample_posteriors_[i] -= total_LL;
  return total_LL;
}

void update_streaming_log_sum_exp(double log_val, double& max_val, double& total){
  if (log_val <= max_val)
    total += fasterexp(log_val - max_val);
  else {
    total  *= fasterexp(max_val-log_val);
    total  += 1.0;
    max_val = log_val;
  }
}

double finish_streaming_log_sum_exp(double max_val, double total){
  return max_val + fasterlog(total);
}

bool EMAnalyzer::run_EM(double error_rate, double& LL){
  // For now, these are fixed. But we may want to estimate them as well
  double log_correct = log(1.0-error_rate);
  double log_error   = log(error_rate);

  int num_iter = 1;
  LL = -DBL_MAX;

  // Initialize sample fractions
  init_log_sample_priors();

  while (num_iter <= MAX_EM_ITER){
    // E-Step: Compute each alignment's sample posteriors and the total LL of the new configuration
    double new_LL = 0;
    for (unsigned int i = 0; i < alignments_.size(); i++)
      new_LL += alignments_[i].recalc_log_sample_posteriors(log_sample_priors_, log_correct, log_error);

    // Occasionally the log-likelihood isn't monotonically increasing b/c of the pseudocounts
    // used when recalculating the sample fractions. Let's return true anyways
    assert(new_LL <= TOLERANCE);
    if (new_LL < LL+TOLERANCE){
      LL = new_LL;
      return true;
    }

    // M-Step: Reestimate the sample fractions
    // For each sample, use a read ownsership pseudocount of exp(LOG_READ_PRIOR_PC)
    // Use streaming log-sum-exp trick to avoid underflow
    std::vector<double> log_priors_max(num_samples_, LOG_READ_PRIOR_PC);
    std::vector<double> log_priors_total(num_samples_, 1.0);
    for (unsigned int i = 0; i < alignments_.size(); i++){
      double* log_read_sample_posteriors = alignments_[i].get_log_sample_posteriors();
      for (unsigned int j = 0; j < num_samples_; j++)
	update_streaming_log_sum_exp(log_read_sample_posteriors[j], log_priors_max[j], log_priors_total[j]);
    }
    for (unsigned int i = 0; i < num_samples_; i++)
      log_sample_priors_[i] = finish_streaming_log_sum_exp(log_priors_max[i], log_priors_total[i]);

    // Check for convergence
    double abs_change  = new_LL - LL;
    double frac_change = -(new_LL - LL)/LL;
    if (abs_change < ABS_LL_CONVERGE && frac_change < FRAC_LL_CONVERGE){
      LL = new_LL;
      return true;
    }
    
    LL = new_LL;
    num_iter++;
  }
  return false;
}

bool EMAnalyzer::analyze(BamTools::BamMultiReader& reader, double error_rate){
  // Load and store all of the relevant alignment information
  BamTools::BamAlignment alignment;
  std::vector<SNP> snps;
  while (reader.GetNextAlignmentCore(alignment)){
    if (!alignment.IsMapped() || alignment.Position == 0 || alignment.CigarData.size() == 0 || alignment.Length == 0)
      continue;

    printErrorAndDie("Extract coords from alignment");
    std::string chrom;
    int32_t start      = alignment.Position;
    int32_t end        = alignment.GetEndPosition();
    
    if(!snp_vcf_.set_region(chrom, start, end))
      printErrorAndDie("Failed to set the region to chromosome " + chrom + " in the SNP VCF. Please check the SNP VCF and rerun the analysis");

    VCF::Variant snp_variant;
    while (snp_vcf_.get_next_variant(snp_variant)){
      if (!snp_variant.is_biallelic_snp())
	continue;

      // Check if missing any genotypes or only one allele present in genotypes
      bool missing = false;
      std::vector<int> gts_a(num_samples_, -1), gts_b(num_samples_, -1);
      std::vector<int> allele_counts(2, 0);
      for (unsigned int sample_index = 0; sample_index < num_samples_; sample_index++){
	if (snp_variant.sample_call_missing(sample_index)){
	  missing = true;
	  break;
	}
	snp_variant.get_genotype(sample_index, gts_a[sample_index], gts_b[sample_index]);
	allele_counts[gts_a[sample_index]] += 1;
	allele_counts[gts_b[sample_index]] += 1;
      }
      if (missing)
	continue;
      if (allele_counts[0] == 0 || allele_counts[1] == 0)
	continue;


      //SNP(char ref, char alt, std::vector<int>& gt_a, std::vector<int>& gt_b){


      printErrorAndDie("Implement alignment SNP extraction...");
    }
  }

  // Run the EM procedure
  double LL;
  return run_EM(error_rate, LL);
}
