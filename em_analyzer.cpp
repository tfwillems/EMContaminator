#include <assert.h>

#include <cfloat>

#include "em_analyzer.h"
#include "mathops.h"

double AlignmentEntry::recalc_log_sample_posteriors(double* log_sample_priors, double log_correct, double log_error){
  std::memcpy(log_sample_posteriors_, log_sample_priors, num_samples_*sizeof(double));
  if (snps_.empty())
    return 0;
  for (unsigned int i = 0; i < snps_.size(); i++){
    //snps_[i].print_base_log_likelihoods(bases_[i], qualities_[i], log_correct, log_error);
    for (unsigned int j = 0; j < num_samples_; j++)
      log_sample_posteriors_[j] += snps_[i].get_base_log_likelihood(bases_[i], qualities_[i], j, log_correct, log_error);
  }

  double total_LL = log_sum_exp(log_sample_posteriors_, log_sample_posteriors_+num_samples_);
  for (int i = 0; i < num_samples_; i++)
    log_sample_posteriors_[i] -= total_LL;
  assert(total_LL < TOLERANCE);
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
    for (unsigned int i = 0; i < aln_entries_.size(); i++)
      new_LL += aln_entries_[i].recalc_log_sample_posteriors(log_sample_priors_, log_correct, log_error);

    std::cerr << "ITER #" << num_iter << " , LL =" << new_LL << std::endl;
    print_sample_priors(std::cerr);
    std::cerr << std::endl << std::endl;

    // Occasionally the log-likelihood isn't monotonically increasing b/c of the pseudocounts
    // used when recalculating the sample fractions. Let's return true anyways
    assert(new_LL <= TOLERANCE);
    if (new_LL < LL+TOLERANCE){
      LL = new_LL;
      return true;
    }

    // M-Step: Reestimate the sample fractions
    // For each sample, use a read ownership pseudocount of exp(LOG_READ_PRIOR_PC)
    // Use streaming log-sum-exp trick to avoid underflow
    std::vector<double> log_priors_max(num_samples_, LOG_READ_PRIOR_PC);
    std::vector<double> log_priors_total(num_samples_, 1.0);
    for (unsigned int i = 0; i < aln_entries_.size(); i++){
      double* log_read_sample_posteriors = aln_entries_[i].get_log_sample_posteriors();
      for (unsigned int j = 0; j < num_samples_; j++)
	update_streaming_log_sum_exp(log_read_sample_posteriors[j], log_priors_max[j], log_priors_total[j]);
    }
    for (unsigned int i = 0; i < num_samples_; i++)
      log_sample_priors_[i] = finish_streaming_log_sum_exp(log_priors_max[i], log_priors_total[i]);

    // Normalize the read ownership counts to obtain posteriors
    double log_total_count = log_sum_exp(log_sample_priors_, log_sample_priors_+num_samples_);
    for (unsigned int i = 0; i < num_samples_; i++)
      log_sample_priors_[i] -= log_total_count;

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

bool use_chrom(std::string& chrom){
  if (chrom.find("_") != std::string::npos)
    return false;
  if (chrom.find("X") != std::string::npos)
    return false;
  if (chrom.find("Y") != std::string::npos)
    return false;
  if (chrom.find("M") != std::string::npos)
    return false;
  return true;
}

bool EMAnalyzer::analyze(BamTools::BamMultiReader& reader, double error_rate){
  int32_t snp_match_count = 0, snp_mismatch_count = 0;
  BamTools::RefVector ref_vector = reader.GetReferenceData();

  // Load and store all of the relevant alignment information
  int32_t aln_count = 0, aln_with_snps_count = 0;
  BamTools::BamAlignment alignment;
  while (reader.GetNextAlignment(alignment)){
    aln_count++;
    if (aln_count % 1000 == 0)
      std::cerr << aln_count << " " << aln_with_snps_count << std::endl;
    if (!alignment.IsMapped() || alignment.Position == 0 || alignment.CigarData.size() == 0 || alignment.Length == 0)
      continue;

    std::string chrom  = ref_vector[alignment.RefID].RefName;
    int32_t start      = alignment.Position;
    int32_t end        = alignment.GetEndPosition()-1;
    if (!use_chrom(chrom))
      continue;
    if(!snp_vcf_.set_region(chrom, start+1, end+1))
      printErrorAndDie("Failed to set the region to chromosome " + chrom + " in the SNP VCF. Please check the SNP VCF and rerun the analysis");

    std::vector<SNP> snps;
    std::vector<int32_t> snp_positions;
    VCF::Variant snp_variant;
    while (snp_vcf_.get_next_variant(snp_variant)){
      snp_positions.push_back(snp_variant.get_position());
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

      char ref = snp_variant.get_allele(0)[0];
      char alt = snp_variant.get_allele(1)[0];
      snps.push_back(SNP(snp_variant.get_position()-1, ref, alt, gts_a, gts_b));
    }

    // Extract the bases and quality scores stored at each SNP
    // SNPs without a matching base in the alignment will have a '-' character inserted
    std::vector<char> bases, quals;
    std::vector<int> dists;
    extract_bases_and_qualities(alignment, snps, bases, quals, dists);

    // Construct a new alignment entry
    aln_entries_.push_back(AlignmentEntry(num_samples_));
    for (unsigned int i = 0; i < snps.size(); i++){
      if (bases[i] != '-' && dists[i] >= DIST_BASE_FROM_INDEL){
	aln_entries_.back().add_snp(snps[i], bases[i], quals[i]);
	if (bases[i] == snps[i].ref() || bases[i] == snps[i].alt())
	  snp_match_count++;
	else
	  snp_mismatch_count++;
      }
    }

    if (aln_entries_.back().num_snps() > 0)
      aln_with_snps_count++;
    else
      aln_entries_.pop_back();
    // TO DO: Discard entries with 0 informative SNPs?
  }

  std::cout << "SNP match    count = " << snp_match_count    << "\n"
	    << "SNP mismatch count = " << snp_mismatch_count << "\n" << "\n";

  // Run the EM procedure and return whether or not it succeeded
  double LL;
  return run_EM(error_rate, LL);
}
