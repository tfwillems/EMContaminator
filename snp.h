#ifndef SNP_H_
#define SNP_H_

#include <assert.h>
#include <stdint.h>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"

class SNP {
 private:
  int32_t pos_;
  char ref_, alt_;
  int num_samples_;
  int* dosage_; // Number of non-reference SNPs for each sample
  
 public:
  SNP(int32_t pos, char ref, char alt, std::vector<char>& gt_a, std::vector<char>& gt_b){
    assert(gt_a.size() == gt_b.size());
    pos_         = pos;
    ref_         = ref;
    alt_         = alt;
    num_samples_ = gt_a.size();
    
    dosage_ = new int[num_samples_];
    for (unsigned int i = 0; i < gt_a.size(); i++){
      assert(gt_a[i] == ref_ || gt_a[i] == alt_);
      assert(gt_b[i] == ref_ || gt_b[i] == alt_);
      dosage_[i] = (gt_a[i] == ref_ ? 0 : 1) + (gt_b[i] == ref_ ? 0 : 1);
    }
  }

  SNP(int32_t pos, char ref, char alt, std::vector<int>& gt_a, std::vector<int>& gt_b){
    assert(gt_a.size() == gt_b.size());
    pos_         = pos;
    ref_         = ref;
    alt_         = alt;
    num_samples_ = gt_a.size();

    dosage_ = new int[num_samples_];
    for (unsigned int i = 0; i < gt_a.size(); i++){
      assert(gt_a[i] == 0 || gt_a[i] == 1);
      assert(gt_b[i] == 0 || gt_b[i] == 1);
      dosage_[i] = (gt_a[i] + gt_b[i]);
    }
  }

  SNP(const SNP& snp){
    pos_         = snp.pos_;
    ref_         = snp.ref_;
    alt_         = snp.alt_;
    num_samples_ = snp.num_samples_;
    dosage_      = new int[num_samples_];
    for (unsigned int i = 0; i < num_samples_; i++)
      dosage_[i] = snp.dosage_[i];
  }

  int32_t pos() const { return pos_; }
  char    ref() const { return ref_; }
  char    alt() const { return alt_; }
  
  ~SNP(){
    delete [] dosage_;
  }
  
  double get_base_log_likelihood(char base, char quality, int sample_index, double log_correct, double log_error);

  void print_base_log_likelihoods(char base, char quality, double log_correct, double log_error);
};


class SNPSorter {
 public:
  bool operator() (const SNP& snp_a, const SNP& snp_b){
    return snp_a.pos() < snp_b.pos();
  }
};


void extract_bases_and_qualities(BamTools::BamAlignment& aln, std::vector<SNP>& snps,
				 std::vector<char>& bases,  std::vector<char>& quals, std::vector<int>& dists);
#endif
