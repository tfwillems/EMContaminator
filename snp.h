#ifndef SNP_H_
#define SNP_H_

#include <assert.h>
#include <vector>

class SNP {
 private:
  char ref_, alt_;
  int num_samples_;
  int* dosage_; // Number of non-reference SNPs for each sample
  
 public:
  SNP(char ref, char alt, std::vector<char>& gt_a, std::vector<char>& gt_b){
    assert(gt_a.size() == gt_b.size());
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

  SNP(char ref, char alt, std::vector<int>& gt_a, std::vector<int>& gt_b){
    assert(gt_a.size() == gt_b.size());
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

  char ref() const { return ref_; }
  char alt() const { return alt_; }
  
  ~SNP(){
    delete [] dosage_;
  }
  
  double get_base_log_likelihood(char base, char quality, int sample_index, double log_correct, double log_error);
};

#endif
