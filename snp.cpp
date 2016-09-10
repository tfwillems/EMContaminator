#include "snp.h"

#include "error.h"
#include "mathops.h"

double SNP::get_base_log_likelihood(char base, char quality, int sample_index, double log_correct, double log_error){
  switch(dosage_[sample_index]){
  case 0:
    return (base == ref_ ? log_correct : log_error);
    break;
  case 1:
    return ((base == ref_ || base == alt_) ? LOG_ONE_HALF : log_error);
    break;
  case 2:
    return (base == ref_ ? log_error : log_correct);
    break;
  default:
    printErrorAndDie("Invalid SNP dosage");
    break;
  }
}
