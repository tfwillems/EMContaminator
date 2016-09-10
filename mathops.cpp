#include <algorithm>
#include <assert.h>
#include <math.h>

#include "mathops.h"

const double LOG_ONE_HALF  = log(0.5);
const double TOLERANCE     = 1e-10;
const double LOG_E_BASE_10 = 0.4342944819;

double sum(double* begin, double* end){
  double total = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += *iter;
  return total;
}

double sum(std::vector<double>& vals){
  double total = 0.0;
  for (auto iter = vals.begin(); iter != vals.end(); iter++)
    total += *iter;
  return total;
}

double log_sum_exp(double* begin, double* end){
  double max_val = *std::max_element(begin, end);
  double total   = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}

double log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2)
    return log_v1 + log(1 + exp(log_v2-log_v1));
  else
    return log_v2 + log(1 + exp(log_v1-log_v2));
}

double log_sum_exp(std::vector<double>& log_vals){
  double max_val = *std::max_element(log_vals.begin(), log_vals.end());
  double total   = 0;
  for (auto iter = log_vals.begin(); iter != log_vals.end(); iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}
