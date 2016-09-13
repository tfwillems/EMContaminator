#include <assert.h>
#include <algorithm>
#include <iostream>

#include "snp.h"
#include "error.h"
#include "mathops.h"

void SNP::print_base_log_likelihoods(char base, char quality, double log_correct, double log_error){
  int gt;
  if (base == ref_)
    gt = 0;
  else if (base == alt_)
    gt = 1;
  else
    gt = -1;

  std::cerr << "BASE LLs:" << "\n";
  for (unsigned int i = 0; i < num_samples_; i++)
    std::cerr << " " << i << " " << dosage_[i] << " " << gt << " " << get_base_log_likelihood(base, quality, i, log_correct, log_error) << "\n";
  std::cerr << std::endl;
}

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

void printCigarString(BamTools::BamAlignment& aln, std::ostream& out){
  for (unsigned int i = 0; i < aln.CigarData.size(); ++i)
    out << aln.CigarData[i].Length << aln.CigarData[i].Type;
}

inline bool less_than(int32_t v1, int32_t v2)   { return v1 < v2; }
inline bool greater_than(int32_t v1, int32_t v2){ return v1 > v2; }

template <typename StringIterator, typename SNPListIterator, typename CigarOpListIterator>
void extract_bases_and_qualities_helper(BamTools::BamAlignment& aln,    bool forward,
					StringIterator     base_iter,   StringIterator      qual_iter,
					SNPListIterator    snp_iter,    SNPListIterator     snp_end_iter,
					CigarOpListIterator cigar_iter, CigarOpListIterator cigar_end_iter,
					std::vector<char>& bases, std::vector<char>& quals, std::vector<int>& dists){
  int32_t pos        = (forward ? aln.Position : aln.GetEndPosition()-1);
  int32_t dist       = 0;
  const int32_t sign = (forward ? 1 : -1);
  bool (*comparator)(int32_t, int32_t) = (forward ? &less_than : &greater_than);

  while (snp_iter != snp_end_iter && cigar_iter != cigar_end_iter){
    switch(cigar_iter->Type){
    case 'M': case '=': case 'X':
      if (comparator(snp_iter->pos(), pos + sign*cigar_iter->Length)){
	bases.push_back(*(base_iter + sign*(snp_iter->pos() - pos)));
	quals.push_back(*(qual_iter + sign*(snp_iter->pos() - pos)));
	dists.push_back(sign*(snp_iter->pos() - pos) + dist);
	snp_iter++;
      }
      else {
	dist      += cigar_iter->Length;
	pos       += sign*cigar_iter->Length;
	base_iter += cigar_iter->Length;
	qual_iter += cigar_iter->Length;
	cigar_iter++;
      }
      break;
    case 'D':
      if (comparator(snp_iter->pos(), pos + sign*cigar_iter->Length)){
	bases.push_back('-');
	quals.push_back('-');
	dists.push_back(-1);
	snp_iter++;
      }
      else {
	dist = 0;
	pos += sign*cigar_iter->Length;
	cigar_iter++;
      }
      break;
    case 'I':
      dist       = 0;
      base_iter += cigar_iter->Length;
      qual_iter += cigar_iter->Length;
      cigar_iter++;
      break;
    case 'S':
      if (comparator(snp_iter->pos(), pos)){
	bases.push_back('-'); // Ignore bases in soft clips by marking them as dummy deletions
	quals.push_back('-');
	dists.push_back(-1);
	snp_iter++;
      }
      else {
	dist       = 0;
	base_iter += cigar_iter->Length;
	qual_iter += cigar_iter->Length;
	cigar_iter++;
      }
      break;
    case 'H':
      dist = 0;
      cigar_iter++;
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered");
      break;
    }
  }
}




void extract_bases_and_qualities(BamTools::BamAlignment& aln, std::vector<SNP>& snps,
                                 std::vector<char>& bases,    std::vector<char>& quals, std::vector<int>& dists){
  assert(bases.size() == 0 && quals.size() == 0 && dists.size() == 0);
  assert(aln.CigarData.size() > 0);
  assert(std::is_sorted(snps.begin(), snps.end(), SNPSorter()));

  std::vector<char> fw_bases, fw_quals, bw_bases, bw_quals;
  std::vector<int>  fw_dists, bw_dists;
  extract_bases_and_qualities_helper(aln, true,  aln.QueryBases.begin(),  aln.Qualities.begin(),  snps.begin(),  snps.end(),  aln.CigarData.begin(),  aln.CigarData.end(),  fw_bases, fw_quals, fw_dists);
  extract_bases_and_qualities_helper(aln, false, aln.QueryBases.rbegin(), aln.Qualities.rbegin(), snps.rbegin(), snps.rend(), aln.CigarData.rbegin(), aln.CigarData.rend(), bw_bases, bw_quals, bw_dists);

  assert(fw_bases.size() == snps.size() && fw_quals.size() == snps.size() && fw_dists.size() == snps.size());
  assert(bw_bases.size() == snps.size() && bw_quals.size() == snps.size() && bw_dists.size() == snps.size());

  // SNP information is extracted in the reverse direction, so reverse it to make it match the forward information
  std::reverse(bw_bases.begin(), bw_bases.end());
  std::reverse(bw_quals.begin(), bw_quals.end());
  std::reverse(bw_dists.begin(), bw_dists.end());

  for (unsigned int i = 0; i < snps.size(); i++){
    assert(fw_bases[i] == bw_bases[i] && fw_quals[i] == bw_quals[i]);
    bases.push_back(fw_bases[i]);
    quals.push_back(fw_quals[i]);
    dists.push_back(std::min(fw_dists[i], bw_dists[i]));
  }
}
