#include <assert.h>
#include <iostream>

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

void printCigarString(BamTools::BamAlignment& aln, std::ostream& out){
  for (unsigned int i = 0; i < aln.CigarData.size(); ++i)
    out << aln.CigarData[i].Length << aln.CigarData[i].Type;
}

void extract_bases_and_qualities(BamTools::BamAlignment& aln, std::vector<SNP>& snps,
				 std::vector<char>& bases,    std::vector<char>& quals){
  assert(bases.size() == 0 && quals.size() == 0);
  assert(aln.CigarData.size() > 0);
  assert(std::is_sorted(snps.begin(), snps.end(), SNPSorter()));

  int32_t pos = aln.Position;
  unsigned int snp_index = 0, cigar_index = 0, base_index = 0;
  while (snp_index < snps.size() && cigar_index < aln.CigarData.size()){
    switch(aln.CigarData[cigar_index].Type){
    case 'M': case '=': case 'X':
      if (snps[snp_index].pos() < pos + aln.CigarData[cigar_index].Length){
	bases.push_back(aln.QueryBases.at(snps[snp_index].pos() - pos + base_index));
	quals.push_back(aln.Qualities.at(snps[snp_index].pos() - pos + base_index));
	snp_index++;
      }
      else {
	pos += aln.CigarData[cigar_index].Length;
	base_index += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'D':
      if (snps[snp_index].pos() < pos + aln.CigarData[cigar_index].Length){
	bases.push_back('-');
	quals.push_back('-');
	snp_index++;
      }
      else {
	pos += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'I':
      base_index += aln.CigarData[cigar_index].Length;
      cigar_index++;
      break;
    case 'S':
      if (snps[snp_index].pos() < pos){
	bases.push_back('-'); // Ignore bases in soft clips by marking them as dummy deletions
	quals.push_back('-');
	snp_index++;
      }
      else {
	base_index += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'H':
      cigar_index++;
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered");
      break;
    }
  }
  assert(bases.size() == snps.size() && snp_index == snps.size());
}
