#include <climits>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"

#include "em_analyzer.h"
#include "error.h"
#include "stringops.h"
#include "vcf_reader.h"
#include "version.h"

bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

void print_usage(){
  std::cerr << "Usage: EMContaminator --bam <reads.bam> --vcf <snps.vcf.gz>" << "\n"
	    << "\t" << "--bam   <reads.bam>    "  << "\t" << "BAM file containing reads from sequencing run" << "\n"
	    << "\t" << "--vcf   <snps.vcf.gz>  "  << "\t" << "Bgzipped VCF file containing SNP genotypes for a population of individuals" << "\n\n"
	    << "Optional arguments:" << "\n"
	    << "\t" << "--help                 "  << "\t" << "Print this help message and exit" << "\n"
	    << "\t" << "--version              "  << "\t" << "Print the version and exit" << "\n" << "\n";
}
  
void parse_command_line_args(int argc, char** argv, std::string& bam_file, std::string& snp_vcf_file){
  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage();
    exit(0);
  }
  
  int print_help = 0, print_version  = 0;
  static struct option long_options[] = {
    {"bam",            required_argument, 0, 'b'},
    {"h",              no_argument, &print_help, 1},
    {"help",           no_argument, &print_help, 1},
    {"vcf",           required_argument, 0, 'v'}, 
    {"version",        no_argument, &print_version, 1},
    {0, 0, 0, 0}
  };

  std::string filename;
  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "b:v:", long_options, &option_index);
    if (c == -1)
      break;

    if (optarg != NULL){
      std::string val(optarg);
      if (string_starts_with(val, "--"))
	printErrorAndDie("Argument to option --" + std::string(long_options[option_index].name) + " cannot begin with \"--\"\n\tBad argument: " + val);
    }

    switch(c){
    case 0:
      break;
    case 'b':
      bam_file = std::string(optarg);
      break;
    case 'v':
      snp_vcf_file = std::string(optarg);
      break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }

  if (optind < argc) {
    std::stringstream msg;
    msg << "Did not recognize the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./HipSTR --help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_version == 1){
    std::cerr << "EMContaminator version " << VERSION << std::endl;
    exit(0);
  }
  if (print_help){
    print_usage();
    exit(0);
  }
}

int main(int argc, char** argv){
  double total_time = clock();

  // Parameters to consider
  bool init_randomly = false;  // If true, randomly initialize contamination fractions. Otherwise, use a uniform prior
  double error_rate  = 0.01;   // Error rate of sequencing data

  std::stringstream full_command_ss;
  full_command_ss << "EMContaminator-" << VERSION;
  for (int i = 1; i < argc; i++)
    full_command_ss << " " << argv[i];
  std::string full_command = full_command_ss.str();

  std::string bam_file = "", snp_vcf_file = "";
  parse_command_line_args(argc, argv, bam_file, snp_vcf_file);
  if (bam_file.empty())
    printErrorAndDie("You must specify the --bam option");
  if (snp_vcf_file.empty())
    printErrorAndDie("You must specify the --vcf option");

  // Open all BAM files
  BamTools::BamMultiReader reader;
  std::vector<std::string> bam_files;
  bam_files.push_back(bam_file);
  if (!reader.Open(bam_files)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM files");
  }

  if (!string_ends_with(snp_vcf_file, ".gz"))
    printErrorAndDie("SNP VCF file must be bgzipped (and end in .gz)");
  
  // Check that the VCF exists
  if (!file_exists(snp_vcf_file))
    printErrorAndDie("SNP VCF file " + snp_vcf_file + " does not exist. Please ensure that the path provided to --snp-vcf is valid");

  // Check that tabix index exists
  if (!file_exists(snp_vcf_file + ".tbi"))
    printErrorAndDie("No .tbi index found for the SNP VCF file. Please index using tabix and rerun HipSTR");

  // Run analysis
  EMAnalyzer em_analyzer(snp_vcf_file, init_randomly);
  if (em_analyzer.analyze(reader, error_rate)){
    std::cout << "Sample fraction estimates" << "\n";
    em_analyzer.print_sample_priors(std::cout);
  }
  else
    std::cout << "Sample fraction estimation failed to converge" << "\n";

  reader.Close();
  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  return 0;  
}
