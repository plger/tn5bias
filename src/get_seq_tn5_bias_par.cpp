// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
#include <map>
#include <string>
#include <cmath>

std::string reverse_complement(const std::string& seq) {
  std::unordered_map<char, char> complement = {
    {'A', 'T'}, {'T', 'A'},
    {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'},
    {'c', 'g'}, {'g', 'c'},
    {'N', 'N'}, {'n', 'n'}
  };
  
  std::string rev_comp;
  rev_comp.reserve(seq.size());
  
  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    char base = *it;
    rev_comp += complement.count(base) ? complement[base] : 'N';
  }
  
  return rev_comp;
}

std::unordered_map<std::string, double> build_bias_map(DataFrame bias_table){
  // Build bias map
  std::unordered_map<std::string, double> bias_map;
  CharacterVector kmers = bias_table["kmer"];
  NumericVector biases = bias_table["bias"];
  for (int i = 0; i < kmers.size(); ++i) {
    bias_map[as<std::string>(kmers[i])] = biases[i];
  }
  return bias_map;
}

double get_pos_bias(std::size_t i, const std::string& genome, 
                    const std::unordered_map<std::string, double>& bias_map,
                    long unsigned int genome_length){
  std::string kmer = genome.substr(i - 5, 10);
  double bias1 = bias_map.count(kmer) ? bias_map.at(kmer) : 1;
  if (i < (genome_length-9)) {
    kmer = reverse_complement(genome.substr(i + 4, 10));
    double bias2 = bias_map.count(kmer) ? bias_map.at(kmer) : 1;
    bias1 = std::sqrt(bias1 * bias2);
  }
  return bias1;
}


// [[Rcpp::export]]
NumericVector get_seq_tn5_bias(std::string genome, DataFrame bias_table) {
  
  long unsigned int genome_length = genome.size();
  std::unordered_map<std::string, double> bias_map = build_bias_map(bias_table);
  
  NumericVector bias(genome_length);
  
  for(std::size_t i = 5; i < (genome_length-5); ++i){
    bias[i] = get_pos_bias(i, genome, bias_map, genome_length);
  }
  
  return bias;
}


// Worker class
struct BiasWorker : public Worker {
  const std::string& genome;
  const std::unordered_map<std::string, double>& bias_map;
  long unsigned int genome_length;
  RVector<double> bias;
  
  BiasWorker(const std::string& genome,
             const std::unordered_map<std::string, double>& bias_map,
             long unsigned int genome_length,
             NumericVector bias)
    : genome(genome), bias_map(bias_map), genome_length(genome_length), bias(bias) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      bias[i] = get_pos_bias(i, genome, bias_map, genome_length);
    }
  }
};

// [[Rcpp::export]]
NumericVector get_seq_tn5_bias_par(std::string genome, DataFrame bias_table) {
  
  long unsigned int genome_length = genome.size();
  std::unordered_map<std::string, double> bias_map = build_bias_map(bias_table);
  
  NumericVector bias(genome_length);

  // Run parallel worker
  BiasWorker worker(genome, bias_map, genome_length, bias);
  parallelFor(5, genome_length - 5, worker);
  
  return bias;
}
