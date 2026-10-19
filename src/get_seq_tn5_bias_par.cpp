#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

inline int base_to_int(char c) {
  switch(c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1; 
  }
}

// Follows SELMA's bias_exp_cleavage_ATAC: for a cut at (0-based) position g,
// the bias is the geometric mean of two 10-mers, because the two cleavage
// sites of a Tn5 dimer are 9 bp apart:
//  + strand: kmer [g-5, g+4] and the reverse complement of kmer [g+4, g+13]
//  - strand: reverse complement of kmer [g-4, g+5] and kmer [g-13, g-4]
// The - strand cut at g+8 has the same bias as the + strand cut at g.
struct BiasWorker : public Worker {
  const std::string& genome;
  const std::vector<double>& bias_vec;
  const std::vector<int>& rc_map;
  bool shift;
  bool complement;
  RVector<double> output;
  
  BiasWorker(const std::string& genome, const std::vector<double>& bias_vec, 
             const std::vector<int>& rc_map, bool shift, bool complement, NumericVector output)
    : genome(genome), bias_vec(bias_vec), rc_map(rc_map), 
      shift(shift), complement(complement), output(output) {}
  
  // bias of the kmer starting at `start` (1 if it contains non-ACGT), 
  // optionally of its reverse complement
  double kmer_bias(int start, bool rc) const {
    int k_idx = 0;
    for (int j = 0; j < 10; ++j) {
      int val = base_to_int(genome[start + j]);
      if (val == -1) return 1.0;
      k_idx = (k_idx << 2) | val;
    }
    return rc ? bias_vec[rc_map[k_idx]] : bias_vec[k_idx];
  }
  
  double plus_cut(int g) const {
    return std::sqrt(kmer_bias(g - 5, false) * kmer_bias(g + 4, true));
  }
  
  double minus_cut(int g) const {
    return std::sqrt(kmer_bias(g - 4, true) * kmer_bias(g - 13, false));
  }
  
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int g = (int)i;
      if (shift) {
        // i is the insertion site, i.e. 4 bp after the + strand cut
        output[i] = plus_cut(g - 4);
      } else {
        output[i] = complement ? minus_cut(g) : plus_cut(g);
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector get_seq_tn5_bias_par(std::string genome, NumericVector bias_vec, 
                                   IntegerVector rc_map, bool shift, bool complement) {
  int n = genome.size();
  NumericVector output(n, 1.0);
  std::vector<double> std_bias = as<std::vector<double>>(bias_vec);
  std::vector<int> std_rc = as<std::vector<int>>(rc_map);
  
  BiasWorker worker(genome, std_bias, std_rc, shift, complement, output);
  
  // Boundary safety: the widest window is 17 bases around the position
  // (shifted: i-4-5 .. i-4+13; minus strand: i-13 .. i+5).
  const int margin = 18;
  if (n > 2 * margin) parallelFor(margin, n - margin, worker);
  
  return output;
}
