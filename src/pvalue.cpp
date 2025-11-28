#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector pvalue(SEXP a, SEXP b) {
  NumericVector observed(a);
  NumericVector permuted(b);
  NumericVector pvalues(observed.length());
  
  int j = 0;
  for (int i=0; i<observed.length(); i++) {
    while(j<permuted.length() && permuted[j]>=observed[i]) {
      j++;
    }
    pvalues[i] = double(j) / double(permuted.length());
  }

  return pvalues;
}
