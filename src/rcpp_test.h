#ifndef test_H_INCLUDED
#define test_H_INCLUDED

#include <Rcpp.h>
using namespace Rcpp;

// SEXP _BioEcoPar_ratio(SEXP numSEXP, SEXP denumSEXP);
SEXP _BioEcoPar_stripe_ecoCDF(SEXP out_EcoDCFSEXP, SEXP FlistSEXP, int ind_t);
SEXP _BioEcoPar_init_dataframe(int nrow, int ncol,  SEXP dimnamesSEXP, SEXP DimCstSEXP);

#endif // test_H_INCLUDED

