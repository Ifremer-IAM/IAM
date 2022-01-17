#include <Rcpp.h>
#include "rcpp_test.h"

using namespace Rcpp;

// NumericVector ratio(NumericVector num, NumericVector denum){
//     NumericVector res = num / denum;
//     return(res);
// }

// SEXP _BioEcoPar_ratio(SEXP numSEXP, SEXP denumSEXP) {
// BEGIN_RCPP
//     Rcpp::RObject rcpp_result_gen;
//     Rcpp::RNGScope rcpp_rngScope_gen;
//     Rcpp::traits::input_parameter< NumericVector >::type num(numSEXP);
//     Rcpp::traits::input_parameter< NumericVector >::type denum(denumSEXP);
//     rcpp_result_gen = Rcpp::wrap(ratio(num, denum));
//     return rcpp_result_gen;
// END_RCPP
// }

List stripe_ecoCDF(List EcoDCF, List Flist, int ind_t){

    CharacterVector dimCstF = Rf_install("DimCst");

    NumericVector nbv = Flist["nbv_f"]; 
    NumericVector K = Flist["K_f"];
    NumericVector inv = Flist["inv_f"]; 

    DataFrame GVLav (EcoDCF["GVLav_f_out"]);     NumericVector GVLav_t = GVLav[ind_t];
    DataFrame gp (EcoDCF["gp_f_out"]);           NumericVector gp_t = gp[ind_t];
    DataFrame ncf (EcoDCF["ncf_f_out"]);         NumericVector ncf_t = ncf[ind_t];
    DataFrame np (EcoDCF["np_f_out"]);           NumericVector np_t = np[ind_t];

    DataFrame ratio_GVL_K(EcoDCF["ratio_GVL_K_f_out"]);
    DataFrame ratio_gp_K(EcoDCF["ratio_gp_K_f_out"]);
    DataFrame RoFTA(EcoDCF["RoFTA_f_out"]);
    DataFrame ROI(EcoDCF["ROI_f_out"]);
    DataFrame ratio_np_K(EcoDCF["ratio_np_K_f_out"]);
    
    //-- 43. ratio_GVL_K_f
    ratio_GVL_K[ind_t] = GVLav_t / K;
    //-- 44. ratio_gp_K_f
    ratio_gp_K[ind_t] = gp_t / K;
    //-- 45. RoFTA_f
    RoFTA[ind_t] = ncf_t / K;
    //-- 46. ROI_f
    ROI[ind_t] = (gp_t - inv) / inv;
    //-- 47. ratio_np_K_f
    ratio_np_K[ind_t] = np_t / K;

    List res = List::create(
        Named("nothing") = K,
        _["ratio_GVL_K"] = ratio_GVL_K, _["ratio_gp_K"] = ratio_gp_K, _["RoFTA"] = RoFTA, _["ROI"] = ROI, _["ratio_np_K"] = ratio_np_K
    );
    return(res);
}

SEXP _BioEcoPar_stripe_ecoCDF(SEXP out_EcoDCFSEXP, SEXP FlistSEXP, int ind_t) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type out_EcoDCF(out_EcoDCFSEXP);
    Rcpp::traits::input_parameter< List >::type Flist(FlistSEXP);
    rcpp_result_gen = Rcpp::wrap(stripe_ecoCDF(out_EcoDCF, Flist, ind_t));
    return rcpp_result_gen;
END_RCPP
}



DataFrame init_dataframe(int nrow, int ncol, List dimnames, IntegerVector DimCst){

    List result(ncol);
    NumericVector empty = rep(R_NaN, nrow);
    CharacterVector cnms = dimnames[1], rnms = dimnames[0];

    for(int i = 0; i < ncol; i++){
        result[i] = empty;
    }

    DataFrame res(result);
    res.attr("DimCst") = DimCst;
    res.attr("names") = cnms;
    res.attr("row.names") = rnms;
    res.attr("DimCst") = DimCst;
    return(res);
}

SEXP _BioEcoPar_init_dataframe(int nrow, int ncol,  SEXP dimnamesSEXP, SEXP DimCstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    // Rcpp::traits::input_parameter< NumericVector >::type nrow(nrowSEXP);
    // Rcpp::traits::input_parameter< NumericVector >::type ncol(ncolSEXP);
    Rcpp::traits::input_parameter< List >::type dimnames(dimnamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type DimCst(DimCstSEXP);
    rcpp_result_gen = Rcpp::wrap(init_dataframe(nrow, ncol, dimnames, DimCst));
    return rcpp_result_gen;
END_RCPP
}


