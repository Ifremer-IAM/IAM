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
    NumericVector effort1_f = Flist["effort1_f"]; NumericVector effort2_f = Flist["effort2_f"];
    // NumericVector cshr = Flist["cshr_f"]; 
    // NumericVector cnb = Flist["cnb_f"]; 
    // NumericVector persc = Flist["persc_f"]; 
    // NumericVector eec = Flist["eec_f"]; 
    // NumericVector mwh = Flist["mwh_f"]; 
    // // NumericVector rep = Flist["rep_f"]; 
    // // NumericVector gc = Flist["gc_f"]; 
    // // NumericVector fixc = Flist["fixc_f"]; 
    // NumericVector FTE = Flist["FTE_f"]; 
    // NumericVector dep = Flist["dep_f"]; 
    // NumericVector ic = Flist["ic_f"]; 

    NumericMatrix cnb_out = EcoDCF["cnb_f_out"];
    NumericMatrix GVLav = EcoDCF["GVLav_f_out"];
    NumericMatrix gp = EcoDCF["gp_f_out"];
    NumericMatrix ncf = EcoDCF["ncf_f_out"];
    NumericMatrix np = EcoDCF["np_f_out"];

    NumericMatrix npmargin = EcoDCF["npmargin_f_out"];
    NumericMatrix prof = EcoDCF["prof_f_out"];
    NumericMatrix npmargin_trend = EcoDCF["npmargin_trend_f_out"];
    NumericMatrix ssTot = EcoDCF["ssTot_f_out"];
    NumericMatrix ratio_GVL_K = EcoDCF["ratio_GVL_K_f_out"];
    NumericMatrix ratio_gp_K = EcoDCF["ratio_gp_K_f_out"];
    NumericMatrix RoFTA = EcoDCF["RoFTA_f_out"];
    NumericMatrix ROI = EcoDCF["ROI_f_out"];
    NumericMatrix ratio_np_K = EcoDCF["ratio_np_K_f_out"];
    NumericMatrix ratio_GVL_cnb_ue = EcoDCF["ratio_GVL_cnb_ue_f_out"]; 
    
    //-- 32. npmargin_f
    npmargin( _ , ind_t) = np( _ , ind_t) / GVLav( _ , ind_t);
    //-- 33. prof_f
    NumericVector prof_t = npmargin( _ , ind_t);
    prof_t[prof_t < 0] = -1; prof_t[prof_t > 0.1] = 1; 
    prof_t[(prof_t >= 0) & (prof_t < 0.1)] = 0;
    prof( _ , ind_t) = prof_t;
    //-- 34. npmargin_trend_f
    NumericVector npmargin_trend_t = npmargin( _ , ind_t);
    if(ind_t >= 5){
        NumericMatrix sub_npmargin = npmargin( _ , Range(ind_t-5, ind_t -1));
        npmargin_trend_t = npmargin_trend_t / (0.2 * rowSums(sub_npmargin));

        npmargin_trend_t[npmargin_trend_t < -0.05] = -1; npmargin_trend_t[npmargin_trend_t > 0.05] = 1; 
        npmargin_trend_t[(npmargin_trend_t >= -0.05) & (npmargin_trend_t <= 0.05)] = 0; 
    } else {
        npmargin_trend_t = rep(-1, npmargin_trend.nrow()); // this is nbF
    }
    npmargin_trend( _ , ind_t) = npmargin_trend_t;
    //-- 35. ssTot_f
    ssTot( _ , ind_t) = gp( _ , ind_t) * nbv;
    //-- 43. ratio_GVL_K_f
    ratio_GVL_K( _ , ind_t) = GVLav( _ , ind_t) / K;
    //-- 44. ratio_gp_K_f
    ratio_gp_K( _ , ind_t) = gp( _ , ind_t) / K;
    //-- 45. RoFTA_f
    RoFTA( _ , ind_t) = ncf( _ , ind_t) / K;
    //-- 46. ROI_f
    ROI( _ , ind_t) = (gp( _ , ind_t) - inv) / inv;
    //-- 47. ratio_np_K_f
    ratio_np_K( _ , ind_t) = np( _ , ind_t) / K;
    //-- 48. ratio_GVL_cnb_ue_f
    ratio_GVL_cnb_ue( _ , ind_t) = GVLav( _ , ind_t) / (cnb_out( _ , ind_t) * effort1_f * effort2_f);

    List res = List::create(
        Named("nothing") = K, 
        _["ssTot"] = ssTot,
        _["ratio_GVL_K"] = ratio_GVL_K, _["ratio_gp_K"] = ratio_gp_K, _["RoFTA"] = RoFTA, _["ROI"] = ROI, _["ratio_np_K"] = ratio_np_K,
         _["ratio_GVL_cnb_ue"] = ratio_GVL_cnb_ue
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


NumericMatrix init_matrix(int nrow, int ncol, List dimnames, IntegerVector DimCst){

    NumericVector empty = rep(R_NaN, nrow*ncol);
    CharacterVector cnms = dimnames[1], rnms = dimnames[0];

    NumericMatrix res( nrow , ncol , empty.begin() );
    rownames(res) = rnms;
    colnames(res) = cnms;
    res.attr("DimCst") = DimCst;
    return(res);
}

SEXP _BioEcoPar_init_matrix(int nrow, int ncol,  SEXP dimnamesSEXP, SEXP DimCstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type dimnames(dimnamesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type DimCst(DimCstSEXP);
    rcpp_result_gen = Rcpp::wrap(init_matrix(nrow, ncol, dimnames, DimCst));
    return rcpp_result_gen;
END_RCPP
}


IntegerVector init_DimCst(int nbF, int nbM, int nbI, int nbT){
    // IntegerVector x = {nbF, nbM, nbI, nbT};
    IntegerVector x = IntegerVector::create(nbF, nbM, nbI, nbT);
    return(x);
}

SEXP _BioEcoPar_init_DimCst(int nbF, int nbM, int nbI, int nbT) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(init_DimCst(nbF, nbM, nbI, nbT));
    return rcpp_result_gen;
END_RCPP
}

IntegerVector n_DimCst(IntegerVector DimCst){
    IntegerVector res = DimCst[DimCst > 0];
    return(res);
}

SEXP _BioEcoPar_n_DimCst(SEXP DimCstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type DimCst(DimCstSEXP);
    rcpp_result_gen = Rcpp::wrap(n_DimCst(DimCst));
    return rcpp_result_gen;
END_RCPP
}
