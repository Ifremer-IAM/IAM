// #include <stdlib.h>
// #include <stdio.h>
// #include <time.h>
// #include <vector>
// #include <math.h>
// #include <string>
// #include <sstream>
// #include <iostream>
// #include <fstream>
// #include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
//#include <Rcpp.h>

#include "BioEcoPar.h" // Class is defined in this file.
// #include "Modules.h" // Contient tout les modules

//using namespace Rcpp;
// using namespace std;



//---------------------------------
//
// Module de marche de quotas
//
//---------------------------------

extern "C" {

void BioEcoPar::QuotaMarket(SEXP list, SEXP pQuotaIni, SEXP pQuotaMin, SEXP pQuotaMax, double lambdaQ, double sdmax, double ftol, int itmax, SEXP paramBehav, int ind_t, int persCalc ) //ind_t>0
{

// ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.QuotaMarket.txt", ios::out | ios::trunc);
//ofstream fichier;
//if (ind_t ==1) fichier.open ("C:\\Users\\fbriton\\Dropbox\\These\\IAM_Dvt\\test.QuotaMarket.txt");
//fichier << "Start" << endl;
 //time_t my_time;

    SEXP listTemp, eVarCopy, alphaBhv, pQuota ,nam_eQuota,nam_eQuota_dyn,
        dimCstF,dimCstFM, nDimT,nDimT2, eFACTf, eFACTfm;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_effSup = REAL(effSupMat); // specified in objArgs@arguments$Gestion$effSup = relates to effort1*effort2
//fichier << "QM0.1"  << endl;
    double *r_out_effort1_f = REAL(NBDSF); // efforts from output for past values
    double *r_out_effort2_f = REAL(EFF2F);
    double *r_out_effort1_f_m = REAL(NBDSFM);
    double *r_out_effort2_f_m = REAL(EFF2FM);
//fichier << "QM0.2"  << endl;

    double *g_effort1FM = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f_m")); // efforts from list to update to then send to other modules
    double *g_effort1F = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f"));
//fichier << "QM0.3"  << endl;
    double* r_L_f_m_e;

    SEXP rtbs_f_m_out, rtbs_f_out, ccw_f_out, rep_f, gc_f, fixc_f, dep_f, ic_f, GVLtot_f_m_out, cshrT_f_m_out,GVLtot_f_out;
    double *r_rtbs_f_m_out, *r_rtbs_f_out, *r_ccw_f_out, *r_rep_f, *r_gc_f, *r_fixc_f, *r_dep_f, *r_ic_f, *r_GVLtot_f_m_out, *r_GVLtot_f_m_e_ref, *r_P_f_m_e,*r_cshrT_f_m_out,*r_GVLtot_f_out;
    PROTECT(rep_f = getListElement(getListElement(listTemp, "Fleet"), "rep_f"));
    r_rep_f = REAL(rep_f);
    PROTECT(gc_f = getListElement(getListElement(listTemp, "Fleet"), "gc_f"));
    r_gc_f = REAL(gc_f);
    PROTECT(fixc_f = getListElement(getListElement(listTemp, "Fleet"), "fixc_f"));
    r_fixc_f = REAL(fixc_f);
    PROTECT(dep_f = getListElement(getListElement(listTemp, "Fleet"), "dep_f"));
    r_dep_f = REAL(dep_f);
    PROTECT(ic_f = getListElement(getListElement(listTemp, "Fleet"), "ic_f"));
    r_ic_f = REAL(ic_f);
    PROTECT(GVLtot_f_m_out = VECTOR_ELT(out_EcoDCF, 6));
    r_GVLtot_f_m_out = REAL(GVLtot_f_m_out);
    PROTECT(rtbs_f_m_out = VECTOR_ELT(out_EcoDCF, 15));
    r_rtbs_f_m_out = REAL(rtbs_f_m_out);
    PROTECT(ccw_f_out = VECTOR_ELT(out_EcoDCF, 29));
    r_ccw_f_out = REAL(ccw_f_out);
    PROTECT(cshrT_f_m_out = VECTOR_ELT(out_EcoDCF, 18));
    r_cshrT_f_m_out = REAL(cshrT_f_m_out);
    PROTECT(rtbs_f_out = VECTOR_ELT(out_EcoDCF, 16));
    r_rtbs_f_out = REAL(rtbs_f_out);

    SEXP dc_rep_f ,dc_gc_f,dc_fixc_f,dc_dep_f;
    int *dim_rep_f,*dim_gc_f,*dim_fixc_f,*dim_dep_f;
    PROTECT(dc_rep_f = iDim(INTEGER(getAttrib(rep_f, install("DimCst")))));
    PROTECT(dc_gc_f = iDim(INTEGER(getAttrib(gc_f, install("DimCst")))));
    PROTECT(dc_fixc_f = iDim(INTEGER(getAttrib(fixc_f, install("DimCst")))));
    PROTECT(dc_dep_f = iDim(INTEGER(getAttrib(dep_f, install("DimCst")))));
    dim_rep_f = INTEGER(dc_rep_f);
    dim_gc_f = INTEGER(dc_gc_f);
    dim_fixc_f = INTEGER(dc_fixc_f);
    dim_dep_f = INTEGER(dc_dep_f);

    SEXP PQuot_temp, diffLQ;
    double *r_PQuot_temp, *r_diffLQ;


//fichier << "QM0.4"  << endl;

    SEXP ProfUE_f_m, ProfUE_f_m_ctr, ProfUE_f, ExpProfUE_f, Profmin_f;
    PROTECT(ProfUE_f_m = NEW_NUMERIC(nbF*nbMe));
    double *r_ProfUE_f_m = REAL(ProfUE_f_m);
    PROTECT(ProfUE_f_m_ctr = NEW_NUMERIC(nbF*nbMe));
    double *r_ProfUE_f_m_ctr = REAL(ProfUE_f_m_ctr);
    PROTECT(ProfUE_f = NEW_NUMERIC(nbF));
    double *r_ProfUE_f = REAL(ProfUE_f);
    double *r_allocEff_f_m = REAL(out_allocEff_fm);
    PROTECT(ExpProfUE_f = NEW_NUMERIC(nbF));
    double *r_ExpProfUE_f = REAL(ExpProfUE_f);
//fichier << "QM0.4"  << endl;
    double *r_GoFish_f = REAL(intermGoFish);
//    PrintValue(intermGoFish);
    PROTECT(Profmin_f = NEW_NUMERIC(nbF));
    double *r_Profmin_f = REAL(Profmin_f);
    double ratio_p;
    double *rans_multPrice;
    double GVLtot_temp_f_m;
//fichier << "QM0.5"  << endl;
    int ind_t_last, ind_t_price;

    PROTECT(alphaBhv = getListElement(paramBehav, "ALPHA"));
    double *r_alphaBhv = REAL(alphaBhv);
//fichier << "QM0.6"  << endl;

    PROTECT(dimCstF = allocVector(INTSXP, 4));
    PROTECT(dimCstFM = allocVector(INTSXP, 4));
    int *dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
    int *dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;
    PROTECT(eFACTf = iDim(dCF));
    PROTECT(eFACTfm = iDim(dCFM));
    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    //fichier << "dCFM[0]: " << dCFM[0] << "; dCFM[1]: " << dCFM[1] << "; dCFM[2]: " << dCFM[2] << "; dCFM[3]: " << dCFM[3] << endl;
    //fichier << "eF_fm[0]: " << eF_fm[0] << "; eF_fm[1]: " << eF_fm[1] << "; eF_fm[2]: " << eF_fm[2] << "; eF_fm[3]: " << eF_fm[3] << endl;

    PROTECT(nDimT = allocVector(INTSXP,4));
    int *ndT = INTEGER(nDimT); ndT[0] = 0; ndT[1] = 0; ndT[2] = 0; ndT[3] = nbT;
    PROTECT(nDimT2 = allocVector(INTSXP,2));
    int *ndT2 = INTEGER(nDimT2); ndT2[0] = 0; ndT2[1] = nbT;

    double *r_pQuota;
    double* r_pQuotaIni;
    SEXP ans_multPrice;
    int ind_ePrice;
    SEXP Lmodel,TACmodel;
    double *r_Lmodel, *r_TACmodel  ;
    SEXP reducelambda;
    PROTECT(reducelambda = NEW_LOGICAL(nbEQuotaMarket_dyn));
    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) LOGICAL(reducelambda)[int_eQuota_dyn] = false;

    double lambdaQ_iter;

    //---------------------------------------------
    // 1 - Tatonnement of quota market
    //---------------------------------------------

    // Initialize quota price
//    fichier << "QM0.7"  << endl;
    //fichier << "nbEQuotaMarket:" << nbEQuotaMarket  << endl;

    // year used to calculate ProfUE_fm and base fishing decisions:
    ind_t_last=0 ; // initial year
    //if you want it to be the last year the vessel has been active, uncomment lines 14774-14783

    // year used to retrieve fish prices and calculate GVL (mostly to account for dynamics in fish price while using ind_t_last=0
    if (ind_t == 1){ind_t_price = 0;} else {ind_t_price = 1;} // to capture market dynamics
    //ind_t_price = ind_t-1; // previous year (less stable than previous option)
//    fichier << "ind_t_price: " << ind_t_price << endl;

    for (int int_eQuota = 0 ; int_eQuota  < nbEQuotaMarket ; int_eQuota++) { // initialize quota prices
            ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListQM,int_eQuota))); // index of species in MultPrice
            PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
            PROTECT(pQuota = VECTOR_ELT(out_PQuot_et,int_eQuota));
            r_pQuota = REAL(pQuota);
            r_pQuotaIni = REAL(VECTOR_ELT(pQuotaIni,int_eQuota));
            r_pQuota[ind_t] = r_pQuotaIni[ind_t] * REAL(ans_multPrice)[ind_t_price]; // adjust quota price for change in fish price in previous year (only value when entering the loop, but to keep ratio quota price/fish price stable)
            UNPROTECT(2);
    }
//                PrintValue(sppListAll);
//                PrintValue(multPrice);
//                PrintValue(sppListQM);
//                PrintValue(out_PQuot_et);

//    fichier << "QM0.8"  << endl;

    bool toinit;
    bool keep;
    bool GoOn = true;

    int itQ=0;
    double min_sum_diffLQ, sum_diffLQ, nb_overTAC, max_diffLQ, min_max_diffLQ, max_nb_overTAC;
    int min_itQ;
    //fichier << "itmax" << itmax  << endl;

    //Calculate L_fme for quoted species
    SEXP L_f_m_e, L_f_m_e_i, ans_L_f_m_e;
    PROTECT(L_f_m_e = allocVector(VECSXP, nbEQuotaMarket));
    setAttrib(L_f_m_e, R_NamesSymbol, sppListQM);

     for (int ind_eQ = 0 ; ind_eQ< nbEQuotaMarket ; ind_eQ++) {
        PROTECT(nam_eQuota = STRING_ELT(sppListQM,ind_eQ));

        if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))){ // espece dyn
                PROTECT(L_f_m_e_i = getListElement(out_L_efmit, CHAR(nam_eQuota)));
                PROTECT(ans_L_f_m_e = aggregObj(L_f_m_e_i,dimCstFM));
        } else{ // espece stat
                PROTECT(ans_L_f_m_e = getListElement(out_Lstat, CHAR(nam_eQuota)));
        }

        SET_VECTOR_ELT(L_f_m_e, ind_eQ, ans_L_f_m_e);

        UNPROTECT(2);
        if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))) UNPROTECT(1);
     }
     //PrintValue(L_f_m_e);


    while (GoOn & (itQ<itmax)){
//            fichier << "itQ: " << itQ << endl;
            //Rprintf("itQ = %i \n", itQ);

//            fichier << "QM1.1"  << endl;
            // Initialization ProfUE_f_m
            for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = 0.0;

                }
                r_ProfUE_f[ind_f] = 0.0;
                r_ExpProfUE_f[ind_f] = 0.0;
            }
            //fichier << "QM1.2"  << endl;
           //Rprintf("QM1.2\n");


        // Calculate metier profitability for each fleet and metier
       // my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                //if (ind_f==5) fichier << "f: " << CHAR(STRING_ELT(fleetList,ind_f)) << endl;
                toinit = true;

                // retrieve last ind_t when active
//                keep = true;
//                ind_t_last = ind_t-1;
//                while(keep){
//                    if (r_GoFish_f[ind_f+nbF*ind_t_last]>0){
//                            keep = false;
//                    } else{
//                    ind_t_last = ind_t_last -1;}
//
//                }

                //fichier << "QM1.3"  << endl;
                //Rprintf("QM1.3\n");
//                if ((itQ==0) & (ind_t<5)) Rprintf("f: %d, ind_t_last: %d \n", ind_f, ind_t_last);

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//                    if (ind_f==6) fichier << "ind_m: " << ind_m << endl;

                     // Adjust GVL by species to account for market (fish price) dynamics
                        for (int e = 0 ; e < nbE+nbEstat ; e++) {
                                if ((nbE>0) & (e<nbE)) {
                                    r_GVLtot_f_m_e_ref = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
                                    ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e)));
                                    //Rprintf("e = %s, ind_ePrice = %d \n", CHAR(STRING_ELT(sppList,e)),ind_ePrice);
                                    PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
                                    rans_multPrice = REAL(ans_multPrice);

                                } if ((nbEstat>0) & (e>=nbE)) {
                                    r_GVLtot_f_m_e_ref = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));
                                    ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e-nbE)));
                                    //Rprintf("e = %s, ind_ePrice = %d \n", CHAR(STRING_ELT(sppListStat,e-nbE)),ind_ePrice);
                                    PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
                                    rans_multPrice = REAL(ans_multPrice);
                                    }
//                                    fichier << "QM1.4"  << endl;

                        ratio_p = rans_multPrice[ind_t_price] / rans_multPrice[ind_t_last];

                        if (!ISNA(r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]* ratio_p)){
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] +
                                r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * ratio_p;} //adjustment for variations in fish price

                        UNPROTECT(1);

//                        if (ind_f==6) fichier << "e = " <<  e <<
//                            "; GVL init = " << r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]  <<
//                            "; ratio = " << ratio_p <<
//                            "; GVL actual = " << r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * ratio_p <<
//                            "; ProfUE_f_m = " << r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
                        }
                        GVLtot_temp_f_m = r_ProfUE_f_m[ind_f + nbF*ind_m]; // temporary save to use to calculate crew share if persCalc=5

                        // Deduce variable costs : RTBS
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                            (r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] - r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]);
//                       if (ind_f==6) fichier << "Deduce var costs : r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;

                        //Deduce costs
                         // 1- crew costs
                        if ((persCalc == 1) | (persCalc == 2) | (persCalc == 3) | (persCalc == 4)) { // crew share RTBS
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] *
                                        (1 - r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]);
                                        }else if (persCalc == 5){ // crew share GVL
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                        (r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] /
                                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *
                                        GVLtot_temp_f_m);
                        } else { // fixed wages
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                                          r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]  *
                                                          (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) /
                                                          (r_out_effort1_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] * r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]);
                        }

//                        if (ind_f==6) {
//                                if ((persCalc == 1) | (persCalc == 2) | (persCalc == 3) | (persCalc == 4)) {
//                                fichier << " cshrT = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " ; rtbs = " <<  r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " crew share = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl; } else if (persCalc == 5){
//                                    fichier << " cshrT ref = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " ; GVL ref = " <<  r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " crew share = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                "; GVL  =  " <<  GVLtot_temp_f_m <<
//                                    "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                                    } else{
//                                fichier << " crew costs_f =  " <<  r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] <<
//                                "; ratio Effort" << (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) /
//                                                          (r_out_effort1_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] * r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]) <<
//                                "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                                }}

                        // 2- Quota costs

                        for (int ind_eQ = 0 ; ind_eQ< nbEQuotaMarket ; ind_eQ++) { //deduce quota expenses for dynamic markets
                                PROTECT(nam_eQuota = STRING_ELT(sppListQM,ind_eQ));

                                r_L_f_m_e = REAL(getListElement(L_f_m_e, CHAR(nam_eQuota)));

//                                if (ind_f==6) fichier << "nam_eQuota: " << CHAR(nam_eQuota) << endl;


                                if (!ISNA(r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]])){

                                    PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota)));
                                    r_pQuota = REAL(pQuota);

                                    r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m]  -
                                                                r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *1000 * r_pQuota[ind_t];
//                                    if (ind_f==6){
//                                            fichier << "index: " << ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3] << endl;
//                                            fichier << "r_L_f_m_e: " << r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *1000 << endl;
//                                            fichier << "r_pQuota: " << r_pQuota[ind_t] << endl;
//                                            fichier << "Deduce quota costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;}


                                    UNPROTECT(1);
                               }
                               UNPROTECT(1);

                        }


                        // Divide by effort
                        if(r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]>0){
                                r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] /
                                (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) ;
                        }

//                        if (ind_f==6){
//                            fichier << "Effort_f_m = " << r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] << endl;
//                            fichier << "Divide by effort: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                       }



                        // 3- fixed costs
//                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
//                                                          (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
//                                                           r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
//                                                           r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]) /
//                                                          g_effSup[ind_f + nbF*ind_t];

                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                                          (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
                                                           r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
                                                           r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]+
                                                           r_dep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]]+
                                                           r_ic_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]]) /
                                                          g_effSup[ind_f + nbF*ind_t];

//                        if (ind_f==6) fichier << " fixed costs_f =  " <<  r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
//                                                                            r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
//                                                                            r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]+
//                                                                            r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t_last*dim_dep_f[3]]+
//                                                                            r_ic_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] <<
//                                "; Eff sup = " << g_effSup[ind_f + nbF*ind_t] <<
//                                "; Deduce fixed costs per UE: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;




                        if (!ISNA(r_ProfUE_f_m[ind_f + nbF*ind_m])){
                            if (toinit){
                                r_Profmin_f[ind_f] = r_ProfUE_f_m[ind_f + nbF*ind_m];
                                toinit=false;
                            }
                            if(r_ProfUE_f_m[ind_f + nbF*ind_m] < r_Profmin_f[ind_f]) r_Profmin_f[ind_f] =  r_ProfUE_f_m[ind_f + nbF*ind_m];
                        }
                }
        }
//        fichier << "QM1.4"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
        //Rprintf("QM1.4\n");

        // Center Profitabilities from the minimal value so that they are all positive
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m]  - r_Profmin_f[ind_f];
                       if(!ISNA(r_ProfUE_f_m_ctr[ind_f + nbF*ind_m])) r_ProfUE_f[ind_f] = r_ProfUE_f[ind_f] + r_ProfUE_f_m_ctr[ind_f + nbF*ind_m];
            }
            if (r_ProfUE_f[ind_f] ==0.0) { //in case all metiers have a profitability equal to the minimum one
                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       if(r_ProfUE_f_m_ctr[ind_f + nbF*ind_m]==0.0) {
                            r_ProfUE_f_m_ctr[ind_f + nbF*ind_m]=1.0;
                            r_ProfUE_f[ind_f] = r_ProfUE_f[ind_f] + r_ProfUE_f_m_ctr[ind_f + nbF*ind_m];
                       }
                }
            }
        }
//        fichier << "QM2"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
        //Rprintf("QM2\n");

        // Update effort allocation: (weighted average between profitability and habit)
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_alphaBhv[ind_f + nbF*ind_t] * r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] / r_ProfUE_f[ind_f] +
                                                           (1-r_alphaBhv[ind_f + nbF*ind_t]) * (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) /
                                                           (r_out_effort1_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]  * r_out_effort2_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]);
                        if (ISNA(r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;
//                        if(ind_f ==6) fichier << "m= " << ind_m << "; attract = " << r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] / r_ProfUE_f[ind_f] << "; trad = " << (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) /
//                                                          (r_out_effort1_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]  * r_out_effort2_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]) << "; allocEff = " << r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] << endl;
                        if (!ISNA(r_ProfUE_f_m[ind_f + nbF*ind_m])) r_ExpProfUE_f[ind_f] = r_ExpProfUE_f[ind_f] + r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * r_ProfUE_f_m[ind_f + nbF*ind_m];
//                        if (ind_f ==6){ fichier << "r_allocEff_f_m = " << r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] << "; r_ProfUE_f_m = " << r_ProfUE_f_m[ind_f + nbF*ind_m] << "; r_ExpProfUE_f = " << r_ExpProfUE_f[ind_f] << endl;}

                }

          //Decision to go fishing or not based on expected profitability
                if(r_ExpProfUE_f[ind_f]>=0.0) {
                        r_GoFish_f[ind_f+nbF*ind_t] = 1.0;
                } else {
                    r_GoFish_f[ind_f+nbF*ind_t] = 0.0;
                    }
        }
//        fichier << "QM3"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
         //Rprintf("QM3\n");

        // Update effort
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

                       g_effort1FM[ind_f + nbF*ind_m] = r_GoFish_f[ind_f+nbF*ind_t] * g_effSup[ind_f + nbF*ind_t] * r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] / r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]];
//                      if(ind_f ==6) fichier << "m= " << ind_m << "; Eff = " << g_effort1FM[ind_f + nbF*ind_m] << endl;
                }
                g_effort1F[ind_f] = r_GoFish_f[ind_f+nbF*ind_t] * g_effSup[ind_f + nbF*ind_t] / r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];
//                if(ind_f ==6) fichier << "; Eff_f = " << g_effort1F[ind_f] << endl;

        }
        //PrintValue(getListElement(getListElement(listTemp, "Fleet"), "effort1_f"));
        //PrintValue(getListElement(getListElement(listTemp, "Fleet"), "effort1_f_m"));


//       fichier << "QM4.1"  << endl;
        // Rprintf("QM4.1\n");

        Mortalite(listTemp, ind_t, eVarCopy);
//        fichier << "QM4.2"  << endl;
        // Rprintf("QM4.2\n");


        DynamicPop(listTemp, ind_t, eVarCopy,false);
//        fichier << "QM4.3"  << endl;
        // Rprintf("QM4.3\n");

        CatchDL(listTemp, ind_t, eVarCopy, 0);
//        fichier << "QM4.4"  << endl;
        // Rprintf("QM4.4\n");
       // my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;

        // Adjust quota prices for all species
        GoOn = false;
        //lambdaQ_iter = lambdaQ * ( 1 - ((double)rand() / (double)RAND_MAX) * pow(((double)itQ/(double)itmax),2.0) * sdmax);
        lambdaQ_iter = lambdaQ * ( 1 - pow(((double)itQ/(double)itmax),2.0) * sdmax); //increase lambda with itQ to avoid oscillating between the same values towards the end of convergence (optionnal, you can also keep it constant when commenting this line)

        //fichier << "LambdaQ = "<< lambdaQ <<  "; itQ = "  << itQ  <<  "; Lambda iter = "<< lambdaQ_iter <<endl;
        sum_diffLQ = 0.0;
        max_diffLQ=0.0;
        nb_overTAC = 0.0;

        for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {

                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                //fichier << "nam_eQuota_dyn: " << CHAR(nam_eQuota_dyn) << endl;
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));

                if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(Lmodel = aggregObj(getListElement(out_L_eit, CHAR(nam_eQuota_dyn)),nDimT));
                                    //PROTECT(Lmodel = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),nDimT));
                                } else{ // espece stat
                                    PROTECT(Lmodel = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),nDimT));
                                }

                //PrintValue(Lmodel);
                r_Lmodel = REAL(Lmodel);
//                fichier << "QM5.1"  << endl;
                PROTECT(TACmodel = getListElement(TAC, CHAR(nam_eQuota_dyn)));
                //PROTECT(TACmodel = aggregObj(getListElement(TACbyF, CHAR(nam_eQuota_dyn)),nDimT)); //ne fonctionne plus si quota pas detenu par navires modelises
                r_TACmodel = REAL(TACmodel);
                r_pQuota = REAL(pQuota);

//                fichier << "QM5.2"  << endl;

                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                r_PQuot_temp = REAL(PQuot_temp);
                PROTECT(diffLQ = getListElement(out_diffLQ, CHAR(nam_eQuota_dyn)));
                r_diffLQ = REAL(diffLQ);
//                fichier << "QM5.3"  << endl;

                r_diffLQ[itQ + itmax*ind_t] =  (r_Lmodel[ind_t] - r_TACmodel[ind_t]) / r_TACmodel[ind_t]; //save value for current iteration to check algorithm convergence
                sum_diffLQ = sum_diffLQ + fabs(r_diffLQ[itQ + itmax*ind_t]);

//                fichier << "L_eit = " << r_Lmodel[ind_t] << "; L_efmit = " <<  REAL(aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),nDimT))[ind_t] <<
//                 "; TAC = " << r_TACmodel[ind_t] << "; diff = " << r_diffLQ[itQ + itmax*ind_t] << endl;

                if (int_eQuota_dyn ==0){
                      max_diffLQ = fabs(r_diffLQ[itQ + itmax*ind_t]);
                } else if (fabs(r_diffLQ[itQ + itmax*ind_t]) > max_diffLQ) max_diffLQ = fabs(r_diffLQ[itQ + itmax*ind_t]);

                if (r_diffLQ[itQ + itmax*ind_t]>0) nb_overTAC = nb_overTAC+1;
                r_PQuot_temp[itQ + itmax*ind_t] = r_pQuota[ind_t]; //save value for current iteration to check algorithm convergence

//                fichier << "QM5"  << endl;
                    if ((!LOGICAL(reducelambda)[int_eQuota_dyn]) & (r_diffLQ[itQ + itmax*ind_t]*r_diffLQ[itQ-1 + itmax*ind_t] < 0)) LOGICAL(reducelambda)[int_eQuota_dyn] = true;
                    //Rprintf("Reduce lambda at iter %i: %d \n",itQ,LOGICAL(reducelambda)[int_eQuota_dyn]);


                    if(!LOGICAL(reducelambda)[int_eQuota_dyn]){
                        if ((lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t]) > 0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 1.2);
                        } else if ((lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t]) < -0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 0.8);
                        } else {r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * (1 + lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t])); }
                    } else {
                        if ((lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t]) > 0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 1.2);
                        } else if ((lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t]) < -0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 0.8);
                        } else {r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * (1 + lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t])); }}


                    if((fabs(r_diffLQ[itQ + itmax*ind_t]) > ftol) | (itQ < floor(itmax/2))) GoOn=true;

                    //fichier << "itQ: " << itQ << ", diffLQ = " << r_diffLQ[itQ + itmax*ind_t] << "Quota price = "<< r_pQuota[ind_t] << "GoOn =" << GoOn << endl;
                UNPROTECT(6);
        }
//        fichier << "QM5"  << endl;
//        if (itQ == 0){ //floor(itmax/2)
//            min_sum_diffLQ = sum_diffLQ;
//            min_itQ = itQ;
//            max_nb_overTAC = nb_overTAC;
//        } else if (itQ > 0){
//            if(nb_overTAC > max_nb_overTAC) {
//                min_sum_diffLQ = sum_diffLQ;
//                min_itQ = itQ;
//                max_nb_overTAC = nb_overTAC;
//            } else if ((nb_overTAC == max_nb_overTAC) & (sum_diffLQ < min_sum_diffLQ)){
//                min_sum_diffLQ = sum_diffLQ;
//                min_itQ = itQ;
//            }

//fichier << "itQ = "  << itQ << ", max_diffLQ = " << max_diffLQ << ", min_max_diffLQ = " << min_max_diffLQ << endl;
//Save best iteration (the one with the minimum max_diffLQ (=max difference between landing and TAC across stocks)):
        if (itQ == 0){ //floor(itmax/2)
            min_itQ = itQ;
            min_max_diffLQ = fabs(max_diffLQ);
            min_sum_diffLQ = sum_diffLQ;
            } else if (fabs(max_diffLQ) < min_max_diffLQ) {
                //Rprintf("Cas 1 , old diffLQ = %f; new diffLQ = %f; min itQ = %d \n",min_max_diffLQ,fabs(max_diffLQ) ,itQ);
                min_itQ = itQ;
                min_max_diffLQ = fabs(max_diffLQ);
                min_sum_diffLQ = sum_diffLQ;
            } else if ((fabs(max_diffLQ) == min_max_diffLQ) & (sum_diffLQ < min_sum_diffLQ)){
                //Rprintf("Cas 2 , diffLQ = %f; old sum_diffLQ = %f; new sum_diffLQ = %f; min itQ = %d \n",fabs(max_diffLQ),min_sum_diffLQ,sum_diffLQ ,itQ);
                min_itQ = itQ;
                min_max_diffLQ = fabs(max_diffLQ);
                min_sum_diffLQ = sum_diffLQ;
            }
//fichier << "min_max_diffLQ = " << min_max_diffLQ << endl;

         //rerun the best iteration which will be saved at index itmax-1 (last iteration)
        if (((GoOn == false) & (itQ<itmax-1)) | (itQ == itmax-2)){
            GoOn = true;
            itQ = itmax-2;

            for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                REAL(pQuota)[ind_t] = REAL(PQuot_temp)[min_itQ + itmax*ind_t]; // retrieve quota price for the best iteration min_itQ

                UNPROTECT(3);
            }
        }

        // for the last iteration, save the final price
        if (itQ == itmax-1){
            for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                REAL(pQuota)[ind_t] = REAL(PQuot_temp)[itQ + itmax*ind_t];

                UNPROTECT(3);
            }
        }

        itQ = itQ+1;


    }

    //---------------------------------------------
    // 2 - Final adjustment quota prices to account for current stock status and market
    //---------------------------------------------
    // Run Market and EcoDCF modules
    double max_ratio_GOS_QuotaExp;
    SEXP GOS_f,QuotaExp_f;
    double *r_GOS_f, *r_QuotaExp_f;

    Marche(listTemp, ind_t);

    GoOn=TRUE;
    itQ=0;

    while(GoOn && (itQ<10)){
//            fichier << "itQ = " << itQ << endl;

        EcoDCF(listTemp, ind_t, EcoIndCopy[4], drCopy);

        PROTECT(GOS_f = VECTOR_ELT(out_EcoDCF, 35));
        r_GOS_f = REAL(GOS_f);
        PROTECT(QuotaExp_f = VECTOR_ELT(out_EcoDCF, 59));
        r_QuotaExp_f = REAL(QuotaExp_f);

        max_ratio_GOS_QuotaExp = 1.0;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                if ((r_GoFish_f[ind_f+nbF*ind_t]==1) && (r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] > max_ratio_GOS_QuotaExp)){
                    max_ratio_GOS_QuotaExp = r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]  ;
                }
//                fichier << "GoFish = " << r_GoFish_f[ind_f+nbF*ind_t] <<
//                "; GOS = " <<  r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]<<
//                "; QuotaExp = " << r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] <<
//                "; ratio = " << r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]<<
//                "; max ratio = " << max_ratio_GOS_QuotaExp << endl;
        }

        if (max_ratio_GOS_QuotaExp == 1.0){ //OK
                GoOn=FALSE;
        } else{// Adjust Quota price by the ratio
                for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                        PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
//                        fichier << "Pquot before = " << REAL(pQuota)[ind_t] << endl;
                        REAL(pQuota)[ind_t] = REAL(pQuota)[ind_t]/(1.05 * max_ratio_GOS_QuotaExp);
//                        fichier << " PQuota after = " << REAL(pQuota)[ind_t] << endl;

                        UNPROTECT(2);
                    }
        }
//        fichier << "GoOn = " << GoOn << endl;
        itQ ++ ;
        UNPROTECT(2);
    }



    //---------------------------------------------
    // 3 - Trade quotas
    //---------------------------------------------

    // Rank fleet in order of profitability
    SEXP rank_Prof_f, ExpProfUE_f_copy;
    PROTECT(rank_Prof_f = NEW_INTEGER(nbF));
    int* r_rank_Prof_f = INTEGER(rank_Prof_f);
    ExpProfUE_f_copy = duplicate(ExpProfUE_f);
    double* r_ExpProfUE_f_copy = REAL(ExpProfUE_f_copy);

    double maxTemp = r_ExpProfUE_f_copy[0];
    int ind_maxTemp ;

    for (int ind_rank = 0 ; ind_rank < nbF ; ind_rank++){
        toinit=true;
        ind_maxTemp = 0;
        while ((ind_maxTemp <nbF) & toinit){
            if (!ISNA(r_ExpProfUE_f_copy[ind_maxTemp])){
                toinit=false;
                maxTemp = r_ExpProfUE_f_copy[ind_maxTemp];
            } else{ind_maxTemp++;}
        }

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            if ((!ISNA(r_ExpProfUE_f_copy[ind_f])) & (r_ExpProfUE_f_copy[ind_f] >= maxTemp)){
                maxTemp = r_ExpProfUE_f_copy[ind_f] ;
                ind_maxTemp = ind_f;
            }
        }
        r_ExpProfUE_f_copy[ind_maxTemp] = NA_REAL;
        r_rank_Prof_f[ind_rank] = ind_maxTemp; //index of the fleet of rank ind_rank
    }
    //PrintValue(ExpProfUE_f);
    //PrintValue(rank_Prof_f);

//    fichier << "QM6"  << endl;


    if (false){ //Old version

    SEXP demandQ_f_e, offerQ_f_e, L_ef,TACbyF_e,TACbyF_ex, QuotaTrade_fe ;
    PROTECT(demandQ_f_e = NEW_NUMERIC(nbF));
    PROTECT(offerQ_f_e = NEW_NUMERIC(nbF));
    double *r_L_ef;
    double *r_demandQ_f_e = REAL(demandQ_f_e);
    double *r_offerQ_f_e = REAL(offerQ_f_e);
    double *r_TACbyF_e;
    double *r_TACbyF_ex;
    double *r_QuotaTrade_fe;
    int ind_buyer, ind_seller;

    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
        //fichier << "nam_eQuota_dyn: " << CHAR(nam_eQuota_dyn) << endl;
        if (!isNull(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(L_ef = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),dimCstF));
                } else{ // espece stat
                                    PROTECT(L_ef = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),dimCstF));
                                }
        r_L_ef = REAL(L_ef);
        PROTECT(TACbyF_e = getListElement(TACbyF, CHAR(nam_eQuota_dyn)));
        r_TACbyF_e = REAL(TACbyF_e);
        PROTECT(TACbyF_ex = duplicate(TACbyF_e));
        r_TACbyF_ex = REAL(TACbyF_ex);

        //fichier << "QM7.1"  << endl;
//        if (ind_t==4){
//                Rprintf("sp = %s Expected catches: \n",CHAR(nam_eQuota_dyn));
//                PrintValue(L_ef);
//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Expected catches:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_L_ef[ind_f+nbF*ind_t]  << endl;
//                }
//        }

        // Calculate net demand for quota by fleet = expected catches - holdings
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            r_demandQ_f_e[ind_f] = r_L_ef[ind_f+nbF*ind_t] - r_TACbyF_e[ind_f+nbF*ind_t] ;
            r_offerQ_f_e[ind_f] = - r_demandQ_f_e[ind_f] ;
        }
        //fichier << "QM7.2"  << endl;
        //Rprintf("TAC before trade\n");
        //PrintValue(TACbyF_e);
        //Rprintf("Expected catches\n");
        //PrintValue(L_ef);
        //Rprintf("Demand before trade\n");
        //PrintValue(demandQ_f_e);

        //Trade quotas
        int rank_buyer = 0;
        ind_buyer = r_rank_Prof_f[rank_buyer];
        int rank_seller = nbF-1;
        ind_seller = r_rank_Prof_f[rank_seller];
        double trade;



                while ((rank_buyer < nbF) &&  (rank_seller >=0)){
                        //fichier << "rank buyer = " << rank_buyer  << "; ind buyer = " << ind_buyer << endl;
                    while ((rank_seller >=0) && (r_demandQ_f_e[ind_buyer] > 0)){
                        //fichier << "rank seller = " << rank_seller << "; ind seller = " << ind_seller << endl;

                        if(r_offerQ_f_e[ind_seller] > 0){
                            trade = fmin2(r_demandQ_f_e[ind_buyer],r_offerQ_f_e[ind_seller]);
                            //fichier << "demand = " << r_demandQ_f_e[ind_buyer] << ", offer = " << r_offerQ_f_e[ind_seller] << ", trade = " << trade << endl;
                            r_demandQ_f_e[ind_buyer] = r_demandQ_f_e[ind_buyer] - trade;
                            r_offerQ_f_e[ind_seller] = r_offerQ_f_e[ind_seller] - trade;

                            r_TACbyF_e[ind_buyer+nbF*ind_t] = r_TACbyF_e[ind_buyer+nbF*ind_t] + trade;
                            r_TACbyF_e[ind_seller+nbF*ind_t] = r_TACbyF_e[ind_seller+nbF*ind_t] - trade;
                        } else {
                            rank_seller --;
                            ind_seller = r_rank_Prof_f[rank_seller];
                        }
                    }
                    rank_buyer ++;
                    ind_buyer = r_rank_Prof_f[rank_buyer];

                }



//        fichier << "QM7.4"  << endl;
        //Rprintf("Demand after trade\n");
        //PrintValue(demandQ_f_e);
        //Rprintf("Offer after trade\n");
        //PrintValue(offerQ_f_e);
        //Rprintf("TAC after trade\n");
        //PrintValue(TACbyF_e);
//
//        if (ind_t==4){
//                Rprintf("sp = %s , Quota after trade: \n",CHAR(nam_eQuota_dyn));
//                PrintValue(TACbyF_e);
//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Quota after trade:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_TACbyF_e[ind_f+nbF*ind_t]  << endl;
//                }
//        }

        // REcord quota trades
         PROTECT(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota_dyn)));
         r_QuotaTrade_fe= REAL(QuotaTrade_fe);

         for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                r_QuotaTrade_fe[ind_f+ nbF*ind_t] = r_TACbyF_e[ind_f+ nbF*ind_t] - r_TACbyF_ex[ind_f+ nbF*ind_t];
         }
         //PrintValue(QuotaTrade_fe);

        UNPROTECT(5);
    }
    UNPROTECT(2);
    } else{ // New

    SEXP demandQ_f_e, offerQ_f_e, L_ef,TACbyF_e, QuotaTrade_fe,Qholdings_e ;
    PROTECT(demandQ_f_e = NEW_NUMERIC(nbF+1));
    PROTECT(offerQ_f_e = NEW_NUMERIC(nbF+1));
    double *r_L_ef;
    double *r_demandQ_f_e = REAL(demandQ_f_e);
    double *r_offerQ_f_e = REAL(offerQ_f_e);
    double *r_TACbyF_e, *r_Qholdings;
    double *r_QuotaTrade_fe;
    int ind_buyer, ind_seller;

    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
        if (!isNull(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(L_ef = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),dimCstF));
                } else{ // espece stat
                                    PROTECT(L_ef = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),dimCstF));
                                }
        r_L_ef = REAL(L_ef);
        PROTECT(TACbyF_e = getListElement(TACbyF, CHAR(nam_eQuota_dyn)));
        r_TACbyF_e = REAL(TACbyF_e);
        PROTECT(Qholdings_e = getListElement(Qholdings, CHAR(nam_eQuota_dyn)));
        r_Qholdings = REAL(Qholdings_e);

//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Expected catches:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_L_ef[ind_f+nbF*ind_t]  << endl;
//                }


        // Calculate net demand for quota by fleet = expected catches - holdings
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            r_demandQ_f_e[ind_f] = r_L_ef[ind_f+nbF*ind_t] - r_Qholdings[ind_f+(nbF+1)*ind_t] ;
            r_offerQ_f_e[ind_f] = - r_demandQ_f_e[ind_f] ;
        }

        // And for external quota holders
        r_demandQ_f_e[nbF] = 0.0;
        r_offerQ_f_e[nbF] = r_Qholdings[nbF+(nbF+1)*ind_t];

        //Trade quotas
        int rank_buyer = 0;
        ind_buyer = r_rank_Prof_f[rank_buyer];
        int rank_seller = nbF; // start with external holders
        ind_seller = nbF;
        double trade;



                while ((rank_buyer < nbF) &&  (rank_seller >=0)){
                        //fichier << "rank buyer = " << rank_buyer  << "; ind buyer = " << ind_buyer << endl;
                    while ((rank_seller >=0) && (r_demandQ_f_e[ind_buyer] > 0)){
                        //fichier << "rank seller = " << rank_seller << "; ind seller = " << ind_seller << endl;

                        if(r_offerQ_f_e[ind_seller] > 0){
                            trade = fmin2(r_demandQ_f_e[ind_buyer],r_offerQ_f_e[ind_seller]);
                            //fichier << "demand = " << r_demandQ_f_e[ind_buyer] << ", offer = " << r_offerQ_f_e[ind_seller] << ", trade = " << trade << endl;
                            r_demandQ_f_e[ind_buyer] = r_demandQ_f_e[ind_buyer] - trade;
                            r_offerQ_f_e[ind_seller] = r_offerQ_f_e[ind_seller] - trade;

                            r_TACbyF_e[ind_buyer+nbF*ind_t] = r_TACbyF_e[ind_buyer+nbF*ind_t] + trade;

                            if(ind_seller<nbF)// Only for vessels, not external holders
                            r_TACbyF_e[ind_seller+nbF*ind_t] = r_TACbyF_e[ind_seller+nbF*ind_t] - trade;
                        } else {
                            rank_seller --;
                            ind_seller = r_rank_Prof_f[rank_seller];
                        }
                    }
                    rank_buyer ++;
                    ind_buyer = r_rank_Prof_f[rank_buyer];

                }

//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Quota after trade:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_TACbyF_e[ind_f+nbF*ind_t]  << endl;
//                }


        // REcord quota trades
         PROTECT(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota_dyn)));
         r_QuotaTrade_fe= REAL(QuotaTrade_fe);

         for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                r_QuotaTrade_fe[ind_f+ nbF*ind_t] = r_TACbyF_e[ind_f+ nbF*ind_t] - r_Qholdings[ind_f+ (nbF+1)*ind_t];
         }


//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Traded quota:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_QuotaTrade_fe[ind_f+nbF*ind_t]  << endl;
//                }


        UNPROTECT(5);
    }
    UNPROTECT(2);
    }
//    fichier << "QM7.5"  << endl;



    UNPROTECT(31);
//  fichier.close();
}
}

