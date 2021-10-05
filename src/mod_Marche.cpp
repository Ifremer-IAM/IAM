#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
//#include <Rcpp.h>

#include "BioEcoPar.h" // Class is defined in this file.
// #include "Modules.h" // Contient tout les modules

//using namespace Rcpp;
using namespace std;

//
//------------------------------------------
// Module 'March�'
//------------------------------------------
//
extern "C" {

void BioEcoPar::Marche(SEXP list, int ind_t)
{

//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.Marche.txt", ios::out | ios::trunc);

    SEXP    elmt, intC, v_P_fmce, v_icat, v_L_efmit, dimCst_P_fmce, dimCst_L_efmit, dimCst_L_efmct, Dim_L_efmct, //Dim_P_fmce,
            ans_L_efmct = R_NilValue, ans_P_fmce = R_NilValue, dimnames_Lc = R_NilValue, dimnames_P = R_NilValue,
            rnames_Esp, cFACTc, cFACTi, cFACTp, cFACTpini, cFACTpStat, cFACTpStatini,
            ans_DD_efmc = R_NilValue, ans_LD_efmc = R_NilValue, v_DD_efmit, v_LD_efmit;
    SEXP    v_P_eStat, dimCst_P_eStat, dimCst_P_eStatR, ans_P_eStat = R_NilValue, dimnames_Pstat = R_NilValue;

    int *dim_P_fmce, *dim_P_fmcet, *dim_L_efmit, *dim_icat, *dim_L_efmct, *dim_P_eStat, *dim_P_eStat_t, *dim_P_eStat_tR, *dimLc;//, *dimP;

    int nbI, nbC;

    double *rans_L_efmct, *r_L_efmit, *r_P_fmce, *r_icat, *r_Pstat, *r_P_fmceIni, *r_PstatIni,
            *rans_DD_efmc, *rans_LD_efmc, *r_DD_efmit, *r_LD_efmit;
    double mult_p; //multiplicateur de prix (1 par defaut sinon fonction de relation prix quantite
    SEXP MarketList= R_NilValue, v_ep= R_NilValue, v_beta_pp= R_NilValue;
    double *r_beta_pp = &NA_REAL;
    int *r_ep = &NA_INTEGER;
    int ind_p;
    SEXP ans_MultPrice;


////Rprintf("CCC1");

//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 5));////Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));////Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 7));////Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 61));////Rprintf("Mort20.5\n");



    if (ind_t==0) {

        PROTECT(rnames_Esp = getAttrib(out_L_efmit,R_NamesSymbol));
        setAttrib(out_L_efmct, R_NamesSymbol, rnames_Esp);
        setAttrib(out_DD_efmc, R_NamesSymbol, rnames_Esp);
        setAttrib(out_LD_efmc, R_NamesSymbol, rnames_Esp);
        setAttrib(out_P_t, R_NamesSymbol, rnames_Esp);
        if (nbEstat>0) setAttrib(out_Pstat, R_NamesSymbol, sppListStat);

    }



if(pUpdate){


    PROTECT(MarketList = getListElement(list,"Market"));
    PROTECT(v_ep = getListElement(MarketList,"ep")); //PrintValue(v_ep);
    PROTECT(v_beta_pp = getListElement(MarketList,"beta_pp")); //PrintValue(v_beta_pp);

    r_beta_pp = REAL(v_beta_pp);
    r_ep = INTEGER(v_ep);

    //---------
    // Calcul debarquements produits
    //---------

    SEXP rnames_p = R_NilValue, ans_L_pt = R_NilValue, dimCstL_pt, dimNamL_pt;
    int *int_dimCstL_pt;
    double L_p;
    double *rans_L_pt;

    if (ind_t==0) {
        PROTECT(rnames_p = allocVector(STRSXP, nbP));
        setAttrib(out_L_pt, R_NamesSymbol, rnames_p);
    }

    for (int p = 0; p < nbP; p++) {

            if (ind_t==0) {
            PROTECT(ans_L_pt = NEW_NUMERIC(nbT));

            PROTECT(dimCstL_pt = allocVector(INTSXP, 1));
            int_dimCstL_pt = INTEGER(dimCstL_pt);  int_dimCstL_pt[0] = nbT;
            setAttrib(ans_L_pt, R_DimSymbol, dimCstL_pt);

            PROTECT( dimNamL_pt = allocVector(VECSXP,1));
            SET_VECTOR_ELT(dimNamL_pt, 0, times);
            setAttrib(ans_L_pt, R_DimNamesSymbol, dimNamL_pt);

            rans_L_pt = REAL(ans_L_pt);
            } else {
                rans_L_pt = REAL(VECTOR_ELT(out_L_pt, p));
                    }

            //Colonne de la matrice ep : 1 si espece correspond au meme produit, 0 sinon
            L_p = 0.0;
            for (int k = 0; k < nbEall; k++) {
                    L_p = L_p + r_ep[p*nbEall + k] * REAL(VECTOR_ELT(out_L_et, k))[ind_t];
                    //fichier << "p = " << CHAR(STRING_ELT(pList,p)) << ", k = " << CHAR(STRING_ELT(sppListAll,k)) << ", r_ep = " << r_ep[p*nbEall + k] << ", L_p = " << L_p << endl;
            }

            rans_L_pt [ind_t] = L_p; //fichier << "p = " << CHAR(STRING_ELT(pList,p)) << ", rans_L_pt = " << rans_L_pt [ind_t] << endl;

            if (ind_t==0) {
                SET_VECTOR_ELT(out_L_pt, p, ans_L_pt);
                SET_STRING_ELT(rnames_p, p, STRING_ELT(pList,p));
                UNPROTECT(3);
                }

        }
    }

if (nbE>0) {

    for (int e = 0 ; e < nbE ; e++) {
//Rprintf("M2\n");fichier << "M2" << endl;
//fichier << "e = " << CHAR(STRING_ELT(sppList,e)) << endl;

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        intC = getListElement(elmt, "modC");
        nbC = length(intC);

        PROTECT(v_P_fmce = getListElement(elmt, "P_fmce"));
        PROTECT(v_icat = getListElement(elmt, "icat"));
        PROTECT(v_L_efmit = getListElement(out_L_efmit, CHAR(STRING_ELT(sppList,e))));
        PROTECT(v_DD_efmit = getListElement(out_DD_efmi, CHAR(STRING_ELT(sppList,e))));
        PROTECT(v_LD_efmit = getListElement(out_LD_efmi, CHAR(STRING_ELT(sppList,e))));

        PROTECT(dimCst_P_fmce = getAttrib(v_P_fmce, install("DimCst")));
        PROTECT(dimCst_L_efmit = getAttrib(v_L_efmit, install("DimCst")));
//Rprintf("M3\n");fichier << "M3" << endl;

        //tests sur les dimensions :
        dim_P_fmce = INTEGER(dimCst_P_fmce);////Rprintf("AAA1");
        if (((dim_P_fmce[0]!=0) & (dim_P_fmce[0]!=nbF)) | ((dim_P_fmce[1]!=0) & (dim_P_fmce[1]!=nbMe)) |
            ((dim_P_fmce[2]!=0) & (dim_P_fmce[2]!=nbC)) | ((dim_P_fmce[3]!=0) & (dim_P_fmce[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in P_fmce element. Check .ini biological parameters files !!\n");
        }

        dim_L_efmit = INTEGER(dimCst_L_efmit);////Rprintf("AAA2");
        if (((dim_L_efmit[0]!=0) & (dim_L_efmit[0]!=nbF)) | ((dim_L_efmit[1]!=0) & (dim_L_efmit[1]!=nbM)) |
            ((dim_L_efmit[2]!=0) & (dim_L_efmit[2]!=nbI)) | ((dim_L_efmit[3]!=0) & (dim_L_efmit[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in L_efmit element. Check .ini biological parameters files !!\n");
        }

        dim_icat = INTEGER(getAttrib(v_icat, R_DimSymbol));////Rprintf("AAA3");
        if ((dim_icat[0]!=nbI) & (dim_icat[1]!=nbC))
        {
            error("Non_homogeneous dimensions in icat element. Check .ini biological parameters files !!\n");
        }
//Rprintf("M4\n");fichier << "M4" << endl;

//        ---------
//         Ajustement prix avec relation prix/quantite
//        ---------
        mult_p = 1.0;
        if(pUpdate){
//Rprintf("M5\n");fichier << "M5" << endl;
                int ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e))); //indice dans matrice ep
                //fichier << "e = " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

                double *r_L_pt;

                //Trouver l'indice du produit correspondant
                for (ind_p = 0 ; ind_p < nbP ; ind_p++) {if (r_ep[ind_p*nbEall + ind_e] == 1) break;}

                //fichier << "p = " << CHAR(STRING_ELT(pList,ind_p)) << endl;

                for ( int k = 0; k < nbP; k++) {

                    r_L_pt = REAL(VECTOR_ELT(out_L_pt, k));
                    if (r_L_pt[0]!=0) mult_p = mult_p * pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) ;
                    //fichier << "k = " << CHAR(STRING_ELT(pList,k)) << ", ratio L = " << r_L_pt[ind_t] / r_L_pt[0] << ", beta = " << r_beta_pp[k*nbP + ind_p] << ", pow = "<<  pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) << ", mult = " << mult_p << endl;

                }

                PROTECT(ans_MultPrice = VECTOR_ELT(multPrice, ind_e));
                REAL(ans_MultPrice)[ind_t] = mult_p;
                UNPROTECT(1);

        }




        //---------
        // calcul de L_efmct
        //---------

        PROTECT(dimCst_L_efmct = allocVector(INTSXP, 4));
        dim_L_efmct = INTEGER(dimCst_L_efmct);
        dim_L_efmct[0] = dim_L_efmit[0] ; dim_L_efmct[1] = dim_L_efmit[1] ; dim_L_efmct[2] = nbC; dim_L_efmct[3] = dim_L_efmit[3];

        int count = 0, prod = 1, count2 = 0, count3 = 0;
        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                count++;
                prod = prod * dim_L_efmct[k];
            }

        }

        PROTECT(Dim_L_efmct = allocVector(INTSXP, count));
        dimLc = INTEGER(Dim_L_efmct);

//Rprintf("M5\n");fichier << "M5" << endl;
        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                dimLc[count2] = dim_L_efmct[k];
                count2++;
            }

        }

//--------------------------

        PROTECT(dimCst_P_fmce = allocVector(INTSXP, 4));
        dim_P_fmcet = INTEGER(dimCst_P_fmce);
        dim_P_fmcet[0] = nbF ; dim_P_fmcet[1] = nbMe ; dim_P_fmcet[2] = nbC; dim_P_fmcet[3] = nbT;



if (ind_t==0){
//Rprintf("M6\n");fichier << "M6" << endl;
        //on cr�e le tableau r�sultat pour l'esp�ce en question
        PROTECT(ans_L_efmct = NEW_NUMERIC(prod));
        setAttrib(ans_L_efmct, R_DimSymbol, Dim_L_efmct);

        PROTECT(ans_DD_efmc = NEW_NUMERIC(prod));
        setAttrib(ans_DD_efmc, R_DimSymbol, Dim_L_efmct);

        PROTECT(ans_LD_efmc = NEW_NUMERIC(prod));
        setAttrib(ans_LD_efmc, R_DimSymbol, Dim_L_efmct);

        PROTECT(dimnames_Lc = allocVector(VECSXP,count));
        if (dim_L_efmct[0]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, fleetList) ; count3++;}
        if (dim_L_efmct[1]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, metierList) ; count3++;}
        if (dim_L_efmct[2]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, intC) ; count3++;}
        if (dim_L_efmct[3]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, times) ; count3++;}

        rans_L_efmct = REAL(ans_L_efmct);
        rans_DD_efmc = REAL(ans_DD_efmc);
        rans_LD_efmc = REAL(ans_LD_efmc);

        PROTECT(ans_P_fmce = NEW_NUMERIC(nbF*nbMe*nbC*nbT));
        setAttrib(ans_P_fmce, R_DimSymbol, dimCst_P_fmce);

        PROTECT(dimnames_P = allocVector(VECSXP,4));
        SET_VECTOR_ELT(dimnames_P, 0, fleetList);
        SET_VECTOR_ELT(dimnames_P, 1, metierListEco);
        SET_VECTOR_ELT(dimnames_P, 2, intC);
        SET_VECTOR_ELT(dimnames_P, 3, times);

        r_P_fmce = REAL(ans_P_fmce);

//Rprintf("M7\n");fichier << "M7" << endl;
} else {

        rans_L_efmct = REAL(VECTOR_ELT(out_L_efmct, e));
        rans_DD_efmc = REAL(VECTOR_ELT(out_DD_efmc, e));
        rans_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
        r_P_fmce = REAL(VECTOR_ELT(out_P_t, e));

//Rprintf("M8\n");fichier << "M8" << endl;
}

        r_L_efmit = REAL(v_L_efmit);
        r_DD_efmit = REAL(v_DD_efmit);
        r_LD_efmit = REAL(v_LD_efmit);
        r_icat = REAL(v_icat);

//Rprintf("M9\n");fichier << "M9" << endl;
        //facteurs des indices
        PROTECT(cFACTc = iDim(dim_L_efmct));
        PROTECT(cFACTi = iDim(dim_L_efmit));
        PROTECT(cFACTp = iDim(dim_P_fmcet));
        PROTECT(cFACTpini = iDim(dim_P_fmce));

        int *fact_Cc = INTEGER(cFACTc);
        int *fact_Ci = INTEGER(cFACTi);
        int *fact_P = INTEGER(cFACTp);
        int *fact_Pini = INTEGER(cFACTpini);

        r_P_fmceIni = REAL(v_P_fmce);

        //�quation n�1 : conversion �ge/catg�gorie

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
        for (int ind_c = 0 ; ind_c < nbC ; ind_c++) {

            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                if (ind_i ==0) {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

           rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_LD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

           rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_DD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                } else {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

            rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_LD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

            rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_DD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                }
            }
        // valable car nbM = nbMe (M=Me)
            //r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] =
            //    r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + ind_t*fact_Pini[3]];

            r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] =
                r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + 0*fact_Pini[3]] * mult_p; //Florence 05/2019

                //if ((ind_f==0) & (ind_m==0)) fichier << "r_P_fmceIni = " << r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + 0*fact_Pini[3]] << ", mult_p = " << mult_p << ", r_P_fmce = " << r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] << endl;

        }


//Rprintf("M10\n");fichier << "M10" << endl;

if (ind_t==0) {

        setAttrib(ans_L_efmct, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_L_efmct, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_L_efmct, e, ans_L_efmct);

        setAttrib(ans_DD_efmc, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_DD_efmc, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_DD_efmc, e, ans_DD_efmc);

        setAttrib(ans_LD_efmc, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_LD_efmc, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_LD_efmc, e, ans_LD_efmc);

      //  SET_VECTOR_ELT(out_L_efmct2, e, ans_L_efmct);

        setAttrib(ans_P_fmce, R_DimNamesSymbol, dimnames_P);
        setAttrib(ans_P_fmce, install("DimCst"), dimCst_P_fmce);

        SET_VECTOR_ELT(out_P_t, e, ans_P_fmce);

}

//Rprintf("M11\n");fichier << "M11" << endl;
if (ind_t==0) UNPROTECT(6);
UNPROTECT(15);

}
}


if (nbEstat>0) {
    for (int e = 0 ; e < nbEstat ; e++) {
//Rprintf("M2\n");fichier << "M2" << endl;
        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e))));

//        ---------
//         Ajustement prix avec relation prix/quantite
//        ---------
        mult_p = 1.0;
        if(pUpdate){
//Rprintf("M5\n");fichier << "M5" << endl;
                int ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e))); //indice dans matrice ep
                //fichier << "e = " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

                double *r_L_pt;

                //Trouver l'indice du produit correspondant
                for (ind_p = 0 ; ind_p < nbP ; ind_p++) {if (r_ep[ind_p*nbEall + ind_e] == 1) break;}

                //fichier << "p = " << CHAR(STRING_ELT(pList,ind_p)) << endl;

                for ( int k = 0; k < nbP; k++) {

                    r_L_pt = REAL(VECTOR_ELT(out_L_pt, k));
                    if (r_L_pt[0]!=0) mult_p = mult_p * pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) ;
                    //fichier << "k = " << CHAR(STRING_ELT(pList,k)) << ", ratio L = " << r_L_pt[ind_t] / r_L_pt[0] << ", beta = " << r_beta_pp[k*nbP + ind_p] << ", pow = "<<  pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) << ", mult = " << mult_p << endl;

                }
            PROTECT(ans_MultPrice = VECTOR_ELT(multPrice, ind_e));
            REAL(ans_MultPrice)[ind_t] = mult_p;
            UNPROTECT(1);


        }


        PROTECT(v_P_eStat = getListElement(elmt, "P_fme"));
        PROTECT(dimCst_P_eStat = getAttrib(v_P_eStat, install("DimCst")));

        //tests sur les dimensions :
        dim_P_eStat = INTEGER(dimCst_P_eStat);////Rprintf("AAA1");
        if (((dim_P_eStat[0]!=0) & (dim_P_eStat[0]!=nbF)) | ((dim_P_eStat[1]!=0) & (dim_P_eStat[1]!=nbMe)) |
            (dim_P_eStat[2]!=0) | ((dim_P_eStat[3]!=0) & (dim_P_eStat[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in P_fme element. Check .ini biological parameters files !!\n");
        }



        PROTECT(dimCst_P_eStat = allocVector(INTSXP, 4));
        dim_P_eStat_t = INTEGER(dimCst_P_eStat);
        dim_P_eStat_t[0] = nbF ; dim_P_eStat_t[1] = nbMe ; dim_P_eStat_t[2] = 0; dim_P_eStat_t[3] = nbT;
        PROTECT(dimCst_P_eStatR = allocVector(INTSXP, 3));
        dim_P_eStat_tR = INTEGER(dimCst_P_eStatR);
        dim_P_eStat_tR[0] = nbF ; dim_P_eStat_tR[1] = nbMe ; dim_P_eStat_tR[2] = nbT;


if (ind_t==0){
//Rprintf("M6\n");fichier << "M6" << endl;

        PROTECT(ans_P_eStat = NEW_NUMERIC(nbF*nbMe*nbT));
        setAttrib(ans_P_eStat, R_DimSymbol, dimCst_P_eStatR);

        PROTECT(dimnames_Pstat = allocVector(VECSXP,3));
        SET_VECTOR_ELT(dimnames_Pstat, 0, fleetList);
        SET_VECTOR_ELT(dimnames_Pstat, 1, metierListEco);
        SET_VECTOR_ELT(dimnames_Pstat, 2, times);

        r_Pstat = REAL(ans_P_eStat);

//Rprintf("M7\n");fichier << "M7" << endl;
} else {

        r_Pstat = REAL(VECTOR_ELT(out_Pstat, e));

//Rprintf("M8\n");fichier << "M8" << endl;
}
//Rprintf("M9\n");fichier << "M9" << endl;
        //facteurs des indices
        PROTECT(cFACTpStat = iDim(dim_P_eStat_t));
        PROTECT(cFACTpStatini = iDim(dim_P_eStat));

        int *fact_Pstat = INTEGER(cFACTpStat);
        int *fact_Pstatini = INTEGER(cFACTpStatini);

        r_PstatIni = REAL(v_P_eStat);

        //�quation n�1 : conversion �ge/catg�gorie

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

        // valable car nbM = nbMe (M=Me)
           // r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] =
             //   r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + ind_t*fact_Pstatini[3]];

            r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] =
                r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + 0*fact_Pstatini[3]] * mult_p;// Florence 05/2019


                //if ((ind_f==0) & (ind_m==0)) fichier << "r_PstatIni = " << r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + 0*fact_Pstatini[3]] << ", mult_p = " << mult_p << ", r_Pstat = " << r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] << endl;


        }

//Rprintf("M10\n");fichier << "M10" << endl;

if (ind_t==0) {

        setAttrib(ans_P_eStat, R_DimNamesSymbol, dimnames_Pstat);
        setAttrib(ans_P_eStat, install("DimCst"), dimCst_P_eStat);

        SET_VECTOR_ELT(out_Pstat, e, ans_P_eStat);

}

//Rprintf("M11\n");fichier << "M11" << endl;
if (ind_t==0) UNPROTECT(2);
UNPROTECT(7);

}

}

if (pUpdate){
        if(ind_t == 0) UNPROTECT(1);
    UNPROTECT(3);
}
if (ind_t==0) UNPROTECT(1);


//fichier.close();

}}