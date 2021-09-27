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

//------------------------------------------
// Module 'Dynamique de population' : 'out_F_fmi' = output de la fonction 'Mortalite' MAJ 27/09/2010 ajout de l'output SSB_et
//------------------------------------------

extern "C" {

void BioEcoPar::DynamicPop(SEXP list, int ind_t, SEXP EVAR, bool Reality) //Reality : si True, on appelle DynamicPop pour une vraie projection, sinon uniquement utilis� pour estimer TAC comme dans WG
                                                                           // ~ si True, arbitrage RecParamList VS MeanRec_Ftarg en faveur du premier (ie sinon, en faveur du deuxi�me)
{

////Rprintf("G0");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 5));//Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 6));//Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 7));//Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 61));//Rprintf("Mort20.5\n");

//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.DynamicPop.txt", ios::out | ios::trunc);

if (dUpdate) {
//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.DynamicPop.txt", ios::out | ios::trunc);
//fichier << "dUpdate = " << dUpdate << endl;
//Rprintf("dUpdate = %f \n" ,dUpdate);

    SEXP    elmt, dFACT1, dFACT2, dFACT3, dFACT4, dFACT5, dFACT6, dFACT7, dFACT8, dFACT9, dFACT10,
            dimCst1, dimCst2, dimCst3, dimCst4, Dim1, Dim2, Dim3, Dim4,
            dimCst_Fr_efmit=R_NilValue, dimCst_M_ei=R_NilValue, dimCst_w_ei=R_NilValue,
            dimCst_N_ei0=R_NilValue, dimCst_N_e0t=R_NilValue, dimCst_mat_ei=R_NilValue,
            intAge, v_Fr_efmit=R_NilValue, v_F_efmit=R_NilValue, v_M_ei=R_NilValue, v_w_ei=R_NilValue, v_N_ei0=R_NilValue, v_N_e0t=R_NilValue,
            v_mat_ei=R_NilValue, v_Fbar=R_NilValue,
            v_r=R_NilValue, v_K=R_NilValue, v_n=R_NilValue, v_B=R_NilValue,
            v_Fr_efmit_G1=R_NilValue, v_Fr_efmit_G2=R_NilValue, v_F_efmit_G1=R_NilValue,  v_F_efmit_G2=R_NilValue,
            v_M_ei_G1=R_NilValue, v_M_ei_G2=R_NilValue,v_w_ei_G1=R_NilValue, v_w_ei_G2=R_NilValue,  v_Fbar_G1=R_NilValue, v_Fbar_G2=R_NilValue,
            v_N_ei0_G1=R_NilValue, v_N_ei0_G2=R_NilValue, v_N_e0t_G1=R_NilValue,v_N_e0t_G2=R_NilValue, v_mat_ei_G1=R_NilValue, v_mat_ei_G2=R_NilValue;

    SEXP v_Fr_efmit_S1M1 = R_NilValue, v_Fr_efmit_S1M2 = R_NilValue, v_Fr_efmit_S1M3 = R_NilValue, v_Fr_efmit_S1M4 = R_NilValue,
         v_Fr_efmit_S2M1 = R_NilValue, v_Fr_efmit_S2M2 = R_NilValue, v_Fr_efmit_S2M3 = R_NilValue, v_Fr_efmit_S2M4 = R_NilValue,
         v_Fr_efmit_S3M1 = R_NilValue, v_Fr_efmit_S3M2 = R_NilValue, v_Fr_efmit_S3M3 = R_NilValue, v_Fr_efmit_S3M4 = R_NilValue,
         v_Fr_efmit_S4M1 = R_NilValue, v_Fr_efmit_S4M2 = R_NilValue, v_Fr_efmit_S4M3 = R_NilValue, v_Fr_efmit_S4M4 = R_NilValue,
         v_F_efmit_S1M1 = R_NilValue, v_F_efmit_S1M2 = R_NilValue, v_F_efmit_S1M3 = R_NilValue, v_F_efmit_S1M4 = R_NilValue,
         v_F_efmit_S2M1 = R_NilValue, v_F_efmit_S2M2 = R_NilValue, v_F_efmit_S2M3 = R_NilValue, v_F_efmit_S2M4 = R_NilValue,
         v_F_efmit_S3M1 = R_NilValue, v_F_efmit_S3M2 = R_NilValue, v_F_efmit_S3M3 = R_NilValue, v_F_efmit_S3M4 = R_NilValue,
         v_F_efmit_S4M1 = R_NilValue, v_F_efmit_S4M2 = R_NilValue, v_F_efmit_S4M3 = R_NilValue, v_F_efmit_S4M4 = R_NilValue,
         v_N_ei0_S1M1 = R_NilValue, v_N_ei0_S1M2 = R_NilValue, v_N_ei0_S1M3 = R_NilValue, v_N_ei0_S1M4 = R_NilValue,
         v_N_e0t_S1M1 = R_NilValue, v_N_e0t_S2M2 = R_NilValue, v_N_e0t_S3M3 = R_NilValue, v_N_e0t_S4M4 = R_NilValue,
         v_iniNt0q_S1M1 = R_NilValue, v_iniNt0q_S1M2 = R_NilValue, v_iniNt0q_S1M3 = R_NilValue, v_iniNt0q_S1M4 = R_NilValue,
         v_iniNt0q_S2M1 = R_NilValue, v_iniNt0q_S2M2 = R_NilValue, v_iniNt0q_S2M3 = R_NilValue, v_iniNt0q_S2M4 = R_NilValue,
         v_iniNt0q_S3M1 = R_NilValue, v_iniNt0q_S3M2 = R_NilValue, v_iniNt0q_S3M3 = R_NilValue, v_iniNt0q_S3M4 = R_NilValue,
         v_iniNt0q_S4M1 = R_NilValue, v_iniNt0q_S4M2 = R_NilValue, v_iniNt0q_S4M3 = R_NilValue, v_iniNt0q_S4M4 = R_NilValue,
         v_matwt_M1 = R_NilValue, v_matwt_M2 = R_NilValue, v_matwt_M3 = R_NilValue, v_matwt_M4 = R_NilValue;

    SEXP dimnames1=R_NilValue, dimnames2=R_NilValue, dimnames3=R_NilValue, dimnames4=R_NilValue, rnames_Esp=R_NilValue;

    int *dim_Fr_efmit, *dim_M_ei, *dim_w_ei, *dim_N_ei0, *dim_N_e0t, *dim_mat_ei, *dimC1, *dimC2=0, *dimC3, *dimC4, *dim1, *dim2, *dim3, *dim4, *fact4_D=0;
    int nbI;

    double *rans_Z_eit=&NA_REAL, *rans_N_eit=&NA_REAL, *rans_B_et, *rans_Fbar_et=&NA_REAL, *rans_SSB_et, *r_Fr_efmit=&NA_REAL,
          *r_F_efmit=&NA_REAL, *r_Fbar=&NA_REAL, *r_M_ei=&NA_REAL, *r_w_ei=&NA_REAL, *r_N_ei0=&NA_REAL, *r_N_e0t=&NA_REAL, *r_mat_ei=&NA_REAL,
          *rans_Z_eit_G1=&NA_REAL,*rans_Z_eit_G2=&NA_REAL, *rans_N_eit_G1=&NA_REAL, *rans_N_eit_G2=&NA_REAL,
          *r_Fr_efmit_G1=&NA_REAL, *r_Fr_efmit_G2=&NA_REAL,
          *r_Fbar_G1=&NA_REAL, *r_Fbar_G2=&NA_REAL, *r_N_ei0_G1=&NA_REAL, *r_N_ei0_G2=&NA_REAL,
          *r_M_ei_G1=&NA_REAL, *r_M_ei_G2=&NA_REAL, *r_w_ei_G1=&NA_REAL, *r_w_ei_G2=&NA_REAL,
          *r_mat_ei_G1=&NA_REAL, *r_mat_ei_G2=&NA_REAL, *r_r=&NA_REAL, *r_K=&NA_REAL, *r_n=&NA_REAL, *r_B=&NA_REAL;

    double delta_r = 0.0, delta_K = 0.0;

    double *rans_Z_eit_S1M1=&NA_REAL, *rans_Z_eit_S1M2=&NA_REAL, *rans_Z_eit_S1M3=&NA_REAL, *rans_Z_eit_S1M4=&NA_REAL,
           *rans_Z_eit_S2M1=&NA_REAL, *rans_Z_eit_S2M2=&NA_REAL, *rans_Z_eit_S2M3=&NA_REAL, *rans_Z_eit_S2M4=&NA_REAL,
           *rans_Z_eit_S3M1=&NA_REAL, *rans_Z_eit_S3M2=&NA_REAL, *rans_Z_eit_S3M3=&NA_REAL, *rans_Z_eit_S3M4=&NA_REAL,
           *rans_Z_eit_S4M1=&NA_REAL, *rans_Z_eit_S4M2=&NA_REAL, *rans_Z_eit_S4M3=&NA_REAL, *rans_Z_eit_S4M4=&NA_REAL,
           *rans_N_eit_S1M1=&NA_REAL, *rans_N_eit_S1M2=&NA_REAL, *rans_N_eit_S1M3=&NA_REAL, *rans_N_eit_S1M4=&NA_REAL,
           *rans_N_eit_S2M1=&NA_REAL, *rans_N_eit_S2M2=&NA_REAL, *rans_N_eit_S2M3=&NA_REAL, *rans_N_eit_S2M4=&NA_REAL,
           *rans_N_eit_S3M1=&NA_REAL, *rans_N_eit_S3M2=&NA_REAL, *rans_N_eit_S3M3=&NA_REAL, *rans_N_eit_S3M4=&NA_REAL,
           *rans_N_eit_S4M1=&NA_REAL, *rans_N_eit_S4M2=&NA_REAL, *rans_N_eit_S4M3=&NA_REAL, *rans_N_eit_S4M4=&NA_REAL,
           *r_Fr_efmit_S1M1=&NA_REAL, *r_Fr_efmit_S1M2=&NA_REAL, *r_Fr_efmit_S1M3=&NA_REAL, *r_Fr_efmit_S1M4=&NA_REAL,
           *r_Fr_efmit_S2M1=&NA_REAL, *r_Fr_efmit_S2M2=&NA_REAL, *r_Fr_efmit_S2M3=&NA_REAL, *r_Fr_efmit_S2M4=&NA_REAL,
           *r_Fr_efmit_S3M1=&NA_REAL, *r_Fr_efmit_S3M2=&NA_REAL, *r_Fr_efmit_S3M3=&NA_REAL, *r_Fr_efmit_S3M4=&NA_REAL,
           *r_Fr_efmit_S4M1=&NA_REAL, *r_Fr_efmit_S4M2=&NA_REAL, *r_Fr_efmit_S4M3=&NA_REAL, *r_Fr_efmit_S4M4=&NA_REAL,
           *r_F_efmit_S1M1=&NA_REAL, *r_F_efmit_S1M2=&NA_REAL, *r_F_efmit_S1M3=&NA_REAL, *r_F_efmit_S1M4=&NA_REAL,
           *r_F_efmit_S2M1=&NA_REAL, *r_F_efmit_S2M2=&NA_REAL, *r_F_efmit_S2M3=&NA_REAL, *r_F_efmit_S2M4=&NA_REAL,
           *r_F_efmit_S3M1=&NA_REAL, *r_F_efmit_S3M2=&NA_REAL, *r_F_efmit_S3M3=&NA_REAL, *r_F_efmit_S3M4=&NA_REAL,
           *r_F_efmit_S4M1=&NA_REAL, *r_F_efmit_S4M2=&NA_REAL, *r_F_efmit_S4M3=&NA_REAL, *r_F_efmit_S4M4=&NA_REAL,
//           *r_N_e0t_S1M1=&NA_REAL, *r_N_e0t_S2M2=&NA_REAL, *r_N_e0t_S3M3=&NA_REAL, *r_N_e0t_S4M4=&NA_REAL,
//           *r_N_ei0_S1M1=&NA_REAL, *r_N_ei0_S1M2=&NA_REAL, *r_N_ei0_S1M3=&NA_REAL, *r_N_ei0_S1M4=&NA_REAL,
           *r_iniNt0q_S1M1=&NA_REAL, *r_iniNt0q_S1M2=&NA_REAL, *r_iniNt0q_S1M3=&NA_REAL, *r_iniNt0q_S1M4=&NA_REAL,
           *r_iniNt0q_S2M1=&NA_REAL, *r_iniNt0q_S2M2=&NA_REAL, *r_iniNt0q_S2M3=&NA_REAL, *r_iniNt0q_S2M4=&NA_REAL,
           *r_iniNt0q_S3M1=&NA_REAL, *r_iniNt0q_S3M2=&NA_REAL, *r_iniNt0q_S3M3=&NA_REAL, *r_iniNt0q_S3M4=&NA_REAL,
           *r_iniNt0q_S4M1=&NA_REAL, *r_iniNt0q_S4M2=&NA_REAL, *r_iniNt0q_S4M3=&NA_REAL, *r_iniNt0q_S4M4=&NA_REAL,
           *r_matwt_M1=&NA_REAL, *r_matwt_M2=&NA_REAL, *r_matwt_M3=&NA_REAL, *r_matwt_M4=&NA_REAL;


    SEXP ans_Z_eit=R_NilValue, ans_N_eit=R_NilValue, ans_Fbar_et=R_NilValue, ans_B_et=R_NilValue, ans_SSB_et=R_NilValue,
         ans_Z_eit_G1=R_NilValue, ans_Z_eit_G2=R_NilValue, ans_N_eit_G1=R_NilValue, ans_N_eit_G2=R_NilValue;


    SEXP ans_Z_eit_S1M1=R_NilValue, ans_Z_eit_S1M2=R_NilValue, ans_Z_eit_S1M3=R_NilValue, ans_Z_eit_S1M4=R_NilValue,
         ans_Z_eit_S2M1=R_NilValue, ans_Z_eit_S2M2=R_NilValue, ans_Z_eit_S2M3=R_NilValue, ans_Z_eit_S2M4=R_NilValue,
         ans_Z_eit_S3M1=R_NilValue, ans_Z_eit_S3M2=R_NilValue, ans_Z_eit_S3M3=R_NilValue, ans_Z_eit_S3M4=R_NilValue,
         ans_Z_eit_S4M1=R_NilValue, ans_Z_eit_S4M2=R_NilValue, ans_Z_eit_S4M3=R_NilValue, ans_Z_eit_S4M4=R_NilValue,
         ans_N_eit_S1M1=R_NilValue, ans_N_eit_S1M2=R_NilValue, ans_N_eit_S1M3=R_NilValue, ans_N_eit_S1M4=R_NilValue,
         ans_N_eit_S2M1=R_NilValue, ans_N_eit_S2M2=R_NilValue, ans_N_eit_S2M3=R_NilValue, ans_N_eit_S2M4=R_NilValue,
         ans_N_eit_S3M1=R_NilValue, ans_N_eit_S3M2=R_NilValue, ans_N_eit_S3M3=R_NilValue, ans_N_eit_S3M4=R_NilValue,
         ans_N_eit_S4M1=R_NilValue, ans_N_eit_S4M2=R_NilValue, ans_N_eit_S4M3=R_NilValue, ans_N_eit_S4M4=R_NilValue;

    double fmax, sumWt;

if (ind_t==0) {


    //� t=0, pr�paration des outputs

    PROTECT(rnames_Esp = allocVector(STRSXP, nbE));

    setAttrib(out_Z_eit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Fbar_et, R_NamesSymbol, rnames_Esp);
    setAttrib(out_B_et, R_NamesSymbol, rnames_Esp);
    setAttrib(out_SSB_et, R_NamesSymbol, rnames_Esp);

    setAttrib(out_Z_eit_G1, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_G1, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_G2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_G2, R_NamesSymbol, rnames_Esp);

    setAttrib(out_Z_eit_S1M1, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S1M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S1M3, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S1M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S2M1, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S2M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S2M3, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S2M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S3M1, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S3M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S3M3, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S3M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S4M1, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S4M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Z_eit_S4M3, R_NamesSymbol, rnames_Esp); setAttrib(out_Z_eit_S4M4, R_NamesSymbol, rnames_Esp);

    setAttrib(out_N_eit_S1M1, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S1M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S1M3, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S1M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S2M1, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S2M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S2M3, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S2M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S3M1, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S3M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S3M3, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S3M4, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S4M1, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S4M2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_N_eit_S4M3, R_NamesSymbol, rnames_Esp); setAttrib(out_N_eit_S4M4, R_NamesSymbol, rnames_Esp);

}


for (int e = 0 ; e < nbE ; e++) { //Rprintf("G1one");fichier << "G1one" << endl;

                                    double *Ztemp = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,e))));
//fichier << "G2one" << endl;
                                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));//Rprintf("g1");fichier << "g1" << endl;
                                    PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));//Rprintf("g2");fichier << "g2" << endl;

                                    nbI = length(getListElement(elmt, "modI"));//Rprintf("g3");fichier << "g3" << endl;



                             if ((Qvec[e]==1) & (Svec[e]==0)) {
//PrintValue(out_Fr_fmi_S1M1);
                                    PROTECT(v_Fr_efmit = getListElement(out_Fr_fmi, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit = getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));

                                    PROTECT(v_Fr_efmit_S1M1 = getListElement(out_Fr_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));//Rprintf("g7");fichier << "g7" << endl;
                                    PROTECT(v_Fr_efmit_S1M2 = getListElement(out_Fr_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));//Rprintf("g8");fichier << "g8" << endl;
                                    PROTECT(v_Fr_efmit_S1M3 = getListElement(out_Fr_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S1M4 = getListElement(out_Fr_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S2M1 = getListElement(out_Fr_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S2M2 = getListElement(out_Fr_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S2M3 = getListElement(out_Fr_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S2M4 = getListElement(out_Fr_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S3M1 = getListElement(out_Fr_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S3M2 = getListElement(out_Fr_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S3M3 = getListElement(out_Fr_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S3M4 = getListElement(out_Fr_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S4M1 = getListElement(out_Fr_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S4M2 = getListElement(out_Fr_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S4M3 = getListElement(out_Fr_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_S4M4 = getListElement(out_Fr_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));
//Rprintf("G1.1");fichier << "G1.1" << endl;
                                    PROTECT(v_F_efmit_S1M1 = getListElement(out_F_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S1M2 = getListElement(out_F_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S1M3 = getListElement(out_F_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S1M4 = getListElement(out_F_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S2M1 = getListElement(out_F_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S2M2 = getListElement(out_F_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S2M3 = getListElement(out_F_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S2M4 = getListElement(out_F_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S3M1 = getListElement(out_F_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S3M2 = getListElement(out_F_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S3M3 = getListElement(out_F_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S3M4 = getListElement(out_F_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S4M1 = getListElement(out_F_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S4M2 = getListElement(out_F_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S4M3 = getListElement(out_F_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_S4M4 = getListElement(out_F_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));
//Rprintf("G1.2");fichier << "G1.2" << endl;
                                    PROTECT(v_N_ei0_S1M1 = getListElement(elmt, "Nt0_S1M1"));
                                    PROTECT(v_N_ei0_S1M2 = getListElement(elmt, "Nt0_S1M2"));
                                    PROTECT(v_N_ei0_S1M3 = getListElement(elmt, "Nt0_S1M3"));
                                    PROTECT(v_N_ei0_S1M4 = getListElement(elmt, "Nt0_S1M4"));

                                    PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                                    PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                                    PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                                    PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                                    PROTECT(v_iniNt0q_S1M1 = getListElement(elmt, "iniNt0q_S1M1"));
                                    PROTECT(v_iniNt0q_S1M2 = getListElement(elmt, "iniNt0q_S1M2"));
                                    PROTECT(v_iniNt0q_S1M3 = getListElement(elmt, "iniNt0q_S1M3"));
                                    PROTECT(v_iniNt0q_S1M4 = getListElement(elmt, "iniNt0q_S1M4"));
                                    PROTECT(v_iniNt0q_S2M1 = getListElement(elmt, "iniNt0q_S2M1"));
                                    PROTECT(v_iniNt0q_S2M2 = getListElement(elmt, "iniNt0q_S2M2"));
                                    PROTECT(v_iniNt0q_S2M3 = getListElement(elmt, "iniNt0q_S2M3"));
                                    PROTECT(v_iniNt0q_S2M4 = getListElement(elmt, "iniNt0q_S2M4"));
                                    PROTECT(v_iniNt0q_S3M1 = getListElement(elmt, "iniNt0q_S3M1"));
                                    PROTECT(v_iniNt0q_S3M2 = getListElement(elmt, "iniNt0q_S3M2"));
                                    PROTECT(v_iniNt0q_S3M3 = getListElement(elmt, "iniNt0q_S3M3"));
                                    PROTECT(v_iniNt0q_S3M4 = getListElement(elmt, "iniNt0q_S3M4"));
                                    PROTECT(v_iniNt0q_S4M1 = getListElement(elmt, "iniNt0q_S4M1"));
                                    PROTECT(v_iniNt0q_S4M2 = getListElement(elmt, "iniNt0q_S4M2"));
                                    PROTECT(v_iniNt0q_S4M3 = getListElement(elmt, "iniNt0q_S4M3"));
                                    PROTECT(v_iniNt0q_S4M4 = getListElement(elmt, "iniNt0q_S4M4"));

                                    PROTECT(v_matwt_M1 = getListElement(elmt, "matwt_M1"));
                                    PROTECT(v_matwt_M2 = getListElement(elmt, "matwt_M2"));
                                    PROTECT(v_matwt_M3 = getListElement(elmt, "matwt_M3"));
                                    PROTECT(v_matwt_M4 = getListElement(elmt, "matwt_M4"));

                                    PROTECT(v_M_ei = getListElement(elmt, "M_i"));//Rprintf("g4");fichier << "g4" << endl;
                                    PROTECT(v_w_ei = getListElement(elmt, "wStock_i"));//Rprintf("g5");fichier << "g5" << endl;
                                    PROTECT(v_mat_ei = getListElement(elmt, "mat_i"));//Rprintf("g6");fichier << "g6" << endl;
                                    PROTECT(v_Fbar = getListElement(elmt, "Fbar"));

                                    PROTECT(dimCst_Fr_efmit = getAttrib(v_Fr_efmit_S1M1, install("DimCst")));
                                    PROTECT(dimCst_N_ei0 = getAttrib(v_N_ei0_S1M1, install("DimCst")));
                                    PROTECT(dimCst_N_e0t = getAttrib(v_N_e0t_S1M1, install("DimCst")));
                                    PROTECT(dimCst_M_ei = getAttrib(v_M_ei, install("DimCst")));
                                    PROTECT(dimCst_w_ei = getAttrib(v_w_ei, install("DimCst")));
                                    PROTECT(dimCst_mat_ei = getAttrib(v_mat_ei, install("DimCst")));

                                    PROTECT(v_r = getListElement(elmt, "r"));//Rprintf("g4");   fichier << "g4" << endl;           //added SPiCT 19/07/2016
                                    PROTECT(v_K = getListElement(elmt, "K"));//Rprintf("g5");    fichier << "g5" << endl;          //
                                    PROTECT(v_n = getListElement(elmt, "n"));                              //
                                    PROTECT(v_B = getListElement(elmt, "B_i"));



                             }  else if ((Qvec[e]==0) & (Svec[e]==0)) {


                                    PROTECT(v_r = getListElement(elmt, "r"));//Rprintf("g4");   fichier << "g4" << endl;           //added SPiCT 19/07/2016
                                    PROTECT(v_K = getListElement(elmt, "K"));//Rprintf("g5");    fichier << "g5" << endl;          //
                                    PROTECT(v_n = getListElement(elmt, "n"));                              //
                                    PROTECT(v_B = getListElement(elmt, "B_i"));
                                    PROTECT(v_M_ei = getListElement(elmt, "M_i")); //PrintValue(v_M_ei);
                                    PROTECT(v_w_ei = getListElement(elmt, "wStock_i")); //PrintValue(v_w_ei);//Rprintf("g5");fichier << "g5" << endl;
                                    PROTECT(v_mat_ei = getListElement(elmt, "mat_i")); //PrintValue(v_mat_ei);//Rprintf("g6");fichier << "g6" << endl;
                                    PROTECT(v_Fbar = getListElement(elmt, "Fbar"));

                                    PROTECT(v_Fr_efmit = getListElement(out_Fr_fmi, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit = getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_N_ei0 = getListElement(elmt, "N_it0"));
                                    PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
//Rprintf("G1.3");fichier << "G1.3" << endl;
                                    PROTECT(dimCst_Fr_efmit = getAttrib(v_Fr_efmit, install("DimCst")));
                                    PROTECT(dimCst_N_ei0 = getAttrib(v_N_ei0, install("DimCst")));
                                    PROTECT(dimCst_N_e0t = getAttrib(v_N_e0t, install("DimCst")));
                                    PROTECT(dimCst_M_ei = getAttrib(v_M_ei, install("DimCst")));
                                    PROTECT(dimCst_w_ei = getAttrib(v_w_ei, install("DimCst")));
                                    PROTECT(dimCst_mat_ei = getAttrib(v_mat_ei, install("DimCst")));



                             } else if ((Qvec[e]==0) & (Svec[e]==1)){
                                    PROTECT(v_Fr_efmit_G1 = getListElement(out_Fr_fmi_G1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_Fr_efmit_G2 = getListElement(out_Fr_fmi_G2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_G1 = getListElement(out_F_fmi_G1, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit_G2 = getListElement(out_F_fmi_G2, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_N_ei0_G1 = getListElement(elmt, "N_it0_G1"));
                                    PROTECT(v_N_ei0_G2 = getListElement(elmt, "N_it0_G2"));
                                    PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));
                                    PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2"));
//Rprintf("G1.3");fichier << "G1.3" << endl;

                                    PROTECT(v_M_ei_G1 = getListElement(elmt, "M_i_G1"));//Rprintf("g4");fichier << "g4" << endl;
                                    PROTECT(v_M_ei_G2 = getListElement(elmt, "M_i_G2"));
                                    PROTECT(v_w_ei_G1 = getListElement(elmt, "wStock_i_G1"));
                                    PROTECT(v_w_ei_G2 = getListElement(elmt, "wStock_i_G2"));//Rprintf("g5");fichier << "g5" << endl;
                                    PROTECT(v_mat_ei_G1 = getListElement(elmt, "mat_i_G1"));
                                    PROTECT(v_mat_ei_G2 = getListElement(elmt, "mat_i_G2"));//Rprintf("g6");fichier << "g6" << endl;
                                    PROTECT(v_B = getListElement(elmt, "B_i"));
                                    PROTECT(v_Fbar_G1 = getListElement(elmt, "Fbar_G1"));
                                    PROTECT(v_Fbar_G2 = getListElement(elmt, "Fbar_G2"));

                                    PROTECT(dimCst_Fr_efmit = getAttrib(v_Fr_efmit_G1, install("DimCst")));
                                    PROTECT(dimCst_N_ei0 = getAttrib(v_N_ei0_G1, install("DimCst")));
                                    PROTECT(dimCst_N_e0t = getAttrib(v_N_e0t_G1, install("DimCst")));
                                    PROTECT(dimCst_M_ei = getAttrib(v_M_ei_G1, install("DimCst")));
                                    PROTECT(dimCst_w_ei = getAttrib(v_w_ei_G1, install("DimCst")));
                                    PROTECT(dimCst_mat_ei = getAttrib(v_mat_ei_G1, install("DimCst")));

                             }




//Rprintf("G1.4");fichier << "G1.4" << endl;
                                    dim_M_ei = INTEGER(dimCst_M_ei); //fichier << "dim_M_ei =" << dim_M_ei << endl;
                                    dim_w_ei = INTEGER(dimCst_w_ei); //fichier << dim_w_ei << endl;
                                    dim_mat_ei = INTEGER(dimCst_mat_ei);
                                    dim_N_ei0 = INTEGER(dimCst_N_ei0);
                                    dim_N_e0t = INTEGER(dimCst_N_e0t);

                            if (nbI>1) {
                                    //tests sur les dimensions (pas pour SPiCT) :

                                    if ((dim_M_ei[0]!=0) | (dim_M_ei[1]!=0) |
                                        ((dim_M_ei[2]!=0) & (dim_M_ei[2]!=nbI)) | ((dim_M_ei[3]!=0) & (dim_M_ei[3]!=nbT))) //on laisse une ouverture pour un indice temporel
                                    {
                                        error("Non_homogeneous dimensions in M_ei element. Check .ini biological parameters files !!\n");
                                    }


                                    if ((dim_w_ei[0]!=0) | (dim_w_ei[1]!=0) |
                                        ((dim_w_ei[2]!=0) & (dim_w_ei[2]!=nbI)) | (dim_w_ei[3]!=0))
                                    {
                                        error("Non_homogeneous dimensions in w_ei element. Check .ini biological parameters files !!\n");
                                    }


                                    if ((dim_mat_ei[0]!=0) | (dim_mat_ei[1]!=0) |
                                        ((dim_mat_ei[2]!=0) & (dim_mat_ei[2]!=nbI)) | (dim_mat_ei[3]!=0))
                                    {
                                        error("Non_homogeneous dimensions in mat_ei element. Check .ini biological parameters files !!\n");
                                    }


                                    if ((dim_N_ei0[0]!=0) | (dim_N_ei0[1]!=0) |
                                        ((dim_N_ei0[2]!=0) & (dim_N_ei0[2]!=nbI))) // | (dim_N_ei0[3]!=0)) --> peu importe, on ne prendra de toute fa�on que la donn�e � t0
                                    {
                                        error("Non_homogeneous dimensions in N_ei0 element. Check .ini biological parameters files !!\n");
                                    }


                                    if ((dim_N_e0t[0]!=0) | (dim_N_e0t[1]!=0) |
                                        (dim_N_e0t[2]!=0) | ((dim_N_e0t[3]!=0) & (dim_N_e0t[3]!=nbT)))
                                    {
                                        error("Non_homogeneous dimensions in N_e0t element. Check .ini biological parameters files !!\n");
                                    }

                            }

                                    dim_Fr_efmit = INTEGER(dimCst_Fr_efmit);
                                    if (((dim_Fr_efmit[0]!=0) & (dim_Fr_efmit[0]!=nbF)) | ((dim_Fr_efmit[1]!=0) & (dim_Fr_efmit[1]!=nbM)) |
                                        ((dim_Fr_efmit[2]!=0) & (dim_Fr_efmit[2]!=nbI)) | ((dim_Fr_efmit[3]!=0) & (dim_Fr_efmit[3]!=nbT)))
                                    {
                                        error("Non_homogeneous dimensions in Fr_efmit element. Check .ini biological parameters files !!\n");
                                    }

                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // initialisation de out_Z_eit
                                    //---------
//Rprintf("G1.5");fichier << "G1.5" << endl;
                                    //on d�termine l'attribut Dimension de Z_eit
                                    PROTECT(dimCst1 = allocVector(INTSXP, 4));
                                    dimC1 = INTEGER(dimCst1);
                                    dimC1[0] = 0 ; dimC1[1] = 0 ; dimC1[2] = imax2(dim_M_ei[2] , dim_Fr_efmit[2]);
                                    dimC1[3] = imax2(dim_M_ei[3] , dim_Fr_efmit[3]);
                                    //Rprintf("Fr %i %i %i %i",dim_Fr_efmit[0],dim_Fr_efmit[1],dim_Fr_efmit[2],dim_Fr_efmit[3]);
                                    //Rprintf("M %i %i %i %i",dim_M_ei[0],dim_M_ei[1],dim_M_ei[2],dim_M_ei[3]);
                                    //Rprintf("C1 %i %i %i %i",dimC1[0],dimC1[1],dimC1[2],dimC1[3]);

                                    int count = 0, prod = 1, count2 = 0, count3 = 0, count4 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC1[k]>0) {
                                            count++;
                                            prod = prod * dimC1[k];
                                        }
                                    }
//Rprintf("G1.6");fichier << "G1.6" << endl;
                                    PROTECT(Dim1 = allocVector(INTSXP, count));//Rprintf("G1.61");
                                    dim1 = INTEGER(Dim1);//Rprintf("G1.62");

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC1[k]>0) {
                                            dim1[count2] = dimC1[k];
                                            count2++;
                                            }
                                    }
//Rprintf("G1.63");//Rprintf("%i ",prod);fichier << "G1.63" << endl;
                            if (ind_t==0){
                                    //on cr�e le tableau r�sultat pour l'esp�ce en question

//Rprintf("G1.7");fichier << "G1.7" << endl;
                                    if ((Qvec[e]==1) & (Svec[e]==0)) {
                                    PROTECT(ans_Z_eit = NEW_NUMERIC(prod));//Rprintf("G1.64");fichier << "G1.64" << endl;
                                    setAttrib(ans_Z_eit, R_DimSymbol, Dim1);//Rprintf("G1.65");fichier << "G1.65" << endl;

                                    PROTECT(ans_Z_eit_S1M1 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S1M1, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S1M2 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S1M2, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S1M3 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S1M3, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S1M4 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S1M4, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S2M1 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S2M1, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S2M2 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S2M2, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S2M3 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S2M3, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S2M4 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S2M4, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S3M1 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S3M1, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S3M2 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S3M2, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S3M3 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S3M3, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S3M4 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S3M4, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S4M1 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S4M1, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S4M2 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S4M2, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S4M3 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S4M3, R_DimSymbol, Dim1);
                                    PROTECT(ans_Z_eit_S4M4 = NEW_NUMERIC(prod)); setAttrib(ans_Z_eit_S4M4, R_DimSymbol, Dim1);

                                    } else if ((Qvec[e]==0) & (Svec[e]==0)){
                                    PROTECT(ans_Z_eit = NEW_NUMERIC(prod));//Rprintf("G1.64");fichier << "G1.64" << endl;
                                    setAttrib(ans_Z_eit, R_DimSymbol, Dim1);//Rprintf("G1.65");fichier << "G1.65" << endl;

                                    } else if ((Qvec[e]==0) & (Svec[e]==1)){
                                    PROTECT(ans_Z_eit_G1 = NEW_NUMERIC(prod));//Rprintf("G1.64");fichier << "G1.64" << endl;
                                    setAttrib(ans_Z_eit_G1, R_DimSymbol, Dim1);//Rprintf("G1.65");fichier << "G1.65" << endl;
                                    PROTECT(ans_Z_eit_G2 = NEW_NUMERIC(prod));//Rprintf("G1.64");fichier << "G1.64" << endl;
                                    setAttrib(ans_Z_eit_G2, R_DimSymbol, Dim1);

                                    }
//Rprintf("G1.8");fichier << "G1.8" << endl;
                                    PROTECT(dimnames1 = allocVector(VECSXP,count));
                                    if (dimC1[0]>0) {SET_VECTOR_ELT(dimnames1, count3, fleetList) ; count3++;}
                                    if (dimC1[1]>0) {SET_VECTOR_ELT(dimnames1, count3, metierList) ; count3++;}
                                    if (dimC1[2]>0) {SET_VECTOR_ELT(dimnames1, count3, intAge) ; count3++;}
                                    if (dimC1[3]>0) {SET_VECTOR_ELT(dimnames1, count3, times) ; count3++;}



                                    if ((Qvec[e]==1) & (Svec[e]==0)) {
                                    rans_Z_eit = REAL(ans_Z_eit);

                                    rans_Z_eit_S1M1 = REAL(ans_Z_eit_S1M1); rans_Z_eit_S1M2 = REAL(ans_Z_eit_S1M2);
                                    rans_Z_eit_S1M3 = REAL(ans_Z_eit_S1M3); rans_Z_eit_S1M4 = REAL(ans_Z_eit_S1M4);
                                    rans_Z_eit_S2M1 = REAL(ans_Z_eit_S2M1); rans_Z_eit_S2M2 = REAL(ans_Z_eit_S2M2);
                                    rans_Z_eit_S2M3 = REAL(ans_Z_eit_S2M3); rans_Z_eit_S2M4 = REAL(ans_Z_eit_S2M4);
                                    rans_Z_eit_S3M1 = REAL(ans_Z_eit_S3M1); rans_Z_eit_S3M2 = REAL(ans_Z_eit_S3M2);
                                    rans_Z_eit_S3M3 = REAL(ans_Z_eit_S3M3); rans_Z_eit_S3M4 = REAL(ans_Z_eit_S3M4);
                                    rans_Z_eit_S4M1 = REAL(ans_Z_eit_S4M1); rans_Z_eit_S4M2 = REAL(ans_Z_eit_S4M2);
                                    rans_Z_eit_S4M3 = REAL(ans_Z_eit_S4M3); rans_Z_eit_S4M4 = REAL(ans_Z_eit_S4M4);
//Rprintf("G1.9");fichier << "G1.9" << endl;
                                    } else if ((Qvec[e]==0) & (Svec[e]==0)){
                                        rans_Z_eit = REAL(ans_Z_eit);
                                    } else if ((Qvec[e]==0) & (Svec[e]==1)){
                                        rans_Z_eit_G1 = REAL(ans_Z_eit_G1);
                                        rans_Z_eit_G2 = REAL(ans_Z_eit_G2);
                                    }

                            } else {

                                    if ((Qvec[e]==1) & (Svec[e]==0)) {
                                    rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit, e));

                                    rans_Z_eit_S1M1 = REAL(VECTOR_ELT(out_Z_eit_S1M1, e)); rans_Z_eit_S1M2 = REAL(VECTOR_ELT(out_Z_eit_S1M2, e));
                                    rans_Z_eit_S1M3 = REAL(VECTOR_ELT(out_Z_eit_S1M3, e)); rans_Z_eit_S1M4 = REAL(VECTOR_ELT(out_Z_eit_S1M4, e));
                                    rans_Z_eit_S2M1 = REAL(VECTOR_ELT(out_Z_eit_S2M1, e)); rans_Z_eit_S2M2 = REAL(VECTOR_ELT(out_Z_eit_S2M2, e));
                                    rans_Z_eit_S2M3 = REAL(VECTOR_ELT(out_Z_eit_S2M3, e)); rans_Z_eit_S2M4 = REAL(VECTOR_ELT(out_Z_eit_S2M4, e));
                                    rans_Z_eit_S3M1 = REAL(VECTOR_ELT(out_Z_eit_S3M1, e)); rans_Z_eit_S3M2 = REAL(VECTOR_ELT(out_Z_eit_S3M2, e));
                                    rans_Z_eit_S3M3 = REAL(VECTOR_ELT(out_Z_eit_S3M3, e)); rans_Z_eit_S3M4 = REAL(VECTOR_ELT(out_Z_eit_S3M4, e));
                                    rans_Z_eit_S4M1 = REAL(VECTOR_ELT(out_Z_eit_S4M1, e)); rans_Z_eit_S4M2 = REAL(VECTOR_ELT(out_Z_eit_S4M2, e));
                                    rans_Z_eit_S4M3 = REAL(VECTOR_ELT(out_Z_eit_S4M3, e)); rans_Z_eit_S4M4 = REAL(VECTOR_ELT(out_Z_eit_S4M4, e));

                                    } else if ((Qvec[e]==0) & (Svec[e]==0)){
                                        rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit, e));

                                    } else if ((Qvec[e]==0) & (Svec[e]==1)){
                                        rans_Z_eit_G1 = REAL(VECTOR_ELT(out_Z_eit_G1, e));
                                        rans_Z_eit_G2 = REAL(VECTOR_ELT(out_Z_eit_G2, e));
                                    }

                            }
//Rprintf("G1.10");fichier << "G1.10" << endl;

sumWt = 0.0; fmax = 0.0;

    if ((Qvec[e]==0) & (Svec[e]==0)) { // Age-based or global
//Rprintf("Start_Annual\n");fichier << "Start_Annual\n" << endl;
//Rprintf("G1.10.1");fichier << "G1.10.1" << endl;
                                    r_Fr_efmit = REAL(v_Fr_efmit);
                                    r_F_efmit = REAL(v_F_efmit);
                                    r_M_ei = REAL(v_M_ei);
                                    r_Fbar = REAL(v_Fbar);
//Rprintf("G1.10.2");fichier << "G1.10.2" << endl;
                                    //facteurs des indices
                                    PROTECT(dFACT1 = iDim(dimC1));
                                    PROTECT(dFACT2 = iDim(dim_Fr_efmit));
                                    PROTECT(dFACT3 = iDim(dim_M_ei));

                                    int *fact1_D = INTEGER(dFACT1);
                                    int *fact2_D = INTEGER(dFACT2);
                                    int *fact3_D = INTEGER(dFACT3);
//Rprintf("G1.10.3");fichier << "G1.10.3" << endl;
                                    double *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));

//Rprintf("G1.10.4");fichier << "G1.10.4" << endl;
                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0, tempCap = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        }

                                    if (Zoptim_use & (e==eTemp)) {

                                     rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                        Zoptim[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
//Rprintf("G1.11");fichier << "G1.11" << endl;

                                    } else {

                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i[ind_i + nbI*ind_t];

                                           //on initialise aussi Ztemp (attention : index� � partir de 1)

                                           //if (e==eTemp) Ztemp[ind_i+1] = rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        }

                                        Ztemp[ind_i] = rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                    }

                                    //on en profite pour calculer Fbar

                                    fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                                    sumWt = sumWt + r_Fbar[ind_i];

                                    }


                            if (ind_t==0) {

                                    setAttrib(ans_Z_eit, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit, install("DimCst"), dimCst1);

                                    SET_VECTOR_ELT(out_Z_eit, e, ans_Z_eit);

                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 9, dimCst_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 10, v_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 11, v_M_ei);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 12, dFACT1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 13, dFACT2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 14, dFACT3);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 59, v_Fbar);

//Rprintf("G1.12");fichier << "G1.12" << endl;

                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de N_eit
                                    //---------

                                    //on d�termine l'attribut Dimension de N_eit
                                    PROTECT(dimCst2 = allocVector(INTSXP, 4));
                                    dimC2 = INTEGER(dimCst2);

                                    dimC2[0] = 0 ; dimC2[1] = 0 ; dimC2[2] = imax2(dimC1[2] , dim_N_ei0[2]) ; dimC2[3] = nbT;

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            count++;
                                            prod = prod * dimC2[k];
                                        }
                                    }

                                    PROTECT(Dim2 = allocVector(INTSXP, count));
                                    dim2 = INTEGER(Dim2);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            dim2[count2] = dimC2[k];
                                            count2++;
                                            }
                                    }

//Rprintf("G1.13");fichier << "G1.13" << endl;
                            if (ind_t==0) {

                                    //on cr�e le tableau r�sultat pour l'esp�ce en question
                                    PROTECT(ans_N_eit = NEW_NUMERIC(prod));
                                    setAttrib(ans_N_eit, R_DimSymbol, Dim2);

                                    PROTECT(dimnames2 = allocVector(VECSXP,count));
                                    if (dimC2[0]>0) {SET_VECTOR_ELT(dimnames2, count3, fleetList) ; count3++;}
                                    if (dimC2[1]>0) {SET_VECTOR_ELT(dimnames2, count3, metierList) ; count3++;}
                                    if (dimC2[2]>0) {SET_VECTOR_ELT(dimnames2, count3, intAge) ; count3++;}
                                    if (dimC2[3]>0) {SET_VECTOR_ELT(dimnames2, count3, times) ; count3++;}

                                    rans_N_eit = REAL(ans_N_eit);

                            } else {

                                    rans_N_eit = REAL(VECTOR_ELT(out_N_eit, e));

                            }

                                    r_N_ei0 = REAL(v_N_ei0);

                                    r_N_e0t = REAL(v_N_e0t);


                                    //facteurs des indices
                                    PROTECT(dFACT4 = iDim(dimC2));
                                    PROTECT(dFACT5 = iDim(dim_N_ei0));
                                    PROTECT(dFACT6 = iDim(dim_N_e0t));

                                    fact4_D = INTEGER(dFACT4);
                                    int *fact5_D = INTEGER(dFACT5);
                                    int *fact6_D = INTEGER(dFACT6);


                                    //ajout 24/04/2018 pour prise en compte for�age recrutement XSA
                                    if ((!isNull(getListElement(recList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & (nbI>1) & Reality) {
                                       SRInd[e]=0;     //on n'utilise alors pas l'information de l'objet Args
                                       r_N_e0t[ind_t*fact6_D[3]] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[ind_t];
                                       if (ind_t<(nbT-1)) r_N_e0t[ind_t*fact6_D[3]+1] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[ind_t+1];
                                    }



                                    //�quation

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                            if (ind_i == 0) { //recrutement

                                                if ((SRInd[e]==1) & (ind_t>0)) {

                                                    if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                            REAL(VECTOR_ELT(out_SRmod,e))[ind_t];

                                                    } else {

                                                     if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial d�fini

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                                    }}
//Rprintf("G1.14");fichier << "G1.14" << endl;

                                                } else {

                                                    if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial d�fini

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                                    }
                                                }

                                            } else {

                                                if (ind_t == 0) {

                                                    rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                                } else {

                                                    if (ind_i == (nbI-1)) {  //groupe d'�ge +

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                    }
                                                }
                                            }
                                    }


                            if (ind_t==0) {

                                    setAttrib(ans_N_eit, R_DimNamesSymbol, dimnames2);
                                    setAttrib(ans_N_eit, install("DimCst"), dimCst2);

                                    SET_VECTOR_ELT(out_N_eit, e, ans_N_eit);

                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 15, v_N_ei0);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 16, v_N_e0t);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 17, dFACT4);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 18, dFACT5);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 19, dFACT6);

//Rprintf("End_Annual\n");fichier << "End_Annual\n" << endl;

    } else if ((Qvec[e]==0) & (Svec[e]==1)) { // Age and sex-based
//Rprintf("Start_Sex based\n");fichier << "Start_Sex based\n" << endl;



//Rprintf("G1.15");fichier << "G1.15" << endl;
                                    r_Fr_efmit_G1 = REAL(v_Fr_efmit_G1);
                                    r_Fr_efmit_G2 = REAL(v_Fr_efmit_G2);

//Rprintf("G1.16");fichier << "G1.16" << endl;
                                    //r_F_efmit_G1 = REAL(v_F_efmit_G1);
                                    //r_F_efmit_G2 = REAL(v_F_efmit_G2);

//Rprintf("G1.17");fichier << "G1.17" << endl;
                                    r_M_ei_G1 = REAL(v_M_ei_G1);
                                    r_M_ei_G2 = REAL(v_M_ei_G2);
                                    r_Fbar_G1 = REAL(v_Fbar_G1);
                                    r_Fbar_G2 = REAL(v_Fbar_G2);
//Rprintf("G1.18");fichier << "G1.18" << endl;
                                    //facteurs des indices
                                    PROTECT(dFACT1 = iDim(dimC1));
                                    PROTECT(dFACT2 = iDim(dim_Fr_efmit));
                                    PROTECT(dFACT3 = iDim(dim_M_ei));

                                    int *fact1_D = INTEGER(dFACT1);
                                    int *fact2_D = INTEGER(dFACT2);
                                    int *fact3_D = INTEGER(dFACT3);

//Rprintf("G1.19");fichier << "G1.19" << endl;
                                    //double *r_Foth_i_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 224));
                                    //double *r_Foth_i_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 225));

//Rprintf("G1.20");fichier << "G1.20" << endl;
                                    double *r_Froth_i_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 226));
                                    double *r_Froth_i_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 227));

//Rprintf("G1.21");fichier << "G1.21" << endl;
                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_G1[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : G1 -> 0

                                        //if (e==eTemp) Ztemp[ind_i+1+(0*nbI)] = rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(0*nbI)] = rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_G2[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : G2 -> 1

                                        //if (e==eTemp) Ztemp[ind_i+1+(1*nbI)] = rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(1*nbI)] = rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                    }


                              if (ind_t==0) {

                                    setAttrib(ans_Z_eit_G1, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_G2, R_DimNamesSymbol, dimnames1);

                                    setAttrib(ans_Z_eit_G1, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_G2, install("DimCst"), dimCst1);

                                    SET_VECTOR_ELT(out_Z_eit_G1, e, ans_Z_eit_G1);
                                    SET_VECTOR_ELT(out_Z_eit_G2, e, ans_Z_eit_G2);
                                    ////Rprintf("AAAhhh!!!");//PrintValue(out_Z_eit_G1);
                              }


                                    //---------
                                    // calcul de N_eit
                                    //---------

                                    //on d�termine l'attribut Dimension de N_eit
                                    PROTECT(dimCst2 = allocVector(INTSXP, 4));
                                    dimC2 = INTEGER(dimCst2);

                                    dimC2[0] = 0 ; dimC2[1] = 0 ; dimC2[2] = imax2(dimC1[2] , dim_N_ei0[2]) ; dimC2[3] = nbT;

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            count++;
                                            prod = prod * dimC2[k];
                                        }
                                    }

                                    //fichier << "prod = " << prod << endl;

                                    PROTECT(Dim2 = allocVector(INTSXP, count));
                                    dim2 = INTEGER(Dim2);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            dim2[count2] = dimC2[k];
                                            count2++;
                                            }
                                    }


                            if (ind_t==0) {

                                    //on cr�e le tableau r�sultat pour l'esp�ce en question

                                    PROTECT(dimnames2 = allocVector(VECSXP,count));
                                    if (dimC2[0]>0) {SET_VECTOR_ELT(dimnames2, count3, fleetList) ; count3++;}
                                    if (dimC2[1]>0) {SET_VECTOR_ELT(dimnames2, count3, metierList) ; count3++;}
                                    if (dimC2[2]>0) {SET_VECTOR_ELT(dimnames2, count3, intAge) ; count3++;}
                                    if (dimC2[3]>0) {SET_VECTOR_ELT(dimnames2, count3, times) ; count3++;}

                                    PROTECT(ans_N_eit_G1 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_G1, R_DimSymbol, Dim2); rans_N_eit_G1 = REAL(ans_N_eit_G1);
                                    PROTECT(ans_N_eit_G2 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_G2, R_DimSymbol, Dim2); rans_N_eit_G2 = REAL(ans_N_eit_G2);

                            } else {


                                    rans_N_eit_G1 = REAL(VECTOR_ELT(out_N_eit_G1, e));
                                    rans_N_eit_G2 = REAL(VECTOR_ELT(out_N_eit_G2, e));
                            }

                            r_N_ei0_G1 = REAL(v_N_ei0_G1);
                            r_N_ei0_G2 = REAL(v_N_ei0_G2);


                                    //facteurs des indices
                                    PROTECT(dFACT4 = iDim(dimC2));
                                    PROTECT(dFACT5 = iDim(dim_N_ei0));
                                    PROTECT(dFACT6 = iDim(dim_N_e0t));

                                    fact4_D = INTEGER(dFACT4);
                                    //int *fact5_D = INTEGER(dFACT5);
                                    //int *fact6_D = INTEGER(dFACT6);

                                    //�quation

                            if (ind_t==0) {

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                 rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_N_ei0_G1[ind_i];
                                 rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_N_ei0_G2[ind_i];
                                }

                            }


                            if (ind_t==0) {


                             setAttrib(ans_N_eit_G1, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_G1, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_G2, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_G2, install("DimCst"), dimCst2);

                             SET_VECTOR_ELT(out_N_eit_G1, e, ans_N_eit_G1);
                             SET_VECTOR_ELT(out_N_eit_G2, e, ans_N_eit_G2);

                            }

// on peut d�sormais �valuer F, Z et N au niveau annuel et global --> pertinent?
//            for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++)
//            for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++)
//            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
//
//                r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =
//
//                            (r_F_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_F_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_F_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_F_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_F_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));
//
//
//                r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =
//
//                            (r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +
//
//                            (r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                            r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
//                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
//                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
//                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));
//
//            }




//for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
//
//     r_Froth_i[ind_i+ind_t*nbI] =
//      (r_Froth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
//      r_Froth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
//      r_Froth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
//      r_Froth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +
//
//      (r_Froth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
//      r_Froth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
//      r_Froth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
//      r_Froth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +
//
//      (r_Froth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
//      r_Froth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
//      r_Froth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
//      r_Froth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +
//
//      (r_Froth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
//      r_Froth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
//      r_Froth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
//      r_Froth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));
//
//
//    r_Foth_i[ind_i+ind_t*nbI] =
//      (r_Foth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
//      r_Foth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
//      r_Foth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
//      r_Foth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +
//
//      (r_Foth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
//      r_Foth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
//      r_Foth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
//      r_Foth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +
//
//      (r_Foth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
//      r_Foth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
//      r_Foth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
//      r_Foth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +
//
//      (r_Foth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
//      r_Foth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
//      r_Foth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
//      r_Foth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
//       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
//           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));
//
//}



                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double tempG1 = 0.0;
                                        double tempG2 = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            tempG1 = tempG1 +  r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        if (!ISNA(r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            tempG2 = tempG2 +  r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        }


                                        if (ISNA(r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            tempG1 + r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_G1[ind_i + nbI*ind_t];

                                        }

                                        if (ISNA(r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            tempG2 + r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_G2[ind_i + nbI*ind_t];

                                        }


                                    //on en profite pour calculer Fbar


                                    fmax = fmax + (tempG1 + r_Froth_i_G1[ind_i + nbI*ind_t]) * r_Fbar_G1[ind_i] * rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] /
                                                    (rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] + rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] ) +
                                                  (tempG2 + r_Froth_i_G2[ind_i + nbI*ind_t]) * r_Fbar_G2[ind_i] * rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] /
                                                    (rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] + rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] );
                                    sumWt = sumWt + (r_Fbar_G1[ind_i]  + r_Fbar_G2[ind_i])/2 ;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; F_i_G1 = " << tempG1 + r_Froth_i_G1[ind_i + nbI*ind_t] << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; F_i_G2 = " << tempG2 + r_Froth_i_G2[ind_i + nbI*ind_t] << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; N_i_G1 = " << rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; N_i_G2 = " << rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]<< endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; Ntot_i = " << rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]+rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]<< endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; Fbar_G1 = " << r_Fbar_G1[ind_i] << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; Fbar_G2 = " << r_Fbar_G2[ind_i] << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; fmax = " << fmax << endl;
                                    //fichier << "e = " << CHAR(STRING_ELT(sppList,e)) <<"; Age = " << ind_i << "; sumW = " << sumWt << endl;

                                    }




       if (ind_t==0) {
         setAttrib(ans_Z_eit_G1, R_DimNamesSymbol, dimnames1);
         setAttrib(ans_Z_eit_G1, install("DimCst"), dimCst1);
         SET_VECTOR_ELT(out_Z_eit_G1, e, ans_Z_eit_G1);

         setAttrib(ans_Z_eit_G2, R_DimNamesSymbol, dimnames1);
         setAttrib(ans_Z_eit_G2, install("DimCst"), dimCst1);
         SET_VECTOR_ELT(out_Z_eit_G2, e, ans_Z_eit_G2);

         setAttrib(ans_N_eit_G1, R_DimNamesSymbol, dimnames2);
         setAttrib(ans_N_eit_G1, install("DimCst"), dimCst2);
         SET_VECTOR_ELT(out_N_eit_G1, e, ans_N_eit_G1);

         setAttrib(ans_N_eit_G2, R_DimNamesSymbol, dimnames2);
         setAttrib(ans_N_eit_G2, install("DimCst"), dimCst2);
         SET_VECTOR_ELT(out_N_eit_G2, e, ans_N_eit_G2);
       }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 9, dimCst_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 238, v_Fr_efmit_G1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 239, v_Fr_efmit_G2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 236, v_M_ei_G1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 237, v_M_ei_G2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 12, dFACT1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 13, dFACT2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 14, dFACT3);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 234, v_Fbar_G1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 235, v_Fbar_G2);

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 17, dFACT4);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 18, dFACT5);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 19, dFACT6);

//Rprintf("End_Sex-based\n");fichier << "End_Sex-based\n" << endl;

    } else if ((Qvec[e]==1) & (Svec[e]==0)) {
//Rprintf("Start_Quarterly\n");fichier << "Start_Quarterly\n" << endl;

//Rprintf("G1.14.1");fichier << "G1.14.1" << endl;
                                    r_Fr_efmit = REAL(v_Fr_efmit);
                                    r_F_efmit = REAL(v_F_efmit);//Rprintf("G1.14.2");
//Rprintf("G1.15");fichier << "G1.15" << endl;
                                    r_Fr_efmit_S1M1 = REAL(v_Fr_efmit_S1M1);  r_Fr_efmit_S1M2 = REAL(v_Fr_efmit_S1M2);
                                    r_Fr_efmit_S1M3 = REAL(v_Fr_efmit_S1M3);  r_Fr_efmit_S1M4 = REAL(v_Fr_efmit_S1M4);
                                    r_Fr_efmit_S2M1 = REAL(v_Fr_efmit_S2M1);  r_Fr_efmit_S2M2 = REAL(v_Fr_efmit_S2M2);
                                    r_Fr_efmit_S2M3 = REAL(v_Fr_efmit_S2M3);  r_Fr_efmit_S2M4 = REAL(v_Fr_efmit_S2M4);
                                    r_Fr_efmit_S3M1 = REAL(v_Fr_efmit_S3M1);  r_Fr_efmit_S3M2 = REAL(v_Fr_efmit_S3M2);
                                    r_Fr_efmit_S3M3 = REAL(v_Fr_efmit_S3M3);  r_Fr_efmit_S3M4 = REAL(v_Fr_efmit_S3M4);
                                    r_Fr_efmit_S4M1 = REAL(v_Fr_efmit_S4M1);  r_Fr_efmit_S4M2 = REAL(v_Fr_efmit_S4M2);
                                    r_Fr_efmit_S4M3 = REAL(v_Fr_efmit_S4M3);  r_Fr_efmit_S4M4 = REAL(v_Fr_efmit_S4M4);
//Rprintf("G1.16");fichier << "G1.16" << endl;
                                    r_F_efmit_S1M1 = REAL(v_F_efmit_S1M1);  r_F_efmit_S1M2 = REAL(v_F_efmit_S1M2);
                                    r_F_efmit_S1M3 = REAL(v_F_efmit_S1M3);  r_F_efmit_S1M4 = REAL(v_F_efmit_S1M4);
                                    r_F_efmit_S2M1 = REAL(v_F_efmit_S2M1);  r_F_efmit_S2M2 = REAL(v_F_efmit_S2M2);
                                    r_F_efmit_S2M3 = REAL(v_F_efmit_S2M3);  r_F_efmit_S2M4 = REAL(v_F_efmit_S2M4);
                                    r_F_efmit_S3M1 = REAL(v_F_efmit_S3M1);  r_F_efmit_S3M2 = REAL(v_F_efmit_S3M2);
                                    r_F_efmit_S3M3 = REAL(v_F_efmit_S3M3);  r_F_efmit_S3M4 = REAL(v_F_efmit_S3M4);
                                    r_F_efmit_S4M1 = REAL(v_F_efmit_S4M1);  r_F_efmit_S4M2 = REAL(v_F_efmit_S4M2);
                                    r_F_efmit_S4M3 = REAL(v_F_efmit_S4M3);  r_F_efmit_S4M4 = REAL(v_F_efmit_S4M4);
//Rprintf("G1.17");fichier << "G1.17" << endl;
                                    r_M_ei = REAL(v_M_ei);
                                    r_Fbar = REAL(v_Fbar);
//Rprintf("G1.18");fichier << "G1.18" << endl;
                                    //facteurs des indices
                                    PROTECT(dFACT1 = iDim(dimC1));
                                    PROTECT(dFACT2 = iDim(dim_Fr_efmit));
                                    PROTECT(dFACT3 = iDim(dim_M_ei));

                                    int *fact1_D = INTEGER(dFACT1);
                                    int *fact2_D = INTEGER(dFACT2);
                                    int *fact3_D = INTEGER(dFACT3);

                                    double *r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44)); //Rprintf("Dans EVAR (l.6145), Fothi = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                                    double *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));
//Rprintf("G1.19");fichier << "G1.19" << endl;
                                    double *r_Foth_i_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 116));
                                    double *r_Foth_i_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 117));
                                    double *r_Foth_i_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 118));
                                    double *r_Foth_i_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 119));
                                    double *r_Foth_i_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 120));
                                    double *r_Foth_i_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 121));
                                    double *r_Foth_i_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 122));
                                    double *r_Foth_i_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 123));
                                    double *r_Foth_i_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 124));
                                    double *r_Foth_i_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 125));
                                    double *r_Foth_i_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 126));
                                    double *r_Foth_i_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 127));
                                    double *r_Foth_i_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 128));
                                    double *r_Foth_i_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 129));
                                    double *r_Foth_i_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 130));
                                    double *r_Foth_i_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 131));
//Rprintf("G1.20");fichier << "G1.20" << endl;
                                    double *r_Froth_i_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 132));
                                    double *r_Froth_i_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 133));
                                    double *r_Froth_i_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 134));
                                    double *r_Froth_i_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 135));
                                    double *r_Froth_i_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 136));
                                    double *r_Froth_i_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 137));
                                    double *r_Froth_i_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 138));
                                    double *r_Froth_i_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 139));
                                    double *r_Froth_i_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 140));
                                    double *r_Froth_i_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 141));
                                    double *r_Froth_i_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 142));
                                    double *r_Froth_i_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 143));
                                    double *r_Froth_i_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 144));
                                    double *r_Froth_i_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 145));
                                    double *r_Froth_i_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 146));
                                    double *r_Froth_i_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 147));
//Rprintf("G1.21");fichier << "G1.21" << endl;
                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M1[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S1M1 -> 0

                                        //if (e==eTemp) Ztemp[ind_i+1+(0*nbI)] = rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(0*nbI)] = rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M2[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S1M2 -> 1

                                        //if (e==eTemp) Ztemp[ind_i+1+(1*nbI)] = rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(1*nbI)] = rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M3[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S1M3 -> 2

                                        //if (e==eTemp) Ztemp[ind_i+1+(2*nbI)] = rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(2*nbI)] = rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M4[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S1M4 -> 3

                                        //if (e==eTemp) Ztemp[ind_i+1+(3*nbI)] = rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(3*nbI)] = rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M1[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S2M1 -> 4

                                        //if (e==eTemp) Ztemp[ind_i+1+(4*nbI)] = rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(4*nbI)] = rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M2[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S2M2 -> 5

                                        //if (e==eTemp) Ztemp[ind_i+1+(5*nbI)] = rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(5*nbI)] = rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M3[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S2M3 -> 6

                                        //if (e==eTemp) Ztemp[ind_i+1+(6*nbI)] = rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(6*nbI)] = rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M4[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S2M4 -> 7

                                        //if (e==eTemp) Ztemp[ind_i+1+(7*nbI)] = rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(7*nbI)] = rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M1[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S3M1 -> 8

                                        //if (e==eTemp) Ztemp[ind_i+1+(8*nbI)] = rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(8*nbI)] = rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M2[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S3M2 -> 9

                                        //if (e==eTemp) Ztemp[ind_i+1+(9*nbI)] = rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(9*nbI)] = rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M3[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S3M3 -> 10

                                        //if (e==eTemp) Ztemp[ind_i+1+(10*nbI)] = rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(10*nbI)] = rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M4[ind_i + nbI*ind_t];
                                          if (ind_i==0) rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S3M4 -> 11

                                        //if (e==eTemp) Ztemp[ind_i+1+(11*nbI)] = rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(11*nbI)] = rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M1[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S4M1 -> 12

                                        //if (e==eTemp) Ztemp[ind_i+1+(12*nbI)] = rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(12*nbI)] = rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M2[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S4M2 -> 13

                                        //if (e==eTemp) Ztemp[ind_i+1+(13*nbI)] = rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(13*nbI)] = rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M3[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S4M3 -> 14

                                        //if (e==eTemp) Ztemp[ind_i+1+(14*nbI)] = rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(14*nbI)] = rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M4[ind_i + nbI*ind_t];
                                        }

                                        //on initialise aussi Ztemp (attention : index� � partir de 1) : S4M4 -> 15

                                        //if (e==eTemp) Ztemp[ind_i+1+(15*nbI)] = rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        Ztemp[ind_i+(15*nbI)] = rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];

                                    }


                              if (ind_t==0) {

                                    setAttrib(ans_Z_eit_S1M1, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S1M2, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S1M3, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S1M4, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S2M1, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S2M2, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S2M3, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S2M4, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S3M1, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S3M2, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S3M3, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S3M4, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S4M1, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S4M2, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit_S4M3, R_DimNamesSymbol, dimnames1); setAttrib(ans_Z_eit_S4M4, R_DimNamesSymbol, dimnames1);

                                    setAttrib(ans_Z_eit_S1M1, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S1M2, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S1M3, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S1M4, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S2M1, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S2M2, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S2M3, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S2M4, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S3M1, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S3M2, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S3M3, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S3M4, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S4M1, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S4M2, install("DimCst"), dimCst1);
                                    setAttrib(ans_Z_eit_S4M3, install("DimCst"), dimCst1); setAttrib(ans_Z_eit_S4M4, install("DimCst"), dimCst1);

                                    SET_VECTOR_ELT(out_Z_eit_S1M1, e, ans_Z_eit_S1M1); SET_VECTOR_ELT(out_Z_eit_S1M2, e, ans_Z_eit_S1M2);
                                    SET_VECTOR_ELT(out_Z_eit_S1M3, e, ans_Z_eit_S1M3); SET_VECTOR_ELT(out_Z_eit_S1M4, e, ans_Z_eit_S1M4);
                                    SET_VECTOR_ELT(out_Z_eit_S2M1, e, ans_Z_eit_S2M1); SET_VECTOR_ELT(out_Z_eit_S2M2, e, ans_Z_eit_S2M2);
                                    SET_VECTOR_ELT(out_Z_eit_S2M3, e, ans_Z_eit_S2M3); SET_VECTOR_ELT(out_Z_eit_S2M4, e, ans_Z_eit_S2M4);
                                    SET_VECTOR_ELT(out_Z_eit_S3M1, e, ans_Z_eit_S3M1); SET_VECTOR_ELT(out_Z_eit_S3M2, e, ans_Z_eit_S3M2);
                                    SET_VECTOR_ELT(out_Z_eit_S3M3, e, ans_Z_eit_S3M3); SET_VECTOR_ELT(out_Z_eit_S3M4, e, ans_Z_eit_S3M4);
                                    SET_VECTOR_ELT(out_Z_eit_S4M1, e, ans_Z_eit_S4M1); SET_VECTOR_ELT(out_Z_eit_S4M2, e, ans_Z_eit_S4M2);
                                    SET_VECTOR_ELT(out_Z_eit_S4M3, e, ans_Z_eit_S4M3); SET_VECTOR_ELT(out_Z_eit_S4M4, e, ans_Z_eit_S4M4);
                                    ////Rprintf("AAAhhh!!!");//PrintValue(out_Z_eit_S1M1);
                              }






                                    //---------
                                    // calcul de N_eit
                                    //---------

                                    //on d�termine l'attribut Dimension de N_eit
                                    PROTECT(dimCst2 = allocVector(INTSXP, 4));
                                    dimC2 = INTEGER(dimCst2);

                                    dimC2[0] = 0 ; dimC2[1] = 0 ; dimC2[2] = imax2(dimC1[2] , dim_N_ei0[2]) ; dimC2[3] = nbT;

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            count++;
                                            prod = prod * dimC2[k];
                                        }
                                    }

                                    PROTECT(Dim2 = allocVector(INTSXP, count));
                                    dim2 = INTEGER(Dim2);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            dim2[count2] = dimC2[k];
                                            count2++;
                                            }
                                    }


                            if (ind_t==0) {

                                    //on cr�e le tableau r�sultat pour l'esp�ce en question
                                    PROTECT(ans_N_eit = NEW_NUMERIC(prod));
                                    setAttrib(ans_N_eit, R_DimSymbol, Dim2);

                                    PROTECT(dimnames2 = allocVector(VECSXP,count));
                                    if (dimC2[0]>0) {SET_VECTOR_ELT(dimnames2, count3, fleetList) ; count3++;}
                                    if (dimC2[1]>0) {SET_VECTOR_ELT(dimnames2, count3, metierList) ; count3++;}
                                    if (dimC2[2]>0) {SET_VECTOR_ELT(dimnames2, count3, intAge) ; count3++;}
                                    if (dimC2[3]>0) {SET_VECTOR_ELT(dimnames2, count3, times) ; count3++;}

                                    rans_N_eit = REAL(ans_N_eit);

                                    PROTECT(ans_N_eit_S1M1 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S1M1, R_DimSymbol, Dim2); rans_N_eit_S1M1 = REAL(ans_N_eit_S1M1);
                                    PROTECT(ans_N_eit_S1M2 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S1M2, R_DimSymbol, Dim2); rans_N_eit_S1M2 = REAL(ans_N_eit_S1M2);
                                    PROTECT(ans_N_eit_S1M3 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S1M3, R_DimSymbol, Dim2); rans_N_eit_S1M3 = REAL(ans_N_eit_S1M3);
                                    PROTECT(ans_N_eit_S1M4 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S1M4, R_DimSymbol, Dim2); rans_N_eit_S1M4 = REAL(ans_N_eit_S1M4);
                                    PROTECT(ans_N_eit_S2M1 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S2M1, R_DimSymbol, Dim2); rans_N_eit_S2M1 = REAL(ans_N_eit_S2M1);
                                    PROTECT(ans_N_eit_S2M2 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S2M2, R_DimSymbol, Dim2); rans_N_eit_S2M2 = REAL(ans_N_eit_S2M2);
                                    PROTECT(ans_N_eit_S2M3 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S2M3, R_DimSymbol, Dim2); rans_N_eit_S2M3 = REAL(ans_N_eit_S2M3);
                                    PROTECT(ans_N_eit_S2M4 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S2M4, R_DimSymbol, Dim2); rans_N_eit_S2M4 = REAL(ans_N_eit_S2M4);
                                    PROTECT(ans_N_eit_S3M1 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S3M1, R_DimSymbol, Dim2); rans_N_eit_S3M1 = REAL(ans_N_eit_S3M1);
                                    PROTECT(ans_N_eit_S3M2 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S3M2, R_DimSymbol, Dim2); rans_N_eit_S3M2 = REAL(ans_N_eit_S3M2);
                                    PROTECT(ans_N_eit_S3M3 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S3M3, R_DimSymbol, Dim2); rans_N_eit_S3M3 = REAL(ans_N_eit_S3M3);
                                    PROTECT(ans_N_eit_S3M4 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S3M4, R_DimSymbol, Dim2); rans_N_eit_S3M4 = REAL(ans_N_eit_S3M4);
                                    PROTECT(ans_N_eit_S4M1 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S4M1, R_DimSymbol, Dim2); rans_N_eit_S4M1 = REAL(ans_N_eit_S4M1);
                                    PROTECT(ans_N_eit_S4M2 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S4M2, R_DimSymbol, Dim2); rans_N_eit_S4M2 = REAL(ans_N_eit_S4M2);
                                    PROTECT(ans_N_eit_S4M3 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S4M3, R_DimSymbol, Dim2); rans_N_eit_S4M3 = REAL(ans_N_eit_S4M3);
                                    PROTECT(ans_N_eit_S4M4 = NEW_NUMERIC(prod)); setAttrib(ans_N_eit_S4M4, R_DimSymbol, Dim2); rans_N_eit_S4M4 = REAL(ans_N_eit_S4M4);

                            } else {

                                    rans_N_eit = REAL(VECTOR_ELT(out_N_eit, e));

                                    rans_N_eit_S1M1 = REAL(VECTOR_ELT(out_N_eit_S1M1, e)); rans_N_eit_S1M2 = REAL(VECTOR_ELT(out_N_eit_S1M2, e));
                                    rans_N_eit_S1M3 = REAL(VECTOR_ELT(out_N_eit_S1M3, e)); rans_N_eit_S1M4 = REAL(VECTOR_ELT(out_N_eit_S1M4, e));
                                    rans_N_eit_S2M1 = REAL(VECTOR_ELT(out_N_eit_S2M1, e)); rans_N_eit_S2M2 = REAL(VECTOR_ELT(out_N_eit_S2M2, e));
                                    rans_N_eit_S2M3 = REAL(VECTOR_ELT(out_N_eit_S2M3, e)); rans_N_eit_S2M4 = REAL(VECTOR_ELT(out_N_eit_S2M4, e));
                                    rans_N_eit_S3M1 = REAL(VECTOR_ELT(out_N_eit_S3M1, e)); rans_N_eit_S3M2 = REAL(VECTOR_ELT(out_N_eit_S3M2, e));
                                    rans_N_eit_S3M3 = REAL(VECTOR_ELT(out_N_eit_S3M3, e)); rans_N_eit_S3M4 = REAL(VECTOR_ELT(out_N_eit_S3M4, e));
                                    rans_N_eit_S4M1 = REAL(VECTOR_ELT(out_N_eit_S4M1, e)); rans_N_eit_S4M2 = REAL(VECTOR_ELT(out_N_eit_S4M2, e));
                                    rans_N_eit_S4M3 = REAL(VECTOR_ELT(out_N_eit_S4M3, e)); rans_N_eit_S4M4 = REAL(VECTOR_ELT(out_N_eit_S4M4, e));

                            }

//   r_N_ei0_S1M1 = REAL(v_N_ei0_S1M1); r_N_ei0_S1M2 = REAL(v_N_ei0_S1M2); r_N_ei0_S1M3 = REAL(v_N_ei0_S1M3); r_N_ei0_S1M4 = REAL(v_N_ei0_S1M4);
//   r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1); r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2); r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3); r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);
   r_iniNt0q_S1M1 = REAL(v_iniNt0q_S1M1); r_iniNt0q_S1M2 = REAL(v_iniNt0q_S1M2); r_iniNt0q_S1M3 = REAL(v_iniNt0q_S1M3); r_iniNt0q_S1M4 = REAL(v_iniNt0q_S1M4);
   r_iniNt0q_S2M1 = REAL(v_iniNt0q_S2M1); r_iniNt0q_S2M2 = REAL(v_iniNt0q_S2M2); r_iniNt0q_S2M3 = REAL(v_iniNt0q_S2M3); r_iniNt0q_S2M4 = REAL(v_iniNt0q_S2M4);
   r_iniNt0q_S3M1 = REAL(v_iniNt0q_S3M1); r_iniNt0q_S3M2 = REAL(v_iniNt0q_S3M2); r_iniNt0q_S3M3 = REAL(v_iniNt0q_S3M3); r_iniNt0q_S3M4 = REAL(v_iniNt0q_S3M4);
   r_iniNt0q_S4M1 = REAL(v_iniNt0q_S4M1); r_iniNt0q_S4M2 = REAL(v_iniNt0q_S4M2); r_iniNt0q_S4M3 = REAL(v_iniNt0q_S4M3); r_iniNt0q_S4M4 = REAL(v_iniNt0q_S4M4);


                                    //facteurs des indices
                                    PROTECT(dFACT4 = iDim(dimC2));
                                    PROTECT(dFACT5 = iDim(dim_N_ei0));
                                    PROTECT(dFACT6 = iDim(dim_N_e0t));

                                    fact4_D = INTEGER(dFACT4);
                                    //int *fact5_D = INTEGER(dFACT5);
                                    //int *fact6_D = INTEGER(dFACT6);

                                    //�quation

                            if (ind_t==0) {

                                        //S1
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                 rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S1M1[ind_i];
                                 rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S1M2[ind_i];
                                 rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S1M3[ind_i];
                                 rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S1M4[ind_i];
                                 rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S2M1[ind_i];
                                 rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S2M2[ind_i];
                                 rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S2M3[ind_i];
                                 rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S2M4[ind_i];
                                 rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S3M1[ind_i];
                                 rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S3M2[ind_i];
                                 rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S3M3[ind_i];
                                 rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S3M4[ind_i];
                                 rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S4M1[ind_i];
                                 rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S4M2[ind_i];
                                 rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S4M3[ind_i];
                                 rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = r_iniNt0q_S4M4[ind_i];

                                }

                            }


                            if (ind_t==0) {

                             setAttrib(ans_N_eit, R_DimNamesSymbol, dimnames2);
                             setAttrib(ans_N_eit, install("DimCst"), dimCst2);

                             setAttrib(ans_N_eit_S1M1, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S1M1, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S1M2, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S1M2, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S1M3, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S1M3, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S1M4, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S1M4, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S2M1, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S2M1, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S2M2, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S2M2, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S2M3, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S2M3, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S2M4, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S2M4, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S3M1, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S3M1, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S3M2, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S3M2, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S3M3, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S3M3, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S3M4, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S3M4, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S4M1, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S4M1, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S4M2, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S4M2, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S4M3, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S4M3, install("DimCst"), dimCst2);
                             setAttrib(ans_N_eit_S4M4, R_DimNamesSymbol, dimnames2); setAttrib(ans_N_eit_S4M4, install("DimCst"), dimCst2);

                             SET_VECTOR_ELT(out_N_eit, e, ans_N_eit);

                             SET_VECTOR_ELT(out_N_eit_S1M1, e, ans_N_eit_S1M1); SET_VECTOR_ELT(out_N_eit_S1M2, e, ans_N_eit_S1M2);
                             SET_VECTOR_ELT(out_N_eit_S1M3, e, ans_N_eit_S1M3); SET_VECTOR_ELT(out_N_eit_S1M4, e, ans_N_eit_S1M4);
                             SET_VECTOR_ELT(out_N_eit_S2M1, e, ans_N_eit_S2M1); SET_VECTOR_ELT(out_N_eit_S2M2, e, ans_N_eit_S2M2);
                             SET_VECTOR_ELT(out_N_eit_S2M3, e, ans_N_eit_S2M3); SET_VECTOR_ELT(out_N_eit_S2M4, e, ans_N_eit_S2M4);
                             SET_VECTOR_ELT(out_N_eit_S3M1, e, ans_N_eit_S3M1); SET_VECTOR_ELT(out_N_eit_S3M2, e, ans_N_eit_S3M2);
                             SET_VECTOR_ELT(out_N_eit_S3M3, e, ans_N_eit_S3M3); SET_VECTOR_ELT(out_N_eit_S3M4, e, ans_N_eit_S3M4);
                             SET_VECTOR_ELT(out_N_eit_S4M1, e, ans_N_eit_S4M1); SET_VECTOR_ELT(out_N_eit_S4M2, e, ans_N_eit_S4M2);
                             SET_VECTOR_ELT(out_N_eit_S4M3, e, ans_N_eit_S4M3); SET_VECTOR_ELT(out_N_eit_S4M4, e, ans_N_eit_S4M4);

                            }

// on peut d�sormais �valuer F, Z et N au niveau annuel et global
            for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++)
            for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++)
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =

                            (r_F_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));


                r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =

                            (r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));

            }




for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

     r_Froth_i[ind_i+ind_t*nbI] =
      (r_Froth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
      r_Froth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
      r_Froth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
      r_Froth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
      r_Froth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
      r_Froth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
      r_Froth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
      r_Froth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
      r_Froth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
      r_Froth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
      r_Froth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
      r_Froth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
      r_Froth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));


    r_Foth_i[ind_i+ind_t*nbI] =
      (r_Foth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
      r_Foth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
      r_Foth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
      r_Foth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
      r_Foth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
      r_Foth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
      r_Foth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
      r_Foth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
      r_Foth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
      r_Foth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
      r_Foth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
      r_Foth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
      r_Foth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));

}



                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        }


                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i[ind_i + nbI*ind_t];

                                        }


                                    //on en profite pour calculer Fbar

                                    fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                                    sumWt = sumWt + r_Fbar[ind_i];

                                    // et remplir N_eit (effectifs � la saison 1)
                                    rans_N_eit[ind_i + nbI*ind_t] = rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
                                         rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI];

                                    }




       if (ind_t==0) {
         setAttrib(ans_Z_eit, R_DimNamesSymbol, dimnames1);
         setAttrib(ans_Z_eit, install("DimCst"), dimCst1);

         SET_VECTOR_ELT(out_Z_eit, e, ans_Z_eit);

         setAttrib(ans_N_eit, R_DimNamesSymbol, dimnames2);
         setAttrib(ans_N_eit, install("DimCst"), dimCst2);

         SET_VECTOR_ELT(out_N_eit, e, ans_N_eit);
       }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 9, dimCst_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 10, v_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 11, v_M_ei);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 12, dFACT1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 13, dFACT2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 14, dFACT3);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 59, v_Fbar);

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 17, dFACT4);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 18, dFACT5);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 19, dFACT6);

//Rprintf("End_Quarterly\n");fichier << "End_Quarterly\n" << endl;
    }

//--------------------------------------------------------------------------------------



                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de SSB_et
                                    //---------
//Rprintf("Calcul_SSB\n");fichier << "Calcul_SSB" << endl;
                                    //on d�termine l'attribut Dimension de SSB_et
                                    PROTECT(dimCst4 = allocVector(INTSXP, 4));
                                    dimC4 = INTEGER(dimCst4);

                                    dimC4[0] = 0 ; dimC4[1] = 0 ; dimC4[2] = 0 ; dimC4[3] = dimC2[3];

                                    count = 0 ; prod = 1 ; count2 = 0 ; count4 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC4[k]>0) {
                                            count++;
                                            prod = prod * dimC4[k];
                                        }
                                    }


                                    PROTECT(Dim4 = allocVector(INTSXP, count));
                                    dim4 = INTEGER(Dim4);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC4[k]>0) {
                                            dim4[count2] = dimC4[k];
                                            count2++;
                                            }
                                    }

                            if (ind_t==0) {
                                    //on cr�e le tableau r�sultat pour l'esp�ce en question (on en profite pour faire de m�me avec Fbar --> m�me dimension)
                                    PROTECT(ans_SSB_et = NEW_NUMERIC(prod));
                                    PROTECT(ans_Fbar_et = NEW_NUMERIC(prod));

                                    if (count>0) { //valable seulement si SSB_et n'est pas seulement un scalaire

                                    setAttrib(ans_SSB_et, R_DimSymbol, Dim4);
                                    setAttrib(ans_Fbar_et, R_DimSymbol, Dim4);

                                    }

                                    PROTECT(dimnames4 = allocVector(VECSXP,count));
                                    if (dimC4[0]>0) {SET_VECTOR_ELT(dimnames4, count4, fleetList) ; count4++;}
                                    if (dimC4[1]>0) {SET_VECTOR_ELT(dimnames4, count4, metierList) ; count4++;}
                                    if (dimC4[2]>0) {SET_VECTOR_ELT(dimnames4, count4, intAge) ; count4++;}
                                    if (dimC4[3]>0) {SET_VECTOR_ELT(dimnames4, count4, times) ; count4++;}


                                    rans_SSB_et = REAL(ans_SSB_et);
                                    rans_Fbar_et = REAL(ans_Fbar_et);

                            } else {


                                    rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et, e));
                                    rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et, e));

                            }



                                    //facteurs des indices
                                    PROTECT(dFACT9 = iDim(dimC4));
                                    PROTECT(dFACT8 = iDim(dim_w_ei));
                                    PROTECT(dFACT10 = iDim(dim_mat_ei));

                                    int *fact9_D = INTEGER(dFACT9);
                                    int *fact8_D = INTEGER(dFACT8);
                                    int *fact10_D = INTEGER(dFACT10);

                                    //�quation
                        if ((Qvec[e]==0) & (Svec[e]==0))  {
                                    r_w_ei = REAL(v_w_ei);
                                    r_mat_ei = REAL(v_mat_ei);

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                                temp = temp +
                                                 rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                                 r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                            rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                            rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;

                                        }
                        } else if ((Qvec[e]==0) & (Svec[e]==1))  {
                                    r_w_ei_G1 = REAL(v_w_ei_G1);
                                    r_w_ei_G2 = REAL(v_w_ei_G2);
                                    r_mat_ei_G1 = REAL(v_mat_ei_G1);
                                    r_mat_ei_G2 = REAL(v_mat_ei_G2);

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                                temp = temp +
                                                 rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_mat_ei_G1[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                                 r_w_ei_G1[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000 +

                                                 rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_mat_ei_G2[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                                 r_w_ei_G2[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                            rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                            rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;

                                        }
                        } else if ((Qvec[e]==1) & (Svec[e]==0)){

                        r_matwt_M1 = REAL(v_matwt_M1) ; r_matwt_M2 = REAL(v_matwt_M2) ; r_matwt_M3 = REAL(v_matwt_M3) ; r_matwt_M4 = REAL(v_matwt_M4) ;
                        double tempSSB = 0.0;
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                            tempSSB = tempSSB + r_matwt_M1[ind_i]*rans_N_eit_S1M1[ind_i+ind_t*nbI] + r_matwt_M2[ind_i]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
                                  r_matwt_M3[ind_i]*rans_N_eit_S1M3[ind_i+ind_t*nbI] + r_matwt_M4[ind_i]*rans_N_eit_S1M4[ind_i+ind_t*nbI];

                            rans_SSB_et[ind_t] = tempSSB ;
                            rans_Fbar_et[ind_t] = fmax/sumWt;

                        }

                            if(ind_t==0) {

                                    if (count>0) setAttrib(ans_SSB_et, R_DimNamesSymbol, dimnames4);
                                    setAttrib(ans_SSB_et, install("DimCst"), dimCst4);

                                    SET_VECTOR_ELT(out_SSB_et, e, ans_SSB_et);

                                    if (count>0) setAttrib(ans_Fbar_et, R_DimNamesSymbol, dimnames4);
                                    setAttrib(ans_Fbar_et, install("DimCst"), dimCst4);

                                    SET_VECTOR_ELT(out_Fbar_et, e, ans_Fbar_et);

                                    SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 48, dFACT9);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 49, dFACT10);
                                if ((Qvec[e]==0) & (Svec[e]==0)) SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 47, v_mat_ei);




                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de B_et
                                    //---------
//Rprintf("Calcul_B\n");fichier << "Calcul_B" << endl;
//Rprintf("K700\n");fichier << "K700" << endl;
                                    //on d�termine l'attribut Dimension de B_et
                                    PROTECT(dimCst3 = allocVector(INTSXP, 4));
                                    dimC3 = INTEGER(dimCst3);

                                    dimC3[0] = 0 ; dimC3[1] = 0 ; dimC3[2] = 0 ; dimC3[3] = dimC2[3];

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC3[k]>0) {
                                            count++;
                                            prod = prod * dimC3[k];
                                        }
                                    }

                                    PROTECT(Dim3 = allocVector(INTSXP, count));
                                    dim3 = INTEGER(Dim3);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC3[k]>0) {
                                            dim3[count2] = dimC3[k];
                                            count2++;
                                            }
                                    }
//Rprintf("K70\n");fichier << "K70" << endl;
                            if (ind_t==0) {
                                    //on cr�e le tableau r�sultat pour l'esp�ce en question
                                    PROTECT(ans_B_et = NEW_NUMERIC(prod));

                                    if (count>0) { //valable seulement si B_et n'est pas seulement un scalaire

                                    setAttrib(ans_B_et, R_DimSymbol, Dim3);

                                    }

                                    PROTECT(dimnames3 = allocVector(VECSXP,count));
                                    if (dimC3[0]>0) {SET_VECTOR_ELT(dimnames3, count3, fleetList) ; count3++;}
                                    if (dimC3[1]>0) {SET_VECTOR_ELT(dimnames3, count3, metierList) ; count3++;}
                                    if (dimC3[2]>0) {SET_VECTOR_ELT(dimnames3, count3, intAge) ; count3++;}
                                    if (dimC3[3]>0) {SET_VECTOR_ELT(dimnames3, count3, times) ; count3++;}

                                    rans_B_et = REAL(ans_B_et);

                            } else {

                                    rans_B_et = REAL(VECTOR_ELT(out_B_et, e));

                            }

                                    //facteurs des indices
                                    PROTECT(dFACT7 = iDim(dimC3));

                                    int *fact7_D = INTEGER(dFACT7);

                            if (Svec[e]==0){
                                    r_r = REAL(v_r);
                                    r_K = REAL(v_K);
                                    r_n = REAL(v_n);
                                    r_B = REAL(v_B);
//Rprintf("K71\n");fichier << "K71" << endl;
                                    //�quation

                                    if (nbI>1) {

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                                temp = temp +
                                                 rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                            rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                        }

                                    } else {
                                        if ((!isNull(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e)))))){ // si aleatoire sur les parametres K et r
                                            delta_r = REAL(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e))))[ind_t];
                                            delta_K = REAL(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e))))[ind_t+nbT];
                                        }
                                        //Rprintf("t = ",ind_t,"; sp = ",CHAR(STRING_ELT(sppList,e)), "delta_r = ",  delta_r,"\n");
                                        //Rprintf("t = ",ind_t,"; sp = ",CHAR(STRING_ELT(sppList,e)), "delta_K = ",  delta_K,"\n");

                                        rans_B_et[ind_t] = r_B[0]; //biomasse initiale SPiCT

                                        SEXP ans_interm = R_NilValue ;
                                        PROTECT(ans_interm = NEW_NUMERIC(1*nbT)); //PROTECT(ans_interm = NEW_NUMERIC(16*nbT));
                                        double *rans_interm = REAL(ans_interm); for (int yy = 0 ; yy < 1*nbT ; yy++) rans_interm[yy] = NA_REAL; //for (int yy = 0 ; yy < 16*nbT ; yy++) rans_interm[yy] = NA_REAL; //initialisation
                                        double *r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44)); //Rprintf("Dans EVAR (l.7217), Fothi = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                                        double *r_F_efmit = REAL(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                                        double temp = 0.0;
                                        int *dim_F_efmit;
                                        SEXP dimCst_F_efmit, cFACT2;
                                        PROTECT(dimCst_F_efmit = getAttrib(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))), install("DimCst")));
                                        dim_F_efmit = INTEGER(dimCst_F_efmit);
                                        PROTECT(cFACT2 = iDim(dim_F_efmit));
                                        int *fact2_C = INTEGER(cFACT2);
                                        for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {
                                            if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]]))
                                            temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]];
                                        }
                                        //il faut aussi remplir les biomasses par 16�me de temps
                                        rans_interm[0] = r_B[0];
                                        // Spict equation
                                        //rans_interm[1] = rans_interm[0] + (r_r[0]/(r_n[0]-1))*rans_interm[0]*(1-pow(rans_interm[0]/r_K[0],r_n[0]-1)) -
                                         //                                        (temp + r_Foth_i[0 + ind_t*1])*rans_interm[0];

                                        // Fox equation
                                        rans_interm[1] = rans_interm[0] * (1 + (r_r[0]+delta_r) * log ( (r_K[0]+delta_K) / rans_interm[0]) - (temp + r_Foth_i[ind_t])) ;
                                        //Rprintf("ransinterm[0] = %f, ransinterm[1] = %f, r = %f, K = %f, F = %f",rans_interm[0],rans_interm[1], r_r[0], r_K[0], temp + r_Foth_i[ind_t]);
                                        //for (int ib = 1 ; ib < 17 ; ib++) rans_interm[ib] = rans_interm[ib-1] + (r_r[0]/(r_n[0]-1))*rans_interm[ib-1]*(1-pow(rans_interm[ib-1]/r_K[0],r_n[0]-1))/16 -
                                        //                                         (temp + r_Foth_i[0 + ind_t*1])*rans_interm[ib-1]/16;
                                        SET_VECTOR_ELT(intermBIOMspict, e, ans_interm);
                                        UNPROTECT(3);

                                    }
                            } else {

//Rprintf("K71\n");fichier << "K71" << endl;
                                    //�quation

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                                temp = temp +
                                                 (rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_w_ei_G1[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] +
                                                 rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_w_ei_G2[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] )/ 1000;

                                            rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                        }

                            }
//Rprintf("K72\n");fichier << "K72" << endl;

                            if(ind_t==0) {

                                    if (count>0) setAttrib(ans_B_et, R_DimNamesSymbol, dimnames3);
                                    setAttrib(ans_B_et, install("DimCst"), dimCst3);

                                    SET_VECTOR_ELT(out_B_et, e, ans_B_et);
                                    SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));
                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 20, v_w_ei);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 21, dFACT7);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 22, dFACT8);
//Rprintf("K73\n");fichier << "K73" << endl;
    UNPROTECT(12);
    if ((Qvec[e]==1) & (Svec[e]==0)) UNPROTECT(76+8);
    if ((Qvec[e]==0) & (Svec[e]==0)) UNPROTECT(18+8);
    if ((Qvec[e]==0) & (Svec[e]==1)) UNPROTECT(23+8);

    if (ind_t==0){
        UNPROTECT(6);
        if ((Qvec[e]==1) & (Svec[e]==0)) UNPROTECT(17+18);
        if ((Qvec[e]==0) & (Svec[e]==0)) UNPROTECT(1+2);
        if ((Qvec[e]==0) & (Svec[e]==1)) UNPROTECT(2+3);

    }
}
//PrintValue(out_Fbar_et);
dUpdate = false;
UNPROTECT(1);
//fichier.close();

} else {

//fichier << "dUpdate = " << dUpdate << endl;
//Rprintf("dUpdate = %f \n" ,dUpdate);

for (int e = 0 ; e < nbE ; e++) {
//Rprintf("G1");fichier << "G1" << endl;

                    SEXP elmt;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                    int nbI = length(VECTOR_ELT(namDC, e));
                    SEXP v_N_e0t, v_N_e0t_S1M1, v_N_e0t_S2M2, v_N_e0t_S3M3, v_N_e0t_S4M4, v_N_e0t_G1, v_N_e0t_G2,
                         v_N_ei0_S1M1, v_N_ei0_S1M2, v_N_ei0_S1M3, v_N_ei0_S1M4,
                         v_matwt_M1,v_matwt_M2,v_matwt_M3,v_matwt_M4;


                    double  *rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et,e));//Rprintf("G31");
                    double  *rans_B_et = REAL(VECTOR_ELT(out_B_et,e));//Rprintf("G32");
                    double  *rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et,e));//Rprintf("G33");


                    //double  *r_B = REAL(getListElement(elmt, "B_i"));
                    //double  *r_Ytot = REAL(VECTOR_ELT(out_Y_eit,e)); //un seul �ge si SPiCT
                    int *dim_Fr_efmit = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 9)),
                        *fact1_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 12)),
                        *fact2_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 13)),
                        *fact3_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 14)),
                        *fact4_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 17)),
                        *fact5_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 18)),
                        *fact6_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 19)),
                        *fact7_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 21)),
                        *fact8_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 22)),
                        *fact9_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 48)),
                        *fact10_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 49));

                    double fmax = 0.0, sumWt = 0.0;


//Rprintf("G7\n");fichier << "G7" << endl;

    if ((Qvec[e]==0) & (Svec[e]==0)) {
    //Rprintf("Start_Annual\n");fichier << "Start_Annual" << endl;

            PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
            double  *r_Fr_efmit = REAL(VECTOR_ELT(out_Fr_fmi, e));//Rprintf("G34");
            double  *r_F_efmit = REAL(VECTOR_ELT(out_F_fmi, e));
            double  *r_M_ei = REAL(getListElement(elmt, "M_i"));
            double  *r_N_ei0 = REAL(getListElement(elmt, "N_it0"));

            double  *r_N_e0t = REAL(v_N_e0t);
            double  *r_w_ei = REAL(getListElement(elmt, "wStock_i"));////PrintValue(getListElement(elmt, "wStock_i"));
            double  *r_mat_ei = REAL(getListElement(elmt, "mat_i"));

            double  *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));
            double  *r_Fbar = REAL(getListElement(elmt, "Fbar"));

            double  *r_r = REAL(getListElement(elmt, "r"));
            double  *r_K = REAL(getListElement(elmt, "K"));
            double  *r_n = REAL(getListElement(elmt, "n"));

            double delta_r = 0.0;
            double delta_K = 0.0;

            double  *rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit,e));
            double  *rans_N_eit = REAL(VECTOR_ELT(out_N_eit,e));


        //Rprintf("avant T %i\n",ind_t);PrintValue(v_N_e0t);PrintValue(STRING_ELT(sppList,e));PrintValue(getListElement(recList,CHAR(STRING_ELT(sppList,e))));
        //ajout 24/04/2018 pour prise en compte for�age recrutement XSA
        if ((!isNull(getListElement(recList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & (nbI>1) & Reality) { //seulement applicable � t=2
            SRInd[e]=0;
            r_N_e0t[ind_t*fact6_D[3]] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[ind_t]; //t
            if (ind_t<(nbT-1)) r_N_e0t[ind_t*fact6_D[3] + 1] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[ind_t+1]; //t+1
        }
        //Rprintf("apr�s T %i\n",ind_t);PrintValue(v_N_e0t);

           //�quation n�1 : out_Z_eit

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0, tempCap = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                }

                            if (Zoptim_use & (e==eTemp)) {

                               rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    Zoptim[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];


                            } else {

                                if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i[ind_i + ind_t*nbI];
                                }
                            }
                            //on en profite pour calculer Fbar

                            fmax = fmax + (temp + r_Froth_i[ind_i + ind_t*nbI])*r_Fbar[ind_i];
                            sumWt = sumWt + r_Fbar[ind_i];

                            }
//Rprintf("G17");
                    //�quation n�2 : out_N_eit

                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ind_i == 0) {

                                        if ((SRInd[e]==1) & (ind_t>0)) {

                                            if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                    REAL(VECTOR_ELT(out_SRmod,e))[ind_t];

                                            } else {

                                             if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }}




                                        } else {

                                            if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }
                                        }

                                    } else {

                                        if (ind_t == 0) {

                                            rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                        } else {

                                            if (ind_i == (nbI-1)) {  //groupe d'�ge +

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            }
                                        }
                                    }
                            }

//Rprintf("apr�s T %i\n",ind_t);PrintValue(VECTOR_ELT(out_N_eit,e));


                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                    rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;
                                }

     //ajout 01/06/2018 : recrutement al�toire sur la base de recParamList
//Rprintf("G8\n");fichier << "G8" << endl;

        if ((!isNull(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & (nbI>1) & Reality) { //seulement applicable � t>0 et pour une dynamique XSA: ici Reality==TRUE donc forcage avec recParamList remplace la valeur de recList
            double *param = REAL(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param"));
            int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type")); //Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type"));
            int del = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"));

            if ((!ISNA(param[ind_t])) & (ind_t>=del)) {
                double recr = 0.0;

                if (typeSR[ind_t]==1){ // Hockey Stick
                    if ((1/param[ind_t + 1*nbT])>rans_SSB_et[ind_t - del]) {
                        recr = param[ind_t + 0*nbT] * rans_SSB_et[ind_t - del] * param[ind_t + 2*nbT];
                    } else{
                        recr = param[ind_t + 0*nbT] * param[ind_t + 2*nbT] / param[ind_t + 1*nbT];
                    }
                } else if (typeSR[ind_t]==2){ // Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta])
                    recr = (4*param[ind_t + 0*nbT] * param[ind_t + 1*nbT] * rans_SSB_et[ind_t - del]) /
                    (param[ind_t + 2*nbT]*(1-param[ind_t + 0*nbT]) + rans_SSB_et[ind_t - del]*(5*param[ind_t + 0*nbT]-1)) *
                    param[ind_t + 3*nbT] * param[ind_t + 4*nbT];
                }

                //Rprintf("Recruitment in reality = %f\n",recr);fichier << "e = " << CHAR(STRING_ELT(sppList,e)) << "; Recruitment in reality = " << recr << endl;
                r_N_e0t[ind_t] = recr;
                rans_N_eit[0*fact4_D[2] + ind_t*fact4_D[3]] = recr;

            }
        }

    //ajout 01/06/2018 : --------------------------------------------------

//Rprintf("G9\n");fichier << "G9" << endl;
                                //biomasse

                             if (nbI>1) {

                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                }

                             } else { //SPiCT
                                    if ((!isNull(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e)))))){ // si aleatoire sur les parametres K et r
                                        delta_r = REAL(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e))))[ind_t];
                                        delta_K = REAL(getListElement(ParamSPMlist,CHAR(STRING_ELT(sppList,e))))[ind_t+nbT];
                                    }
                                    //Rprintf("t = %i; sp = %s; delta_r = %f \n",ind_t,CHAR(STRING_ELT(sppList,e)), delta_r);
                                    //Rprintf("t = %i; sp = %s; delta_K = %f \n",ind_t,CHAR(STRING_ELT(sppList,e)), delta_K);


                                    double *Bspict = REAL(VECTOR_ELT(intermBIOMspict, e));
                                    double *r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                                    double *r_F_efmit = REAL(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                                    double temp = 0.0;
                                    SEXP dimCst_F_efmit, cFACT2;
                                    PROTECT(dimCst_F_efmit = getAttrib(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))), install("DimCst")));
                                    int *dim_F_efmit = INTEGER(dimCst_F_efmit);
                                    PROTECT(cFACT2 = iDim(dim_F_efmit));
                                    int *fact2_C = INTEGER(cFACT2);
                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {
                                        if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]];
                                    }

                                    rans_B_et[ind_t] =  Bspict[1*ind_t];
                                    if (ind_t <(nbT-1)) {
                                      //Spict equation
                                      //Bspict[ind_t+1] = Bspict[ind_t] + (r_r[0]/(r_n[0]-1))*Bspict[ind_t]*(1-pow(Bspict[ind_t]/r_K[0],r_n[0]-1)) -
                                       //                                          (temp + r_Foth_i[0 + ind_t*1])*Bspict[ind_t];

                                      // Fox equation
                                      Bspict[ind_t+1] = Bspict[ind_t] * (1 + (r_r[0]+delta_r) * log( (r_K[0]+delta_K) / Bspict[ind_t]) - (temp + r_Foth_i[ind_t]));

                                    }
                                    //Rprintf(" Bspict[ind_t+1] = %f, Bspict[ind_t] = %f, r = %f, K = %f, F = %f",Bspict[ind_t+1],Bspict[ind_t], r_r[0], r_K[0], temp + r_Foth_i[ind_t]);
                                    UNPROTECT(2);

                             }

//Rprintf("End_Annual\n");fichier << "End_Annual" << endl;

    } else if ((Qvec[e]==0) & (Svec[e]==1)) {
    //Rprintf("Start_Sex-based\n");fichier << "Start_Sex-based" << endl;

            PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));
            PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2"));
            double  *r_Fr_efmit_G1 = REAL(VECTOR_ELT(out_Fr_fmi_G1, e));//Rprintf("G34");
            double  *r_Fr_efmit_G2 = REAL(VECTOR_ELT(out_Fr_fmi_G2, e));
            double  *r_F_efmit_G1 = REAL(VECTOR_ELT(out_F_fmi_G1, e));
            double  *r_F_efmit_G2 = REAL(VECTOR_ELT(out_F_fmi_G2, e));
            double  *r_M_ei_G1 = REAL(getListElement(elmt, "M_i_G1"));
            double  *r_M_ei_G2  = REAL(getListElement(elmt, "M_i_G2"));

            double  *r_N_e0t_G1 = REAL(v_N_e0t_G1);
            double  *r_N_e0t_G2 = REAL(v_N_e0t_G2);
            double  *r_N_ei0_G1 = REAL(getListElement(elmt, "N_it0_G1"));
            double  *r_N_ei0_G2 = REAL(getListElement(elmt, "N_it0_G2"));
            double  *r_w_ei_G1 = REAL(getListElement(elmt, "wStock_i_G1"));////PrintValue(getListElement(elmt, "wStock_i"));
            double  *r_w_ei_G2 = REAL(getListElement(elmt, "wStock_i_G2"));
            double  *r_mat_ei_G1 = REAL(getListElement(elmt, "mat_i_G1"));
            double  *r_mat_ei_G2 = REAL(getListElement(elmt, "mat_i_G2"));

            double  *r_Froth_i_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 226));
            double  *r_Froth_i_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 227));
            double  *r_Fbar_G1 = REAL(getListElement(elmt, "Fbar_G1"));
            double  *r_Fbar_G2 = REAL(getListElement(elmt, "Fbar_G2"));

            double  *rans_Z_eit_G1 = REAL(VECTOR_ELT(out_Z_eit_G1,e));
            double  *rans_Z_eit_G2 = REAL(VECTOR_ELT(out_Z_eit_G2,e));
            double  *rans_N_eit_G1 = REAL(VECTOR_ELT(out_N_eit_G1,e));
            double  *rans_N_eit_G2 = REAL(VECTOR_ELT(out_N_eit_G2,e));



        //ajout 24/04/2018 pour prise en compte for�age recrutement
        if ((!isNull(getListElement(recList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & (nbI>1) & Reality) { //seulement applicable � t=2
            SRInd[e]=0;
            r_N_e0t_G1[ind_t*fact6_D[3]] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[0 + 2*ind_t]; //fichier << "Rec.1: e = " << CHAR(STRING_ELT(sppList,e)) << ", r_N_e0t_G1 = " << r_N_e0t_G1[ind_t*fact6_D[3]] << endl;
            r_N_e0t_G2[ind_t*fact6_D[3]] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[1 + 2*ind_t];//t
            if (ind_t<(nbT-1)) {
                    r_N_e0t_G1[ind_t*fact6_D[3] + 1] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[0 + 2*(ind_t+1)]; //t+1
                    r_N_e0t_G2[ind_t*fact6_D[3] + 1] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[1 + 2*(ind_t+1)];
            }
        }

           //�quation n�1 : out_Z_eit

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double tempG1 = 0.0, tempCapG1 = 0.0;
                                double tempG2 = 0.0, tempCapG2 = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                if (!ISNA(r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempG1 = tempG1 +  r_Fr_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCapG1 = tempCapG1 +  r_F_efmit_G1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                if (!ISNA(r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempG2 = tempG2 +  r_Fr_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCapG2 = tempCapG2 +  r_F_efmit_G2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                }

                                if (ISNA(r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit_G1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    tempG1 + r_M_ei_G1[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i_G1[ind_i + ind_t*nbI];
                                }

                                if (ISNA(r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit_G2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    tempG2 + r_M_ei_G2[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i_G2[ind_i + ind_t*nbI];
                                }

                            //on en profite pour calculer Fbar

                            fmax = fmax + (tempG1 + r_Froth_i_G1[ind_i + nbI*ind_t]) * r_Fbar_G1[ind_i] * rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] /
                                                    (rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] + rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] ) +
                                                  (tempG2 + r_Froth_i_G2[ind_i + nbI*ind_t]) * r_Fbar_G2[ind_i] * rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] /
                                                    (rans_N_eit_G1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] + rans_N_eit_G2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] );
                                    sumWt = sumWt + (r_Fbar_G1[ind_i] + r_Fbar_G2[ind_i])/2;
                            }

//Rprintf("G17");
                    //�quation n�2 : out_N_eit

                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ind_i == 0) {

                                        if ((SRInd[e]==1) & (ind_t>0)) {

                                            // a faire




                                        } else {

                                            if (ISNA(r_N_e0t_G1[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0_G1[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];
                                                  //fichier << "Rec.2.1: e = " << CHAR(STRING_ELT(sppList,e)) << ", rans_N_eit_G1 = " << rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] << endl;


                                            } else {

                                                rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t_G1[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];
                                                  //fichier << "Rec.2.2: e = " << CHAR(STRING_ELT(sppList,e)) << ", rans_N_eit_G1 = " << rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] << endl;

                                            }

                                            if (ISNA(r_N_e0t_G2[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0_G2[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                            } else {

                                                rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t_G2[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }
                                        }

                                    } else {

                                        if (ind_t == 0) {

                                            rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0_G1[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                              rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0_G2[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                        } else {

                                            if (ind_i == (nbI-1)) {  //groupe d'�ge +

                                                rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G1[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                  rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G1[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                  rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G2[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                  rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G2[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            } else {

                                                rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G1[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                  rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit_G2[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            }
                                        }
                                    }
                            }

//Rprintf("apr�s T %i\n",ind_t);PrintValue(VECTOR_ELT(out_N_eit,e));


                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypoth�se que la dimension �ge est toujours pr�sente
                                        temp = temp +
                                         rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei_G1[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei_G1[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000 +

                                         rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei_G2[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei_G2[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                    rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;
                                }

     //ajout 01/06/2018 : recrutement al�toire sur la base de recParamList
//Rprintf("G8\n");fichier << "G8" << endl;

        if ((!isNull(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & (nbI>1) & Reality) { //seulement applicable � t>0: ici Reality==TRUE donc forcage avec recParamList remplace la valeur de recList

            double *param = REAL(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param"));
            double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"ventil")); //Rprintf("ventil = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"ventil"));
            int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type")); //Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type"));
            int del = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"));

            if ((!ISNA(param[ind_t])) & (ind_t>=del)) {
                double recr = 0.0;

                if (typeSR[ind_t]==1){ // Hockey Stick
                    if ((1/param[ind_t + 1*nbT])>rans_SSB_et[ind_t - del]) {
                        recr = param[ind_t + 0*nbT] * rans_SSB_et[ind_t - del] * param[ind_t + 2*nbT];
                    } else{
                        recr = param[ind_t + 0*nbT] * param[ind_t + 2*nbT] / param[ind_t + 1*nbT];
                    }
                } else if (typeSR[ind_t]==2){ // Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta])
                    recr = (4*param[ind_t + 0*nbT] * param[ind_t + 1*nbT] * rans_SSB_et[ind_t - del]) /
                    (param[ind_t + 2*nbT]*(1-param[ind_t + 0*nbT]) + rans_SSB_et[ind_t - del]*(5*param[ind_t + 0*nbT]-1)) *
                    param[ind_t + 3*nbT] * param[ind_t + 4*nbT];
                }


                //Rprintf("Recruitment in reality = %f\n",recr);fichier << "e = " << CHAR(STRING_ELT(sppList,e)) << "Recruitment in reality = " << recr << endl;
                r_N_e0t_G1[ind_t] = recr*ventil[0];
                rans_N_eit_G1[0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[0]; //fichier << "Rec.3: e = " << CHAR(STRING_ELT(sppList,e)) << ", delay = " << del <<", SSB = " << rans_SSB_et[ind_t - del] << ", r_N_e0t_G1 = " << r_N_e0t_G1[ind_t] << endl;
                //Rprintf("Recruitment in reality G1 = %f\n",recr);fichier << "e = " << CHAR(STRING_ELT(sppList,e)) << "Recruitment in reality G1 = " << rans_N_eit_G1[0*fact4_D[2] + ind_t*fact4_D[3]] << endl;
                r_N_e0t_G2[ind_t] = recr*ventil[1];
                rans_N_eit_G2[0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[1];

            }
        }



    //ajout 01/06/2018 : --------------------------------------------------

//Rprintf("G9\n");fichier << "G9" << endl;
                                //biomasse

                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                        temp = temp +
                                         rans_N_eit_G1[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei_G1[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000 +

                                         rans_N_eit_G2[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei_G2[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000 ;

                                    rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                }

//Rprintf("End_Sex-based\n");fichier << "End_Sex-based" << endl;

    } else if ((Qvec[e]==1) & (Svec[e]==0)) {
//Rprintf("Start_Quarterly\n");fichier << "Start_Quarterly\n" << endl;

                    PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                    PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                    PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                    PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                    PROTECT(v_N_ei0_S1M1 = getListElement(elmt, "Nt0_S1M1"));
                    PROTECT(v_N_ei0_S1M2 = getListElement(elmt, "Nt0_S1M2"));
                    PROTECT(v_N_ei0_S1M3 = getListElement(elmt, "Nt0_S1M3"));
                    PROTECT(v_N_ei0_S1M4 = getListElement(elmt, "Nt0_S1M4")); // ### +4

                    PROTECT(v_matwt_M1 = getListElement(elmt, "matwt_M1"));
                    PROTECT(v_matwt_M2 = getListElement(elmt, "matwt_M2"));
                    PROTECT(v_matwt_M3 = getListElement(elmt, "matwt_M3"));
                    PROTECT(v_matwt_M4 = getListElement(elmt, "matwt_M4"));  // ### +4

                    double  *rans_Z_eit_S1M1 = REAL(VECTOR_ELT(out_Z_eit_S1M1,e)); double  *rans_Z_eit_S1M2 = REAL(VECTOR_ELT(out_Z_eit_S1M2,e));
                    double  *rans_Z_eit_S1M3 = REAL(VECTOR_ELT(out_Z_eit_S1M3,e)); double  *rans_Z_eit_S1M4 = REAL(VECTOR_ELT(out_Z_eit_S1M4,e));
                    double  *rans_Z_eit_S2M1 = REAL(VECTOR_ELT(out_Z_eit_S2M1,e)); double  *rans_Z_eit_S2M2 = REAL(VECTOR_ELT(out_Z_eit_S2M2,e));
                    double  *rans_Z_eit_S2M3 = REAL(VECTOR_ELT(out_Z_eit_S2M3,e)); double  *rans_Z_eit_S2M4 = REAL(VECTOR_ELT(out_Z_eit_S2M4,e));
                    double  *rans_Z_eit_S3M1 = REAL(VECTOR_ELT(out_Z_eit_S3M1,e)); double  *rans_Z_eit_S3M2 = REAL(VECTOR_ELT(out_Z_eit_S3M2,e));
                    double  *rans_Z_eit_S3M3 = REAL(VECTOR_ELT(out_Z_eit_S3M3,e)); double  *rans_Z_eit_S3M4 = REAL(VECTOR_ELT(out_Z_eit_S3M4,e));
                    double  *rans_Z_eit_S4M1 = REAL(VECTOR_ELT(out_Z_eit_S4M1,e)); double  *rans_Z_eit_S4M2 = REAL(VECTOR_ELT(out_Z_eit_S4M2,e));
                    double  *rans_Z_eit_S4M3 = REAL(VECTOR_ELT(out_Z_eit_S4M3,e)); double  *rans_Z_eit_S4M4 = REAL(VECTOR_ELT(out_Z_eit_S4M4,e));

                    double  *rans_N_eit_S1M1 = REAL(VECTOR_ELT(out_N_eit_S1M1,e)); double  *rans_N_eit_S1M2 = REAL(VECTOR_ELT(out_N_eit_S1M2,e));
                    double  *rans_N_eit_S1M3 = REAL(VECTOR_ELT(out_N_eit_S1M3,e)); double  *rans_N_eit_S1M4 = REAL(VECTOR_ELT(out_N_eit_S1M4,e));
                    double  *rans_N_eit_S2M1 = REAL(VECTOR_ELT(out_N_eit_S2M1,e)); double  *rans_N_eit_S2M2 = REAL(VECTOR_ELT(out_N_eit_S2M2,e));
                    double  *rans_N_eit_S2M3 = REAL(VECTOR_ELT(out_N_eit_S2M3,e)); double  *rans_N_eit_S2M4 = REAL(VECTOR_ELT(out_N_eit_S2M4,e));
                    double  *rans_N_eit_S3M1 = REAL(VECTOR_ELT(out_N_eit_S3M1,e)); double  *rans_N_eit_S3M2 = REAL(VECTOR_ELT(out_N_eit_S3M2,e));
                    double  *rans_N_eit_S3M3 = REAL(VECTOR_ELT(out_N_eit_S3M3,e)); double  *rans_N_eit_S3M4 = REAL(VECTOR_ELT(out_N_eit_S3M4,e));
                    double  *rans_N_eit_S4M1 = REAL(VECTOR_ELT(out_N_eit_S4M1,e)); double  *rans_N_eit_S4M2 = REAL(VECTOR_ELT(out_N_eit_S4M2,e));
                    double  *rans_N_eit_S4M3 = REAL(VECTOR_ELT(out_N_eit_S4M3,e)); double  *rans_N_eit_S4M4 = REAL(VECTOR_ELT(out_N_eit_S4M4,e));

                    double  *r_Fr_efmit_S1M1 = REAL(VECTOR_ELT(out_Fr_fmi_S1M1,e)); double  *r_Fr_efmit_S1M2 = REAL(VECTOR_ELT(out_Fr_fmi_S1M2,e));
                    double  *r_Fr_efmit_S1M3 = REAL(VECTOR_ELT(out_Fr_fmi_S1M3,e)); double  *r_Fr_efmit_S1M4 = REAL(VECTOR_ELT(out_Fr_fmi_S1M4,e));
                    double  *r_Fr_efmit_S2M1 = REAL(VECTOR_ELT(out_Fr_fmi_S2M1,e)); double  *r_Fr_efmit_S2M2 = REAL(VECTOR_ELT(out_Fr_fmi_S2M2,e));
                    double  *r_Fr_efmit_S2M3 = REAL(VECTOR_ELT(out_Fr_fmi_S2M3,e)); double  *r_Fr_efmit_S2M4 = REAL(VECTOR_ELT(out_Fr_fmi_S2M4,e));
                    double  *r_Fr_efmit_S3M1 = REAL(VECTOR_ELT(out_Fr_fmi_S3M1,e)); double  *r_Fr_efmit_S3M2 = REAL(VECTOR_ELT(out_Fr_fmi_S3M2,e));
                    double  *r_Fr_efmit_S3M3 = REAL(VECTOR_ELT(out_Fr_fmi_S3M3,e)); double  *r_Fr_efmit_S3M4 = REAL(VECTOR_ELT(out_Fr_fmi_S3M4,e));
                    double  *r_Fr_efmit_S4M1 = REAL(VECTOR_ELT(out_Fr_fmi_S4M1,e)); double  *r_Fr_efmit_S4M2 = REAL(VECTOR_ELT(out_Fr_fmi_S4M2,e));
                    double  *r_Fr_efmit_S4M3 = REAL(VECTOR_ELT(out_Fr_fmi_S4M3,e)); double  *r_Fr_efmit_S4M4 = REAL(VECTOR_ELT(out_Fr_fmi_S4M4,e));

                    double  *r_F_efmit_S1M1 = REAL(VECTOR_ELT(out_F_fmi_S1M1,e)); double  *r_F_efmit_S1M2 = REAL(VECTOR_ELT(out_F_fmi_S1M2,e));
                    double  *r_F_efmit_S1M3 = REAL(VECTOR_ELT(out_F_fmi_S1M3,e)); double  *r_F_efmit_S1M4 = REAL(VECTOR_ELT(out_F_fmi_S1M4,e));
                    double  *r_F_efmit_S2M1 = REAL(VECTOR_ELT(out_F_fmi_S2M1,e)); double  *r_F_efmit_S2M2 = REAL(VECTOR_ELT(out_F_fmi_S2M2,e));
                    double  *r_F_efmit_S2M3 = REAL(VECTOR_ELT(out_F_fmi_S2M3,e)); double  *r_F_efmit_S2M4 = REAL(VECTOR_ELT(out_F_fmi_S2M4,e));
                    double  *r_F_efmit_S3M1 = REAL(VECTOR_ELT(out_F_fmi_S3M1,e)); double  *r_F_efmit_S3M2 = REAL(VECTOR_ELT(out_F_fmi_S3M2,e));
                    double  *r_F_efmit_S3M3 = REAL(VECTOR_ELT(out_F_fmi_S3M3,e)); double  *r_F_efmit_S3M4 = REAL(VECTOR_ELT(out_F_fmi_S3M4,e));
                    double  *r_F_efmit_S4M1 = REAL(VECTOR_ELT(out_F_fmi_S4M1,e)); double  *r_F_efmit_S4M2 = REAL(VECTOR_ELT(out_F_fmi_S4M2,e));
                    double  *r_F_efmit_S4M3 = REAL(VECTOR_ELT(out_F_fmi_S4M3,e)); double  *r_F_efmit_S4M4 = REAL(VECTOR_ELT(out_F_fmi_S4M4,e));

                    double  *r_N_ei0_S1M1 = REAL(v_N_ei0_S1M1); double  *r_N_ei0_S1M2 = REAL(v_N_ei0_S1M2);
                    double  *r_N_ei0_S1M3 = REAL(v_N_ei0_S1M3); double  *r_N_ei0_S1M4 = REAL(v_N_ei0_S1M4);

                    double  *r_matwt_M1 = REAL(v_matwt_M1); double  *r_matwt_M2 = REAL(v_matwt_M2);
                    double  *r_matwt_M3 = REAL(v_matwt_M3); double  *r_matwt_M4 = REAL(v_matwt_M4);

                    double  *r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1); double  *r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                    double  *r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3); double  *r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                    double  *r_M_ei = REAL(getListElement(elmt, "M_i"));
                    double  *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));
                    double  *r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                    double  *r_w_ei = REAL(getListElement(elmt, "wStock_i"));
                    double  *r_Fbar = REAL(getListElement(elmt, "Fbar"));

                    double  *rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit,e));
                    double  *rans_N_eit = REAL(VECTOR_ELT(out_N_eit,e));

                    double  *r_Fr_efmit = REAL(VECTOR_ELT(out_Fr_fmi, e));//Rprintf("G34");
                    double  *r_F_efmit = REAL(VECTOR_ELT(out_F_fmi, e));


        //ajout 24/04/2018 pour prise en compte for�age recrutement XSA
        if ((!isNull(getListElement(recList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>1) & (nbI>1) & Reality) { //seulement applicable � t=3
            r_N_e0t_S1M1[0] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[0 + 4*ind_t];
            r_N_e0t_S2M2[0] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[1 + 4*ind_t];
            r_N_e0t_S3M3[0] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[2 + 4*ind_t];
            r_N_e0t_S4M4[0] = REAL(getListElement(recList,CHAR(STRING_ELT(sppList,e))))[3 + 4*ind_t];
        }


                    double  *r_Foth_i_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 116));
                    double  *r_Foth_i_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 117));
                    double  *r_Foth_i_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 118));
                    double  *r_Foth_i_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 119));
                    double  *r_Foth_i_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 120));
                    double  *r_Foth_i_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 121));
                    double  *r_Foth_i_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 122));
                    double  *r_Foth_i_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 123));
                    double  *r_Foth_i_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 124));
                    double  *r_Foth_i_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 125));
                    double  *r_Foth_i_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 126));
                    double  *r_Foth_i_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 127));
                    double  *r_Foth_i_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 128));
                    double  *r_Foth_i_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 129));
                    double  *r_Foth_i_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 130));
                    double  *r_Foth_i_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 131));

                    double  *r_Froth_i_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 132));
                    double  *r_Froth_i_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 133));
                    double  *r_Froth_i_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 134));
                    double  *r_Froth_i_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 135));
                    double  *r_Froth_i_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 136));
                    double  *r_Froth_i_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 137));
                    double  *r_Froth_i_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 138));
                    double  *r_Froth_i_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 139));
                    double  *r_Froth_i_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 140));
                    double  *r_Froth_i_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 141));
                    double  *r_Froth_i_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 142));
                    double  *r_Froth_i_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 143));
                    double  *r_Froth_i_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 144));
                    double  *r_Froth_i_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 145));
                    double  *r_Froth_i_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 146));
                    double  *r_Froth_i_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 147));


                                    //Z
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M1[ind_i + nbI*ind_t];
                                        }

                                   if (ind_i==0) {rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M2[ind_i + nbI*ind_t];
                                        }
                                   }


                                   if (ind_i==0) {rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M3[ind_i + nbI*ind_t];
                                        }
                                   }

                                   if (ind_i==0) {rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S1M4[ind_i + nbI*ind_t];
                                        }
                                   }


                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M1[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M2[ind_i + nbI*ind_t];
                                        }

                                   if (ind_i==0) {rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M3[ind_i + nbI*ind_t];
                                        }
                                   }

                                   if (ind_i==0) {rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S2M4[ind_i + nbI*ind_t];
                                        }
                                   }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M1[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M2[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M3[ind_i + nbI*ind_t];
                                        }

                                   if (ind_i==0) {rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = 0.0;

                                   } else {

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S3M4[ind_i + nbI*ind_t];
                                        }
                                   }


                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M1[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M2[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M3[ind_i + nbI*ind_t];
                                        }

                                        temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {
                                          if (!ISNA(r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                              temp = temp +  r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        }
                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i_S4M4[ind_i + nbI*ind_t];
                                        }

                                    }

//fichier << "mm" << endl;
                    if (ZoptSS3) {

                        double *Ztemp = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,e))));

                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                           rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(0*nbI)];
                           rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(1*nbI)];
                           rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(2*nbI)];
                           rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(3*nbI)];

                           rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(4*nbI)];
                           rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(5*nbI)];
                           rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(6*nbI)];
                           rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(7*nbI)];

                           rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(8*nbI)];
                           rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(9*nbI)];
                           rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(10*nbI)];
                           rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(11*nbI)];

                           rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(12*nbI)];
                           rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(13*nbI)];
                           rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(14*nbI)];
                           rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = Ztemp[ind_i+(15*nbI)];

                        }

                    }

//fichier << "mm2" << endl;

                                    //---------
                                    // calcul de N_eit
                                    //---------
//Rprintf("Calcul.N_eit\n"); fichier << "G2\n" << endl;
if (ind_t==1) {


                                        //S1
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                            if (ind_i == 0) { //recrutement

                                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t_S1M1[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;
                                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;
                                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;

                                            } else {

                                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                 r_N_ei0_S1M1[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                 r_N_ei0_S1M2[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                 r_N_ei0_S1M3[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                 r_N_ei0_S1M4[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];

                                            }
                                        }

                                        //S2
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S2M2[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                           } else {

                                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }

                                        //S3
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S3M3[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                           } else {

                                              rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }

                                        //S4
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S4M4[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];

                                           } else {

                                              rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }

}




if (ind_t>1) {

                                        //S1
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                            if (ind_i == 0) { //recrutement

                                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t_S1M1[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;
                                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;
                                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] = 0.0;

                                            } else {

                                                if (ind_i==(nbI-1)) {

                                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4) +
                                                rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4) +
                                                rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4) +
                                                rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4) +
                                                rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);

                                                } else {

                                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M1[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M2[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M3[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);
                                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]]*
                                                exp(-rans_Z_eit_S4M4[0*fact1_D[0] + 0*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]/4);

                                                }
                                            }
                                        }

                                        //S2
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S2M2[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                           } else {

                                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S1M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }

                                        //S3
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S3M3[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];
                                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                           } else {

                                              rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S2M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }

                                        //S4
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                           if (ind_i == 0) { //recrutement

                                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4) +
                                                  r_N_e0t_S4M4[0*fact5_D[0] + 0*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]];

                                           } else {

                                              rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);
                                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                                  exp(-rans_Z_eit_S3M4[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]/4);

                                            }

                                        }


}




// on peut d�sormais �valuer F, Z et N au niveau annuel et global
            for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++)
            for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++)
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =

                            (r_F_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_F_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_F_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));


                r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]] =

                            (r_Fr_efmit_S1M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S1M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S1M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S2M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S2M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S2M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S3M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S3M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S3M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]])) +

                            (r_Fr_efmit_S4M1[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M2[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M3[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                            r_Fr_efmit_S4M4[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]*
                               rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]) /
                               (4*(rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] +
                                   rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]]));

            }




for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

     r_Froth_i[ind_i+ind_t*nbI] =
      (r_Froth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
      r_Froth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
      r_Froth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
      r_Froth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
      r_Froth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
      r_Froth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
      r_Froth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
      r_Froth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
      r_Froth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
      r_Froth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +

      (r_Froth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
      r_Froth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
      r_Froth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
      r_Froth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));


    r_Foth_i[ind_i+ind_t*nbI] =
      (r_Foth_i_S1M1[ind_i+ind_t*nbI]*rans_N_eit_S1M1[ind_i+ind_t*nbI] +
      r_Foth_i_S1M2[ind_i+ind_t*nbI]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
      r_Foth_i_S1M3[ind_i+ind_t*nbI]*rans_N_eit_S1M3[ind_i+ind_t*nbI] +
      r_Foth_i_S1M4[ind_i+ind_t*nbI]*rans_N_eit_S1M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
           rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S2M1[ind_i+ind_t*nbI]*rans_N_eit_S2M1[ind_i+ind_t*nbI] +
      r_Foth_i_S2M2[ind_i+ind_t*nbI]*rans_N_eit_S2M2[ind_i+ind_t*nbI] +
      r_Foth_i_S2M3[ind_i+ind_t*nbI]*rans_N_eit_S2M3[ind_i+ind_t*nbI] +
      r_Foth_i_S2M4[ind_i+ind_t*nbI]*rans_N_eit_S2M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S2M1[ind_i+ind_t*nbI] + rans_N_eit_S2M2[ind_i+ind_t*nbI] +
           rans_N_eit_S2M3[ind_i+ind_t*nbI] + rans_N_eit_S2M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S3M1[ind_i+ind_t*nbI]*rans_N_eit_S3M1[ind_i+ind_t*nbI] +
      r_Foth_i_S3M2[ind_i+ind_t*nbI]*rans_N_eit_S3M2[ind_i+ind_t*nbI] +
      r_Foth_i_S3M3[ind_i+ind_t*nbI]*rans_N_eit_S3M3[ind_i+ind_t*nbI] +
      r_Foth_i_S3M4[ind_i+ind_t*nbI]*rans_N_eit_S3M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S3M1[ind_i+ind_t*nbI] + rans_N_eit_S3M2[ind_i+ind_t*nbI] +
           rans_N_eit_S3M3[ind_i+ind_t*nbI] + rans_N_eit_S3M4[ind_i+ind_t*nbI])) +

      (r_Foth_i_S4M1[ind_i+ind_t*nbI]*rans_N_eit_S4M1[ind_i+ind_t*nbI] +
      r_Foth_i_S4M2[ind_i+ind_t*nbI]*rans_N_eit_S4M2[ind_i+ind_t*nbI] +
      r_Foth_i_S4M3[ind_i+ind_t*nbI]*rans_N_eit_S4M3[ind_i+ind_t*nbI] +
      r_Foth_i_S4M4[ind_i+ind_t*nbI]*rans_N_eit_S4M4[ind_i+ind_t*nbI]) /
       (4*(rans_N_eit_S4M1[ind_i+ind_t*nbI] + rans_N_eit_S4M2[ind_i+ind_t*nbI] +
           rans_N_eit_S4M3[ind_i+ind_t*nbI] + rans_N_eit_S4M4[ind_i+ind_t*nbI]));

}


                                double fmax = 0.0, sumWt = 0.0;

                                    //�quation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        }


                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i[ind_i + nbI*ind_t];

                                        }


                                    //on en profite pour calculer Fbar

                                    fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                                    sumWt = sumWt + r_Fbar[ind_i];

                                    // et remplir N_eit (effectifs � la saison 1)
                                    rans_N_eit[ind_i + nbI*ind_t] = rans_N_eit_S1M1[ind_i+ind_t*nbI] + rans_N_eit_S1M2[ind_i+ind_t*nbI] +
                                         rans_N_eit_S1M3[ind_i+ind_t*nbI] + rans_N_eit_S1M4[ind_i+ind_t*nbI];

                                    }


                        double tempSSB = 0.0;
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                            tempSSB = tempSSB + r_matwt_M1[ind_i]*rans_N_eit_S1M1[ind_i+ind_t*nbI] + r_matwt_M2[ind_i]*rans_N_eit_S1M2[ind_i+ind_t*nbI] +
                                  r_matwt_M3[ind_i]*rans_N_eit_S1M3[ind_i+ind_t*nbI] + r_matwt_M4[ind_i]*rans_N_eit_S1M4[ind_i+ind_t*nbI];
                        }


                        rans_SSB_et[ind_t] = tempSSB ;


    //ajout 01/06/2018 : recrutement al�toire sur la base de recParamList

        if ((!isNull(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))))) & (ind_t>0) & Reality) {//seulement applicable � t>0 et pour une dynamique SS3
            double *param = REAL(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"param"));
            double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"ventil")); //Rprintf("ventil = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"ventil"));
            int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type")); //Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"type"));
            int del = INTEGER(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(STRING_ELT(sppList,e))),"delay"));

            if ((!ISNA(param[ind_t])) & (ind_t>=del)) {
                double recr = 0.0;

                if (typeSR[ind_t]==1){ // Hockey Stick
                    if ((1/param[ind_t + 1*nbT])>rans_SSB_et[ind_t - del]) {
                        recr = param[ind_t + 0*nbT] * rans_SSB_et[ind_t - del] * param[ind_t + 2*nbT];
                    } else{
                        recr = param[ind_t + 0*nbT] * param[ind_t + 2*nbT] / param[ind_t + 1*nbT];
                    }
                } else if (typeSR[ind_t]==2){ // Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta])
                    recr = (4*param[ind_t + 0*nbT] * param[ind_t + 1*nbT] * rans_SSB_et[ind_t - del]) /
                    (param[ind_t + 2*nbT]*(1-param[ind_t + 0*nbT]) + rans_SSB_et[ind_t - del]*(5*param[ind_t + 0*nbT]-1)) *
                    param[ind_t + 3*nbT] * param[ind_t + 4*nbT];
                }

                        //S1
                rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[0];

                        //S2
                rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                         rans_N_eit_S1M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                         exp(-rans_Z_eit_S1M1[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[1];

                        //S3
                rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                          rans_N_eit_S2M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                          exp(-rans_Z_eit_S2M1[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                          rans_N_eit_S2M2[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                          exp(-rans_Z_eit_S2M2[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[2];

                         //S4
                rans_N_eit_S4M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                           rans_N_eit_S3M1[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                           exp(-rans_Z_eit_S3M1[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S4M2[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                           rans_N_eit_S3M2[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                           exp(-rans_Z_eit_S3M2[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S4M3[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] =
                           rans_N_eit_S3M3[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]]*
                           exp(-rans_Z_eit_S3M3[0*fact1_D[0] + 0*fact1_D[1] + 0*fact1_D[2] + ind_t*fact1_D[3]]/4);
                rans_N_eit_S4M4[0*fact4_D[0] + 0*fact4_D[1] + 0*fact4_D[2] + ind_t*fact4_D[3]] = recr*ventil[3];


                rans_N_eit[0 + nbI*ind_t] = rans_N_eit_S1M1[0+ind_t*nbI] + rans_N_eit_S1M2[0+ind_t*nbI] +
                           rans_N_eit_S1M3[0+ind_t*nbI] + rans_N_eit_S1M4[0+ind_t*nbI];


            }
        }

    //ajout 01/06/2018 : --------------------------------------------------


                      rans_Fbar_et[ind_t] = fmax/sumWt;



                                for (int ind_f = 0 ; ind_f < 1 ; ind_f++)
                                for (int ind_m = 0 ; ind_m < 1 ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                }



//Rprintf("End_Quarterly\n");fichier << "End_Quarterly\n" << endl;
    }
//Rprintf("K10\n");fichier << "K10" << endl;
                    UNPROTECT(1);//Rprintf("K11\n");
                    if ((Qvec[e]==1) & (Svec[e]==0)) UNPROTECT(4+4+4);
                    if ((Qvec[e]==0) & (Svec[e]==0)) UNPROTECT(1);
                    if ((Qvec[e]==0) & (Svec[e]==1)) UNPROTECT(2);

}}

//Rprintf("End\n");fichier << "End" << endl;
//fichier.close();
}
}

