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

//------------------------------------------
// Module 'Captures, rejets et d�barquements'
//------------------------------------------

extern "C" {

void BioEcoPar::CatchDL(SEXP list, int ind_t, SEXP EVAR, int VERBOSE = 0)
{

//        string str1b, str2b, str3b;
//        str1b = "debugHKEcatch";//str1b = "C:\\Users\\mmerzere\\Desktop\\test3\\debugHKEcatch_V";//
//        str3b = "_V";
//        str2b = ".txt";
//
//        std::stringstream ssi, ssy;
//        ssi << ind_t;
//        ssy << EcoIndCopy[0];
//        str1b = str1b + ssi.str() + str3b + ssy.str() + str2b;
//
//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.Catch_DL.txt" , ios::out | ios::trunc);

//
//       fichier << "D�but" << endl;
//


SEXP Flist;
PROTECT(Flist = getListElement(list, "Fleet")); //unp

double *rans_Yothsue_fm = REAL(getListElement(Flist, "Yothsue_f_m"));
double *reff1 = REAL(getListElement(Flist, "effort1_f_m"));
double *reff2 = REAL(getListElement(Flist, "effort2_f_m"));
double *rnbv = REAL(getListElement(Flist, "nbv_f_m"));


if (cUpdate) {
if(VERBOSE){Rprintf(" pUpdate : ");}

    SEXP    elmt, dimCst, Dim, dimCst_F_efmit, dimCst_N_eit, dimCst_Z_eit, dimCst_wL_ei, dimCst_wD_ei, dimCst_d_efmit,
            dimCst_d_eStat, dimCst_LPUE_eStat, dimCst_eStat, v_d_eStat, v_LPUE_eStat, v_B_et,
            intAge, v_F_efmit=R_NilValue, v_N_eit=R_NilValue, v_Z_eit=R_NilValue, v_wL_ei=R_NilValue, v_wD_ei=R_NilValue,
             v_d_efmit=R_NilValue, v_doth_eit=R_NilValue, dimCst2, Dim2,
             v_F_efmit_G1=R_NilValue, v_N_eit_G1=R_NilValue, v_Z_eit_G1=R_NilValue, v_wL_ei_G1=R_NilValue, v_wD_ei_G1=R_NilValue,
              v_d_efmit_G1=R_NilValue, v_doth_eit_G1=R_NilValue,
             v_F_efmit_G2=R_NilValue, v_N_eit_G2=R_NilValue, v_Z_eit_G2=R_NilValue, v_wL_ei_G2=R_NilValue, v_wD_ei_G2=R_NilValue,
              v_d_efmit_G2=R_NilValue, v_doth_eit_G2=R_NilValue,
            cFACT1, cFACT2, cFACT3, cFACT4, cFACT5, cFACT6, cFACT7,
            dimYtot, dimCstYtot, dimNamYtot,
            dimCstOQ_ft, dimCstOQ_t, dimnames_oqD_eft, dimnames_oqD_et, dimCstL_et, dimNamL_et, nDim;

    SEXP ans_C_efmit=R_NilValue, ans_Y_efmit=R_NilValue, ans_D_efmit=R_NilValue, ans_L_efmit=R_NilValue,
         dimnames=R_NilValue, rnames_Esp=R_NilValue, ans_C_eit=R_NilValue, ans_Y_eit=R_NilValue, ans_L_eit=R_NilValue, dimnames2=R_NilValue,
         ans_Ystat=R_NilValue, ans_Lstat=R_NilValue, ans_Dstat=R_NilValue, dimnames_eStat=R_NilValue, rnames_eStat=R_NilValue,
         ans_DD_efmit=R_NilValue, ans_LD_efmit=R_NilValue,
         ans_C_eit_G1=R_NilValue, ans_C_eit_G2=R_NilValue,
         ans_C_efmit_G1=R_NilValue,
         ans_C_efmit_G2=R_NilValue,
         ans_statDD=R_NilValue, ans_statLD=R_NilValue, ans_statLDst=R_NilValue, ans_statLDor=R_NilValue,
         ans_oqD_eft=R_NilValue, ans_oqD_et=R_NilValue, ans_oqDstat=R_NilValue,
         rnames_eAll=R_NilValue, ans_L_et=R_NilValue;

    int *dim_F_efmit, *dim_N_eit, *dim_Z_eit, *dim_wL_ei, *dim_wD_ei, *dim_d_efmit, *dimC, *dim, *dim2, *dimcst2,
            *dim_d_eStat, *dim_LPUE_eStat, *dim_eStat, *int_dimYtot, *int_dimCstYtot, *dimOQ_ft, *dimOQ_t, *nd, *int_dimCstL_et;
    int nbI, ind_e;

    double *rans_C_efmit=&NA_REAL, *rans_Y_efmit=&NA_REAL, *rans_D_efmit=&NA_REAL, *rans_L_efmit=&NA_REAL, *r_F_efmit=&NA_REAL, *r_N_eit=&NA_REAL,
           *r_Z_eit=&NA_REAL, *r_wL_ei=&NA_REAL, *r_wD_ei=&NA_REAL, *r_d_efmit=&NA_REAL, *r_B_et=&NA_REAL,
           *rans_C_efmit_G1=&NA_REAL, *rans_C_eit_G1=&NA_REAL, *r_F_efmit_G1=&NA_REAL, *r_N_eit_G1=&NA_REAL,
           *rans_C_efmit_G2=&NA_REAL, *rans_C_eit_G2=&NA_REAL, *r_F_efmit_G2=&NA_REAL, *r_N_eit_G2=&NA_REAL,
           *r_Z_eit_G1=&NA_REAL, *r_wL_ei_G1=&NA_REAL, *r_wD_ei_G1=&NA_REAL, *r_d_efmit_G1=&NA_REAL,
           *r_Z_eit_G2=&NA_REAL, *r_wL_ei_G2=&NA_REAL, *r_wD_ei_G2=&NA_REAL, *r_d_efmit_G2=&NA_REAL,
            *rans_C_eit=&NA_REAL, *rans_Y_eit=&NA_REAL,*rans_L_eit=&NA_REAL, *rans_Ystat=&NA_REAL, *rans_Lstat=&NA_REAL, *rans_Dstat=&NA_REAL,
            *rans_Ytot_fm=&NA_REAL, *rans_DD_efmit=&NA_REAL,
            *rans_LD_efmit=&NA_REAL, *rans_statDD=&NA_REAL, *rans_statLD=&NA_REAL, *rans_statLDst=&NA_REAL, *rans_statLDor=&NA_REAL, *doth_eit=&NA_REAL,
            *doth_eit_G1=&NA_REAL,*doth_eit_G2=&NA_REAL,*rans_oqD_eft=&NA_REAL, *rans_oqD_et=&NA_REAL, *rans_oqDstat=&NA_REAL,
            *r_Foth_i=&NA_REAL, *r_Foth_i_G1=&NA_REAL, *r_Foth_i_G2=&NA_REAL, *rans_L_et;

if (ind_t==0) {

    PROTECT(rnames_Esp = allocVector(STRSXP, nbE));
    setAttrib(out_C_efmit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_C_eit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_C_efmit_G1, R_NamesSymbol, rnames_Esp);
    setAttrib(out_C_eit_G1, R_NamesSymbol, rnames_Esp);
    setAttrib(out_C_efmit_G2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_C_eit_G2, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Y_eit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_L_eit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_Y_efmit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_D_efmit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_L_efmit, R_NamesSymbol, rnames_Esp);
    setAttrib(out_DD_efmi, R_NamesSymbol, rnames_Esp);
    setAttrib(out_LD_efmi, R_NamesSymbol, rnames_Esp);
    setAttrib(out_oqD_eft, R_NamesSymbol, rnames_Esp);
    setAttrib(out_oqD_et, R_NamesSymbol, rnames_Esp);

    PROTECT(rnames_eAll = allocVector(STRSXP, nbEall));
    setAttrib(out_L_et, R_NamesSymbol, rnames_eAll);

    PROTECT(rnames_eStat = allocVector(STRSXP, nbEstat)); //+1 t0
    setAttrib(out_Ystat, R_NamesSymbol, rnames_eStat);
    setAttrib(out_Lstat, R_NamesSymbol, rnames_eStat);
    setAttrib(out_Dstat, R_NamesSymbol, rnames_eStat);
    setAttrib(out_statDD_efm, R_NamesSymbol, rnames_eStat);
    setAttrib(out_statLD_efm, R_NamesSymbol, rnames_eStat);
    setAttrib(out_statLDst_efm, R_NamesSymbol, rnames_eStat);
    setAttrib(out_statLDor_efm, R_NamesSymbol, rnames_eStat);
    setAttrib(out_oqDstat, R_NamesSymbol, rnames_eStat);

}


PROTECT(dimCstYtot = allocVector(INTSXP, 4)); //unp
int_dimCstYtot = INTEGER(dimCstYtot); int_dimCstYtot[0] = nbF; int_dimCstYtot[1] =nbM; int_dimCstYtot[2] = 0; int_dimCstYtot[3] = nbT;
PROTECT(dimYtot = allocVector(INTSXP, 3)); //unp
int_dimYtot = INTEGER(dimYtot); int_dimYtot[0] = nbF; int_dimYtot[1] = nbM; int_dimYtot[2] = nbT;
PROTECT(dimNamYtot = allocVector(VECSXP,3)); //unp
  SET_VECTOR_ELT(dimNamYtot, 0, fleetList);
  SET_VECTOR_ELT(dimNamYtot, 1, metierList);
  SET_VECTOR_ELT(dimNamYtot, 2, times);
setAttrib(out_Ytot_fm, R_DimSymbol, dimYtot);
setAttrib(out_Ytot_fm, R_DimNamesSymbol, dimNamYtot);
setAttrib(out_Ytot_fm, install("DimCst"), dimCstYtot);

rans_Ytot_fm = REAL(out_Ytot_fm);

//on l'initialise avec Yothsue_fm


for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
  rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Yothsue_fm[ind_f + nbF*ind_m]*reff1[ind_f + nbF*ind_m]*reff2[ind_f + nbF*ind_m]*rnbv[ind_f + nbF*ind_m];

if(VERBOSE){Rprintf("Dyna sp");}
if (nbE>0) {
//Rprintf("H1\n");fichier << "H1" << endl;

    for (int e = 0 ; e < nbE ; e++) {
//Rprintf("H2\n");fichier << "H2" << endl;

                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                            PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                            nbI = length(getListElement(elmt, "modI"));

                            if (Svec[e]==0){
                                PROTECT(v_wL_ei = getListElement(elmt, "wL_i"));
                                PROTECT(v_wD_ei = getListElement(elmt, "wD_i"));
                                PROTECT(v_d_efmit = getListElement(elmt, "d_i"));
                                PROTECT(v_doth_eit = getListElement(elmt, "doth_i"));

                                PROTECT(v_F_efmit = getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_Z_eit = getListElement(out_Z_eit , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_N_eit = getListElement(out_N_eit , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_B_et = getListElement(out_B_et , CHAR(STRING_ELT(sppList,e))));

                                PROTECT(dimCst_wL_ei = getAttrib(v_wL_ei, install("DimCst")));
                                PROTECT(dimCst_wD_ei = getAttrib(v_wD_ei, install("DimCst")));
                                PROTECT(dimCst_d_efmit = getAttrib(v_d_efmit, install("DimCst")));
                                PROTECT(dimCst_F_efmit = getAttrib(v_F_efmit, install("DimCst")));
                                PROTECT(dimCst_N_eit = getAttrib(v_N_eit, install("DimCst")));
                                PROTECT(dimCst_Z_eit = getAttrib(v_Z_eit, install("DimCst")));

                            } else {
                                PROTECT(v_wL_ei_G1 = getListElement(elmt, "wL_i_G1")); //PrintValue(v_wL_ei_G1);
                                PROTECT(v_wL_ei_G2 = getListElement(elmt, "wL_i_G2"));
                                PROTECT(v_wD_ei_G1 = getListElement(elmt, "wD_i_G1"));
                                PROTECT(v_wD_ei_G2 = getListElement(elmt, "wD_i_G2"));
                                PROTECT(v_d_efmit_G1 = getListElement(elmt, "d_i_G1"));
                                PROTECT(v_d_efmit_G2 = getListElement(elmt, "d_i_G2"));
                                PROTECT(v_doth_eit_G1 = getListElement(elmt, "doth_i_G1"));
                                PROTECT(v_doth_eit_G2 = getListElement(elmt, "doth_i_G2"));

                                PROTECT(v_F_efmit_G1 = getListElement(out_F_fmi_G1, CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_Z_eit_G1 = getListElement(out_Z_eit_G1 , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_N_eit_G1 = getListElement(out_N_eit_G1 , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_F_efmit_G2 = getListElement(out_F_fmi_G2, CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_Z_eit_G2 = getListElement(out_Z_eit_G2 , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(v_N_eit_G2 = getListElement(out_N_eit_G2 , CHAR(STRING_ELT(sppList,e))));

                                PROTECT(v_B_et = getListElement(out_B_et , CHAR(STRING_ELT(sppList,e))));
                                PROTECT(dimCst_wL_ei = getAttrib(v_wL_ei_G1, install("DimCst")));
                                PROTECT(dimCst_wD_ei = getAttrib(v_wD_ei_G1, install("DimCst")));
                                PROTECT(dimCst_d_efmit = getAttrib(v_d_efmit_G1, install("DimCst")));
                                PROTECT(dimCst_F_efmit = getAttrib(v_F_efmit_G1, install("DimCst")));
                                PROTECT(dimCst_N_eit = getAttrib(v_N_eit_G1, install("DimCst")));
                                PROTECT(dimCst_Z_eit = getAttrib(v_Z_eit_G1, install("DimCst")));
                            }

//Rprintf("H3\n");fichier << "H3" << endl;
                            //tests sur les dimensions :
                            dim_d_efmit = INTEGER(dimCst_d_efmit);
                            if (((dim_d_efmit[0]!=0) & (dim_d_efmit[0]!=nbF)) | ((dim_d_efmit[1]!=0) & (dim_d_efmit[1]!=nbM)) |
                                ((dim_d_efmit[2]!=0) & (dim_d_efmit[2]!=nbI)) | ((dim_d_efmit[3]!=0) & (dim_d_efmit[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in d_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_wL_ei = INTEGER(dimCst_wL_ei);
                            if ((dim_wL_ei[0]!=0) | (dim_wL_ei[1]!=0) |
                                ((dim_wL_ei[2]!=0) & (dim_wL_ei[2]!=nbI)) | (dim_wL_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wL_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_wD_ei = INTEGER(dimCst_wD_ei);
                            if ((dim_wD_ei[0]!=0) | (dim_wD_ei[1]!=0) |
                                ((dim_wD_ei[2]!=0) & (dim_wD_ei[2]!=nbI)) | (dim_wD_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wD_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_F_efmit = INTEGER(dimCst_F_efmit);
                            if (((dim_F_efmit[0]!=0) & (dim_F_efmit[0]!=nbF)) | ((dim_F_efmit[1]!=0) & (dim_F_efmit[1]!=nbM)) |
                                ((dim_F_efmit[2]!=0) & (dim_F_efmit[2]!=nbI)) | ((dim_F_efmit[3]!=0) & (dim_F_efmit[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in F_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_N_eit = INTEGER(dimCst_N_eit);
                            if ((dim_N_eit[0]!=0) | (dim_N_eit[1]!=0) |
                                ((dim_N_eit[2]!=0) & (dim_N_eit[2]!=nbI)) | ((dim_N_eit[3]!=0) & (dim_N_eit[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in N_eit element. Check .ini biological parameters files !!\n");
                            }

                            dim_Z_eit = INTEGER(dimCst_Z_eit);
                            if ((dim_Z_eit[0]!=0) | (dim_Z_eit[1]!=0) |
                                ((dim_Z_eit[2]!=0) & (dim_Z_eit[2]!=nbI)) | ((dim_Z_eit[3]!=0) & (dim_Z_eit[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in Z_eit element. Check .ini biological parameters files !!\n");
                            }

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de C_efmit
                            //---------
//Rprintf("Calcul.C_efmit\n");fichier << "Calcul.C_efmit" << endl;

                            //on d�termine l'attribut Dimension de C_efmit
                            PROTECT(dimCst = allocVector(INTSXP, 4));
                            dimC = INTEGER(dimCst);
                            dimC[0] = dim_F_efmit[0] ; dimC[1] = dim_F_efmit[1] ; dimC[2] = dim_F_efmit[2];
                            dimC[3] = imax2(dim_N_eit[3] , dim_F_efmit[3]);

                            int count = 0, prod = 1, count2 = 0, count3 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    count++;
                                    prod = prod * dimC[k];
                                }

                            }

                            PROTECT(Dim = allocVector(INTSXP, count));
                            dim = INTEGER(Dim);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    dim[count2] = dimC[k];
                                    count2++;
                                    }
                            }

                            PROTECT(dimCst2 = allocVector(INTSXP, 4));
                            dimcst2 = INTEGER(dimCst2); dimcst2[0] = 0; dimcst2[1] = 0; dimcst2[2] = nbI; dimcst2[3] = nbT;
                            PROTECT(Dim2 = allocVector(INTSXP, 2));
                            dim2 = INTEGER(Dim2); dim2[0] = nbI; dim2[1] = nbT;
//Rprintf("H4\n");fichier << "H4" << endl;

                    if (ind_t==0){
//Rprintf("H4.1\n");fichier << "H4.1" << endl;
                            //on cr�e le tableau r�sultat pour l'esp�ce en question
                            if (Svec[e]==0){
                                PROTECT(ans_C_efmit = NEW_NUMERIC(prod));
                                setAttrib(ans_C_efmit, R_DimSymbol, Dim);

                                PROTECT(ans_C_eit = NEW_NUMERIC(nbI*nbT));
                                setAttrib(ans_C_eit, R_DimSymbol, Dim2);

                                rans_C_efmit = REAL(ans_C_efmit);
                                rans_C_eit = REAL(ans_C_eit);

                            } else{
                                PROTECT(ans_C_efmit_G1 = NEW_NUMERIC(prod));
                                setAttrib(ans_C_efmit_G1, R_DimSymbol, Dim);
                                PROTECT(ans_C_efmit_G2 = NEW_NUMERIC(prod));
                                setAttrib(ans_C_efmit_G2, R_DimSymbol, Dim);

                                PROTECT(ans_C_eit_G1 = NEW_NUMERIC(nbI*nbT));
                                setAttrib(ans_C_eit_G1, R_DimSymbol, Dim2);
                                PROTECT(ans_C_eit_G2 = NEW_NUMERIC(nbI*nbT));
                                setAttrib(ans_C_eit_G2, R_DimSymbol, Dim2);

                                rans_C_efmit_G1 = REAL(ans_C_efmit_G1);
                                rans_C_efmit_G2 = REAL(ans_C_efmit_G2);
                                rans_C_eit_G1 = REAL(ans_C_eit_G1);
                                rans_C_eit_G2 = REAL(ans_C_eit_G2);
                            }

                            PROTECT(dimnames = allocVector(VECSXP,count));
                            if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                            if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                            if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                            if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}

                            PROTECT(dimnames2 = allocVector(VECSXP,2));
                            SET_VECTOR_ELT(dimnames2, 0, intAge);
                            SET_VECTOR_ELT(dimnames2, 1, times);

//Rprintf("H4.2\n");fichier << "H4.2" << endl;
                    } else {
                        if (Svec[e]==0){

                            rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit, e));
                            rans_C_eit = REAL(VECTOR_ELT(out_C_eit, e));

                        } else {
                            rans_C_efmit_G1 = REAL(VECTOR_ELT(out_C_efmit_G1, e));
                            rans_C_efmit_G2 = REAL(VECTOR_ELT(out_C_efmit_G2, e));
                            rans_C_eit_G1 = REAL(VECTOR_ELT(out_C_eit_G1, e));
                            rans_C_eit_G2 = REAL(VECTOR_ELT(out_C_eit_G2, e));
                        }
//Rprintf("H4.3\n");fichier << "H4.3" << endl;
                    }

                    if (Svec[e]==0){
                       r_F_efmit = REAL(v_F_efmit);
                       r_N_eit = REAL(v_N_eit);
                       r_Z_eit = REAL(v_Z_eit);
                       r_B_et = REAL(v_B_et);
                       r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44)); //Rprintf("Dans EVAR (l.9245), Fothi = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));

                    } else {
                       r_F_efmit_G1 = REAL(v_F_efmit_G1);
                       r_F_efmit_G2 = REAL(v_F_efmit_G2);
                       r_N_eit_G1 = REAL(v_N_eit_G1);
                       r_N_eit_G2 = REAL(v_N_eit_G2);
                       r_Z_eit_G1 = REAL(v_Z_eit_G1);
                       r_Z_eit_G2 = REAL(v_Z_eit_G2);
                       r_B_et = REAL(v_B_et);
                       r_Foth_i_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 224));
                       r_Foth_i_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 225));
                    }



//Rprintf("H4.4\n");fichier << "H4.4" << endl;
                            //facteurs des indices
                            PROTECT(cFACT1 = iDim(dimC));
                            PROTECT(cFACT2 = iDim(dim_F_efmit));
                            PROTECT(cFACT3 = iDim(dim_N_eit));
                            PROTECT(cFACT4 = iDim(dim_Z_eit));

                            int *fact1_C = INTEGER(cFACT1);
                            int *fact2_C = INTEGER(cFACT2);
                            int *fact3_C = INTEGER(cFACT3);
                            int *fact4_C = INTEGER(cFACT4);


//Rprintf("H4.5\n");fichier << "H4.5" << endl;
                            //�quation
                            if ((Qvec[e]==1) & (Svec[e]==0)) {
//Rprintf("H4.6\n");fichier << "H4.6" << endl;
                                double *r_F_fmi_S1M1 = REAL(getListElement(out_F_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M2 = REAL(getListElement(out_F_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M3 = REAL(getListElement(out_F_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M4 = REAL(getListElement(out_F_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M1 = REAL(getListElement(out_F_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M2 = REAL(getListElement(out_F_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M3 = REAL(getListElement(out_F_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M4 = REAL(getListElement(out_F_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M1 = REAL(getListElement(out_F_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M2 = REAL(getListElement(out_F_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M3 = REAL(getListElement(out_F_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M4 = REAL(getListElement(out_F_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M1 = REAL(getListElement(out_F_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M2 = REAL(getListElement(out_F_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M3 = REAL(getListElement(out_F_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M4 = REAL(getListElement(out_F_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));
//Rprintf("H4.7\n");fichier << "H4.7" << endl;//PrintValue(out_Z_eit_S1M1);
                                double *r_Z_eit_S1M1 = REAL(getListElement(out_Z_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.1\n");
                                double *r_Z_eit_S1M2 = REAL(getListElement(out_Z_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.2\n");
                                double *r_Z_eit_S1M3 = REAL(getListElement(out_Z_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.3\n");
                                double *r_Z_eit_S1M4 = REAL(getListElement(out_Z_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.4\n");
                                double *r_Z_eit_S2M1 = REAL(getListElement(out_Z_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M2 = REAL(getListElement(out_Z_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M3 = REAL(getListElement(out_Z_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M4 = REAL(getListElement(out_Z_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M1 = REAL(getListElement(out_Z_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M2 = REAL(getListElement(out_Z_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M3 = REAL(getListElement(out_Z_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M4 = REAL(getListElement(out_Z_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M1 = REAL(getListElement(out_Z_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M2 = REAL(getListElement(out_Z_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M3 = REAL(getListElement(out_Z_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M4 = REAL(getListElement(out_Z_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));
//Rprintf("H4.8\n");fichier << "H4.8" << endl;
                                double *r_N_eit_S1M1 = REAL(getListElement(out_N_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M2 = REAL(getListElement(out_N_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M3 = REAL(getListElement(out_N_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M4 = REAL(getListElement(out_N_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M1 = REAL(getListElement(out_N_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M2 = REAL(getListElement(out_N_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M3 = REAL(getListElement(out_N_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M4 = REAL(getListElement(out_N_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M1 = REAL(getListElement(out_N_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M2 = REAL(getListElement(out_N_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M3 = REAL(getListElement(out_N_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M4 = REAL(getListElement(out_N_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M1 = REAL(getListElement(out_N_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M2 = REAL(getListElement(out_N_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M3 = REAL(getListElement(out_N_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M4 = REAL(getListElement(out_N_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));
//Rprintf("H5aaa\n");fichier << "H5aaa" << endl;
                                double  *r_M_ei = REAL(getListElement(elmt, "M_i"));//Rprintf("H5.1\n");
                                int *fact3_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 14));

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    rans_C_eit[ind_i + ind_t*nbI] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));
//Rprintf("H6\n");fichier << "H6" << endl;

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                      (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));

                                }

                                }


                            } else if ((Qvec[e]==0) & (Svec[e]==0)) {
//Rprintf("H7\n");fichier << "H7" << endl;
                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                      rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    double temp = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }

                                    rans_C_eit[ind_i + ind_t*nbI] =
                                        (temp + r_Foth_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                }

                            } else if ((Qvec[e]==0) & (Svec[e]==1)) {
//Rprintf("H7.1\n");fichier << "H7.1" << endl;
                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                      rans_C_efmit_G1[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit_G1[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                       rans_C_efmit_G2[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit_G2[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];
                                }

//Rprintf("H7.2\n");fichier << "H7.2" << endl;


                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    double tempG1 = 0.0;
                                    double tempG2 = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        tempG1 = tempG1 + r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                        if (!ISNA(r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        tempG2 = tempG2 + r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }


                                    rans_C_eit_G1[ind_i + ind_t*nbI] =
                                        (tempG1 + r_Foth_i_G1[ind_i + ind_t*nbI]) * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                    rans_C_eit_G2[ind_i + ind_t*nbI] =
                                        (tempG2 + r_Foth_i_G2[ind_i + ind_t*nbI]) * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                }

                            }


//Rprintf("H7.3\n");fichier << "H7.3" << endl;
//PrintValue(ans_C_eit_G1);


                    if (ind_t==0){
                            if (Svec[e]==0){
                                setAttrib(ans_C_efmit, R_DimNamesSymbol, dimnames);
                                setAttrib(ans_C_efmit, install("DimCst"), dimCst);

                                setAttrib(ans_C_eit, R_DimNamesSymbol, dimnames2);
                                setAttrib(ans_C_eit, install("DimCst"), dimCst2);

                                SET_VECTOR_ELT(out_C_efmit, e, ans_C_efmit);
                                SET_VECTOR_ELT(out_C_eit, e, ans_C_eit);

                            } else{
                                setAttrib(ans_C_efmit_G1, R_DimNamesSymbol, dimnames);
                                setAttrib(ans_C_efmit_G1, install("DimCst"), dimCst);
                                setAttrib(ans_C_efmit_G2, R_DimNamesSymbol, dimnames);
                                setAttrib(ans_C_efmit_G2, install("DimCst"), dimCst);

                                setAttrib(ans_C_eit_G1, R_DimNamesSymbol, dimnames2);
                                setAttrib(ans_C_eit_G1, install("DimCst"), dimCst2);
                                setAttrib(ans_C_eit_G2, R_DimNamesSymbol, dimnames2);
                                setAttrib(ans_C_eit_G2, install("DimCst"), dimCst2);

                                SET_VECTOR_ELT(out_C_efmit_G1, e, ans_C_efmit_G1);
                                SET_VECTOR_ELT(out_C_eit_G1, e, ans_C_eit_G1);
                                SET_VECTOR_ELT(out_C_efmit_G2, e, ans_C_efmit_G2);
                                SET_VECTOR_ELT(out_C_eit_G2, e, ans_C_eit_G2);
                            }


                    }
//Rprintf("H7.4\n");fichier << "H7.4" << endl;

                    if (Svec[e]==0){

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 23, v_F_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 24, v_N_eit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 25, v_Z_eit);
                    } else {
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 228, v_F_efmit_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 240, v_N_eit_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 242, v_Z_eit_G1);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 229, v_F_efmit_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 241, v_N_eit_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 243, v_Z_eit_G2);
                    }
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 26, cFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 27, cFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 28, cFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 29, cFACT4);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de Y_efmit
                            //---------
//Rprintf("Calcul.Y_efmit\n");fichier << "Calcul.Y_efmit" << endl;

                        //on consid�re les dimensions de C, Y, D et L homog�nes sur tout le module --> pas besoin de les red�finir

                    if (ind_t==0){

                            PROTECT(ans_Y_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_Y_efmit, R_DimSymbol, Dim);

                            rans_Y_efmit = REAL(ans_Y_efmit);

                            PROTECT(ans_Y_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_Y_eit, R_DimSymbol, Dim2);

                            rans_Y_eit = REAL(ans_Y_eit);

                            PROTECT(ans_L_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_L_eit, R_DimSymbol, Dim2);

                            rans_L_eit = REAL(ans_L_eit);


                    } else {

                            rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e));
                            rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e));
                            rans_L_eit = REAL(VECTOR_ELT(out_L_eit,e));

                    }


                            if (Svec[e]==0){
                                    r_wL_ei = REAL(v_wL_ei);
                            } else {
                                r_wL_ei_G1 = REAL(v_wL_ei_G1);
                                r_wL_ei_G2 = REAL(v_wL_ei_G2);
                            }

                                                    //facteurs des indices
                            PROTECT(cFACT5 = iDim(dim_wL_ei));

                            int *fact5_C = INTEGER(cFACT5);

                            //�quation

                    if ((Qvec[e]==0) & (Svec[e]==0)) {

                      if (nbI>1) { //ajout SPiCT

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                  rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                                    rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                            }

                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                  rans_Y_eit[ind_i + ind_t*nbI] =
                                    r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_eit[ind_i + ind_t*nbI] / 1000;

                                  rans_L_eit[ind_i + ind_t*nbI] = NA_REAL;
                           }

                    } else { //SPiCT

                            //double *Bspict = REAL(VECTOR_ELT(intermBIOMspict, e));
                            // on peut sommer avant d'appliquer � F puisque F est suppos� (pour le moment) constant sur l'ensemble de l'ann�e N
                            //double Btemp = 0.0;
                            //for (int ii = ind_t*16 ; ii < (ind_t*16 + 16) ; ii++) Btemp = Btemp + Bspict[ii]/16;

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                                      rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_t*fact1_C[3]] =
                                        r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]] * r_B_et[ind_t*fact3_C[3]];

                            double temp = 0.0;

                            for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                            for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_t*fact2_C[3]];

                            }

                            rans_Y_eit[0 + ind_t*1] = (temp + r_Foth_i[0 + ind_t*1]) * r_B_et[ind_t*fact3_C[3]];
                            rans_L_eit[0 + ind_t*1] = NA_REAL;

                    }
                    } else if ((Qvec[e]==0) & (Svec[e]==1)) {


                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                  rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wL_ei_G1[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit_G1[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                  r_wL_ei_G2[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit_G2[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                                    rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                            }

                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                  rans_Y_eit[ind_i + ind_t*nbI] =
                                    r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_eit_G1[ind_i + ind_t*nbI] / 1000 +

                                    r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_eit_G2[ind_i + ind_t*nbI] / 1000 ;

                                  rans_L_eit[ind_i + ind_t*nbI] = NA_REAL;
                           }
                    }

                    if (ind_t==0) {

                            setAttrib(ans_Y_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_Y_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_Y_efmit, e, ans_Y_efmit);

                            setAttrib(ans_Y_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_Y_eit, install("DimCst"), dimCst2);

                            SET_VECTOR_ELT(out_Y_eit, e, ans_Y_eit);

                            setAttrib(ans_L_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_L_eit, install("DimCst"), dimCst2);

                            SET_VECTOR_ELT(out_L_eit, e, ans_L_eit);


                    }

                    if (Svec[e]==0) {
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 30, v_wL_ei);
                    } else {
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 244, v_wL_ei_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 245, v_wL_ei_G2);
                    }
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 31, cFACT5);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de D_efmit
                            //---------
//Rprintf("Calcul.D_efmit\n");fichier << "Calcul.D_efmit" << endl;

double *r_dd1_efm = REAL(getListElement(elmt, "dd1_f_m_e"));
double *r_dd2_efm = REAL(getListElement(elmt, "dd2_f_m_e"));
double *r_OD_e = REAL(getListElement(elmt, "OD_e"));



                    if (ind_t==0) {


                                PROTECT(ans_D_efmit = NEW_NUMERIC(prod));
                                setAttrib(ans_D_efmit, R_DimSymbol, Dim);

                                rans_D_efmit = REAL(ans_D_efmit);

                                PROTECT(ans_LD_efmit = NEW_NUMERIC(prod));
                                setAttrib(ans_LD_efmit, R_DimSymbol, Dim);

                                rans_LD_efmit = REAL(ans_LD_efmit);

                                PROTECT(ans_DD_efmit = NEW_NUMERIC(prod));
                                setAttrib(ans_DD_efmit, R_DimSymbol, Dim);

                                rans_DD_efmit = REAL(ans_DD_efmit);


                       //if (!(r_OD_e[0]>0.5 & r_OD_e[0]<=(ind_t+1)) & (activeQR!=0)) { //over quota discards sera impl�ment�

                            PROTECT(ans_oqD_eft = NEW_NUMERIC(nbF*nbT));//Rprintf("AA1");
                            PROTECT(ans_oqD_et = NEW_NUMERIC(nbT));

                            PROTECT(dimCstOQ_ft = allocVector(INTSXP, 2));//Rprintf("AA2");
                            dimOQ_ft = INTEGER(dimCstOQ_ft);
                            dimOQ_ft[0] = nbF; dimOQ_ft[1] = nbT;

                            PROTECT(dimCstOQ_t = allocVector(INTSXP, 1));//Rprintf("AA3");
                            dimOQ_t = INTEGER(dimCstOQ_t);
                            dimOQ_t[0] = nbT;

                            setAttrib(ans_oqD_eft, R_DimSymbol, dimCstOQ_ft);//Rprintf("AA31");
                            setAttrib(ans_oqD_et, R_DimSymbol, dimCstOQ_t);//Rprintf("AA32");

                            PROTECT(dimnames_oqD_eft = allocVector(VECSXP,2));//Rprintf("AA4");
                            SET_VECTOR_ELT(dimnames_oqD_eft, 0, fleetList);
                            SET_VECTOR_ELT(dimnames_oqD_eft, 1, times);

                            PROTECT(dimnames_oqD_et = allocVector(VECSXP,1));
                            SET_VECTOR_ELT(dimnames_oqD_et, 0, times);

                            setAttrib(ans_oqD_eft, R_DimNamesSymbol, dimnames_oqD_eft);
                            setAttrib(ans_oqD_et, R_DimNamesSymbol, dimnames_oqD_et); //Rprintf("AA5");

                            rans_oqD_eft = REAL(ans_oqD_eft);
                            for (int tt = 0; tt<nbF*nbT; tt++) rans_oqD_eft[tt] = 0.0;
                            rans_oqD_et = REAL(ans_oqD_et); //Rprintf("BB");
                            for (int tt = 0; tt<nbT; tt++) rans_oqD_et[tt] = 0.0;

                        //}

                    } else {

                            rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e));
                            rans_DD_efmit = REAL(VECTOR_ELT(out_DD_efmi,e));
                            rans_LD_efmit = REAL(VECTOR_ELT(out_LD_efmi,e));

                        //if (!(r_OD_e[0]>0.5 & r_OD_e[0]<=(ind_t+1)) & (activeQR!=0)) {//Rprintf("CC");

                            rans_oqD_eft = REAL(VECTOR_ELT(out_oqD_eft,e));
                            rans_oqD_et = REAL(VECTOR_ELT(out_oqD_et,e));//Rprintf("DD");


                        //}

                    }

               //facteurs des indices
                            PROTECT(cFACT6 = iDim(dim_d_efmit));
                            PROTECT(cFACT7 = iDim(dim_wD_ei));

                            int *fact6_C = INTEGER(cFACT6);
                            int *fact7_C = INTEGER(cFACT7);

                if ((Qvec[e]==0) & (Svec[e]==0)) {

                    r_wD_ei = REAL(v_wD_ei);
                    r_d_efmit = REAL(v_d_efmit);
                    doth_eit = REAL(v_doth_eit);
                            //�quation : 2 mani�res de calculer selon la disponibilit� de wD_i


                    if (all_is_na(v_wD_ei)) { //1�re m�thode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                            }

                            //Loth_eit
if (nbI>1) {
                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit[ind_i]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000;
                            }
} else {

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)

                                rans_L_eit[0 + ind_t*1] =
                                        r_Foth_i[0 + ind_t*1] * r_B_et[ind_t*fact3_C[3]];


                            } else { //pas d'OD

                                rans_L_eit[0 + ind_t*1] =
                                        r_Foth_i[0 + ind_t*1] * r_B_et[ind_t*fact3_C[3]] * (1-doth_eit[0]);

                            }

}

                    } else {                 //2�me m�thode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                            }


                            //Loth_eit : pas d'exemption pour les autres si OD

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit, pas de rejet car pas d'exemption

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit[ind_i] * r_wD_ei[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000;
                            }


                    }



                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd1_efm[ind_f + nbF*ind_m]))) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }

                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                         fmin2( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                        if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd2_efm[ind_f + nbF*ind_m]))){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) ,  1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; // il reste ensuite � integrer Lefmit pour contituer Ltot_i

                    }


                   }  else if ((Qvec[e]==0) & (Svec[e]==1)) {

                    r_wD_ei_G1 = REAL(v_wD_ei_G1);
                    r_wD_ei_G2 = REAL(v_wD_ei_G2);
                    r_d_efmit_G1 = REAL(v_d_efmit_G1);
                    r_d_efmit_G2 = REAL(v_d_efmit_G2);
                    doth_eit_G1 = REAL(v_doth_eit_G1);
                    doth_eit_G2 = REAL(v_doth_eit_G2);
                            //�quation : 2 mani�res de calculer selon la disponibilit� de wD_i


                    if (all_is_na(v_wD_ei_G1) | all_is_na(v_wD_ei_G2)) { //1�re m�thode
//Rprintf("L1\n");
                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                if (ISNA(r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (ISNA(r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    r_wL_ei_G1[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit_G1[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                    r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    r_wL_ei_G2[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit_G2[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 ;
                            }

                            //Loth_eit

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)
//Rprintf("L1.1\n");
                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000 +

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000;


                            } else { //pas d'OD
//Rprintf("L1.2\n");
                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit_G1[ind_i]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000+

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit_G2[ind_i]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]]/1000;
                            }


                    } else {                 //2�me m�thode
//Rprintf("L2\n");
                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (ISNA(r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wD_ei_G1[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_C_efmit_G1[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                    r_wD_ei_G2[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_C_efmit_G2[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                            }


                            //Loth_eit : pas d'exemption pour les autres si OD

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit, pas de rejet car pas d'exemption
//Rprintf("L2.1\n")  ;
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000 +

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;


                            } else { //pas d'OD
//Rprintf("L2.2\n");
                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit_G1[ind_i] * r_wD_ei_G1[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000 +

                                         (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit_G2[ind_i] * r_wD_ei_G2[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000;
                            }


                    }

//Rprintf("e = %s\n",STRING_ELT(sppList,e));
//PrintValue(ans_L_eit);

                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd1_efm[ind_f + nbF*ind_m]))) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }

                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                         fmin2( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                        if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd2_efm[ind_f + nbF*ind_m]))){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) ,  1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; // il reste ensuite � integrer Lefmit pour contituer Ltot_i

                    }

//Rprintf("e = %s\n",STRING_ELT(sppList,e));
//PrintValue(ans_L_eit);

                   }


                    if (ind_t==0) {

                            setAttrib(ans_D_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_D_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_D_efmit, e, ans_D_efmit);

                            setAttrib(ans_LD_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_LD_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_LD_efmi, e, ans_LD_efmit);

                            setAttrib(ans_DD_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_DD_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_DD_efmi, e, ans_DD_efmit);

                            SET_VECTOR_ELT(out_oqD_eft, e, ans_oqD_eft);
                            SET_VECTOR_ELT(out_oqD_et, e, ans_oqD_et);

                    }

                        if (Svec[e]==0){
                            SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 32, v_wD_ei);
                            SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 33, v_d_efmit);
                        } else {
                            // a voir si necessaire de rajouter wD_G1, wD_G2, d_efmit_G1, d_efmit_G2
                        }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 34, cFACT6);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 35, cFACT7);

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de L_efmit
                            //---------

 // Rprintf("Calcul.L_efmit\n");fichier << "Calcul.L_efmit" << endl;
                    if (ind_t==0) {

                            PROTECT(ans_L_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_L_efmit, R_DimSymbol, Dim);

                            rans_L_efmit = REAL(ans_L_efmit);

                    } else {

                            rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e));

                    }

//fichier << "Calcul.L_efmit.1" << endl;

                            //�quation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                  rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                  if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                  rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                    rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; //Ltot constitu�

                            }
//fichier << "Calcul.L_efmit.2" << endl;


                    if (ind_t==0) {

                            setAttrib(ans_L_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_L_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_L_efmit, e, ans_L_efmit);
                            SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));
//Rprintf("K12\n");fichier << "K12" << endl;
                            //UNPROTECT(11);

                            //*if (!(r_OD_e[0]>0.5 & r_OD_e[0]<=(ind_t+1)) & (activeQR!=0)) */UNPROTECT(6); //over quota discards

                    }
//fichier << "Calcul.L_efmit.3" << endl;



                     if ((Qvec[e]==1) & (Svec[e]==0)) {
//Rprintf("H4.6\n");fichier << "H4.6" << endl;

                                double *r_Z_eit_S1M1 = REAL(getListElement(out_Z_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.1\n");fichier << "H4.7.1" << endl;
                                double *r_Z_eit_S1M2 = REAL(getListElement(out_Z_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.2\n");fichier << "H4.7.2" << endl;
                                double *r_Z_eit_S1M3 = REAL(getListElement(out_Z_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.3\n");fichier << "H4.7.3" << endl;
                                double *r_Z_eit_S1M4 = REAL(getListElement(out_Z_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));//Rprintf("H4.7.4\n");fichier << "H4.7.4" << endl;
                                double *r_Z_eit_S2M1 = REAL(getListElement(out_Z_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M2 = REAL(getListElement(out_Z_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M3 = REAL(getListElement(out_Z_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M4 = REAL(getListElement(out_Z_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M1 = REAL(getListElement(out_Z_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M2 = REAL(getListElement(out_Z_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M3 = REAL(getListElement(out_Z_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M4 = REAL(getListElement(out_Z_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M1 = REAL(getListElement(out_Z_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M2 = REAL(getListElement(out_Z_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M3 = REAL(getListElement(out_Z_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M4 = REAL(getListElement(out_Z_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));
//Rprintf("H4.8\n");fichier << "H4.8" << endl;
                                double *r_N_eit_S1M1 = REAL(getListElement(out_N_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M2 = REAL(getListElement(out_N_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M3 = REAL(getListElement(out_N_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M4 = REAL(getListElement(out_N_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M1 = REAL(getListElement(out_N_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M2 = REAL(getListElement(out_N_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M3 = REAL(getListElement(out_N_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M4 = REAL(getListElement(out_N_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M1 = REAL(getListElement(out_N_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M2 = REAL(getListElement(out_N_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M3 = REAL(getListElement(out_N_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M4 = REAL(getListElement(out_N_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M1 = REAL(getListElement(out_N_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M2 = REAL(getListElement(out_N_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M3 = REAL(getListElement(out_N_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M4 = REAL(getListElement(out_N_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));
//Rprintf("H5bbb\n");fichier << "H5bbb" << endl;
                                //double  *r_M_ei = REAL(getListElement(elmt, "M_i"));
//Rprintf("H5bbb2\n");fichier << "H5bbb2" << endl;//PrintValue(out_FRWT_fmi_S1M1);
                                double *r_FRWT_fmi_S1M1 = REAL(getListElement(out_FRWT_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M2 = REAL(getListElement(out_FRWT_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M3 = REAL(getListElement(out_FRWT_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M4 = REAL(getListElement(out_FRWT_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M1 = REAL(getListElement(out_FRWT_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M2 = REAL(getListElement(out_FRWT_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M3 = REAL(getListElement(out_FRWT_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M4 = REAL(getListElement(out_FRWT_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M1 = REAL(getListElement(out_FRWT_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M2 = REAL(getListElement(out_FRWT_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M3 = REAL(getListElement(out_FRWT_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M4 = REAL(getListElement(out_FRWT_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M1 = REAL(getListElement(out_FRWT_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M2 = REAL(getListElement(out_FRWT_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M3 = REAL(getListElement(out_FRWT_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M4 = REAL(getListElement(out_FRWT_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));
//Rprintf("H4.7\n");fichier << "H4.7" << endl;//PrintValue(out_Z_eit_S1M1);
                                double *r_FDWT_fmi_S1M1 = REAL(getListElement(out_FDWT_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M2 = REAL(getListElement(out_FDWT_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M3 = REAL(getListElement(out_FDWT_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M4 = REAL(getListElement(out_FDWT_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M1 = REAL(getListElement(out_FDWT_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M2 = REAL(getListElement(out_FDWT_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M3 = REAL(getListElement(out_FDWT_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M4 = REAL(getListElement(out_FDWT_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M1 = REAL(getListElement(out_FDWT_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M2 = REAL(getListElement(out_FDWT_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M3 = REAL(getListElement(out_FDWT_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M4 = REAL(getListElement(out_FDWT_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M1 = REAL(getListElement(out_FDWT_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M2 = REAL(getListElement(out_FDWT_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M3 = REAL(getListElement(out_FDWT_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M4 = REAL(getListElement(out_FDWT_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));

                                double *r_FRWToth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 176));
                                double *r_FRWToth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 177));
                                double *r_FRWToth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 178));
                                double *r_FRWToth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 179));
                                double *r_FRWToth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 180));
                                double *r_FRWToth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 181));
                                double *r_FRWToth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 182));
                                double *r_FRWToth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 183));
                                double *r_FRWToth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 184));
                                double *r_FRWToth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 185));
                                double *r_FRWToth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 186));
                                double *r_FRWToth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 187));
                                double *r_FRWToth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 188));
                                double *r_FRWToth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 189));
                                double *r_FRWToth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 190));
                                double *r_FRWToth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 191));

                                double *r_FDWToth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 208));
                                double *r_FDWToth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 209));
                                double *r_FDWToth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 210));
                                double *r_FDWToth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 211));
                                double *r_FDWToth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 212));
                                double *r_FDWToth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 213));
                                double *r_FDWToth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 214));
                                double *r_FDWToth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 215));
                                double *r_FDWToth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 216));
                                double *r_FDWToth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 217));
                                double *r_FDWToth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 218));
                                double *r_FDWToth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 219));
                                double *r_FDWToth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 220));
                                double *r_FDWToth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 221));
                                double *r_FDWToth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 222));
                                double *r_FDWToth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 223));


                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                                    rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                      (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));





                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                      (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                                rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];


                                rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                }




                  for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd1_efm[ind_f + nbF*ind_m]))) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2 ( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0 );
                        if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                            rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd2_efm[ind_f + nbF*ind_m]))){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                    fmin2 ( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0 ) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];



                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                     double temp_S1M1 = 0.0, temp_S1M2 = 0.0, temp_S1M3 = 0.0, temp_S1M4 = 0.0,
                                            temp_S2M1 = 0.0, temp_S2M2 = 0.0, temp_S2M3 = 0.0, temp_S2M4 = 0.0,
                                            temp_S3M1 = 0.0, temp_S3M2 = 0.0, temp_S3M3 = 0.0, temp_S3M4 = 0.0,
                                            temp_S4M1 = 0.0, temp_S4M2 = 0.0, temp_S4M3 = 0.0, temp_S4M4 = 0.0,
                                            temp2_S1M1 = 0.0, temp2_S1M2 = 0.0, temp2_S1M3 = 0.0, temp2_S1M4 = 0.0,
                                            temp2_S2M1 = 0.0, temp2_S2M2 = 0.0, temp2_S2M3 = 0.0, temp2_S2M4 = 0.0,
                                            temp2_S3M1 = 0.0, temp2_S3M2 = 0.0, temp2_S3M3 = 0.0, temp2_S3M4 = 0.0,
                                            temp2_S4M1 = 0.0, temp2_S4M2 = 0.0, temp2_S4M3 = 0.0, temp2_S4M4 = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M1 = temp_S1M1 + r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M2 = temp_S1M2 + r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M3 = temp_S1M3 + r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M4 = temp_S1M4 + r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M1 = temp_S2M1 + r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M2 = temp_S2M2 + r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M3 = temp_S2M3 + r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M4 = temp_S2M4 + r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M1 = temp_S3M1 + r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M2 = temp_S3M2 + r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M3 = temp_S3M3 + r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M4 = temp_S3M4 + r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M1 = temp_S4M1 + r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M2 = temp_S4M2 + r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M3 = temp_S4M3 + r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M4 = temp_S4M4 + r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];


                                        if (!ISNA(r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M1 = temp2_S1M1 + r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M2 = temp2_S1M2 + r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M3 = temp2_S1M3 + r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M4 = temp2_S1M4 + r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M1 = temp2_S2M1 + r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M2 = temp2_S2M2 + r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M3 = temp2_S2M3 + r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M4 = temp2_S2M4 + r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M1 = temp2_S3M1 + r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M2 = temp2_S3M2 + r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M3 = temp2_S3M3 + r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M4 = temp2_S3M4 + r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M1 = temp2_S4M1 + r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M2 = temp2_S4M2 + r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M3 = temp2_S4M3 + r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M4 = temp2_S4M4 + r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }

                                    temp_S1M1 = temp_S1M1 + r_FRWToth_it_S1M1[ind_i+ind_t*nbI];
                                    temp_S1M2 = temp_S1M2 + r_FRWToth_it_S1M2[ind_i+ind_t*nbI];
                                    temp_S1M3 = temp_S1M3 + r_FRWToth_it_S1M3[ind_i+ind_t*nbI];
                                    temp_S1M4 = temp_S1M4 + r_FRWToth_it_S1M4[ind_i+ind_t*nbI];
                                    temp_S2M1 = temp_S2M1 + r_FRWToth_it_S2M1[ind_i+ind_t*nbI];
                                    temp_S2M2 = temp_S2M2 + r_FRWToth_it_S2M2[ind_i+ind_t*nbI];
                                    temp_S2M3 = temp_S2M3 + r_FRWToth_it_S2M3[ind_i+ind_t*nbI];
                                    temp_S2M4 = temp_S2M4 + r_FRWToth_it_S2M4[ind_i+ind_t*nbI];
                                    temp_S3M1 = temp_S3M1 + r_FRWToth_it_S3M1[ind_i+ind_t*nbI];
                                    temp_S3M2 = temp_S3M2 + r_FRWToth_it_S3M2[ind_i+ind_t*nbI];
                                    temp_S3M3 = temp_S3M3 + r_FRWToth_it_S3M3[ind_i+ind_t*nbI];
                                    temp_S3M4 = temp_S3M4 + r_FRWToth_it_S3M4[ind_i+ind_t*nbI];
                                    temp_S4M1 = temp_S4M1 + r_FRWToth_it_S4M1[ind_i+ind_t*nbI];
                                    temp_S4M2 = temp_S4M2 + r_FRWToth_it_S4M2[ind_i+ind_t*nbI];
                                    temp_S4M3 = temp_S4M3 + r_FRWToth_it_S4M3[ind_i+ind_t*nbI];
                                    temp_S4M4 = temp_S4M4 + r_FRWToth_it_S4M4[ind_i+ind_t*nbI];

                                    temp2_S1M1 = temp2_S1M1 + r_FDWToth_it_S1M1[ind_i+ind_t*nbI];
                                    temp2_S1M2 = temp2_S1M2 + r_FDWToth_it_S1M2[ind_i+ind_t*nbI];
                                    temp2_S1M3 = temp2_S1M3 + r_FDWToth_it_S1M3[ind_i+ind_t*nbI];
                                    temp2_S1M4 = temp2_S1M4 + r_FDWToth_it_S1M4[ind_i+ind_t*nbI];
                                    temp2_S2M1 = temp2_S2M1 + r_FDWToth_it_S2M1[ind_i+ind_t*nbI];
                                    temp2_S2M2 = temp2_S2M2 + r_FDWToth_it_S2M2[ind_i+ind_t*nbI];
                                    temp2_S2M3 = temp2_S2M3 + r_FDWToth_it_S2M3[ind_i+ind_t*nbI];
                                    temp2_S2M4 = temp2_S2M4 + r_FDWToth_it_S2M4[ind_i+ind_t*nbI];
                                    temp2_S3M1 = temp2_S3M1 + r_FDWToth_it_S3M1[ind_i+ind_t*nbI];
                                    temp2_S3M2 = temp2_S3M2 + r_FDWToth_it_S3M2[ind_i+ind_t*nbI];
                                    temp2_S3M3 = temp2_S3M3 + r_FDWToth_it_S3M3[ind_i+ind_t*nbI];
                                    temp2_S3M4 = temp2_S3M4 + r_FDWToth_it_S3M4[ind_i+ind_t*nbI];
                                    temp2_S4M1 = temp2_S4M1 + r_FDWToth_it_S4M1[ind_i+ind_t*nbI];
                                    temp2_S4M2 = temp2_S4M2 + r_FDWToth_it_S4M2[ind_i+ind_t*nbI];
                                    temp2_S4M3 = temp2_S4M3 + r_FDWToth_it_S4M3[ind_i+ind_t*nbI];
                                    temp2_S4M4 = temp2_S4M4 + r_FDWToth_it_S4M4[ind_i+ind_t*nbI];


                              rans_L_eit[ind_i + ind_t*nbI] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * temp_S1M1 /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * temp_S1M2  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * temp_S1M3  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * temp_S1M4  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * temp_S2M1 /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * temp_S2M2  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * temp_S2M3  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * temp_S2M4  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * temp_S3M1 /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * temp_S3M2  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * temp_S3M3  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * temp_S3M4  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * temp_S4M1 /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * temp_S4M2  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * temp_S4M3  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * temp_S4M4  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                    rans_Y_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * temp2_S1M1 /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * temp2_S1M2  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * temp2_S1M3  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * temp2_S1M4  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * temp2_S2M1 /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * temp2_S2M2  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * temp2_S2M3  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * temp2_S2M4  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * temp2_S3M1 /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * temp2_S3M2  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * temp2_S3M3  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * temp2_S3M4  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * temp2_S4M1 /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * temp2_S4M2  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * temp2_S4M3  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * temp2_S4M4  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                   for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  //une fois Ytot g�n�r� � partir de Ltot (fraction d�barqu�e r�elle), on peut ajouter � Ltot les rejets d�barqu�s
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {
                        if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }



                   if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //et sous OD, on ajoute � Ltot les rejets autres

                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));

                   }

                                }

                            }

//Rprintf("K13\n");fichier << "K13" << endl;


/* insertion over quota management discards pour corriger D et L -> esp�ces dynamiques */

if (!((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & ((activeQR!=0) & (activeQR<=ind_t))) { //pas d'OD appliqu�, et activation du module demand�e


 // on s'occupe d'abord de la partie "autres"

    if (!isNull(getListElement(listQR, CHAR(STRING_ELT(sppList,e)))) & !isNull(getListElement(listQR_f, CHAR(STRING_ELT(sppList,e))))) { //TACs renseign�s aux 2 niveaux

        double *QR = REAL(getListElement(listQR, CHAR(STRING_ELT(sppList,e))));
        double *QR_f = REAL(getListElement(listQR_f, CHAR(STRING_ELT(sppList,e))));

        double QRoth = QR[ind_t];
        bool recal = false;
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) QRoth = QRoth - QR_f[ind_f + nbF*ind_t];

            double Ltot_oth = 0.0, Ytot_oth = 0.0; //, Ytot_othini = 0.0;
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                //double Yothini = rans_Y_eit[ind_i + 0*nbI];
                double Loth = rans_L_eit[ind_i + ind_t*nbI], Yoth = rans_Y_eit[ind_i + ind_t*nbI];
                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){
                //if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                //  Yothini = Yothini - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]];
                if (!ISNA(rans_L_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  Loth = Loth - rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  Yoth = Yoth - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                }
                Ltot_oth = Ltot_oth + Loth; Ytot_oth = Ytot_oth + Yoth; //Ytot_othini = Ytot_othini + Yothini;
            }

            rans_oqD_et[ind_t] = 0.0;

            if (Ltot_oth>QRoth) { //on proc�de � la correction "autres"

                recal = true;
                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    double Doth_i_t = rans_Y_eit[ind_i + ind_t*nbI] - rans_L_eit[ind_i + ind_t*nbI],// Yoth_i_0 = rans_Y_eit[ind_i + 0*nbI],
                           Yoth_i_t = rans_Y_eit[ind_i + ind_t*nbI];
                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++){
                    if (!ISNA(rans_D_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        Doth_i_t = Doth_i_t - rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    //if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                    //    Yoth_i_0 = Yoth_i_0 - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]];
                    if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        Yoth_i_t = Yoth_i_t - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }
                    rans_oqD_et[ind_t] = fmin2((Ltot_oth-QRoth) * (Yoth_i_t - Doth_i_t) / Ltot_oth, Yoth_i_t - Doth_i_t); //fmin2((Ltot_oth-QRoth) * Yoth_i_0 / Ytot_othini, Yoth_i_t - Doth_i_t);
                    Doth_i_t = fmin2( Doth_i_t + (Ltot_oth-QRoth) * (Yoth_i_t - Doth_i_t) / Ltot_oth, Yoth_i_t ); //fmin2( Doth_i_t + (Ltot_oth-QRoth) * Yoth_i_0 / Ytot_othini, Yoth_i_t );

                    if (ISNAN(Doth_i_t)) Doth_i_t = 0.0;
                    if (ISNAN(rans_oqD_et[ind_t])) rans_oqD_et[ind_t] = 0.0;
                    rans_L_eit[ind_i + ind_t*nbI] = Yoth_i_t - Doth_i_t; //on incr�mentera par la suite avec les L recalcul�s

                }
            }




            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  { //on proc�de � la correction "flottilles"


            double sumL = 0.0; //, sumYini = 0.0;

            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  sumL = sumL + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                //if (!ISNA(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                //  sumYini = sumYini + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]];
            }

            rans_oqD_eft[ind_f + nbF*ind_t] = 0.0;

            if (sumL>QR_f[ind_f + nbF*ind_t]) { //on proc�de � la correction sur la flottille detect�e

                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    //if (!ISNA(fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                    //        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                    //        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]])))
                    rans_oqD_eft[ind_f + nbF*ind_t] = rans_oqD_eft[ind_f + nbF*ind_t] +
                      fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) / sumL, //rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                            finite(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) -
                            finite(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]));

                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                       fmin2(
                        finite(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) +
                        (sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) / sumL, //(sumL-QR_f[ind_f + nbF*ind_t]) * rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                        finite(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]));

                    if (ISNAN(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;


                    rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                    if (!recal & !ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] - rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }

                    rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                       rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                       rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                    if (!recal & !ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }


                }
            }

            if (recal) {
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                      if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                      rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                      }
            }

        }
    }
}

    //---------
    // calcul de L_et
    //---------
    PROTECT(nDim = allocVector(INTSXP,4));
    nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e))); //fichier << "e = " << e << ", ind_e = " << ind_e << ", Name in sppList: " << CHAR(STRING_ELT(sppList,e)) << ", Name in sppListAll: " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

    if (ind_t==0) {
            PROTECT(ans_L_et = NEW_NUMERIC(nbT));

            PROTECT(dimCstL_et = allocVector(INTSXP, 1));
            int_dimCstL_et = INTEGER(dimCstL_et);  int_dimCstL_et[0] = nbT;
            setAttrib(ans_L_et, R_DimSymbol, dimCstL_et);

            PROTECT( dimNamL_et = allocVector(VECSXP,1));
            SET_VECTOR_ELT(dimNamL_et, 0, times);
            setAttrib(ans_L_et, R_DimNamesSymbol, dimNamL_et);

            rans_L_et = REAL(ans_L_et);
            } else {
                rans_L_et = REAL(VECTOR_ELT(out_L_et, ind_e));
                    }

            rans_L_et [ind_t] = REAL(aggregObj(ans_L_eit,nDim))[ind_t]; //fichier << "rans_L_et = " << rans_L_et [ind_t] << endl;

    if (ind_t==0) {
            SET_VECTOR_ELT(out_L_et, ind_e, ans_L_et);
            SET_STRING_ELT(rnames_eAll, ind_e, STRING_ELT(sppList,e));
    }

//fichier << "out_L_et = " << REAL(VECTOR_ELT(out_L_et, ind_e))[ind_t] << endl;
/*----------------------------------------------------------------*/



  UNPROTECT(2+4+4+1+2+1);
  if (Svec[e]==0){
    UNPROTECT(14);
  } else {UNPROTECT(21);}

  if (ind_t==0){
    UNPROTECT(2+3+9+1+3);

    if (Svec[e]==0){
        UNPROTECT(2);
    } else {UNPROTECT(4);}
  }
         }

}
//on passe aux esp�ces statiques
if(VERBOSE){Rprintf(" . Static sp");}
if (nbEstat>0) {

//Rprintf("H8\n");fichier << "H8" << endl;
    for (int e = 0 ; e < nbEstat ; e++) {
//Rprintf("H9\n");fichier << "H9" << endl;
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e))));

                            PROTECT(v_LPUE_eStat = getListElement(elmt, "LPUE_f_m_e"));
                            PROTECT(v_d_eStat = getListElement(elmt, "d_f_m_e")); //if (e==2) PrintValue(getListElement(elmt, "d_f_m_e"));

                            PROTECT(dimCst_LPUE_eStat = getAttrib(v_LPUE_eStat, install("DimCst")));
                            PROTECT(dimCst_d_eStat = getAttrib(v_d_eStat, install("DimCst")));

                            //tests sur les dimensions :
                            dim_d_eStat = INTEGER(dimCst_d_eStat);
                            if (((dim_d_eStat[0]!=0) & (dim_d_eStat[0]!=nbF)) | ((dim_d_eStat[1]!=0) & (dim_d_eStat[1]!=nbM)) |
                                (dim_d_eStat[2]!=0) | ((dim_d_eStat[3]!=0) & (dim_d_eStat[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in d_f_m_e element. Check .ini biological parameters files !!\n");
                            }

                            dim_LPUE_eStat = INTEGER(dimCst_LPUE_eStat);
                            if (((dim_LPUE_eStat[0]!=0) & (dim_LPUE_eStat[0]!=nbF)) | ((dim_LPUE_eStat[1]!=0) & (dim_LPUE_eStat[1]!=nbM)) |
                                (dim_LPUE_eStat[2]!=0) | ((dim_LPUE_eStat[3]!=0) & (dim_LPUE_eStat[3]!=nbT)))
                            {
                                error("Non_homogeneous dimensions in LPUE_f_m_e element. Check .ini biological parameters files !!\n");
                            }

//Rprintf("H10\n");fichier << "H10" << endl;
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de Ystat
                            //---------

                            //on d�termine l'attribut Dimension de Ystat
                            PROTECT(dimCst_eStat = allocVector(INTSXP, 4));//Rprintf("H10.1\n");
                            dim_eStat = INTEGER(dimCst_eStat);//Rprintf("H10.2\n");
                            dim_eStat[0] = dim_LPUE_eStat[0] ; dim_eStat[1] = dim_LPUE_eStat[1] ; dim_eStat[2] = 0;//Rprintf("H10.3\n");
                            dim_eStat[3] = nbT;//Rprintf("H10.4\n");

                            //variables d'effort

                    double *r_nbv_f, *r_nbds_f, *r_nbds2_f;//Rprintf("H10.5\n");fichier << "H10.5" << endl;

                    r_nbv_f = REAL(getListElement(Flist, "nbv_f_m"));
                    r_nbds_f = REAL(getListElement(Flist, "effort1_f_m"));
                    r_nbds2_f = REAL(getListElement(Flist, "effort2_f_m"));
//Rprintf("H10.6\n");fichier << "H10.6" << endl;

//                    int *fFactSup1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 50)), //ATTENTION : suppose au moins une esp�ce dynamiquement mod�lis�e
//                        *fFactSup2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 51));

                      int *fFactSup1 = INTEGER(iDim(INTEGER(getAttrib(getListElement(Flist, "nbv_f_m"), install("DimCst"))))),
                          *fFactSup2 = INTEGER(iDim(INTEGER(getAttrib(getListElement(Flist, "effort1_f_m"), install("DimCst")))));


//Rprintf("H10.7\n");fichier << "H10.7" << endl;
double *r_dd1_efm = REAL(getListElement(elmt, "dd1_f_m_e"));
double *r_dd2_efm = REAL(getListElement(elmt, "dd2_f_m_e"));
double *r_OD_e = REAL(getListElement(elmt, "OD_e"));
double *r_dst_efm = REAL(getListElement(elmt, "dst_f_m_e"));


                    if (ind_t==0){

                            //on cr�e le tableau r�sultat pour l'esp�ce en question
//Rprintf("H11\n");fichier << "H11" << endl;
                            PROTECT(Dim = allocVector(INTSXP, 3));
                            dim = INTEGER(Dim); dim[0]=dim_eStat[0]; dim[1]=dim_eStat[1]; dim[2]=dim_eStat[3];

                            PROTECT(ans_Ystat = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_Ystat, R_DimSymbol, Dim);
                            PROTECT(ans_Lstat = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_Lstat, R_DimSymbol, Dim);
                            PROTECT(ans_Dstat = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_Dstat, R_DimSymbol, Dim);
                            PROTECT(ans_statDD = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_statDD, R_DimSymbol, Dim);
                            PROTECT(ans_statLD = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_statLD, R_DimSymbol, Dim);
                            PROTECT(ans_statLDst = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_statLDst, R_DimSymbol, Dim);
                            PROTECT(ans_statLDor = NEW_NUMERIC(dim_eStat[0]*dim_eStat[1]*dim_eStat[3]));
                            setAttrib(ans_statLDor, R_DimSymbol, Dim);

                            PROTECT(dimnames_eStat = allocVector(VECSXP,3));
                            SET_VECTOR_ELT(dimnames_eStat, 0, fleetList);
                            SET_VECTOR_ELT(dimnames_eStat, 1, metierList);
                            SET_VECTOR_ELT(dimnames_eStat, 2, times);

                            rans_Ystat = REAL(ans_Ystat);
                            rans_Lstat = REAL(ans_Lstat);
                            rans_Dstat = REAL(ans_Dstat);
                            rans_statDD = REAL(ans_statDD);
                            rans_statLD = REAL(ans_statLD);
                            rans_statLDst = REAL(ans_statLDst);
                            rans_statLDor = REAL(ans_statLDor);


                         //if (!(r_OD_e[0]>0.5 & r_OD_e[0]<=(ind_t+1)) & (activeQR!=0)) { //over quota discards sera impl�ment�


                            PROTECT(ans_oqDstat = NEW_NUMERIC(nbF*nbT));

                            PROTECT(dimCstOQ_ft = allocVector(INTSXP, 2));
                            dimOQ_ft = INTEGER(dimCstOQ_ft);
                            dimOQ_ft[0] = nbF; dimOQ_ft[1] = nbT;

                            setAttrib(ans_oqDstat, R_DimSymbol, dimCstOQ_ft);

                            PROTECT(dimnames_oqD_eft = allocVector(VECSXP,2));
                            SET_VECTOR_ELT(dimnames_oqD_eft, 0, fleetList);
                            SET_VECTOR_ELT(dimnames_oqD_eft, 1, times);

                            setAttrib(ans_oqDstat, R_DimNamesSymbol, dimnames_oqD_eft);//Rprintf("EE");

                            rans_oqDstat = REAL(ans_oqDstat);// //Rprintf("FF");
                            for (int tt=0; tt<(nbF*nbT); tt++) rans_oqDstat[tt] = 0.0;

                        //}


                    } else {
//Rprintf("H12\n");fichier << "H12" << endl;
                            rans_Ystat = REAL(VECTOR_ELT(out_Ystat, e));
                            rans_Lstat = REAL(VECTOR_ELT(out_Lstat, e));
                            rans_Dstat = REAL(VECTOR_ELT(out_Dstat, e));
                            rans_statDD = REAL(VECTOR_ELT(out_statDD_efm, e));
                            rans_statLD = REAL(VECTOR_ELT(out_statLD_efm, e));
                            rans_statLDst = REAL(VECTOR_ELT(out_statLDst_efm, e));
                            rans_statLDor = REAL(VECTOR_ELT(out_statLDor_efm, e));//Rprintf("GG");

                            rans_oqDstat = REAL(VECTOR_ELT(out_oqDstat, e));//Rprintf("HH");

                    }


                   double *r_LPUE_eStat = REAL(v_LPUE_eStat);
                   double *r_d_eStat = REAL(v_d_eStat); //if (e==2) PrintValue(v_d_eStat);
//Rprintf("H13\n");fichier << "H13" << endl;
                   for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                   for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {

                            rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               r_LPUE_eStat[ind_f + nbF*ind_m + 0*ind_t] *
                               r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + 0*fFactSup1[2] + ind_t*fFactSup1[3]] *
                               r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + 0*fFactSup2[2] + ind_t*fFactSup2[3]]*
                               r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + 0*fFactSup2[2] + ind_t*fFactSup2[3]]/1000.0;

                            if (ISNA(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))
                                  rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                            rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] / (1 - finite(r_d_eStat[ind_f + nbF*ind_m + 0*ind_t]));

                            if (ISNAN(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                            rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]);

                            rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                            rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * r_d_eStat[ind_f + nbF*ind_m + 0*ind_t]);

                            if (ISNA(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                   }


                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & (!ISNA(r_dd1_efm[ind_f + nbF*ind_m]))) {

                           rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                              fmin2( r_dd1_efm[ind_f + nbF*ind_m] * rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] ,
                                      rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] ) ;

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd2_efm[ind_f + nbF*ind_m])){

                                double rYsum = 0.0;
                                if (ind_t==0) {
                                 rYsum = REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m];
                                } else {
                                 rYsum = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                }

                                rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                fmin2( r_dd2_efm[ind_f + nbF*ind_m] * rYsum , rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] );

                        } else {

                             rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t];


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLDst[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * finite(r_dst_efm[ind_f + nbF*ind_m]);


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLDor[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_statLDst[ind_f + nbF*ind_m + nbF*nbMe*ind_t];


/* insertion over quota management discards pour corriger D et L -> esp�ces statiques */

if (!((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & ((activeQR!=0) & (activeQR<=ind_t))) { //pas d'OD appliqu�, et activation du module demand�e

    if (!isNull(getListElement(listQR_f, CHAR(STRING_ELT(sppListStat,e))))) { //TACs renseign�s au niveau flottille

        double *QR_f = REAL(getListElement(listQR_f, CHAR(STRING_ELT(sppListStat,e))));

            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  { //on proc�de � la correction "flottilles"

            double sumL = 0.0;

            for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {
                if (!ISNA(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))
                  sumL = sumL + rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                //if (!ISNA(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0]))
                //  sumYini = sumYini + rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0];
            }

            rans_oqDstat[ind_f + nbF*ind_t] = 0.0;

            if (sumL>QR_f[ind_f + nbF*ind_t]) { //on proc�de � la correction sur la flottille detect�e

                for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                    if (!ISNAN(fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) / sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                                       finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]-rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))))
                    rans_oqDstat[ind_f + nbF*ind_t] = rans_oqDstat[ind_f + nbF*ind_t] +
                                fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) / sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                                       finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]-rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]));

                    rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                       fmin2(finite(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) + (sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])/sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                             finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]));

                    if (ISNAN(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))
                        rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                    rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                    rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                }
            }
        }
    }
}

/*----------------------------------------------------------------*/



//Rprintf("H14\n");fichier << "H14" << endl;
                   if (ind_t==0){

                            setAttrib(ans_Ystat, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_Ystat, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_Ystat, e, ans_Ystat);

                            setAttrib(ans_Lstat, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_Lstat, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_Lstat, e, ans_Lstat);

                            setAttrib(ans_Dstat, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_Dstat, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_Dstat, e, ans_Dstat);

                            SET_STRING_ELT(rnames_eStat, e, STRING_ELT(sppListStat,e));

                            setAttrib(ans_statDD, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_statDD, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_statDD_efm, e, ans_statDD);

                            setAttrib(ans_statLD, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_statLD, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_statLD_efm, e, ans_statLD);

                            setAttrib(ans_statLDst, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_statLDst, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_statLDst_efm, e, ans_statLDst);

                            setAttrib(ans_statLDor, R_DimNamesSymbol, dimnames_eStat);
                            setAttrib(ans_statLDor, install("DimCst"), dimCst_eStat);

                            SET_VECTOR_ELT(out_statLDor_efm, e, ans_statLDor);

                            SET_VECTOR_ELT(out_oqDstat, e, ans_oqDstat);

                    }

    //---------
    //  Calcul L_et
    //---------
    PROTECT(nDim = allocVector(INTSXP,4));
    nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e))); //fichier << "e = " << e << ", ind_e = " << ind_e << ", Name in sppListStat: " << CHAR(STRING_ELT(sppListStat,e)) << ", Name in sppListAll: " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

    if (ind_t==0) {
            PROTECT(ans_L_et = NEW_NUMERIC(nbT));

            PROTECT(dimCstL_et = allocVector(INTSXP, 1));
            int_dimCstL_et = INTEGER(dimCstL_et);  int_dimCstL_et[0] = nbT;
            setAttrib(ans_L_et, R_DimSymbol, dimCstL_et);

            PROTECT( dimNamL_et = allocVector(VECSXP,1));
            SET_VECTOR_ELT(dimNamL_et, 0, times);
            setAttrib(ans_L_et, R_DimNamesSymbol, dimNamL_et);

            rans_L_et = REAL(ans_L_et);
            } else {
                rans_L_et = REAL(VECTOR_ELT(out_L_et, ind_e));
                    }

            rans_L_et [ind_t] = REAL(aggregObj(ans_Ystat,nDim))[ind_t]; //fichier << "rans_L_et = " << rans_L_et[ind_t] << endl;


    if (ind_t==0) {
            SET_VECTOR_ELT(out_L_et, ind_e, ans_L_et);
            SET_STRING_ELT(rnames_eAll, ind_e, STRING_ELT(sppListStat,e));
    }

    //fichier << "out_L_et = " << REAL(VECTOR_ELT(out_L_et, ind_e))[ind_t] << endl;

//Rprintf("H14.1\n");fichier << "H14.1" << endl;
                    if (ind_t==0) {

                      UNPROTECT(9);
                      /*if (!(r_OD_e[0]>0.5 & r_OD_e[0]<=(ind_t+1)) & (activeQR!=0)) */UNPROTECT(3);
                      UNPROTECT(3);

                    }

                    UNPROTECT(6+1);

         }

}


cUpdate = false;



if (ind_t==0) UNPROTECT(3);
UNPROTECT(3);

//PrintValue(out_Fbar_et);

} else {
    if(VERBOSE){Rprintf("cUpdate : ");}
//fichier << "cUpdate = " << cUpdate << endl;
//fichier << "ind_t = " << ind_t << endl;

double *rans_Ytot_fm = REAL(out_Ytot_fm);

for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {
  rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] =
     rans_Yothsue_fm[ind_f + nbF*ind_m]*reff1[ind_f + nbF*ind_m]*reff2[ind_f + nbF*ind_m]*rnbv[ind_f + nbF*ind_m];

  if (ISNA(rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t])) rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = 0.0;
}

if(VERBOSE){Rprintf("Dyna sp");}
if (nbE>0) {

    for (int e = 0 ; e < nbE ; e++) {

//Rprintf("H15\n");fichier << "H15" << endl;
                            int nbI = length(VECTOR_ELT(namDC,e));
                            int ind_e;
                            int *nd;

                            SEXP elmt, nDim, ans_L_eit;
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                            double *rans_C_efmit=&NA_REAL,
                                   *rans_C_eit=&NA_REAL,
                                   *rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e)),
                                   *rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e)),
                                   *rans_L_eit = REAL(VECTOR_ELT(out_L_eit,e)),
                                   *rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e)),
                                   *rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e)),
                                   *rans_DD_efmit = REAL(VECTOR_ELT(out_DD_efmi,e)),
                                   *rans_LD_efmit = REAL(VECTOR_ELT(out_LD_efmi,e)),
                                   *r_F_efmit = &NA_REAL,
                                   *r_Foth_i=&NA_REAL,
                                   *r_Foth_i_G1=&NA_REAL,
                                   *r_Foth_i_G2=&NA_REAL,
                                   *r_N_eit = &NA_REAL,
                                   *r_Z_eit = &NA_REAL,
                                   *rans_L_et = &NA_REAL;

                            double *rans_oqD_eft = REAL(VECTOR_ELT(out_oqD_eft,e));//Rprintf("JJ");
                            double *rans_oqD_et = REAL(VECTOR_ELT(out_oqD_et,e));//Rprintf("KK");
                            double *r_B_et = REAL(VECTOR_ELT(out_B_et,e));
                            SEXP v_wD_ei, v_wD_ei_G1, v_wD_ei_G2 ;

//Rprintf("H15.2\n");fichier << "H15.2" << endl;
                                int     *fact1_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 26)),
                                        *fact2_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 27)),
                                        *fact3_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 28)),
                                        *fact4_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 29)),
                                        *fact5_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 31)),
                                        *fact6_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 34)),
                                        *fact7_C  = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 35));


double *r_dd1_efm = REAL(getListElement(elmt, "dd1_f_m_e"));//Rprintf("H15.3\n");
double *r_dd2_efm = REAL(getListElement(elmt, "dd2_f_m_e"));//Rprintf("H15.4\n");
double *r_OD_e = REAL(getListElement(elmt, "OD_e"));//Rprintf("H15.5\n");fichier << "H15.5" << endl;

                                //�quation n�1
//Rprintf("H16\n");fichier << "H16" << endl;

                            if ((Qvec[e]==1) & (Svec[e]==0)) {

                                PROTECT(v_wD_ei = getListElement(elmt, "wD_i"));
//double Btemp;
                                rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit,e));
                                rans_C_eit = REAL(VECTOR_ELT(out_C_eit,e));


                                double *r_F_fmi_S1M1 = REAL(getListElement(out_F_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M2 = REAL(getListElement(out_F_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M3 = REAL(getListElement(out_F_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S1M4 = REAL(getListElement(out_F_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M1 = REAL(getListElement(out_F_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M2 = REAL(getListElement(out_F_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M3 = REAL(getListElement(out_F_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S2M4 = REAL(getListElement(out_F_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M1 = REAL(getListElement(out_F_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M2 = REAL(getListElement(out_F_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M3 = REAL(getListElement(out_F_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S3M4 = REAL(getListElement(out_F_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M1 = REAL(getListElement(out_F_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M2 = REAL(getListElement(out_F_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M3 = REAL(getListElement(out_F_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_F_fmi_S4M4 = REAL(getListElement(out_F_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));

                                double *r_Z_eit_S1M1 = REAL(getListElement(out_Z_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S1M2 = REAL(getListElement(out_Z_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S1M3 = REAL(getListElement(out_Z_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S1M4 = REAL(getListElement(out_Z_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M1 = REAL(getListElement(out_Z_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M2 = REAL(getListElement(out_Z_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M3 = REAL(getListElement(out_Z_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S2M4 = REAL(getListElement(out_Z_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M1 = REAL(getListElement(out_Z_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M2 = REAL(getListElement(out_Z_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M3 = REAL(getListElement(out_Z_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S3M4 = REAL(getListElement(out_Z_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M1 = REAL(getListElement(out_Z_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M2 = REAL(getListElement(out_Z_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M3 = REAL(getListElement(out_Z_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_Z_eit_S4M4 = REAL(getListElement(out_Z_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));

                                double *r_N_eit_S1M1 = REAL(getListElement(out_N_eit_S1M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M2 = REAL(getListElement(out_N_eit_S1M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M3 = REAL(getListElement(out_N_eit_S1M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S1M4 = REAL(getListElement(out_N_eit_S1M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M1 = REAL(getListElement(out_N_eit_S2M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M2 = REAL(getListElement(out_N_eit_S2M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M3 = REAL(getListElement(out_N_eit_S2M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S2M4 = REAL(getListElement(out_N_eit_S2M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M1 = REAL(getListElement(out_N_eit_S3M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M2 = REAL(getListElement(out_N_eit_S3M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M3 = REAL(getListElement(out_N_eit_S3M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S3M4 = REAL(getListElement(out_N_eit_S3M4 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M1 = REAL(getListElement(out_N_eit_S4M1 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M2 = REAL(getListElement(out_N_eit_S4M2 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M3 = REAL(getListElement(out_N_eit_S4M3 , CHAR(STRING_ELT(sppList,e))));
                                double *r_N_eit_S4M4 = REAL(getListElement(out_N_eit_S4M4 , CHAR(STRING_ELT(sppList,e))));

                                double  *r_M_ei = REAL(getListElement(elmt, "M_i"));
                                int *fact3_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 14));
//Rprintf("H17\n");fichier << "H17" << endl;

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    rans_C_eit[ind_i + ind_t*nbI] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] - r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]]) /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));
//Rprintf("H18\n");fichier << "H18" << endl;

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_F_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));

                                }

                                }

                                //--------------
                                //Rprintf("H4.6\n");fichier << "H4.6" << endl;


                                double *r_FRWT_fmi_S1M1 = REAL(getListElement(out_FRWT_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M2 = REAL(getListElement(out_FRWT_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M3 = REAL(getListElement(out_FRWT_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S1M4 = REAL(getListElement(out_FRWT_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M1 = REAL(getListElement(out_FRWT_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M2 = REAL(getListElement(out_FRWT_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M3 = REAL(getListElement(out_FRWT_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S2M4 = REAL(getListElement(out_FRWT_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M1 = REAL(getListElement(out_FRWT_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M2 = REAL(getListElement(out_FRWT_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M3 = REAL(getListElement(out_FRWT_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S3M4 = REAL(getListElement(out_FRWT_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M1 = REAL(getListElement(out_FRWT_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M2 = REAL(getListElement(out_FRWT_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M3 = REAL(getListElement(out_FRWT_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FRWT_fmi_S4M4 = REAL(getListElement(out_FRWT_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));



//Rprintf("H4.7\n");fichier << "H4.7" << endl;//PrintValue(out_Z_eit_S1M1);
                                double *r_FDWT_fmi_S1M1 = REAL(getListElement(out_FDWT_fmi_S1M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M2 = REAL(getListElement(out_FDWT_fmi_S1M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M3 = REAL(getListElement(out_FDWT_fmi_S1M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S1M4 = REAL(getListElement(out_FDWT_fmi_S1M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M1 = REAL(getListElement(out_FDWT_fmi_S2M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M2 = REAL(getListElement(out_FDWT_fmi_S2M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M3 = REAL(getListElement(out_FDWT_fmi_S2M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S2M4 = REAL(getListElement(out_FDWT_fmi_S2M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M1 = REAL(getListElement(out_FDWT_fmi_S3M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M2 = REAL(getListElement(out_FDWT_fmi_S3M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M3 = REAL(getListElement(out_FDWT_fmi_S3M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S3M4 = REAL(getListElement(out_FDWT_fmi_S3M4, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M1 = REAL(getListElement(out_FDWT_fmi_S4M1, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M2 = REAL(getListElement(out_FDWT_fmi_S4M2, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M3 = REAL(getListElement(out_FDWT_fmi_S4M3, CHAR(STRING_ELT(sppList,e))));
                                double *r_FDWT_fmi_S4M4 = REAL(getListElement(out_FDWT_fmi_S4M4, CHAR(STRING_ELT(sppList,e))));
//Rprintf("H5ddd\n");fichier << "H5ddd" << endl;

                                double *r_FRWToth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 176));
                                double *r_FRWToth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 177));
                                double *r_FRWToth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 178));
                                double *r_FRWToth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 179));
                                double *r_FRWToth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 180));
                                double *r_FRWToth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 181));
                                double *r_FRWToth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 182));
                                double *r_FRWToth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 183));
                                double *r_FRWToth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 184));
                                double *r_FRWToth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 185));
                                double *r_FRWToth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 186));
                                double *r_FRWToth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 187));
                                double *r_FRWToth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 188));
                                double *r_FRWToth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 189));
                                double *r_FRWToth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 190));
                                double *r_FRWToth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 191));
//Rprintf("H5.1\n");fichier << "H5.1" << endl;
                                double *r_FDWToth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 208));
                                double *r_FDWToth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 209));
                                double *r_FDWToth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 210));
                                double *r_FDWToth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 211));
                                double *r_FDWToth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 212));
                                double *r_FDWToth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 213));
                                double *r_FDWToth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 214));
                                double *r_FDWToth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 215));
                                double *r_FDWToth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 216));
                                double *r_FDWToth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 217));
                                double *r_FDWToth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 218));
                                double *r_FDWToth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 219));
                                double *r_FDWToth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 220));
                                double *r_FDWToth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 221));
                                double *r_FDWToth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 222));
                                double *r_FDWToth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 223));
//Rprintf("H5.2\n");fichier << "H5.2" << endl;



                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {


                                    rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                      (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));

//if ((rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]>0) &
//    (rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]<0.00000001)) {
//if (ind_t==1) {
//std::stringstream ff3S1M1;
//ff3S1M1 << r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
//fichier1 << ff3S1M1.str() << endl;
//}


                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                      (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) *
                                       r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                                rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];


                                rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];


                                }






                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd1_efm[ind_f + nbF*ind_m])) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd2_efm[ind_f + nbF*ind_m])){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;
                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];





                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                     double temp_S1M1 = 0.0, temp_S1M2 = 0.0, temp_S1M3 = 0.0, temp_S1M4 = 0.0,
                                            temp_S2M1 = 0.0, temp_S2M2 = 0.0, temp_S2M3 = 0.0, temp_S2M4 = 0.0,
                                            temp_S3M1 = 0.0, temp_S3M2 = 0.0, temp_S3M3 = 0.0, temp_S3M4 = 0.0,
                                            temp_S4M1 = 0.0, temp_S4M2 = 0.0, temp_S4M3 = 0.0, temp_S4M4 = 0.0,
                                            temp2_S1M1 = 0.0, temp2_S1M2 = 0.0, temp2_S1M3 = 0.0, temp2_S1M4 = 0.0,
                                            temp2_S2M1 = 0.0, temp2_S2M2 = 0.0, temp2_S2M3 = 0.0, temp2_S2M4 = 0.0,
                                            temp2_S3M1 = 0.0, temp2_S3M2 = 0.0, temp2_S3M3 = 0.0, temp2_S3M4 = 0.0,
                                            temp2_S4M1 = 0.0, temp2_S4M2 = 0.0, temp2_S4M3 = 0.0, temp2_S4M4 = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M1 = temp_S1M1 + r_FRWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M2 = temp_S1M2 + r_FRWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M3 = temp_S1M3 + r_FRWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S1M4 = temp_S1M4 + r_FRWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M1 = temp_S2M1 + r_FRWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M2 = temp_S2M2 + r_FRWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M3 = temp_S2M3 + r_FRWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S2M4 = temp_S2M4 + r_FRWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M1 = temp_S3M1 + r_FRWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M2 = temp_S3M2 + r_FRWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M3 = temp_S3M3 + r_FRWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S3M4 = temp_S3M4 + r_FRWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M1 = temp_S4M1 + r_FRWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M2 = temp_S4M2 + r_FRWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M3 = temp_S4M3 + r_FRWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp_S4M4 = temp_S4M4 + r_FRWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];


                                        if (!ISNA(r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M1 = temp2_S1M1 + r_FDWT_fmi_S1M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M2 = temp2_S1M2 + r_FDWT_fmi_S1M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M3 = temp2_S1M3 + r_FDWT_fmi_S1M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S1M4 = temp2_S1M4 + r_FDWT_fmi_S1M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M1 = temp2_S2M1 + r_FDWT_fmi_S2M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M2 = temp2_S2M2 + r_FDWT_fmi_S2M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M3 = temp2_S2M3 + r_FDWT_fmi_S2M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S2M4 = temp2_S2M4 + r_FDWT_fmi_S2M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M1 = temp2_S3M1 + r_FDWT_fmi_S3M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M2 = temp2_S3M2 + r_FDWT_fmi_S3M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M3 = temp2_S3M3 + r_FDWT_fmi_S3M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S3M4 = temp2_S3M4 + r_FDWT_fmi_S3M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M1 = temp2_S4M1 + r_FDWT_fmi_S4M1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M2 = temp2_S4M2 + r_FDWT_fmi_S4M2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M3 = temp2_S4M3 + r_FDWT_fmi_S4M3[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];
                                        if (!ISNA(r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp2_S4M4 = temp2_S4M4 + r_FDWT_fmi_S4M4[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];


                                    }

                                    temp_S1M1 = temp_S1M1 + r_FRWToth_it_S1M1[ind_i+ind_t*nbI];
                                    temp_S1M2 = temp_S1M2 + r_FRWToth_it_S1M2[ind_i+ind_t*nbI];
                                    temp_S1M3 = temp_S1M3 + r_FRWToth_it_S1M3[ind_i+ind_t*nbI];
                                    temp_S1M4 = temp_S1M4 + r_FRWToth_it_S1M4[ind_i+ind_t*nbI];
                                    temp_S2M1 = temp_S2M1 + r_FRWToth_it_S2M1[ind_i+ind_t*nbI];
                                    temp_S2M2 = temp_S2M2 + r_FRWToth_it_S2M2[ind_i+ind_t*nbI];
                                    temp_S2M3 = temp_S2M3 + r_FRWToth_it_S2M3[ind_i+ind_t*nbI];
                                    temp_S2M4 = temp_S2M4 + r_FRWToth_it_S2M4[ind_i+ind_t*nbI];
                                    temp_S3M1 = temp_S3M1 + r_FRWToth_it_S3M1[ind_i+ind_t*nbI];
                                    temp_S3M2 = temp_S3M2 + r_FRWToth_it_S3M2[ind_i+ind_t*nbI];
                                    temp_S3M3 = temp_S3M3 + r_FRWToth_it_S3M3[ind_i+ind_t*nbI];
                                    temp_S3M4 = temp_S3M4 + r_FRWToth_it_S3M4[ind_i+ind_t*nbI];
                                    temp_S4M1 = temp_S4M1 + r_FRWToth_it_S4M1[ind_i+ind_t*nbI];
                                    temp_S4M2 = temp_S4M2 + r_FRWToth_it_S4M2[ind_i+ind_t*nbI];
                                    temp_S4M3 = temp_S4M3 + r_FRWToth_it_S4M3[ind_i+ind_t*nbI];
                                    temp_S4M4 = temp_S4M4 + r_FRWToth_it_S4M4[ind_i+ind_t*nbI];

                                    temp2_S1M1 = temp2_S1M1 + r_FDWToth_it_S1M1[ind_i+ind_t*nbI];
                                    temp2_S1M2 = temp2_S1M2 + r_FDWToth_it_S1M2[ind_i+ind_t*nbI];
                                    temp2_S1M3 = temp2_S1M3 + r_FDWToth_it_S1M3[ind_i+ind_t*nbI];
                                    temp2_S1M4 = temp2_S1M4 + r_FDWToth_it_S1M4[ind_i+ind_t*nbI];
                                    temp2_S2M1 = temp2_S2M1 + r_FDWToth_it_S2M1[ind_i+ind_t*nbI];
                                    temp2_S2M2 = temp2_S2M2 + r_FDWToth_it_S2M2[ind_i+ind_t*nbI];
                                    temp2_S2M3 = temp2_S2M3 + r_FDWToth_it_S2M3[ind_i+ind_t*nbI];
                                    temp2_S2M4 = temp2_S2M4 + r_FDWToth_it_S2M4[ind_i+ind_t*nbI];
                                    temp2_S3M1 = temp2_S3M1 + r_FDWToth_it_S3M1[ind_i+ind_t*nbI];
                                    temp2_S3M2 = temp2_S3M2 + r_FDWToth_it_S3M2[ind_i+ind_t*nbI];
                                    temp2_S3M3 = temp2_S3M3 + r_FDWToth_it_S3M3[ind_i+ind_t*nbI];
                                    temp2_S3M4 = temp2_S3M4 + r_FDWToth_it_S3M4[ind_i+ind_t*nbI];
                                    temp2_S4M1 = temp2_S4M1 + r_FDWToth_it_S4M1[ind_i+ind_t*nbI];
                                    temp2_S4M2 = temp2_S4M2 + r_FDWToth_it_S4M2[ind_i+ind_t*nbI];
                                    temp2_S4M3 = temp2_S4M3 + r_FDWToth_it_S4M3[ind_i+ind_t*nbI];
                                    temp2_S4M4 = temp2_S4M4 + r_FDWToth_it_S4M4[ind_i+ind_t*nbI];


                              rans_L_eit[ind_i + ind_t*nbI] =

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * temp_S1M1 /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * temp_S1M2  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * temp_S1M3  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * temp_S1M4  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * temp_S2M1 /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * temp_S2M2  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * temp_S2M3  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * temp_S2M4  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * temp_S3M1 /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * temp_S3M2  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * temp_S3M3  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * temp_S3M4  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * temp_S4M1 /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * temp_S4M2  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * temp_S4M3  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * temp_S4M4  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                            rans_Y_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * temp2_S1M1 /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * temp2_S1M2  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * temp2_S1M3  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * temp2_S1M4  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * temp2_S2M1 /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * temp2_S2M2  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * temp2_S2M3  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * temp2_S2M4  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * temp2_S3M1 /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * temp2_S3M2  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * temp2_S3M3  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * temp2_S3M4  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * temp2_S4M1 /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * temp2_S4M2  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * temp2_S4M3  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * temp2_S4M4  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  //une fois Ytot g�n�r� � partir de Ltot (fraction d�barqu�e r�elle), on peut ajouter � Ltot les rejets d�barqu�s
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {
                        if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }



                   if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //et sous OD, on ajoute � Ltot les rejets autres

                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +

                                      r_N_eit_S1M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S1M1[ind_i + ind_t*nbI] + (r_Z_eit_S1M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M2[ind_i + ind_t*nbI] + (r_Z_eit_S1M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M3[ind_i + ind_t*nbI] + (r_Z_eit_S1M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S1M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S1M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S1M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S1M4[ind_i + ind_t*nbI] + (r_Z_eit_S1M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S2M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S2M1[ind_i + ind_t*nbI] + (r_Z_eit_S2M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M2[ind_i + ind_t*nbI] + (r_Z_eit_S2M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M3[ind_i + ind_t*nbI] + (r_Z_eit_S2M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S2M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S2M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S2M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S2M4[ind_i + ind_t*nbI] + (r_Z_eit_S2M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S3M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S3M1[ind_i + ind_t*nbI] + (r_Z_eit_S3M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M2[ind_i + ind_t*nbI] + (r_Z_eit_S3M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M3[ind_i + ind_t*nbI] + (r_Z_eit_S3M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S3M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S3M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S3M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S3M4[ind_i + ind_t*nbI] + (r_Z_eit_S3M4[ind_i + ind_t*nbI]==0)) +

                                      r_N_eit_S4M1[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M1[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M1[ind_i+ind_t*nbI] /
                                       (r_Z_eit_S4M1[ind_i + ind_t*nbI] + (r_Z_eit_S4M1[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M2[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M2[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M2[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M2[ind_i + ind_t*nbI] + (r_Z_eit_S4M2[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M3[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M3[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M3[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M3[ind_i + ind_t*nbI] + (r_Z_eit_S4M3[ind_i + ind_t*nbI]==0)) +
                                      r_N_eit_S4M4[ind_i + ind_t*nbI] * (1-exp(-r_Z_eit_S4M4[ind_i + ind_t*nbI]/4)) * r_FDWToth_it_S4M4[ind_i+ind_t*nbI]  /
                                       (r_Z_eit_S4M4[ind_i + ind_t*nbI] + (r_Z_eit_S4M4[ind_i + ind_t*nbI]==0));
                                       }
                                       }


                            } else if ((Qvec[e]==0) & (Svec[e]==0)){
                                PROTECT(v_wD_ei = getListElement(elmt, "wD_i"));

//double Btemp;
                                rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit,e));
                                rans_C_eit = REAL(VECTOR_ELT(out_C_eit,e));
                                r_Foth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44)); //Rprintf("Dans EVAR (l.12171), Fothi = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));

                                r_F_efmit = REAL(VECTOR_ELT(out_F_fmi,e));
                                r_N_eit = REAL(VECTOR_ELT(out_N_eit, e));
                                r_Z_eit = REAL(VECTOR_ELT(out_Z_eit, e));



                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                     rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];



                                //�quation

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    double temp = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }

                                    rans_C_eit[ind_i + ind_t*nbI] =
                                        (temp + r_Foth_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                }

                                //------------

                             double *r_wL_ei = REAL(getListElement(elmt, "wL_i")),
                                    *r_wD_ei = REAL(getListElement(elmt, "wD_i")),
                                    *r_d_efmit = REAL(getListElement(elmt, "d_i")),
                                    *doth_eit = REAL(getListElement(elmt, "doth_i"));

if (nbI>1) {
                               //�quation n�2

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;



                                      rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                }


                              for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_Y_eit[ind_i + ind_t*nbI] =
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                        rans_C_eit[ind_i + ind_t*nbI] / 1000;

                                      rans_L_eit[ind_i + ind_t*nbI] = NA_REAL;
                              }
} else {


                           //double *Bspict = REAL(VECTOR_ELT(intermBIOMspict, e));
                            // on peut sommer avant d'appliquer � F puisque F est suppos� (pour le moment) constant sur l'ensemble de l'ann�e N
                           //Btemp = 0.0;
                           //for (int ii = ind_t*16 ; ii < (ind_t*16 + 16) ; ii++) Btemp = Btemp + Bspict[ii]/16;

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                                     rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + 0*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + 0*fact2_C[2] + ind_t*fact2_C[3]] * r_B_et[ind_t*fact3_C[3]];

                                //�quation

                                double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                 if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + 0*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + 0*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                                rans_Y_eit[0 + ind_t*1] = (temp + r_Foth_i[0 + ind_t*1]) * r_B_et[ind_t*fact3_C[3]];//if (nbI==1) {//Rprintf("Yi");PrintValue(out_Y_eit);}

}
                               //�quation n�3

                            if (all_is_na(v_wD_ei)) { // on peut aussi laisser le test SPiCT � l'int�rieur car les deux conditions sont �quivalentes

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                }



if (nbI>1) {
                            //Loth_eit
                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit[ind_i]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;
                            }
} else {
////Rprintf("LtotAvant4");PrintValue(out_L_eit);
                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)

                                rans_L_eit[0 + ind_t*1] = r_Foth_i[0 + ind_t*1] * r_B_et[0*fact5_C[0]  + 0*fact5_C[1] + 0*fact5_C[2] + ind_t*fact5_C[3]];


                            } else { //pas d'OD
////Rprintf("LtotAvant31");PrintValue(out_L_eit);//Rprintf("Foth %f B %f doth %f 1-doth %f",r_Foth_i[0 + ind_t*1],r_B_et[0*fact5_C[0]  + 0*fact5_C[1] + 0*fact5_C[2] + ind_t*fact5_C[3]],doth_eit[0],1-doth_eit[0]);
                                rans_L_eit[0 + ind_t*1] = r_Foth_i[0 + ind_t*1] * r_B_et[ind_t*fact3_C[3]] * (1-doth_eit[0]);
////Rprintf("LtotAvant32");PrintValue(out_L_eit);
                            }
                            ////Rprintf("LtotAvant3");PrintValue(out_L_eit);

}



                            } else {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;


                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                        r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                }


                               //Loth_eit : pas d'exemption pour les autres si OD

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit, pas de rejet car pas d'exemption

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i[ind_i + ind_t*nbI] * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit[ind_i] * r_wD_ei[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000;
                            }



                            }


                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd1_efm[ind_f + nbF*ind_m])) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd2_efm[ind_f + nbF*ind_m])){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }

//if (nbI==1) {//Rprintf("LtotAvant2");PrintValue(out_L_eit);}

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; // il reste ensuite � integrer Lefmit pour contituer Ltot_i

                    }


//if (nbI==1) {//Rprintf("LtotAvant");PrintValue(out_L_eit);}
                               //�quation n�4

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                      if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                        rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; //Ltot constitu�

                                }
//if (nbI==1) {//Rprintf("LtotApres");PrintValue(out_L_eit);}
//---------------

                            } else if ((Qvec[e]==0) & (Svec[e]==1)){
                                PROTECT(v_wD_ei_G1 = getListElement(elmt, "wD_i_G1"));
                                PROTECT(v_wD_ei_G2 = getListElement(elmt, "wD_i_G2"));

                                double  *rans_C_efmit_G1 = REAL(VECTOR_ELT(out_C_efmit_G1,e)),
                                        *rans_C_efmit_G2 = REAL(VECTOR_ELT(out_C_efmit_G2,e)),

                                        *rans_C_eit_G1 = REAL(VECTOR_ELT(out_C_eit_G1,e)),
                                        *rans_C_eit_G2 = REAL(VECTOR_ELT(out_C_eit_G2,e)),

                                        *r_F_efmit_G1 = REAL(VECTOR_ELT(out_F_fmi_G1,e)),
                                        *r_F_efmit_G2 = REAL(VECTOR_ELT(out_F_fmi_G2,e)),
                                        *r_N_eit_G1 = REAL(VECTOR_ELT(out_N_eit_G1, e)),
                                        *r_N_eit_G2 = REAL(VECTOR_ELT(out_N_eit_G2, e)),
                                        *r_Z_eit_G1 = REAL(VECTOR_ELT(out_Z_eit_G1, e)),
                                        *r_Z_eit_G2 = REAL(VECTOR_ELT(out_Z_eit_G2, e));

                                r_Foth_i_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 224));
                                r_Foth_i_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 225));
//Rprintf("H19\n");fichier << "H19" << endl;

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                        rans_C_efmit_G1[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit_G1[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                        rans_C_efmit_G2[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit_G2[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                }

                                //�quation

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    double tempG1 = 0.0;
                                    double tempG2 = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        tempG1 = tempG1 + r_F_efmit_G1[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                        if (!ISNA(r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        tempG2 = tempG2 + r_F_efmit_G2[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }
                                    //fichier << "Age " << ind_i << ", tempG1 = " << tempG1 << endl;
                                    //fichier << "Age " << ind_i << ", FothG1 = " << r_Foth_i_G1[ind_i + ind_t*nbI] << endl;
                                    //fichier << "Age " << ind_i << ", tempG2 = " << tempG2 << endl;
                                    //fichier << "Age " << ind_i << ", FothG2 = " << r_Foth_i_G2[ind_i + ind_t*nbI] << endl;

                                    rans_C_eit_G1[ind_i + ind_t*nbI] =
                                        (tempG1 + r_Foth_i_G1[ind_i + ind_t*nbI]) * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                        //fichier << "Ceit_G1 = " << rans_C_eit_G1[ind_i + ind_t*nbI] << endl;

                                    rans_C_eit_G2[ind_i + ind_t*nbI] =
                                        (tempG2 + r_Foth_i_G2[ind_i + ind_t*nbI]) * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                        //fichier << "Ceit_G2 = " << rans_C_eit_G2[ind_i + ind_t*nbI] << endl;

                                }

                                //-------------------------

                                double *r_wL_ei_G1 = REAL(getListElement(elmt, "wL_i_G1")),
                                        *r_wL_ei_G2 = REAL(getListElement(elmt, "wL_i_G2")),
                                        *r_wD_ei_G1 = REAL(getListElement(elmt, "wD_i_G1")),
                                        *r_wD_ei_G2 = REAL(getListElement(elmt, "wD_i_G2")),
                                        *r_d_efmit_G1 = REAL(getListElement(elmt, "d_i_G1")),
                                        *r_d_efmit_G2 = REAL(getListElement(elmt, "d_i_G2")),
                                        *doth_eit_G1 = REAL(getListElement(elmt, "doth_i_G1")), //Rprintf("II");//PrintValue(out_oqD_eft);
                                        *doth_eit_G2 = REAL(getListElement(elmt, "doth_i_G2"));


                               //�quation n�2

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wL_ei_G1[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit_G1[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                        r_wL_ei_G2[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit_G2[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                                      rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                }


                              for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_Y_eit[ind_i + ind_t*nbI] =
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                        rans_C_eit_G1[ind_i + ind_t*nbI] / 1000 +

                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                        rans_C_eit_G2[ind_i + ind_t*nbI] / 1000;

                                        //fichier << "Yeit = " << rans_Y_eit[ind_i + ind_t*nbI] << endl;


                                      rans_L_eit[ind_i + ind_t*nbI] = NA_REAL;
                                      //fichier << "Leit_ini = " << rans_L_eit[ind_i + ind_t*nbI] << endl;
                              }
                             //�quation n�3

                            if (all_is_na(v_wD_ei_G1) | all_is_na(v_wD_ei_G2)) {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;
                                if (ISNA(r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        r_wL_ei_G1[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit_G1[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                        r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        r_wL_ei_G2[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit_G2[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                }


                            //Loth_eit
                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit (pas d'exemption)

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000 +

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000 ;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit_G1[ind_i]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000 +

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) * (1-doth_eit_G2[ind_i]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;
                            }


                            } else {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (ISNA(r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;


                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wD_ei_G1[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                        r_d_efmit_G1[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_C_efmit_G1[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000 +

                                        r_wD_ei_G2[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                        r_d_efmit_G2[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_C_efmit_G2[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                }


                               //Loth_eit : pas d'exemption pour les autres si OD

                            if ((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) { //OD s'applique, Loth_eit=Yoth_eit, pas de rejet car pas d'exemption

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000 +

                                        (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] / 1000;


                            } else { //pas d'OD

                             for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                rans_L_eit[ind_i + ind_t*nbI] =
                                        (r_Foth_i_G1[ind_i + ind_t*nbI] * r_N_eit_G1[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G1[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei_G1[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit_G1[ind_i] * r_wD_ei_G1[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000 +

                                         (r_Foth_i_G2[ind_i + ind_t*nbI] * r_N_eit_G2[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit_G2[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) *
                                        (r_wL_ei_G2[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] -
                                         doth_eit_G2[ind_i] * r_wD_ei_G2[0*fact7_C[0]  + 0*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]]) / 1000;

                                    //fichier << "Leit_oth = " << rans_L_eit[ind_i + ind_t*nbI] << endl;
                             }
                            }



                            }



                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd1_efm[ind_f + nbF*ind_m])) {


                           double rYsum = 0.0, rDsum = 0.0;
                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                rYsum = rYsum + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                           }
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                           rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd1_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd2_efm[ind_f + nbF*ind_m])){

                                double rYsum = 0.0, rDsum = 0.0;
                                if (ind_t==0) rYsum=REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m]; else rYsum=rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) rDsum = rDsum + rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                              rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] *
                                     fmin2( r_dd2_efm[ind_f + nbF*ind_m] * finite(rYsum / rDsum) , 1.0) ;

                            if (ISNAN(rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;
                        }

                        } else {

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                             rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                        }
                    }
                 }

//if (nbI==1) {Rprintf("LtotAvant2");PrintValue(out_L_eit);}

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                             rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            if (!ISNA(rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                            rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                rans_LD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; // il reste ensuite � integrer Lefmit pour contituer Ltot_i

                    }


//if (nbI==1) {Rprintf("LtotAvant");PrintValue(out_L_eit);}
                               //�quation n�4

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                      rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                                      if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]])){
                                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] +
                                        rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]; //Ltot constitu�

                                        //fichier << "Leit = " << rans_L_eit[ind_i + ind_t*nbI] << endl;
                                      }

                                }
//if (nbI==1) {Rprintf("LtotApres");PrintValue(out_L_eit);}

                            }


/* insertion over quota management discards pour corriger D et L -> esp�ces dynamiques */

if (!((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & ((activeQR!=0) & (activeQR<=ind_t))) { //pas d'OD appliqu�, et activation du module demand�e


 // on s'occupe d'abord de la partie "autres"

    if (!isNull(getListElement(listQR, CHAR(STRING_ELT(sppList,e)))) & !isNull(getListElement(listQR_f, CHAR(STRING_ELT(sppList,e))))) { //TACs renseign�s aux 2 niveaux

        double *QR = REAL(getListElement(listQR, CHAR(STRING_ELT(sppList,e))));
        double *QR_f = REAL(getListElement(listQR_f, CHAR(STRING_ELT(sppList,e))));

        double QRoth = QR[ind_t];
        bool recal = false;
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) QRoth = QRoth - QR_f[ind_f + nbF*ind_t];

            double Ltot_oth = 0.0, Ytot_oth = 0.0; //, Ytot_othini = 0.0;
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                //double Yothini = rans_Y_eit[ind_i + 0*nbI];
                double Loth = rans_L_eit[ind_i + ind_t*nbI], Yoth = rans_Y_eit[ind_i + ind_t*nbI];
                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){
                //if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                //  Yothini = Yothini - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]];
                if (!ISNA(rans_L_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  Loth = Loth - rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  Yoth = Yoth - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                }
                Ltot_oth = Ltot_oth + Loth; Ytot_oth = Ytot_oth + Yoth; //Ytot_othini = Ytot_othini + Yothini;
            }
            //if (e==1 & ind_t==13) //Rprintf("ind_t %i QRoth %f Ltot_oth %f Ytot_oth %f Ytot_othini %f \n",ind_t,QRoth,Ltot_oth,Ytot_oth,Ytot_othini);

            rans_oqD_et[ind_t] = 0.0;

            if (Ltot_oth>QRoth) { //on proc�de � la correction "autres"

                recal = true;
                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    double Doth_i_t = rans_Y_eit[ind_i + ind_t*nbI] - rans_L_eit[ind_i + ind_t*nbI],// Yoth_i_0 = rans_Y_eit[ind_i + 0*nbI],
                           Yoth_i_t = rans_Y_eit[ind_i + ind_t*nbI];//, Loth_i_t = rans_L_eit[ind_i + ind_t*nbI];
                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++){
                    if (!ISNA(rans_D_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        Doth_i_t = Doth_i_t - rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    //if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                    //    Yoth_i_0 = Yoth_i_0 - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + 0*fact1_C[3]];
                    if (!ISNA(rans_Y_efmit [ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        Yoth_i_t = Yoth_i_t - rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }
                    rans_oqD_et[ind_t] = fmin2((Ltot_oth-QRoth) * (Yoth_i_t - Doth_i_t) / Ltot_oth, Yoth_i_t - Doth_i_t); //fmin2((Ltot_oth-QRoth) * Yoth_i_0 / Ytot_othini, Yoth_i_t - Doth_i_t);
                    Doth_i_t = fmin2( Doth_i_t + (Ltot_oth-QRoth) * (Yoth_i_t - Doth_i_t) / Ltot_oth, Yoth_i_t ); //fmin2( Doth_i_t + (Ltot_oth-QRoth) * Yoth_i_0 / Ytot_othini, Yoth_i_t );

                    if (ISNAN(Doth_i_t)) Doth_i_t = 0.0;
                    if (ISNAN(rans_oqD_et[ind_t])) rans_oqD_et[ind_t] = 0.0;
                    rans_L_eit[ind_i + ind_t*nbI] = Yoth_i_t - Doth_i_t; //on incr�mentera par la suite avec les L recalcul�s

                //if (e==1 & ind_t==13) //Rprintf("Yoth_i_t %f Doth_i_t %f rans_L_eit[ind_i + ind_t*nbI] %f \n",Yoth_i_t,Doth_i_t,rans_L_eit[ind_i + ind_t*nbI]);
                }
            }




            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  { //on proc�de � la correction "flottilles"


            double sumL = 0.0;//, sumYini = 0.0;

            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                  sumL = sumL + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                //if (!ISNA(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]]))
                //  sumYini = sumYini + rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]];
            }

            rans_oqD_eft[ind_f + nbF*ind_t] = 0.0;

            if (sumL>QR_f[ind_f + nbF*ind_t]) { //on proc�de � la correction sur la flottille detect�e

                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    //if (!ISNA(fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                    //        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                    //        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]])))
                    rans_oqD_eft[ind_f + nbF*ind_t] = rans_oqD_eft[ind_f + nbF*ind_t] +
                      fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) / sumL, //rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                            finite(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) -
                            finite(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]));

                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                       fmin2(
                        finite(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) +
                        (sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]]) / sumL, //(sumL-QR_f[ind_f + nbF*ind_t]) * rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + 0*fact1_C[3]] / sumYini,
                        finite(rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]));

                    if (ISNAN(rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]))
                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                    rans_DD_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
//if (e==1 & ind_t==13) //Rprintf("1 rans_L_eit %f rans_L_efmit %f rans_D_efmit %f \n",rans_L_eit[ind_i + ind_t*nbI],rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]],rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]);
                    if (!recal & !ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] - rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }
//if (e==1 & ind_t==13) //Rprintf("2 rans_L_eit %f rans_L_efmit %f rans_D_efmit %f \n",rans_L_eit[ind_i + ind_t*nbI],rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]],rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]);
//fichier1 << "STeee" << endl;

                    rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                       rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                       rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
//if (e==1 & ind_t==13) //Rprintf("3 rans_L_eit %f rans_L_efmit %f rans_D_efmit %f \n",rans_L_eit[ind_i + ind_t*nbI],rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]],rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]);

                    if (!recal & !ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                        rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }
//if (e==1 & ind_t==13) //Rprintf("4 rans_L_eit %f rans_L_efmit %f rans_D_efmit %f \n",rans_L_eit[ind_i + ind_t*nbI],rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]],rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]]);


                }
            }

            if (recal) {
                    //if (e==1 & ind_t==13) //Rprintf("recal\n");

                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                    if (!ISNA(rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]])) {
                      rans_L_eit[ind_i + ind_t*nbI] = rans_L_eit[ind_i + ind_t*nbI] + rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                    }
            }

        }
    }
}

/*----------------------------------------------------------------*/
    //---------
    // calcul de L_et
    //---------
    PROTECT(ans_L_eit = VECTOR_ELT(out_L_eit,e));
    PROTECT(nDim = allocVector(INTSXP,4));
    nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e))); //fichier << "e = " << e << ", ind_e = " << ind_e << ", Name in sppList: " << CHAR(STRING_ELT(sppList,e)) << ", Name in sppListAll: " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

    rans_L_et = REAL(VECTOR_ELT(out_L_et, ind_e));
    rans_L_et [ind_t] = REAL(aggregObj(ans_L_eit,nDim))[ind_t]; //fichier << "rans_L_et = " << rans_L_et [ind_t] << "out_L_et = " << REAL(VECTOR_ELT(out_L_et, ind_e))[ind_t] << endl;



                           if (Svec[e]==0){
                                    UNPROTECT(1);
                            } else {UNPROTECT(2);}

                            UNPROTECT(1+2);
//Rprintf("H20\n");fichier << "H20" << endl;
    }
}
//on passe aux esp�ces statiques //??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
if(VERBOSE){Rprintf(" . Static sp");}
if (nbEstat>0) {

//Rprintf("H21\n");fichier << "H21" << endl;
    for (int e = 0 ; e < nbEstat ; e++) {

                    SEXP elmt, nDim, ans_Ystat;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e))));
//Rprintf("H22\n");fichier << "H22" << endl;
                    double *r_LPUE_eStat = REAL(getListElement(elmt, "LPUE_f_m_e"));
                    double *r_d_eStat = REAL(getListElement(elmt, "d_f_m_e"));
                    double *r_dd1_efm = REAL(getListElement(elmt, "dd1_f_m_e"));
                    double *r_dd2_efm = REAL(getListElement(elmt, "dd2_f_m_e"));
                    double *r_OD_e = REAL(getListElement(elmt, "OD_e"));
                    double *r_dst_efm = REAL(getListElement(elmt, "dst_f_m_e"));

                    //variables d'effort

                    double *r_nbv_f = REAL(getListElement(Flist, "nbv_f_m"));
                    double *r_nbds_f = REAL(getListElement(Flist, "effort1_f_m"));
                    double *r_nbds2_f = REAL(getListElement(Flist, "effort2_f_m"));

//                    int *fFactSup1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 50)),   //ATTENTION : suppose au moins une esp�ce dynamiquement mod�lis�e
//                        *fFactSup2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, 0), 51));

                      int *fFactSup1 = INTEGER(iDim(INTEGER(getAttrib(getListElement(Flist, "nbv_f_m"), install("DimCst"))))),
                          *fFactSup2 = INTEGER(iDim(INTEGER(getAttrib(getListElement(Flist, "effort1_f_m"), install("DimCst")))));
                    int ind_e;
                    int *nd;

                    double *rans_Ystat = REAL(VECTOR_ELT(out_Ystat, e));
                    double *rans_Lstat = REAL(VECTOR_ELT(out_Lstat, e));
                    double *rans_Dstat = REAL(VECTOR_ELT(out_Dstat, e));//Rprintf("LL");

                    double *rans_oqDstat = REAL(VECTOR_ELT(out_oqDstat, e));//Rprintf("MM");

                    double *rans_statDD = REAL(VECTOR_ELT(out_statDD_efm, e));
                    double *rans_statLD = REAL(VECTOR_ELT(out_statLD_efm, e));
                    double *rans_statLDst = REAL(VECTOR_ELT(out_statLDst_efm, e));
                    double *rans_statLDor = REAL(VECTOR_ELT(out_statLDor_efm, e));
                    double *rans_L_et;

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {

                            rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               r_LPUE_eStat[ind_f + nbF*ind_m + 0*ind_t] *
                               r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + 0*fFactSup1[2] + ind_t*fFactSup1[3]] *
                               r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + 0*fFactSup2[2] + ind_t*fFactSup2[3]]*
                               r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + 0*fFactSup2[2] + ind_t*fFactSup2[3]]/1000.0;

                            if (ISNA(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                            rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] / (1 - r_d_eStat[ind_f + nbF*ind_m + 0*ind_t]);

                            if (ISNAN(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                            rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) ;

                            rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*ind_t] +
                                        rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                            rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                               rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * r_d_eStat[ind_f + nbF*ind_m + 0*ind_t];

                            if (ISNA(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                   }
                   //PrintValue(out_Lstat);


                 for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                 for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                   if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd1_efm[ind_f + nbF*ind_m])) {

                           rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                           fmin2( r_dd1_efm[ind_f + nbF*ind_m] * rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] ,
                                  rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] );

                  } else {

                        if (((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & !ISNA(r_dd2_efm[ind_f + nbF*ind_m])){

                                double rYsum = 0.0;
                                if (ind_t==0) {
                                 rYsum = REAL(getListElement(Flist, "Lref_f_m"))[ind_f + nbF*ind_m];
                                } else {
                                 rYsum = rans_Ytot_fm[ind_f + nbF*ind_m + nbF*nbM*(ind_t-1)];
                                }

                                rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = fmin2(
                                    r_dd2_efm[ind_f + nbF*ind_m] * rYsum ,
                                    rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] );

                        } else {

                             rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                        }
                    }
                 }


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t];


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLDst[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * r_dst_efm[ind_f + nbF*ind_m];


                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)

                             rans_statLDor[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                                rans_statLD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_statLDst[ind_f + nbF*ind_m + nbF*nbMe*ind_t];



/* insertion over quota management discards pour corriger D et L -> esp�ces statiques */

if (!((r_OD_e[0]>0.5) & (r_OD_e[0]<=(ind_t+1))) & ((activeQR!=0) & (activeQR<=ind_t))) { //pas d'OD appliqu�, et activation du module demand�e

    if (!isNull(getListElement(listQR_f, CHAR(STRING_ELT(sppListStat,e))))) { //TACs renseign�s au niveau flottille

        double *QR_f = REAL(getListElement(listQR_f, CHAR(STRING_ELT(sppListStat,e))));

            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)  { //on proc�de � la correction "flottilles"

            double sumL = 0.0;//, sumYini = 0.0;

            for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {
                if (!ISNA(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))
                  sumL = sumL + rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                //if (!ISNA(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0]))
                //  sumYini = sumYini + rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0];
            }

            rans_oqDstat[ind_f + nbF*ind_t] = 0.0;

            if (sumL>QR_f[ind_f + nbF*ind_t]) { //on proc�de � la correction sur la flottille detect�e

                for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

                    if (!ISNAN(fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) / sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                                       finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]-rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))))
                    rans_oqDstat[ind_f + nbF*ind_t] = rans_oqDstat[ind_f + nbF*ind_t] +
                                fmin2((sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) / sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                                       finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]-rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]));

                    rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] =
                       fmin2(finite(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]) + (sumL-QR_f[ind_f + nbF*ind_t]) * finite(rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t])/sumL, //rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*0] / sumYini,
                             finite(rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]));

                    if (ISNAN(rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t]))
                        rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;

                    rans_statDD[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                    rans_Lstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = rans_Ystat[ind_f + nbF*ind_m + nbF*nbMe*ind_t] - rans_Dstat[ind_f + nbF*ind_m + nbF*nbMe*ind_t];

                }
            }
        }
    }
}

/*----------------------------------------------------------------*/
    //---------
    //  Calcul L_et
    //---------
    PROTECT(ans_Ystat = VECTOR_ELT(out_Ystat, e));
    PROTECT(nDim = allocVector(INTSXP,4));
    nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e))); //fichier << "e = " << e << ", ind_e = " << ind_e << ", Name in sppListStat: " << CHAR(STRING_ELT(sppListStat,e)) << ", Name in sppListAll: " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

    rans_L_et = REAL(VECTOR_ELT(out_L_et, ind_e));
    rans_L_et [ind_t] = REAL(aggregObj(ans_Ystat,nDim))[ind_t]; //fichier << "rans_L_et = " << rans_L_et[ind_t] << "out_L_et = " << REAL(VECTOR_ELT(out_L_et, ind_e))[ind_t] <<  endl;



                   UNPROTECT(1+2);

         }

}

if(VERBOSE){Rprintf(" . ");}


}
////PrintValue(out_Ystat);

UNPROTECT(1);
//Rprintf("End CatchDL\n");fichier << "End" << endl;

//fichier.close();
}}


