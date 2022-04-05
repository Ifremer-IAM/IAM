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
#include "array_fcts.h" // AggregObj function
// #include "Modules.h" // Contient tout les modules

//using namespace Rcpp;
// using namespace std;

//------------------------------------------
// Module 'EstimationTACfromF'
//------------------------------------------


// TODO : change from to int to void
extern "C" {

int BioEcoPar::EstimationTACfromF(int ind_t)
{


//string str1, str2, str3;
//str1 = "testGestion";//"\\home1\\datahome\\fbriton\\AMURE\\Sc_bug_hke\\debugHKE_V";
//str3 = "_V";
//str2 = ".txt";
//
//std::stringstream ss, mp;
//mp << ind_t;
//str1 = str1 + mp.str()+ str3 + ss.str() + str2;
//
//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.EstimationTAC.txt", ios::out | ios::trunc);
//fichier << "D�but" << endl;

//Rprintf("A1\n");
    IND_T = ind_t;
//spQ = spp;

    //double *totFM, *totFM2, *totF, *totF2, *totFF, *totFF2, *tot, *totMod, *totMod2;

    SEXP listTempP, nDim;

    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;

    PROTECT(listTempP = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));


    double *g_effort1FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m"));
    double *g_effort1F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f"));
    double *g_nbTripFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f_m"));
    double *g_nbTripF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f"));

//Rprintf("A2\n");
//for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//    if ((ind_t==1) & (ind_f==0)) {
//
//        std::stringstream ggg1;
//        ggg1 << g_effort1FM[ind_f + nbF*ind_m];
//
//        fichier << "effort_step1T1" << ggg1.str() << endl;
//
//    }
//}



//for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//    if ((ind_t==1) & (ind_f==0)) {
//
//        std::stringstream ggg2;
//        ggg2 << g_effort1FM[ind_f + nbF*ind_m];
//
//        fichier << "effort_step2T1" << ggg2.str() << endl;
//
//    }
//}

//if (false) PrintValue(inpFtarg);
//if (false) {


int nbEtarg = length(getAttrib(inpFtarg, R_NamesSymbol));

int denom=0, denom2=0;

double *r_N_eit_S1M1=&NA_REAL, *r_N_eit_S2M2=&NA_REAL, *r_N_eit_S3M3=&NA_REAL, *r_N_eit_S4M4=&NA_REAL, *rans_N_eit=&NA_REAL, *r_N_eit_G1=&NA_REAL, *r_N_eit_G2=&NA_REAL,
       *r_N_e0t_S1M1=&NA_REAL, *r_N_e0t_S2M2=&NA_REAL, *r_N_e0t_S3M3=&NA_REAL, *r_N_e0t_S4M4=&NA_REAL, *r_N_e0t=&NA_REAL, *r_N_e0t_G1=&NA_REAL, *r_N_e0t_G2=&NA_REAL, *recValues=&NA_REAL, *LTOT=&NA_REAL,
       *TAC_byFleet=&NA_REAL, *TAC_glob=&NA_REAL, *r_W_Ftarg=&NA_REAL, *r_Qholdings;

double *Fothi2=&NA_REAL, *Fothi2_G1=&NA_REAL, *Fothi2_G2=&NA_REAL;

double *Fothi2_S1M1=&NA_REAL,*Fothi2_S1M2=&NA_REAL,*Fothi2_S1M3=&NA_REAL,*Fothi2_S1M4=&NA_REAL,
       *Fothi2_S2M1=&NA_REAL,*Fothi2_S2M2=&NA_REAL,*Fothi2_S2M3=&NA_REAL,*Fothi2_S2M4=&NA_REAL,
       *Fothi2_S3M1=&NA_REAL,*Fothi2_S3M2=&NA_REAL,*Fothi2_S3M3=&NA_REAL,*Fothi2_S3M4=&NA_REAL,
       *Fothi2_S4M1=&NA_REAL,*Fothi2_S4M2=&NA_REAL,*Fothi2_S4M3=&NA_REAL,*Fothi2_S4M4=&NA_REAL,
       *FRWTothi2_S1M1=&NA_REAL,*FRWTothi2_S1M2=&NA_REAL,*FRWTothi2_S1M3=&NA_REAL,*FRWTothi2_S1M4=&NA_REAL,
       *FRWTothi2_S2M1=&NA_REAL,*FRWTothi2_S2M2=&NA_REAL,*FRWTothi2_S2M3=&NA_REAL,*FRWTothi2_S2M4=&NA_REAL,
       *FRWTothi2_S3M1=&NA_REAL,*FRWTothi2_S3M2=&NA_REAL,*FRWTothi2_S3M3=&NA_REAL,*FRWTothi2_S3M4=&NA_REAL,
       *FRWTothi2_S4M1=&NA_REAL,*FRWTothi2_S4M2=&NA_REAL,*FRWTothi2_S4M3=&NA_REAL,*FRWTothi2_S4M4=&NA_REAL,
       *FDWTothi2_S1M1=&NA_REAL,*FDWTothi2_S1M2=&NA_REAL,*FDWTothi2_S1M3=&NA_REAL,*FDWTothi2_S1M4=&NA_REAL,
       *FDWTothi2_S2M1=&NA_REAL,*FDWTothi2_S2M2=&NA_REAL,*FDWTothi2_S2M3=&NA_REAL,*FDWTothi2_S2M4=&NA_REAL,
       *FDWTothi2_S3M1=&NA_REAL,*FDWTothi2_S3M2=&NA_REAL,*FDWTothi2_S3M3=&NA_REAL,*FDWTothi2_S3M4=&NA_REAL,
       *FDWTothi2_S4M1=&NA_REAL,*FDWTothi2_S4M2=&NA_REAL,*FDWTothi2_S4M3=&NA_REAL,*FDWTothi2_S4M4=&NA_REAL;

double newRec=0.0, newRec_Q1=0.0, newRec_Q2=0.0, newRec_Q3=0.0, newRec_Q4=0.0, newRec_G1=0.0, newRec_G2=0.0;
//Rprintf("%i",nbEtarg);


for (int intEspTarg = 0 ; intEspTarg < nbEtarg ; intEspTarg++) {
//Rprintf("A3\n");
    SEXP namVarTarg, elmt, v_N_e0t, v_N_e0t_S1M1, v_N_e0t_S2M2, v_N_e0t_S3M3, v_N_e0t_S4M4, v_N_e0t_G1, v_N_e0t_G2, v_MeanRec_Ftarg, v_W_Ftarg, v_out_L_eit;
    PROTECT(namVarTarg=STRING_ELT(getAttrib(inpFtarg, R_NamesSymbol),intEspTarg));

    //calcul du ratio Ftarg/Fbar
    double r_Ftarg = REAL(getListElement(inpFtarg, CHAR(namVarTarg)))[IND_T];
    double r_Fbar_prev = REAL(getListElement(out_Fbar_et, CHAR(namVarTarg)))[IND_T-1];
    double r_Fbar_init = REAL(getListElement(out_Fbar_et, CHAR(namVarTarg)))[0];
    //Rprintf("Ftarg: %f, Fbar_init: %f, Fbar_prev: %f \n",r_Ftarg, r_Fbar_init,r_Fbar_prev); fichier << "Ftarg:" << r_Ftarg << ", Fbar_prev:" << r_Fbar_prev << ", Fbar_init:" << r_Fbar_init << endl;

    PROTECT(elmt = getListElement(listTempP, CHAR(namVarTarg)));
    int nbI = length(getListElement(elmt, "modI"));

//Rprintf("A4\n");
    //correction des efforts par le ratio pr�c�dent
    for (int indF = 0 ; indF < nbF ; indF++) {

        for (int indM = 0 ; indM<nbMe ; indM++) {

            g_effort1FM[indF + nbF*indM] = g_effort1FM[indF + nbF*indM] * r_Ftarg / r_Fbar_prev;
            g_nbTripFM[indF + nbF*indM] = g_nbTripFM[indF + nbF*indM] * r_Ftarg / r_Fbar_prev;

        }

        //fichier << "Avant T" << ind_t << ": Effort1 Fleet " << indF << " = " << g_effort1F[indF] << "; Nbtrip = " << g_nbTripF[indF]<< endl;

        g_effort1F[indF] = g_effort1F[indF] * r_Ftarg / r_Fbar_prev;
        g_nbTripF[indF] = g_nbTripF[indF] * r_Ftarg / r_Fbar_prev;

        //fichier << "Apres T" << ind_t << ": Effort1 Fleet " << indF << " = " << g_effort1F[indF] << "; Nbtrip = " << g_nbTripF[indF]<< endl;

    }


    //et correction des mortalit�s autres pour les esp�ces dynamiques XSA, Spict et SS3


            int nbi = length(getListElement(getListElement(list, CHAR(namVarTarg)), "modI"));

            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global

                    // Dans eVarCopy
                    Fothi2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44)); //Rprintf("Dans EVARcopy (l.14478), Fothi2 = "); PrintValue(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi age " << ag <<"=" << Fothi2[ag + IND_T*nbi] << endl;
                            Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi age" << ag << "=" << Fothi2[ag + IND_T*nbi] << endl;
                    }

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 44)); //Rprintf("Dans EVARcopy (l.14478), Fothi2 = "); PrintValue(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi age " << ag <<"=" << Fothi2[ag + IND_T*nbi] << endl;
                            Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi age" << ag << "=" << Fothi2[ag + IND_T*nbi] << endl;
                    }

            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)) {//Age and sex-based

                    // Dans EvarCopy
                    Fothi2_G1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 224));
                    Fothi2_G2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 225));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi_G1 age " << ag <<"=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                            Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi_G1 age" << ag << "=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                    }

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2_G1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 224));
                    Fothi2_G2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 225));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi_G1 age " << ag <<"=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                            Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi_G1 age" << ag << "=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                    }
            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly

                    // Dans eVarCopy
                    Fothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 116));
                    Fothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 117));
                    Fothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 118));
                    Fothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 119));
                    Fothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 120));
                    Fothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 121));
                    Fothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 122));
                    Fothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 123));
                    Fothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 124));
                    Fothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 125));
                    Fothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 126));
                    Fothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 127));
                    Fothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 128));
                    Fothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 129));
                    Fothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 130));
                    Fothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 131));

                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FRWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 176));
                    FRWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 177));
                    FRWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 178));
                    FRWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 179));
                    FRWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 180));
                    FRWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 181));
                    FRWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 182));
                    FRWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 183));
                    FRWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 184));
                    FRWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 185));
                    FRWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 186));
                    FRWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 187));
                    FRWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 188));
                    FRWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 189));
                    FRWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 190));
                    FRWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 191));

                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FDWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 208));
                    FDWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 209));
                    FDWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 210));
                    FDWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 211));
                    FDWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 212));
                    FDWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 213));
                    FDWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 214));
                    FDWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 215));
                    FDWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 216));
                    FDWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 217));
                    FDWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 218));
                    FDWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 219));
                    FDWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 220));
                    FDWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 221));
                    FDWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 222));
                    FDWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 223));

                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 116));
                    Fothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 117));
                    Fothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 118));
                    Fothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 119));
                    Fothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 120));
                    Fothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 121));
                    Fothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 122));
                    Fothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 123));
                    Fothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 124));
                    Fothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 125));
                    Fothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 126));
                    Fothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 127));
                    Fothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 128));
                    Fothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 129));
                    Fothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 130));
                    Fothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 131));

                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FRWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 176));
                    FRWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 177));
                    FRWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 178));
                    FRWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 179));
                    FRWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 180));
                    FRWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 181));
                    FRWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 182));
                    FRWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 183));
                    FRWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 184));
                    FRWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 185));
                    FRWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 186));
                    FRWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 187));
                    FRWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 188));
                    FRWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 189));
                    FRWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 190));
                    FRWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 191));

                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FDWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 208));
                    FDWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 209));
                    FDWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 210));
                    FDWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 211));
                    FDWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 212));
                    FDWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 213));
                    FDWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 214));
                    FDWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 215));
                    FDWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 216));
                    FDWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 217));
                    FDWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 218));
                    FDWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 219));
                    FDWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 220));
                    FDWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 221));
                    FDWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 222));
                    FDWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 223));

                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;


            }

    // envoi du module Mortalit�
    //Rprintf("call.Mortalite\n");
   // fichier << "call.Mortalite" << endl;
    Mortalite(listTempP, IND_T, eVarCopy);
    //Rprintf("End.Mortalite\n");
    //fichier << "End.Mortalite" << endl;

    // calcul recrutement
    if ((nbI>1) && !isNull(inpMeanRec_Ftarg) && (!isNull(getListElement(inpMeanRec_Ftarg,CHAR(namVarTarg))))){ // si MeanRec_Ftarg renseigne: forcage selon type 1 (Moyenne sur X dernieres annees) ou 2 (For�age avec valeurs renseignees)

        //fichier << "Recruitment in HCR = from MeanRecFtarg" << endl;

            PROTECT(v_MeanRec_Ftarg = getListElement(inpMeanRec_Ftarg, CHAR(namVarTarg)));
            //PrintValue(v_MeanRec_Ftarg);
        //Rprintf("A5\n");
            if (length(v_MeanRec_Ftarg)==1) {

                denom = INTEGER(v_MeanRec_Ftarg)[0];
                denom2 = denom;

                if (!ISNA(denom)){
                    if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global
        //Rprintf("A6\n");
                      newRec = 0.0;
                      rans_N_eit = REAL(getListElement(out_N_eit,CHAR(namVarTarg)));
                      for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec = newRec + rans_N_eit[(IND_T-index)*nbI] ;
                                                                 }}
                      newRec = newRec / denom2;
                      PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                      r_N_e0t = REAL(v_N_e0t);
                      r_N_e0t[IND_T] = newRec;
                      rans_N_eit[IND_T*nbI] = newRec;
                      UNPROTECT(1);
        //Rprintf("A7\n");
        //fichier << "Type 1 (moyenne) MeanRec: " << rans_N_eit[IND_T*nbI] << endl;

                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)) {//Ageand sex-based
                        r_N_eit_G1 = REAL(getListElement(out_N_eit_G1,CHAR(namVarTarg)));
                        r_N_eit_G2 = REAL(getListElement(out_N_eit_G2,CHAR(namVarTarg)));
                        newRec_G1 = 0.0; newRec_G2 = 0.0;

                        for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec_G1 = newRec_G1 + r_N_eit_G1[(IND_T-index)*nbI] ;
                                                                   newRec_G2 = newRec_G2 + r_N_eit_G2[(IND_T-index)*nbI] ;
                                                                 }}
                        newRec_G1 = newRec_G1 / denom2 ;
                        newRec_G2 = newRec_G2 / denom2 ;

                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2"));

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = newRec_G1;
                        r_N_e0t_G2[IND_T] = newRec_G2;

                        UNPROTECT(2);
       // fichier << "Type 1 (moyenne) MeanRecG1: " << r_N_e0t_G1[IND_T] << "/ MeanRecG2: " << r_N_e0t_G2[IND_T]<< endl;
                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly
       // Rprintf("A8\n");
                      //Quarter 1
                        r_N_eit_S1M1 = REAL(getListElement(out_N_eit_S1M1,CHAR(namVarTarg)));
                        r_N_eit_S2M2 = REAL(getListElement(out_N_eit_S2M2,CHAR(namVarTarg)));
                        r_N_eit_S3M3 = REAL(getListElement(out_N_eit_S3M3,CHAR(namVarTarg)));
                        r_N_eit_S4M4 = REAL(getListElement(out_N_eit_S4M4,CHAR(namVarTarg)));
                        newRec_Q1 = 0.0; newRec_Q2 = 0.0; newRec_Q3 = 0.0; newRec_Q4 = 0.0;
        //Rprintf("A9\n");
                        for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec_Q1 = newRec_Q1 + r_N_eit_S1M1[(IND_T-index)*nbI] ;
                                                                   newRec_Q2 = newRec_Q2 + r_N_eit_S2M2[(IND_T-index)*nbI] ;
                                                                   newRec_Q3 = newRec_Q3 + r_N_eit_S3M3[(IND_T-index)*nbI] ;
                                                                   newRec_Q4 = newRec_Q4 + r_N_eit_S4M4[(IND_T-index)*nbI] ;
                                                                 }}
                        newRec_Q1 = newRec_Q1 / denom2 ;
                        newRec_Q2 = newRec_Q2 / denom2 ;
                        newRec_Q3 = newRec_Q3 / denom2 ;
                        newRec_Q4 = newRec_Q4 / denom2 ;
        //Rprintf("A10\n");
                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = newRec_Q1;
                        r_N_e0t_S2M2[0] = newRec_Q2;
                        r_N_e0t_S3M3[0] = newRec_Q3;
                        r_N_e0t_S4M4[0] = newRec_Q4;

                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...
                        UNPROTECT(4);

                    }

                }
            } else {

              if (length(v_MeanRec_Ftarg)>1) {  //historique XSA ou SS3
        //Rprintf("A11\n");
                recValues = REAL(v_MeanRec_Ftarg);
       // Rprintf("A110\n");
                if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global
        //Rprintf("A12\n");
                      PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                      r_N_e0t = REAL(v_N_e0t);
                      r_N_e0t[IND_T] = recValues[IND_T];
                      UNPROTECT(1);
        //Rprintf("A13\n");
        //fichier << "Type 2 (historique) MeanRec: " << rans_N_eit[IND_T*nbI] << endl;

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){//age and sex-based
        //Rprintf("A14\n");
                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1")); //PrintValue(v_N_e0t_G1);
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2")); //PrintValue(v_N_e0t_G2);

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = recValues[2*IND_T];
                        r_N_e0t_G2[IND_T] = recValues[2*IND_T + 1];
       // Rprintf("A15\n");
       // fichier << "Type 2 (historique) MeanRecG1: " << r_N_e0t_G1[IND_T] << "/ MeanRecG2: " << r_N_e0t_G2[IND_T]<< endl;

                       UNPROTECT(2);

                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly
     //   Rprintf("A14\n");
                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = recValues[4*IND_T];
                        r_N_e0t_S2M2[0] = recValues[4*IND_T + 1];
                        r_N_e0t_S3M3[0] = recValues[4*IND_T + 2];
                        r_N_e0t_S4M4[0] = recValues[4*IND_T + 3];
        //Rprintf("A15\n");
                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...

                       UNPROTECT(4);

                    }

              }

            }

            UNPROTECT(1);

    } else if ((nbI>1) && (!isNull(getListElement(recParamList,CHAR(namVarTarg)))) && (IND_T>0)) { // recrutement = recrutement attendu avec la relation SR (sans incertitude)
//fichier << "Recruitment in HCR = expected from SR" << endl;

            double  *rans_SSB_et = REAL(getListElement(out_SSB_et,CHAR(namVarTarg)));

            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
                double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type")); //Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));

                if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                    double recr = 0.0;

                    if (typeSR[IND_T]==1){ // Hockey Stick
                        if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T - 1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                            recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                        } else{
                            recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                        }
                    } else if (typeSR[IND_T]==2){ // Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                        recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T -1 - del]) /
                        (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T -1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                        param[IND_T + 3*nbT] ;
                    }

                    PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                    r_N_e0t = REAL(v_N_e0t);
                    r_N_e0t[IND_T] = recr; //fichier << "e = " << CHAR(namVarTarg) << " recr = " << recr << ", N_i0t = " << r_N_e0t[IND_T] << endl;
                    UNPROTECT(1);
                }

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
                    double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                    int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));// Rprintf("type = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                    int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));
                    double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil")); //Rprintf("ventil = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil"));

                    if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                        double recr = 0.0;

                        if (typeSR[IND_T]==1){ // Hockey Stick
                            if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T - 1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                                recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                            } else{
                                recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                            }
                        } else if (typeSR[IND_T]==2){ //Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                            recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T -1 - del]) /
                            (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T -1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                            param[IND_T + 3*nbT];
                        }

                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = recr*ventil[0];
                        r_N_e0t_S2M2[0] =  recr*ventil[1];
                        r_N_e0t_S3M3[0] =  recr*ventil[2];
                        r_N_e0t_S4M4[0] =  recr*ventil[3];
                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...

                       UNPROTECT(4);

                    }

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){
                    double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                    int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));// Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                    int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));
                    double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil")); //Rprintf("ventil = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil"));

                    if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                        double recr = 0.0;
                        //fichier << "e = " << CHAR(namVarTarg) << " SSB = " << rans_SSB_et[IND_T - 1 - del] << endl;

                        if (typeSR[IND_T]==1){ // Hockey Stick
                            if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T -1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                                recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                            } else{
                                recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                            }
                        } else if (typeSR[IND_T]==2){ //Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                            //fichier << " Beverton Holt" << endl;
                            recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T - 1 - del]) /
                            (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T - 1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                            param[IND_T + 3*nbT];
                        }

                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));  //PrintValue(v_N_e0t_G1);
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2")); //PrintValue(v_N_e0t_G2);

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = recr*ventil[0]; //fichier << "e = " << CHAR(namVarTarg) << " recr = " << recr << ", N_i0t_G1 = " << r_N_e0t_G1[IND_T] << endl;
                        r_N_e0t_G2[IND_T] = recr*ventil[1];

                         UNPROTECT(2);

                    }
                }

        }



//Rprintf("A16\n");
 //   Rprintf("call.DynamicPop\n");
 //   fichier << "call.DynamicPop" << endl;
    DynamicPop(listTempP, IND_T, eVarCopy, false);
 //   Rprintf("end.DynamicPop\n");
 //   fichier << "end.DynamicPop" << endl;
//Rprintf("A17\n");
  //  Rprintf("call.CatchDL\n");
  //  fichier << "call.CatchDL" << endl;
    CatchDL(listTempP, IND_T, eVarCopy, 0);
 //   Rprintf("end.CatchDL\n");
 //   fichier << "end.CatchDL" << endl;
//Rprintf("A18\n");
    //on peut d�sormais d�duire des d�barquements mod�lis�s les TAC par flottille et totaux

    PROTECT(v_W_Ftarg = getListElement(inpW_Ftarg, CHAR(namVarTarg)));
//Rprintf("A19\n");
    PROTECT(v_out_L_eit = getListElement(out_L_eit, CHAR(namVarTarg)));
//Rprintf("A20\n");
    LTOT = REAL(aggregObj(v_out_L_eit,nDim));
//Rprintf("A21\n");
    //if (IND_T==2) {Rprintf("AA\n"); PrintValue(TACbyF);PrintValue(TAC);}

    TAC_byFleet = REAL(getListElement(TACbyF, CHAR(namVarTarg)));
    TAC_glob = REAL(getListElement(TAC, CHAR(namVarTarg)));
    r_W_Ftarg = REAL(v_W_Ftarg);
    r_Qholdings = REAL(getListElement(Qholdings, CHAR(namVarTarg)));
//Rprintf("A22\n");

    TAC_glob[IND_T] = LTOT[IND_T];
    //fichier << "Ltot t-1: " << LTOT[IND_T-1] << endl;
    //fichier << "Tactot t: " << TAC_glob[IND_T] << endl;

    for (int indF = 0 ; indF < nbF ; indF++) TAC_byFleet[indF + nbF*IND_T] = r_W_Ftarg[indF + (nbF+1)*IND_T] * LTOT[IND_T]; // for use in Gestion F2: only modelled fleets

    for (int indF = 0 ; indF <= nbF ; indF++) r_Qholdings[indF + (nbF+1)*IND_T] = r_W_Ftarg[indF + (nbF+1)*IND_T] * LTOT[IND_T]; // for use in quota trading, contains also external investors
    //PrintValue(getListElement(Qholdings, CHAR(namVarTarg)));
        //re-correction des efforts par l'inverse du ratio pr�c�dent
    for (int indF = 0 ; indF < nbF ; indF++) {

        for (int indM = 0 ; indM<nbMe ; indM++) {

            g_effort1FM[indF + nbF*indM] = g_effort1FM[indF + nbF*indM] * r_Fbar_prev/ r_Ftarg;
            g_nbTripFM[indF + nbF*indM] = g_nbTripFM[indF + nbF*indM] * r_Fbar_prev / r_Ftarg;

        }

        g_effort1F[indF] = g_effort1F[indF] * r_Fbar_prev / r_Ftarg;
        g_nbTripF[indF] = g_nbTripF[indF] * r_Fbar_prev / r_Ftarg;

    }

        //et re-correction des mortalit�s autres pour les esp�ces dynamiques XSA, Spict et SS3

//            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {
//
//                    for (int ag = 0; ag < nbi; ag++) Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){
//                for (int ag = 0; ag < nbi; ag++) {
//                        Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                        Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                }
//
//                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
//
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//            }

     UNPROTECT(4);


}
//UNPROTECT(2);
 //}
     UNPROTECT(3);


   return(0);
   //fichier.close();


  }
}

