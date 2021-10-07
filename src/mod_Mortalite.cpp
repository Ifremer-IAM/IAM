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

// #include "Modules.h" // Contient tout les modules
#include "BioEcoPar.h" // Class is defined in this file.

//using namespace Rcpp;
// using namespace std;

//------------------------------------------
// Module 'Mortalit� par p�che et survie des rejets'
//------------------------------------------

extern "C" {

void BioEcoPar::Mortalite(SEXP list, int ind_t, SEXP EVAR)
{

SEXP Flist;
PROTECT(Flist = getListElement(list, "Fleet"));


//        string str1c, str2c, str3c;
//        str1c = "debugHKEmort_V";//"C:\\Users\\mmerzere\\Desktop\\test3\\debugHKEmort_V";//str1b = "\\home1\\datahome\\fbriton\\AMURE\\Sc_bug_hke\\debugHKEcatch_V";
//        str3c = "_V";
//        str2c = ".txt";
//
//        std::stringstream ssj, ssz;
//        ssj << ind_t;
//        ssz << EcoIndCopy[0];
//        str1c = str1c + ssj.str() + str3c + ssz.str() + str2c;
//
//        ofstream fichier2(str1c.c_str() , ios::out | ios::trunc);
//
//ofstream fichier2("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test_Mortalite.txt", ios::out | ios::trunc);
//fichier2 << "D�but " << endl;

if (fUpdate) {
        //fichier2 << "fUpdate = " << fUpdate << endl;

SEXP    elmt, dimEff,
        dimCst, Dim, dimCst_Sr_e= R_NilValue, dimCst_d_efi= R_NilValue, dimCst_doth_ei= R_NilValue, dimCst_F_efmi= R_NilValue, intAge, //dimCst_Capt_emi, dimCst_Capt_ei,
        v_Sr_e= R_NilValue, v_d_efi= R_NilValue, v_doth_ei= R_NilValue, v_F_efmi = R_NilValue, v_F_efmi2 = R_NilValue, formatEff, dimCstEff, //rDim, v_Capt_emi, v_Capt_ei,
        v_d_efi_G1= R_NilValue, v_d_efi_G2= R_NilValue, v_doth_ei_G1= R_NilValue, v_doth_ei_G2= R_NilValue,
        v_F_efmi_G1 = R_NilValue, v_F_efmi_G2 = R_NilValue, v_F_efmi2_G1 = R_NilValue, v_F_efmi2_G2 = R_NilValue,

        v_F_efmi_S1M1 = R_NilValue, v_F_efmi_S1M2 = R_NilValue, v_F_efmi_S1M3 = R_NilValue, v_F_efmi_S1M4 = R_NilValue,
        v_F_efmi_S2M1 = R_NilValue, v_F_efmi_S2M2 = R_NilValue, v_F_efmi_S2M3 = R_NilValue, v_F_efmi_S2M4 = R_NilValue,
        v_F_efmi_S3M1 = R_NilValue, v_F_efmi_S3M2 = R_NilValue, v_F_efmi_S3M3 = R_NilValue, v_F_efmi_S3M4 = R_NilValue,
        v_F_efmi_S4M1 = R_NilValue, v_F_efmi_S4M2 = R_NilValue, v_F_efmi_S4M3 = R_NilValue, v_F_efmi_S4M4 = R_NilValue,
        v_F_efmi2_S1M1 = R_NilValue, v_F_efmi2_S1M2 = R_NilValue, v_F_efmi2_S1M3 = R_NilValue, v_F_efmi2_S1M4 = R_NilValue,
        v_F_efmi2_S2M1 = R_NilValue, v_F_efmi2_S2M2 = R_NilValue, v_F_efmi2_S2M3 = R_NilValue, v_F_efmi2_S2M4 = R_NilValue,
        v_F_efmi2_S3M1 = R_NilValue, v_F_efmi2_S3M2 = R_NilValue, v_F_efmi2_S3M3 = R_NilValue, v_F_efmi2_S3M4 = R_NilValue,
        v_F_efmi2_S4M1 = R_NilValue, v_F_efmi2_S4M2 = R_NilValue, v_F_efmi2_S4M3 = R_NilValue, v_F_efmi2_S4M4 = R_NilValue,
//        v_Foth_i_S1M1 = R_NilValue, v_Foth_i_S1M2 = R_NilValue, v_Foth_i_S1M3 = R_NilValue, v_Foth_i_S1M4 = R_NilValue,
//        v_Foth_i_S2M1 = R_NilValue, v_Foth_i_S2M2 = R_NilValue, v_Foth_i_S2M3 = R_NilValue, v_Foth_i_S2M4 = R_NilValue,
//        v_Foth_i_S3M1 = R_NilValue, v_Foth_i_S3M2 = R_NilValue, v_Foth_i_S3M3 = R_NilValue, v_Foth_i_S3M4 = R_NilValue,
//        v_Foth_i_S4M1 = R_NilValue, v_Foth_i_S4M2 = R_NilValue, v_Foth_i_S4M3 = R_NilValue, v_Foth_i_S4M4 = R_NilValue,
//        v_Froth_i_S1M1 = R_NilValue, v_Froth_i_S1M2 = R_NilValue, v_Froth_i_S1M3 = R_NilValue, v_Froth_i_S1M4 = R_NilValue,
//        v_Froth_i_S2M1 = R_NilValue, v_Froth_i_S2M2 = R_NilValue, v_Froth_i_S2M3 = R_NilValue, v_Froth_i_S2M4 = R_NilValue,
//        v_Froth_i_S3M1 = R_NilValue, v_Froth_i_S3M2 = R_NilValue, v_Froth_i_S3M3 = R_NilValue, v_Froth_i_S3M4 = R_NilValue,
//        v_Froth_i_S4M1 = R_NilValue, v_Froth_i_S4M2 = R_NilValue, v_Froth_i_S4M3 = R_NilValue, v_Froth_i_S4M4 = R_NilValue,
        Foth_i_S1M1 = R_NilValue, Foth_i_S1M2 = R_NilValue, Foth_i_S1M3 = R_NilValue, Foth_i_S1M4 = R_NilValue,
        Foth_i_S2M1 = R_NilValue, Foth_i_S2M2 = R_NilValue, Foth_i_S2M3 = R_NilValue, Foth_i_S2M4 = R_NilValue,
        Foth_i_S3M1 = R_NilValue, Foth_i_S3M2 = R_NilValue, Foth_i_S3M3 = R_NilValue, Foth_i_S3M4 = R_NilValue,
        Foth_i_S4M1 = R_NilValue, Foth_i_S4M2 = R_NilValue, Foth_i_S4M3 = R_NilValue, Foth_i_S4M4 = R_NilValue,
        Froth_i_S1M1 = R_NilValue, Froth_i_S1M2 = R_NilValue, Froth_i_S1M3 = R_NilValue, Froth_i_S1M4 = R_NilValue,
        Froth_i_S2M1 = R_NilValue, Froth_i_S2M2 = R_NilValue, Froth_i_S2M3 = R_NilValue, Froth_i_S2M4 = R_NilValue,
        Froth_i_S3M1 = R_NilValue, Froth_i_S3M2 = R_NilValue, Froth_i_S3M3 = R_NilValue, Froth_i_S3M4 = R_NilValue,
        Froth_i_S4M1 = R_NilValue, Froth_i_S4M2 = R_NilValue, Froth_i_S4M3 = R_NilValue, Froth_i_S4M4 = R_NilValue,

        v_FRWT_efmi_S1M1 = R_NilValue, v_FRWT_efmi_S1M2 = R_NilValue, v_FRWT_efmi_S1M3 = R_NilValue, v_FRWT_efmi_S1M4 = R_NilValue,
        v_FRWT_efmi_S2M1 = R_NilValue, v_FRWT_efmi_S2M2 = R_NilValue, v_FRWT_efmi_S2M3 = R_NilValue, v_FRWT_efmi_S2M4 = R_NilValue,
        v_FRWT_efmi_S3M1 = R_NilValue, v_FRWT_efmi_S3M2 = R_NilValue, v_FRWT_efmi_S3M3 = R_NilValue, v_FRWT_efmi_S3M4 = R_NilValue,
        v_FRWT_efmi_S4M1 = R_NilValue, v_FRWT_efmi_S4M2 = R_NilValue, v_FRWT_efmi_S4M3 = R_NilValue, v_FRWT_efmi_S4M4 = R_NilValue,
        v_FRWT_efmi2_S1M1 = R_NilValue, v_FRWT_efmi2_S1M2 = R_NilValue, v_FRWT_efmi2_S1M3 = R_NilValue, v_FRWT_efmi2_S1M4 = R_NilValue,
        v_FRWT_efmi2_S2M1 = R_NilValue, v_FRWT_efmi2_S2M2 = R_NilValue, v_FRWT_efmi2_S2M3 = R_NilValue, v_FRWT_efmi2_S2M4 = R_NilValue,
        v_FRWT_efmi2_S3M1 = R_NilValue, v_FRWT_efmi2_S3M2 = R_NilValue, v_FRWT_efmi2_S3M3 = R_NilValue, v_FRWT_efmi2_S3M4 = R_NilValue,
        v_FRWT_efmi2_S4M1 = R_NilValue, v_FRWT_efmi2_S4M2 = R_NilValue, v_FRWT_efmi2_S4M3 = R_NilValue, v_FRWT_efmi2_S4M4 = R_NilValue,
        FRWToth_i_S1M1 = R_NilValue, FRWToth_i_S1M2 = R_NilValue, FRWToth_i_S1M3 = R_NilValue, FRWToth_i_S1M4 = R_NilValue,
        FRWToth_i_S2M1 = R_NilValue, FRWToth_i_S2M2 = R_NilValue, FRWToth_i_S2M3 = R_NilValue, FRWToth_i_S2M4 = R_NilValue,
        FRWToth_i_S3M1 = R_NilValue, FRWToth_i_S3M2 = R_NilValue, FRWToth_i_S3M3 = R_NilValue, FRWToth_i_S3M4 = R_NilValue,
        FRWToth_i_S4M1 = R_NilValue, FRWToth_i_S4M2 = R_NilValue, FRWToth_i_S4M3 = R_NilValue, FRWToth_i_S4M4 = R_NilValue,
//        FRWTroth_i_S1M1 = R_NilValue, FRWTroth_i_S1M2 = R_NilValue, FRWTroth_i_S1M3 = R_NilValue, FRWTroth_i_S1M4 = R_NilValue,
//        FRWTroth_i_S2M1 = R_NilValue, FRWTroth_i_S2M2 = R_NilValue, FRWTroth_i_S2M3 = R_NilValue, FRWTroth_i_S2M4 = R_NilValue,
//        FRWTroth_i_S3M1 = R_NilValue, FRWTroth_i_S3M2 = R_NilValue, FRWTroth_i_S3M3 = R_NilValue, FRWTroth_i_S3M4 = R_NilValue,
//        FRWTroth_i_S4M1 = R_NilValue, FRWTroth_i_S4M2 = R_NilValue, FRWTroth_i_S4M3 = R_NilValue, FRWTroth_i_S4M4 = R_NilValue,

        v_FDWT_efmi_S1M1 = R_NilValue, v_FDWT_efmi_S1M2 = R_NilValue, v_FDWT_efmi_S1M3 = R_NilValue, v_FDWT_efmi_S1M4 = R_NilValue,
        v_FDWT_efmi_S2M1 = R_NilValue, v_FDWT_efmi_S2M2 = R_NilValue, v_FDWT_efmi_S2M3 = R_NilValue, v_FDWT_efmi_S2M4 = R_NilValue,
        v_FDWT_efmi_S3M1 = R_NilValue, v_FDWT_efmi_S3M2 = R_NilValue, v_FDWT_efmi_S3M3 = R_NilValue, v_FDWT_efmi_S3M4 = R_NilValue,
        v_FDWT_efmi_S4M1 = R_NilValue, v_FDWT_efmi_S4M2 = R_NilValue, v_FDWT_efmi_S4M3 = R_NilValue, v_FDWT_efmi_S4M4 = R_NilValue,
        v_FDWT_efmi2_S1M1 = R_NilValue, v_FDWT_efmi2_S1M2 = R_NilValue, v_FDWT_efmi2_S1M3 = R_NilValue, v_FDWT_efmi2_S1M4 = R_NilValue,
        v_FDWT_efmi2_S2M1 = R_NilValue, v_FDWT_efmi2_S2M2 = R_NilValue, v_FDWT_efmi2_S2M3 = R_NilValue, v_FDWT_efmi2_S2M4 = R_NilValue,
        v_FDWT_efmi2_S3M1 = R_NilValue, v_FDWT_efmi2_S3M2 = R_NilValue, v_FDWT_efmi2_S3M3 = R_NilValue, v_FDWT_efmi2_S3M4 = R_NilValue,
        v_FDWT_efmi2_S4M1 = R_NilValue, v_FDWT_efmi2_S4M2 = R_NilValue, v_FDWT_efmi2_S4M3 = R_NilValue, v_FDWT_efmi2_S4M4 = R_NilValue,
        FDWToth_i_S1M1 = R_NilValue, FDWToth_i_S1M2 = R_NilValue, FDWToth_i_S1M3 = R_NilValue, FDWToth_i_S1M4 = R_NilValue,
        FDWToth_i_S2M1 = R_NilValue, FDWToth_i_S2M2 = R_NilValue, FDWToth_i_S2M3 = R_NilValue, FDWToth_i_S2M4 = R_NilValue,
        FDWToth_i_S3M1 = R_NilValue, FDWToth_i_S3M2 = R_NilValue, FDWToth_i_S3M3 = R_NilValue, FDWToth_i_S3M4 = R_NilValue,
        FDWToth_i_S4M1 = R_NilValue, FDWToth_i_S4M2 = R_NilValue, FDWToth_i_S4M3 = R_NilValue, FDWToth_i_S4M4 = R_NilValue,
//        FDWTroth_i_S1M1 = R_NilValue, FDWTroth_i_S1M2 = R_NilValue, FDWTroth_i_S1M3 = R_NilValue, FDWTroth_i_S1M4 = R_NilValue,
//        FDWTroth_i_S2M1 = R_NilValue, FDWTroth_i_S2M2 = R_NilValue, FDWTroth_i_S2M3 = R_NilValue, FDWTroth_i_S2M4 = R_NilValue,
//        FDWTroth_i_S3M1 = R_NilValue, FDWTroth_i_S3M2 = R_NilValue, FDWTroth_i_S3M3 = R_NilValue, FDWTroth_i_S3M4 = R_NilValue,
//        FDWTroth_i_S4M1 = R_NilValue, FDWTroth_i_S4M2 = R_NilValue, FDWTroth_i_S4M3 = R_NilValue, FDWTroth_i_S4M4 = R_NilValue,


        v_iniF_efmi_S1M1 = R_NilValue, v_iniF_efmi_S1M2 = R_NilValue, v_iniF_efmi_S1M3 = R_NilValue, v_iniF_efmi_S1M4 = R_NilValue,
        v_iniF_efmi_S2M1 = R_NilValue, v_iniF_efmi_S2M2 = R_NilValue, v_iniF_efmi_S2M3 = R_NilValue, v_iniF_efmi_S2M4 = R_NilValue,
        v_iniF_efmi_S3M1 = R_NilValue, v_iniF_efmi_S3M2 = R_NilValue, v_iniF_efmi_S3M3 = R_NilValue, v_iniF_efmi_S3M4 = R_NilValue,
        v_iniF_efmi_S4M1 = R_NilValue, v_iniF_efmi_S4M2 = R_NilValue, v_iniF_efmi_S4M3 = R_NilValue, v_iniF_efmi_S4M4 = R_NilValue,
        v_iniF_efmi2_S1M1 = R_NilValue, v_iniF_efmi2_S1M2 = R_NilValue, v_iniF_efmi2_S1M3 = R_NilValue, v_iniF_efmi2_S1M4 = R_NilValue,
        v_iniF_efmi2_S2M1 = R_NilValue, v_iniF_efmi2_S2M2 = R_NilValue, v_iniF_efmi2_S2M3 = R_NilValue, v_iniF_efmi2_S2M4 = R_NilValue,
        v_iniF_efmi2_S3M1 = R_NilValue, v_iniF_efmi2_S3M2 = R_NilValue, v_iniF_efmi2_S3M3 = R_NilValue, v_iniF_efmi2_S3M4 = R_NilValue,
        v_iniF_efmi2_S4M1 = R_NilValue, v_iniF_efmi2_S4M2 = R_NilValue, v_iniF_efmi2_S4M3 = R_NilValue, v_iniF_efmi2_S4M4 = R_NilValue,
//        v_iniFoth_i_S1M1 = R_NilValue, v_iniFoth_i_S1M2 = R_NilValue, v_iniFoth_i_S1M3 = R_NilValue, v_iniFoth_i_S1M4 = R_NilValue,
//        v_iniFoth_i_S2M1 = R_NilValue, v_iniFoth_i_S2M2 = R_NilValue, v_iniFoth_i_S2M3 = R_NilValue, v_iniFoth_i_S2M4 = R_NilValue,
//        v_iniFoth_i_S3M1 = R_NilValue, v_iniFoth_i_S3M2 = R_NilValue, v_iniFoth_i_S3M3 = R_NilValue, v_iniFoth_i_S3M4 = R_NilValue,
//        v_iniFoth_i_S4M1 = R_NilValue, v_iniFoth_i_S4M2 = R_NilValue, v_iniFoth_i_S4M3 = R_NilValue, v_iniFoth_i_S4M4 = R_NilValue,
//        v_iniFroth_i_S1M1 = R_NilValue, v_iniFroth_i_S1M2 = R_NilValue, v_iniFroth_i_S1M3 = R_NilValue, v_iniFroth_i_S1M4 = R_NilValue,
//        v_iniFroth_i_S2M1 = R_NilValue, v_iniFroth_i_S2M2 = R_NilValue, v_iniFroth_i_S2M3 = R_NilValue, v_iniFroth_i_S2M4 = R_NilValue,
//        v_iniFroth_i_S3M1 = R_NilValue, v_iniFroth_i_S3M2 = R_NilValue, v_iniFroth_i_S3M3 = R_NilValue, v_iniFroth_i_S3M4 = R_NilValue,
//        v_iniFroth_i_S4M1 = R_NilValue, v_iniFroth_i_S4M2 = R_NilValue, v_iniFroth_i_S4M3 = R_NilValue, v_iniFroth_i_S4M4 = R_NilValue,
//        iniFoth_i_S1M1 = R_NilValue, iniFoth_i_S1M2 = R_NilValue, iniFoth_i_S1M3 = R_NilValue, iniFoth_i_S1M4 = R_NilValue,
//        iniFoth_i_S2M1 = R_NilValue, iniFoth_i_S2M2 = R_NilValue, iniFoth_i_S2M3 = R_NilValue, iniFoth_i_S2M4 = R_NilValue,
//        iniFoth_i_S3M1 = R_NilValue, iniFoth_i_S3M2 = R_NilValue, iniFoth_i_S3M3 = R_NilValue, iniFoth_i_S3M4 = R_NilValue,
//        iniFoth_i_S4M1 = R_NilValue, iniFoth_i_S4M2 = R_NilValue, iniFoth_i_S4M3 = R_NilValue, iniFoth_i_S4M4 = R_NilValue,
//        iniFroth_i_S1M1 = R_NilValue, iniFroth_i_S1M2 = R_NilValue, iniFroth_i_S1M3 = R_NilValue, iniFroth_i_S1M4 = R_NilValue,
//        iniFroth_i_S2M1 = R_NilValue, iniFroth_i_S2M2 = R_NilValue, iniFroth_i_S2M3 = R_NilValue, iniFroth_i_S2M4 = R_NilValue,
//        iniFroth_i_S3M1 = R_NilValue, iniFroth_i_S3M2 = R_NilValue, iniFroth_i_S3M3 = R_NilValue, iniFroth_i_S3M4 = R_NilValue,
//        iniFroth_i_S4M1 = R_NilValue, iniFroth_i_S4M2 = R_NilValue, iniFroth_i_S4M3 = R_NilValue, iniFroth_i_S4M4 = R_NilValue,

        v_iniFRWT_efmi_S1M1 = R_NilValue, v_iniFRWT_efmi_S1M2 = R_NilValue, v_iniFRWT_efmi_S1M3 = R_NilValue, v_iniFRWT_efmi_S1M4 = R_NilValue,
        v_iniFRWT_efmi_S2M1 = R_NilValue, v_iniFRWT_efmi_S2M2 = R_NilValue, v_iniFRWT_efmi_S2M3 = R_NilValue, v_iniFRWT_efmi_S2M4 = R_NilValue,
        v_iniFRWT_efmi_S3M1 = R_NilValue, v_iniFRWT_efmi_S3M2 = R_NilValue, v_iniFRWT_efmi_S3M3 = R_NilValue, v_iniFRWT_efmi_S3M4 = R_NilValue,
        v_iniFRWT_efmi_S4M1 = R_NilValue, v_iniFRWT_efmi_S4M2 = R_NilValue, v_iniFRWT_efmi_S4M3 = R_NilValue, v_iniFRWT_efmi_S4M4 = R_NilValue,
        v_iniFRWT_efmi2_S1M1 = R_NilValue, v_iniFRWT_efmi2_S1M2 = R_NilValue, v_iniFRWT_efmi2_S1M3 = R_NilValue, v_iniFRWT_efmi2_S1M4 = R_NilValue,
        v_iniFRWT_efmi2_S2M1 = R_NilValue, v_iniFRWT_efmi2_S2M2 = R_NilValue, v_iniFRWT_efmi2_S2M3 = R_NilValue, v_iniFRWT_efmi2_S2M4 = R_NilValue,
        v_iniFRWT_efmi2_S3M1 = R_NilValue, v_iniFRWT_efmi2_S3M2 = R_NilValue, v_iniFRWT_efmi2_S3M3 = R_NilValue, v_iniFRWT_efmi2_S3M4 = R_NilValue,
        v_iniFRWT_efmi2_S4M1 = R_NilValue, v_iniFRWT_efmi2_S4M2 = R_NilValue, v_iniFRWT_efmi2_S4M3 = R_NilValue, v_iniFRWT_efmi2_S4M4 = R_NilValue,
//        iniFRWToth_i_S1M1 = R_NilValue, iniFRWToth_i_S1M2 = R_NilValue, iniFRWToth_i_S1M3 = R_NilValue, iniFRWToth_i_S1M4 = R_NilValue,
//        iniFRWToth_i_S2M1 = R_NilValue, iniFRWToth_i_S2M2 = R_NilValue, iniFRWToth_i_S2M3 = R_NilValue, iniFRWToth_i_S2M4 = R_NilValue,
//        iniFRWToth_i_S3M1 = R_NilValue, iniFRWToth_i_S3M2 = R_NilValue, iniFRWToth_i_S3M3 = R_NilValue, iniFRWToth_i_S3M4 = R_NilValue,
//        iniFRWToth_i_S4M1 = R_NilValue, iniFRWToth_i_S4M2 = R_NilValue, iniFRWToth_i_S4M3 = R_NilValue, iniFRWToth_i_S4M4 = R_NilValue,
//        iniFRWTroth_i_S1M1 = R_NilValue, iniFRWTroth_i_S1M2 = R_NilValue, iniFRWTroth_i_S1M3 = R_NilValue, iniFRWTroth_i_S1M4 = R_NilValue,
//        iniFRWTroth_i_S2M1 = R_NilValue, iniFRWTroth_i_S2M2 = R_NilValue, iniFRWTroth_i_S2M3 = R_NilValue, iniFRWTroth_i_S2M4 = R_NilValue,
//        iniFRWTroth_i_S3M1 = R_NilValue, iniFRWTroth_i_S3M2 = R_NilValue, iniFRWTroth_i_S3M3 = R_NilValue, iniFRWTroth_i_S3M4 = R_NilValue,
//        iniFRWTroth_i_S4M1 = R_NilValue, iniFRWTroth_i_S4M2 = R_NilValue, iniFRWTroth_i_S4M3 = R_NilValue, iniFRWTroth_i_S4M4 = R_NilValue,

        v_iniFDWT_efmi_S1M1 = R_NilValue, v_iniFDWT_efmi_S1M2 = R_NilValue, v_iniFDWT_efmi_S1M3 = R_NilValue, v_iniFDWT_efmi_S1M4 = R_NilValue,
        v_iniFDWT_efmi_S2M1 = R_NilValue, v_iniFDWT_efmi_S2M2 = R_NilValue, v_iniFDWT_efmi_S2M3 = R_NilValue, v_iniFDWT_efmi_S2M4 = R_NilValue,
        v_iniFDWT_efmi_S3M1 = R_NilValue, v_iniFDWT_efmi_S3M2 = R_NilValue, v_iniFDWT_efmi_S3M3 = R_NilValue, v_iniFDWT_efmi_S3M4 = R_NilValue,
        v_iniFDWT_efmi_S4M1 = R_NilValue, v_iniFDWT_efmi_S4M2 = R_NilValue, v_iniFDWT_efmi_S4M3 = R_NilValue, v_iniFDWT_efmi_S4M4 = R_NilValue,
        v_iniFDWT_efmi2_S1M1 = R_NilValue, v_iniFDWT_efmi2_S1M2 = R_NilValue, v_iniFDWT_efmi2_S1M3 = R_NilValue, v_iniFDWT_efmi2_S1M4 = R_NilValue,
        v_iniFDWT_efmi2_S2M1 = R_NilValue, v_iniFDWT_efmi2_S2M2 = R_NilValue, v_iniFDWT_efmi2_S2M3 = R_NilValue, v_iniFDWT_efmi2_S2M4 = R_NilValue,
        v_iniFDWT_efmi2_S3M1 = R_NilValue, v_iniFDWT_efmi2_S3M2 = R_NilValue, v_iniFDWT_efmi2_S3M3 = R_NilValue, v_iniFDWT_efmi2_S3M4 = R_NilValue,
        v_iniFDWT_efmi2_S4M1 = R_NilValue, v_iniFDWT_efmi2_S4M2 = R_NilValue, v_iniFDWT_efmi2_S4M3 = R_NilValue, v_iniFDWT_efmi2_S4M4 = R_NilValue,
//        iniFDWToth_i_S1M1 = R_NilValue, iniFDWToth_i_S1M2 = R_NilValue, iniFDWToth_i_S1M3 = R_NilValue, iniFDWToth_i_S1M4 = R_NilValue,
//        iniFDWToth_i_S2M1 = R_NilValue, iniFDWToth_i_S2M2 = R_NilValue, iniFDWToth_i_S2M3 = R_NilValue, iniFDWToth_i_S2M4 = R_NilValue,
//        iniFDWToth_i_S3M1 = R_NilValue, iniFDWToth_i_S3M2 = R_NilValue, iniFDWToth_i_S3M3 = R_NilValue, iniFDWToth_i_S3M4 = R_NilValue,
//        iniFDWToth_i_S4M1 = R_NilValue, iniFDWToth_i_S4M2 = R_NilValue, iniFDWToth_i_S4M3 = R_NilValue, iniFDWToth_i_S4M4 = R_NilValue,
//        iniFDWTroth_i_S1M1 = R_NilValue, iniFDWTroth_i_S1M2 = R_NilValue, iniFDWTroth_i_S1M3 = R_NilValue, iniFDWTroth_i_S1M4 = R_NilValue,
//        iniFDWTroth_i_S2M1 = R_NilValue, iniFDWTroth_i_S2M2 = R_NilValue, iniFDWTroth_i_S2M3 = R_NilValue, iniFDWTroth_i_S2M4 = R_NilValue,
//        iniFDWTroth_i_S3M1 = R_NilValue, iniFDWTroth_i_S3M2 = R_NilValue, iniFDWTroth_i_S3M3 = R_NilValue, iniFDWTroth_i_S3M4 = R_NilValue,
//        iniFDWTroth_i_S4M1 = R_NilValue, iniFDWTroth_i_S4M2 = R_NilValue, iniFDWTroth_i_S4M3 = R_NilValue, iniFDWTroth_i_S4M4 = R_NilValue,

        v_nbNav_f, v_nbds_f, dim_nbNavCst, dim_nbdsCst, dim_Finput = R_NilValue,// v_Finput, v_fm, v_ventilMoy_f, dim_fmCst, dim_ventilMoyCst,
        fFACT1, fFACT2, fFACT3, fFACT4, fFACT5, fFACT6,
        Foth_i = R_NilValue, Foth_i_G1 = R_NilValue, Foth_i_G2 = R_NilValue,
        Froth_i = R_NilValue, Froth_i_G1 = R_NilValue, Froth_i_G2 = R_NilValue,
        dimI, dimIT, DimIT, fFACTsup1, fFACTsup2; //v_ventil2,

SEXP ans_11 = R_NilValue, ans_11l = R_NilValue, dimnames= R_NilValue, dimnamesIT= R_NilValue, rnames= R_NilValue,
     ans_11_G1 = R_NilValue,ans_11_G2 = R_NilValue, ans_11l_G1 = R_NilValue, ans_11l_G2 = R_NilValue;

SEXP    ans_11_S1M1 = R_NilValue, ans_11_S1M2 = R_NilValue, ans_11_S1M3 = R_NilValue, ans_11_S1M4 = R_NilValue,
        ans_11_S2M1 = R_NilValue, ans_11_S2M2 = R_NilValue, ans_11_S2M3 = R_NilValue, ans_11_S2M4 = R_NilValue,
        ans_11_S3M1 = R_NilValue, ans_11_S3M2 = R_NilValue, ans_11_S3M3 = R_NilValue, ans_11_S3M4 = R_NilValue,
        ans_11_S4M1 = R_NilValue, ans_11_S4M2 = R_NilValue, ans_11_S4M3 = R_NilValue, ans_11_S4M4 = R_NilValue,

        ans_FRWT_S1M1 = R_NilValue, ans_FRWT_S1M2 = R_NilValue, ans_FRWT_S1M3 = R_NilValue, ans_FRWT_S1M4 = R_NilValue,
        ans_FRWT_S2M1 = R_NilValue, ans_FRWT_S2M2 = R_NilValue, ans_FRWT_S2M3 = R_NilValue, ans_FRWT_S2M4 = R_NilValue,
        ans_FRWT_S3M1 = R_NilValue, ans_FRWT_S3M2 = R_NilValue, ans_FRWT_S3M3 = R_NilValue, ans_FRWT_S3M4 = R_NilValue,
        ans_FRWT_S4M1 = R_NilValue, ans_FRWT_S4M2 = R_NilValue, ans_FRWT_S4M3 = R_NilValue, ans_FRWT_S4M4 = R_NilValue,

        ans_FDWT_S1M1 = R_NilValue, ans_FDWT_S1M2 = R_NilValue, ans_FDWT_S1M3 = R_NilValue, ans_FDWT_S1M4 = R_NilValue,
        ans_FDWT_S2M1 = R_NilValue, ans_FDWT_S2M2 = R_NilValue, ans_FDWT_S2M3 = R_NilValue, ans_FDWT_S2M4 = R_NilValue,
        ans_FDWT_S3M1 = R_NilValue, ans_FDWT_S3M2 = R_NilValue, ans_FDWT_S3M3 = R_NilValue, ans_FDWT_S3M4 = R_NilValue,
        ans_FDWT_S4M1 = R_NilValue, ans_FDWT_S4M2 = R_NilValue, ans_FDWT_S4M3 = R_NilValue, ans_FDWT_S4M4 = R_NilValue,

        ans_11l_S1M1 = R_NilValue, ans_11l_S1M2 = R_NilValue, ans_11l_S1M3 = R_NilValue, ans_11l_S1M4 = R_NilValue,
        ans_11l_S2M1 = R_NilValue, ans_11l_S2M2 = R_NilValue, ans_11l_S2M3 = R_NilValue, ans_11l_S2M4 = R_NilValue,
        ans_11l_S3M1 = R_NilValue, ans_11l_S3M2 = R_NilValue, ans_11l_S3M3 = R_NilValue, ans_11l_S3M4 = R_NilValue,
        ans_11l_S4M1 = R_NilValue, ans_11l_S4M2 = R_NilValue, ans_11l_S4M3 = R_NilValue, ans_11l_S4M4 = R_NilValue;

SEXP effort;

//on int�gre la donn�e d'effort (qu'on l'utilise ensuite pour le calcul de la capturabilit�, ou pas)

PROTECT(effort = getListElement(Flist, "effort_f_m_tot"));

PROTECT(dimEff = getAttrib(effort, install("DimCst")));

int *dim_Sr_e, *dim_d_efi, *dim_doth_ei, *dim_F_efmi, *dimC, *dimE, *dimM=0, *dimEffort, *dimF, //*dim_Capt_emi, *dim_Capt_ei
    *dimNav, *dimNbds; //, *dim_fm, *dim_vMoy,*rdim;
int nbI;

double *rans_11=&NA_REAL, *rans_11l=&NA_REAL, *r_Sr_e=&NA_REAL, *r_d_efi=&NA_REAL, *r_doth_ei=&NA_REAL, *r_F_efmi=&NA_REAL, *rEff=&NA_REAL, //*r_nbNav_f,  //*r_fm, *r_ventilMoy_f,
       *rans_11_G1=&NA_REAL, *rans_11_G2=&NA_REAL, *rans_11l_G1=&NA_REAL, *rans_11l_G2=&NA_REAL, *r_d_efi_G1=&NA_REAL, *r_d_efi_G2=&NA_REAL,
       *r_doth_ei_G1=&NA_REAL, *r_doth_ei_G2=&NA_REAL, *r_F_efmi_G1=&NA_REAL,*r_F_efmi_G2=&NA_REAL,

//        *r_F_efmi_S1M1=&NA_REAL, *r_F_efmi_S1M2=&NA_REAL, *r_F_efmi_S1M3=&NA_REAL, *r_F_efmi_S1M4=&NA_REAL,
//        *r_F_efmi_S2M1=&NA_REAL, *r_F_efmi_S2M2=&NA_REAL, *r_F_efmi_S2M3=&NA_REAL, *r_F_efmi_S2M4=&NA_REAL,
//        *r_F_efmi_S3M1=&NA_REAL, *r_F_efmi_S3M2=&NA_REAL, *r_F_efmi_S3M3=&NA_REAL, *r_F_efmi_S3M4=&NA_REAL,
//        *r_F_efmi_S4M1=&NA_REAL, *r_F_efmi_S4M2=&NA_REAL, *r_F_efmi_S4M3=&NA_REAL, *r_F_efmi_S4M4=&NA_REAL,

        *r_iniF_efmi_S1M1=&NA_REAL, *r_iniF_efmi_S1M2=&NA_REAL, *r_iniF_efmi_S1M3=&NA_REAL, *r_iniF_efmi_S1M4=&NA_REAL,
        *r_iniF_efmi_S2M1=&NA_REAL, *r_iniF_efmi_S2M2=&NA_REAL, *r_iniF_efmi_S2M3=&NA_REAL, *r_iniF_efmi_S2M4=&NA_REAL,
        *r_iniF_efmi_S3M1=&NA_REAL, *r_iniF_efmi_S3M2=&NA_REAL, *r_iniF_efmi_S3M3=&NA_REAL, *r_iniF_efmi_S3M4=&NA_REAL,
        *r_iniF_efmi_S4M1=&NA_REAL, *r_iniF_efmi_S4M2=&NA_REAL, *r_iniF_efmi_S4M3=&NA_REAL, *r_iniF_efmi_S4M4=&NA_REAL,

        *r_Foth_i_S1M1=&NA_REAL, *r_Foth_i_S1M2=&NA_REAL, *r_Foth_i_S1M3=&NA_REAL, *r_Foth_i_S1M4=&NA_REAL,
        *r_Foth_i_S2M1=&NA_REAL, *r_Foth_i_S2M2=&NA_REAL, *r_Foth_i_S2M3=&NA_REAL, *r_Foth_i_S2M4=&NA_REAL,
        *r_Foth_i_S3M1=&NA_REAL, *r_Foth_i_S3M2=&NA_REAL, *r_Foth_i_S3M3=&NA_REAL, *r_Foth_i_S3M4=&NA_REAL,
        *r_Foth_i_S4M1=&NA_REAL, *r_Foth_i_S4M2=&NA_REAL, *r_Foth_i_S4M3=&NA_REAL, *r_Foth_i_S4M4=&NA_REAL,
        *r_Froth_i_S1M1=&NA_REAL, *r_Froth_i_S1M2=&NA_REAL, *r_Froth_i_S1M3=&NA_REAL, *r_Froth_i_S1M4=&NA_REAL,
        *r_Froth_i_S2M1=&NA_REAL, *r_Froth_i_S2M2=&NA_REAL, *r_Froth_i_S2M3=&NA_REAL, *r_Froth_i_S2M4=&NA_REAL,
        *r_Froth_i_S3M1=&NA_REAL, *r_Froth_i_S3M2=&NA_REAL, *r_Froth_i_S3M3=&NA_REAL, *r_Froth_i_S3M4=&NA_REAL,
        *r_Froth_i_S4M1=&NA_REAL, *r_Froth_i_S4M2=&NA_REAL, *r_Froth_i_S4M3=&NA_REAL, *r_Froth_i_S4M4=&NA_REAL,
        *rans_11_S1M1=&NA_REAL, *rans_11_S1M2=&NA_REAL, *rans_11_S1M3=&NA_REAL, *rans_11_S1M4=&NA_REAL,
        *rans_11_S2M1=&NA_REAL, *rans_11_S2M2=&NA_REAL, *rans_11_S2M3=&NA_REAL, *rans_11_S2M4=&NA_REAL,
        *rans_11_S3M1=&NA_REAL, *rans_11_S3M2=&NA_REAL, *rans_11_S3M3=&NA_REAL, *rans_11_S3M4=&NA_REAL,
        *rans_11_S4M1=&NA_REAL, *rans_11_S4M2=&NA_REAL, *rans_11_S4M3=&NA_REAL, *rans_11_S4M4=&NA_REAL,
        *rans_11l_S1M1=&NA_REAL, *rans_11l_S1M2=&NA_REAL, *rans_11l_S1M3=&NA_REAL, *rans_11l_S1M4=&NA_REAL,
        *rans_11l_S2M1=&NA_REAL, *rans_11l_S2M2=&NA_REAL, *rans_11l_S2M3=&NA_REAL, *rans_11l_S2M4=&NA_REAL,
        *rans_11l_S3M1=&NA_REAL, *rans_11l_S3M2=&NA_REAL, *rans_11l_S3M3=&NA_REAL, *rans_11l_S3M4=&NA_REAL,
        *rans_11l_S4M1=&NA_REAL, *rans_11l_S4M2=&NA_REAL, *rans_11l_S4M3=&NA_REAL, *rans_11l_S4M4=&NA_REAL,

//        *r_FRWT_efmi_S1M1=&NA_REAL, *r_FRWT_efmi_S1M2=&NA_REAL, *r_FRWT_efmi_S1M3=&NA_REAL, *r_FRWT_efmi_S1M4=&NA_REAL,
//        *r_FRWT_efmi_S2M1=&NA_REAL, *r_FRWT_efmi_S2M2=&NA_REAL, *r_FRWT_efmi_S2M3=&NA_REAL, *r_FRWT_efmi_S2M4=&NA_REAL,
//        *r_FRWT_efmi_S3M1=&NA_REAL, *r_FRWT_efmi_S3M2=&NA_REAL, *r_FRWT_efmi_S3M3=&NA_REAL, *r_FRWT_efmi_S3M4=&NA_REAL,
//        *r_FRWT_efmi_S4M1=&NA_REAL, *r_FRWT_efmi_S4M2=&NA_REAL, *r_FRWT_efmi_S4M3=&NA_REAL, *r_FRWT_efmi_S4M4=&NA_REAL,

        *r_iniFRWT_efmi_S1M1=&NA_REAL, *r_iniFRWT_efmi_S1M2=&NA_REAL, *r_iniFRWT_efmi_S1M3=&NA_REAL, *r_iniFRWT_efmi_S1M4=&NA_REAL,
        *r_iniFRWT_efmi_S2M1=&NA_REAL, *r_iniFRWT_efmi_S2M2=&NA_REAL, *r_iniFRWT_efmi_S2M3=&NA_REAL, *r_iniFRWT_efmi_S2M4=&NA_REAL,
        *r_iniFRWT_efmi_S3M1=&NA_REAL, *r_iniFRWT_efmi_S3M2=&NA_REAL, *r_iniFRWT_efmi_S3M3=&NA_REAL, *r_iniFRWT_efmi_S3M4=&NA_REAL,
        *r_iniFRWT_efmi_S4M1=&NA_REAL, *r_iniFRWT_efmi_S4M2=&NA_REAL, *r_iniFRWT_efmi_S4M3=&NA_REAL, *r_iniFRWT_efmi_S4M4=&NA_REAL,

        *r_FRWToth_i_S1M1=&NA_REAL, *r_FRWToth_i_S1M2=&NA_REAL, *r_FRWToth_i_S1M3=&NA_REAL, *r_FRWToth_i_S1M4=&NA_REAL,
        *r_FRWToth_i_S2M1=&NA_REAL, *r_FRWToth_i_S2M2=&NA_REAL, *r_FRWToth_i_S2M3=&NA_REAL, *r_FRWToth_i_S2M4=&NA_REAL,
        *r_FRWToth_i_S3M1=&NA_REAL, *r_FRWToth_i_S3M2=&NA_REAL, *r_FRWToth_i_S3M3=&NA_REAL, *r_FRWToth_i_S3M4=&NA_REAL,
        *r_FRWToth_i_S4M1=&NA_REAL, *r_FRWToth_i_S4M2=&NA_REAL, *r_FRWToth_i_S4M3=&NA_REAL, *r_FRWToth_i_S4M4=&NA_REAL,
        *rans_FRWT_S1M1=&NA_REAL, *rans_FRWT_S1M2=&NA_REAL, *rans_FRWT_S1M3=&NA_REAL, *rans_FRWT_S1M4=&NA_REAL,
        *rans_FRWT_S2M1=&NA_REAL, *rans_FRWT_S2M2=&NA_REAL, *rans_FRWT_S2M3=&NA_REAL, *rans_FRWT_S2M4=&NA_REAL,
        *rans_FRWT_S3M1=&NA_REAL, *rans_FRWT_S3M2=&NA_REAL, *rans_FRWT_S3M3=&NA_REAL, *rans_FRWT_S3M4=&NA_REAL,
        *rans_FRWT_S4M1=&NA_REAL, *rans_FRWT_S4M2=&NA_REAL, *rans_FRWT_S4M3=&NA_REAL, *rans_FRWT_S4M4=&NA_REAL,

//        *r_FDWT_efmi_S1M1=&NA_REAL, *r_FDWT_efmi_S1M2=&NA_REAL, *r_FDWT_efmi_S1M3=&NA_REAL, *r_FDWT_efmi_S1M4=&NA_REAL,
//        *r_FDWT_efmi_S2M1=&NA_REAL, *r_FDWT_efmi_S2M2=&NA_REAL, *r_FDWT_efmi_S2M3=&NA_REAL, *r_FDWT_efmi_S2M4=&NA_REAL,
//        *r_FDWT_efmi_S3M1=&NA_REAL, *r_FDWT_efmi_S3M2=&NA_REAL, *r_FDWT_efmi_S3M3=&NA_REAL, *r_FDWT_efmi_S3M4=&NA_REAL,
//        *r_FDWT_efmi_S4M1=&NA_REAL, *r_FDWT_efmi_S4M2=&NA_REAL, *r_FDWT_efmi_S4M3=&NA_REAL, *r_FDWT_efmi_S4M4=&NA_REAL,

        *r_iniFDWT_efmi_S1M1=&NA_REAL, *r_iniFDWT_efmi_S1M2=&NA_REAL, *r_iniFDWT_efmi_S1M3=&NA_REAL, *r_iniFDWT_efmi_S1M4=&NA_REAL,
        *r_iniFDWT_efmi_S2M1=&NA_REAL, *r_iniFDWT_efmi_S2M2=&NA_REAL, *r_iniFDWT_efmi_S2M3=&NA_REAL, *r_iniFDWT_efmi_S2M4=&NA_REAL,
        *r_iniFDWT_efmi_S3M1=&NA_REAL, *r_iniFDWT_efmi_S3M2=&NA_REAL, *r_iniFDWT_efmi_S3M3=&NA_REAL, *r_iniFDWT_efmi_S3M4=&NA_REAL,
        *r_iniFDWT_efmi_S4M1=&NA_REAL, *r_iniFDWT_efmi_S4M2=&NA_REAL, *r_iniFDWT_efmi_S4M3=&NA_REAL, *r_iniFDWT_efmi_S4M4=&NA_REAL,

        *r_FDWToth_i_S1M1=&NA_REAL, *r_FDWToth_i_S1M2=&NA_REAL, *r_FDWToth_i_S1M3=&NA_REAL, *r_FDWToth_i_S1M4=&NA_REAL,
        *r_FDWToth_i_S2M1=&NA_REAL, *r_FDWToth_i_S2M2=&NA_REAL, *r_FDWToth_i_S2M3=&NA_REAL, *r_FDWToth_i_S2M4=&NA_REAL,
        *r_FDWToth_i_S3M1=&NA_REAL, *r_FDWToth_i_S3M2=&NA_REAL, *r_FDWToth_i_S3M3=&NA_REAL, *r_FDWToth_i_S3M4=&NA_REAL,
        *r_FDWToth_i_S4M1=&NA_REAL, *r_FDWToth_i_S4M2=&NA_REAL, *r_FDWToth_i_S4M3=&NA_REAL, *r_FDWToth_i_S4M4=&NA_REAL,
        *rans_FDWT_S1M1=&NA_REAL, *rans_FDWT_S1M2=&NA_REAL, *rans_FDWT_S1M3=&NA_REAL, *rans_FDWT_S1M4=&NA_REAL,
        *rans_FDWT_S2M1=&NA_REAL, *rans_FDWT_S2M2=&NA_REAL, *rans_FDWT_S2M3=&NA_REAL, *rans_FDWT_S2M4=&NA_REAL,
        *rans_FDWT_S3M1=&NA_REAL, *rans_FDWT_S3M2=&NA_REAL, *rans_FDWT_S3M3=&NA_REAL, *rans_FDWT_S3M4=&NA_REAL,
        *rans_FDWT_S4M1=&NA_REAL, *rans_FDWT_S4M2=&NA_REAL, *rans_FDWT_S4M3=&NA_REAL, *rans_FDWT_S4M4=&NA_REAL,

        *r_nbds_f, *r_Foth_i=&NA_REAL, *r_Froth_i=&NA_REAL, *r_Foth_i_G1=&NA_REAL, *r_Foth_i_G2=&NA_REAL, *r_Froth_i_G1=&NA_REAL, *r_Froth_i_G2=&NA_REAL;

//pr�paration de l'output
if (ind_t==0) { //Rprintf("Mort1\n");

    PROTECT(rnames = allocVector(STRSXP, nbE));
    setAttrib(out_F_fmi, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_G1, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_G2, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_G1, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_G2, R_NamesSymbol, rnames);

    setAttrib(out_F_fmi_S1M1, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S1M2, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S1M3, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S1M4, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S2M1, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S2M2, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S2M3, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S2M4, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S3M1, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S3M2, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S3M3, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S3M4, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S4M1, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S4M2, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S4M3, R_NamesSymbol, rnames);
    setAttrib(out_F_fmi_S4M4, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S1M1, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S1M2, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S1M3, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S1M4, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S2M1, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S2M2, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S2M3, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S2M4, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S3M1, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S3M2, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S3M3, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S3M4, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S4M1, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S4M2, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S4M3, R_NamesSymbol, rnames);
    setAttrib(out_Fr_fmi_S4M4, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S1M1, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S1M2, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S1M3, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S1M4, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S2M1, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S2M2, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S2M3, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S2M4, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S3M1, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S3M2, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S3M3, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S3M4, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S4M1, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S4M2, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S4M3, R_NamesSymbol, rnames);
    setAttrib(out_FRWT_fmi_S4M4, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S1M1, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S1M2, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S1M3, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S1M4, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S2M1, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S2M2, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S2M3, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S2M4, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S3M1, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S3M2, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S3M3, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S3M4, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S4M1, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S4M2, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S4M3, R_NamesSymbol, rnames);
    setAttrib(out_FDWT_fmi_S4M4, R_NamesSymbol, rnames);

}

for (int e = 0 ; e < nbE ; e++) {
//Rprintf("Mort2\n");
    //---------
    // calcul de Fr_efmit
    //---------

                        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                        PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                        nbI = length(getListElement(elmt, "modI"));


                    //---------------------------------------------------------------------
                    // 1�re �tape : on ventile la mortalit� par les captures si possible
                    //---------------------------------------------------------------------

                    if ((Qvec[e]==1) & (Svec[e]==0)) {
//Rprintf("Mort3.1\n");
                    PROTECT(v_Sr_e = getListElement(elmt, "sr"));
                    PROTECT(v_d_efi = getListElement(elmt, "d_i"));
                    PROTECT(v_doth_ei = getListElement(elmt, "doth_i"));
                    PROTECT(dimCst_Sr_e = getAttrib(v_Sr_e, install("DimCst")));
                    PROTECT(dimCst_d_efi = getAttrib(v_d_efi, install("DimCst")));
                    PROTECT(dimCst_doth_ei = getAttrib(v_doth_ei, install("DimCst")));

                        r_Sr_e = REAL(v_Sr_e);
                        r_d_efi = REAL(v_d_efi);
                        r_doth_ei = REAL(v_doth_ei);


                     PROTECT(v_F_efmi_S1M1 = getListElement(elmt, "Ffmi_S1M1"));

                     //double *ttt = REAL(v_F_efmi_S1M1);

//for (int ind_i = 0 ; ind_i < 1 ; ind_i++)
//for (int ind_m = 0 ; ind_m < 1 ; ind_m++)
//for (int ind_f = 0 ; ind_f < 25 ; ind_f++) {
//
//Rprintf("ind_i %i ind_m %i ind_f %i v_F_efmi_S1M1 intro %f\n",ind_i,ind_m,ind_f,ttt[ind_f + ind_m*nbF + ind_i*nbF*nbMe]*10000000000);
//
//}
                     PROTECT(v_F_efmi_S1M2 = getListElement(elmt, "Ffmi_S1M2"));
                     PROTECT(v_F_efmi_S1M3 = getListElement(elmt, "Ffmi_S1M3"));
                     PROTECT(v_F_efmi_S1M4 = getListElement(elmt, "Ffmi_S1M4"));
                     PROTECT(v_F_efmi_S2M1 = getListElement(elmt, "Ffmi_S2M1"));
                     PROTECT(v_F_efmi_S2M2 = getListElement(elmt, "Ffmi_S2M2"));
                     PROTECT(v_F_efmi_S2M3 = getListElement(elmt, "Ffmi_S2M3"));
                     PROTECT(v_F_efmi_S2M4 = getListElement(elmt, "Ffmi_S2M4"));
                     PROTECT(v_F_efmi_S3M1 = getListElement(elmt, "Ffmi_S3M1"));
                     PROTECT(v_F_efmi_S3M2 = getListElement(elmt, "Ffmi_S3M2"));
                     PROTECT(v_F_efmi_S3M3 = getListElement(elmt, "Ffmi_S3M3"));
                     PROTECT(v_F_efmi_S3M4 = getListElement(elmt, "Ffmi_S3M4"));
                     PROTECT(v_F_efmi_S4M1 = getListElement(elmt, "Ffmi_S4M1"));
                     PROTECT(v_F_efmi_S4M2 = getListElement(elmt, "Ffmi_S4M2"));
                     PROTECT(v_F_efmi_S4M3 = getListElement(elmt, "Ffmi_S4M3"));
                     PROTECT(v_F_efmi_S4M4 = getListElement(elmt, "Ffmi_S4M4"));    //PROTECT + 15

                     PROTECT(v_FRWT_efmi_S1M1 = getListElement(elmt, "FLWfmi_S1M1"));
                     PROTECT(v_FRWT_efmi_S1M2 = getListElement(elmt, "FLWfmi_S1M2"));
                     PROTECT(v_FRWT_efmi_S1M3 = getListElement(elmt, "FLWfmi_S1M3"));
                     PROTECT(v_FRWT_efmi_S1M4 = getListElement(elmt, "FLWfmi_S1M4"));
                     PROTECT(v_FRWT_efmi_S2M1 = getListElement(elmt, "FLWfmi_S2M1"));
                     PROTECT(v_FRWT_efmi_S2M2 = getListElement(elmt, "FLWfmi_S2M2"));
                     PROTECT(v_FRWT_efmi_S2M3 = getListElement(elmt, "FLWfmi_S2M3"));
                     PROTECT(v_FRWT_efmi_S2M4 = getListElement(elmt, "FLWfmi_S2M4"));
                     PROTECT(v_FRWT_efmi_S3M1 = getListElement(elmt, "FLWfmi_S3M1"));
                     PROTECT(v_FRWT_efmi_S3M2 = getListElement(elmt, "FLWfmi_S3M2"));
                     PROTECT(v_FRWT_efmi_S3M3 = getListElement(elmt, "FLWfmi_S3M3"));
                     PROTECT(v_FRWT_efmi_S3M4 = getListElement(elmt, "FLWfmi_S3M4"));
                     PROTECT(v_FRWT_efmi_S4M1 = getListElement(elmt, "FLWfmi_S4M1"));
                     PROTECT(v_FRWT_efmi_S4M2 = getListElement(elmt, "FLWfmi_S4M2"));
                     PROTECT(v_FRWT_efmi_S4M3 = getListElement(elmt, "FLWfmi_S4M3"));
                     PROTECT(v_FRWT_efmi_S4M4 = getListElement(elmt, "FLWfmi_S4M4"));

                     PROTECT(v_FDWT_efmi_S1M1 = getListElement(elmt, "FDWfmi_S1M1"));
                     PROTECT(v_FDWT_efmi_S1M2 = getListElement(elmt, "FDWfmi_S1M2"));
                     PROTECT(v_FDWT_efmi_S1M3 = getListElement(elmt, "FDWfmi_S1M3"));
                     PROTECT(v_FDWT_efmi_S1M4 = getListElement(elmt, "FDWfmi_S1M4"));
                     PROTECT(v_FDWT_efmi_S2M1 = getListElement(elmt, "FDWfmi_S2M1"));
                     PROTECT(v_FDWT_efmi_S2M2 = getListElement(elmt, "FDWfmi_S2M2"));
                     PROTECT(v_FDWT_efmi_S2M3 = getListElement(elmt, "FDWfmi_S2M3"));
                     PROTECT(v_FDWT_efmi_S2M4 = getListElement(elmt, "FDWfmi_S2M4"));
                     PROTECT(v_FDWT_efmi_S3M1 = getListElement(elmt, "FDWfmi_S3M1"));
                     PROTECT(v_FDWT_efmi_S3M2 = getListElement(elmt, "FDWfmi_S3M2"));
                     PROTECT(v_FDWT_efmi_S3M3 = getListElement(elmt, "FDWfmi_S3M3"));
                     PROTECT(v_FDWT_efmi_S3M4 = getListElement(elmt, "FDWfmi_S3M4"));
                     PROTECT(v_FDWT_efmi_S4M1 = getListElement(elmt, "FDWfmi_S4M1"));
                     PROTECT(v_FDWT_efmi_S4M2 = getListElement(elmt, "FDWfmi_S4M2"));
                     PROTECT(v_FDWT_efmi_S4M3 = getListElement(elmt, "FDWfmi_S4M3"));
                     PROTECT(v_FDWT_efmi_S4M4 = getListElement(elmt, "FDWfmi_S4M4"));

                     PROTECT(v_iniF_efmi_S1M1 = getListElement(elmt, "iniFfmi_S1M1"));
                     PROTECT(v_iniF_efmi_S1M2 = getListElement(elmt, "iniFfmi_S1M2"));
                     PROTECT(v_iniF_efmi_S1M3 = getListElement(elmt, "iniFfmi_S1M3"));
                     PROTECT(v_iniF_efmi_S1M4 = getListElement(elmt, "iniFfmi_S1M4"));
                     PROTECT(v_iniF_efmi_S2M1 = getListElement(elmt, "iniFfmi_S2M1"));
                     PROTECT(v_iniF_efmi_S2M2 = getListElement(elmt, "iniFfmi_S2M2"));
                     PROTECT(v_iniF_efmi_S2M3 = getListElement(elmt, "iniFfmi_S2M3"));
                     PROTECT(v_iniF_efmi_S2M4 = getListElement(elmt, "iniFfmi_S2M4"));
                     PROTECT(v_iniF_efmi_S3M1 = getListElement(elmt, "iniFfmi_S3M1"));
                     PROTECT(v_iniF_efmi_S3M2 = getListElement(elmt, "iniFfmi_S3M2"));
                     PROTECT(v_iniF_efmi_S3M3 = getListElement(elmt, "iniFfmi_S3M3"));
                     PROTECT(v_iniF_efmi_S3M4 = getListElement(elmt, "iniFfmi_S3M4"));
                     PROTECT(v_iniF_efmi_S4M1 = getListElement(elmt, "iniFfmi_S4M1"));
                     PROTECT(v_iniF_efmi_S4M2 = getListElement(elmt, "iniFfmi_S4M2"));
                     PROTECT(v_iniF_efmi_S4M3 = getListElement(elmt, "iniFfmi_S4M3"));
                     PROTECT(v_iniF_efmi_S4M4 = getListElement(elmt, "iniFfmi_S4M4"));    //PROTECT + 15

                     PROTECT(v_iniFRWT_efmi_S1M1 = getListElement(elmt, "iniFLWfmi_S1M1"));
                     PROTECT(v_iniFRWT_efmi_S1M2 = getListElement(elmt, "iniFLWfmi_S1M2"));
                     PROTECT(v_iniFRWT_efmi_S1M3 = getListElement(elmt, "iniFLWfmi_S1M3"));
                     PROTECT(v_iniFRWT_efmi_S1M4 = getListElement(elmt, "iniFLWfmi_S1M4"));
                     PROTECT(v_iniFRWT_efmi_S2M1 = getListElement(elmt, "iniFLWfmi_S2M1"));
                     PROTECT(v_iniFRWT_efmi_S2M2 = getListElement(elmt, "iniFLWfmi_S2M2"));
                     PROTECT(v_iniFRWT_efmi_S2M3 = getListElement(elmt, "iniFLWfmi_S2M3"));
                     PROTECT(v_iniFRWT_efmi_S2M4 = getListElement(elmt, "iniFLWfmi_S2M4"));
                     PROTECT(v_iniFRWT_efmi_S3M1 = getListElement(elmt, "iniFLWfmi_S3M1"));
                     PROTECT(v_iniFRWT_efmi_S3M2 = getListElement(elmt, "iniFLWfmi_S3M2"));
                     PROTECT(v_iniFRWT_efmi_S3M3 = getListElement(elmt, "iniFLWfmi_S3M3"));
                     PROTECT(v_iniFRWT_efmi_S3M4 = getListElement(elmt, "iniFLWfmi_S3M4"));
                     PROTECT(v_iniFRWT_efmi_S4M1 = getListElement(elmt, "iniFLWfmi_S4M1"));
                     PROTECT(v_iniFRWT_efmi_S4M2 = getListElement(elmt, "iniFLWfmi_S4M2"));
                     PROTECT(v_iniFRWT_efmi_S4M3 = getListElement(elmt, "iniFLWfmi_S4M3"));
                     PROTECT(v_iniFRWT_efmi_S4M4 = getListElement(elmt, "iniFLWfmi_S4M4"));

                     PROTECT(v_iniFDWT_efmi_S1M1 = getListElement(elmt, "iniFDWfmi_S1M1"));
                     PROTECT(v_iniFDWT_efmi_S1M2 = getListElement(elmt, "iniFDWfmi_S1M2"));
                     PROTECT(v_iniFDWT_efmi_S1M3 = getListElement(elmt, "iniFDWfmi_S1M3"));
                     PROTECT(v_iniFDWT_efmi_S1M4 = getListElement(elmt, "iniFDWfmi_S1M4"));
                     PROTECT(v_iniFDWT_efmi_S2M1 = getListElement(elmt, "iniFDWfmi_S2M1"));
                     PROTECT(v_iniFDWT_efmi_S2M2 = getListElement(elmt, "iniFDWfmi_S2M2"));
                     PROTECT(v_iniFDWT_efmi_S2M3 = getListElement(elmt, "iniFDWfmi_S2M3"));
                     PROTECT(v_iniFDWT_efmi_S2M4 = getListElement(elmt, "iniFDWfmi_S2M4"));
                     PROTECT(v_iniFDWT_efmi_S3M1 = getListElement(elmt, "iniFDWfmi_S3M1"));
                     PROTECT(v_iniFDWT_efmi_S3M2 = getListElement(elmt, "iniFDWfmi_S3M2"));
                     PROTECT(v_iniFDWT_efmi_S3M3 = getListElement(elmt, "iniFDWfmi_S3M3"));
                     PROTECT(v_iniFDWT_efmi_S3M4 = getListElement(elmt, "iniFDWfmi_S3M4"));
                     PROTECT(v_iniFDWT_efmi_S4M1 = getListElement(elmt, "iniFDWfmi_S4M1"));
                     PROTECT(v_iniFDWT_efmi_S4M2 = getListElement(elmt, "iniFDWfmi_S4M2"));
                     PROTECT(v_iniFDWT_efmi_S4M3 = getListElement(elmt, "iniFDWfmi_S4M3"));
                     PROTECT(v_iniFDWT_efmi_S4M4 = getListElement(elmt, "iniFDWfmi_S4M4"));   //### +48

                     PROTECT(dim_Finput = getAttrib(v_F_efmi_S1M1, install("DimCst")));

                    } else if ((Qvec[e]==0) & (Svec[e]==0)){
//Rprintf("Mort3.2\n");
                        PROTECT(v_Sr_e = getListElement(elmt, "sr"));
                        PROTECT(v_d_efi = getListElement(elmt, "d_i"));
                        PROTECT(v_doth_ei = getListElement(elmt, "doth_i"));

                        r_Sr_e = REAL(v_Sr_e);
                        r_d_efi = REAL(v_d_efi);
                        r_doth_ei = REAL(v_doth_ei);

                        PROTECT(dimCst_Sr_e = getAttrib(v_Sr_e, install("DimCst")));
                        PROTECT(dimCst_d_efi = getAttrib(v_d_efi, install("DimCst")));
                        PROTECT(dimCst_doth_ei = getAttrib(v_doth_ei, install("DimCst")));

                     PROTECT(v_F_efmi = getListElement(elmt, "F_fmi"));
                     PROTECT(dim_Finput = getAttrib(v_F_efmi, install("DimCst")));

                    } else if ((Qvec[e]==0) & (Svec[e]==1)){
//Rprintf("Mort3.3\n");
                        PROTECT(v_Sr_e = getListElement(elmt, "sr"));

                        PROTECT(v_d_efi_G1 = getListElement(elmt, "d_i_G1"));
                        PROTECT(v_d_efi_G2 = getListElement(elmt, "d_i_G2"));
                        PROTECT(v_doth_ei_G1 = getListElement(elmt, "doth_i_G1"));
                        PROTECT(v_doth_ei_G2 = getListElement(elmt, "doth_i_G2"));


                        r_Sr_e = REAL(v_Sr_e);
                        r_d_efi_G1 = REAL(v_d_efi_G1);
                        r_doth_ei_G1 = REAL(v_doth_ei_G1);
                        r_d_efi_G2  = REAL(v_d_efi_G2);
                        r_doth_ei_G2 = REAL(v_doth_ei_G2);

                        PROTECT(dimCst_Sr_e = getAttrib(v_Sr_e, install("DimCst")));

                        PROTECT(dimCst_d_efi = getAttrib(v_d_efi_G1, install("DimCst")));
                        PROTECT(dimCst_doth_ei = getAttrib(v_doth_ei_G1, install("DimCst")));

                     PROTECT(v_F_efmi_G1 = getListElement(elmt, "F_fmi_G1"));
                     PROTECT(v_F_efmi_G2 = getListElement(elmt, "F_fmi_G2"));
                     PROTECT(dim_Finput = getAttrib(v_F_efmi_G1, install("DimCst")));

                    }

                    dimF = INTEGER(dim_Finput);


                    PROTECT(v_nbNav_f = getListElement(Flist, "nbv_f_m"));
                    PROTECT(dim_nbNavCst = getAttrib(v_nbNav_f, install("DimCst")));
                    dimNav = INTEGER(dim_nbNavCst);
                    //r_nbNav_f = REAL(v_nbNav_f);

                    PROTECT(v_nbds_f = getListElement(Flist, "effort1_f_m"));
                    PROTECT(dim_nbdsCst = getAttrib(v_nbds_f, install("DimCst")));
                    dimNbds = INTEGER(dim_nbdsCst);
                    r_nbds_f = REAL(v_nbds_f);
                    //r_nbds2_f = REAL(getListElement(Flist, "effort2_f_m"));

                    if (e==0) {

                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            Etemp[ind_f+1] = r_nbds_f[ind_f];

                    }

//Rprintf("Mort4\n");
//Rprintf("Qvec[e] %i\n",Qvec[e]);

                    //on calcule la mortalit� via la capturabilit�
                    if ((Qvec[e]==1) & (Svec[e]==0)) {

                    PROTECT(v_F_efmi2_S1M1 = calcCapturabilite(v_F_efmi_S1M1 , effort));



//fichier2 << "STAAA1S1M1" << endl;
//
//for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
////Rprintf("ind_i %i ind_m %i v_F_efmi_S1M1 %f\n",ind_i,ind_m,REAL(v_F_efmi_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe]*10000000);
////Rprintf("ind_i %i ind_m %i effort %f\n",ind_i,ind_m,REAL(effort)[ind_f + ind_m*nbF]);
////Rprintf("ind_i %i ind_m %i v_F_efmi2_S1M1 %f\n",ind_i,ind_m,REAL(v_F_efmi2_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe]*10000000);
//
//std::stringstream ff4S1M1;
//ff4S1M1 << REAL(v_F_efmi_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe];
//fichier2 << ff4S1M1.str() << endl;
//
//}
//
//fichier2 << "STAAA2S1M1" << endl;
//
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
//std::stringstream ff6S1M1;
//ff6S1M1 << REAL(effort)[ind_f + ind_m*nbF];
//fichier2 << ff6S1M1.str() << endl;
//
//}
//
//
//fichier2 << "STAAA3S1M1" << endl;
//
//for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
//std::stringstream ff5S1M1;
//ff5S1M1 << REAL(v_F_efmi2_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe];
//fichier2 << ff5S1M1.str() << endl;
//
//}
//





                     PROTECT(v_F_efmi2_S1M2 = calcCapturabilite(v_F_efmi_S1M2 , effort));
                     PROTECT(v_F_efmi2_S1M3 = calcCapturabilite(v_F_efmi_S1M3 , effort));
                     PROTECT(v_F_efmi2_S1M4 = calcCapturabilite(v_F_efmi_S1M4 , effort));
                     PROTECT(v_F_efmi2_S2M1 = calcCapturabilite(v_F_efmi_S2M1 , effort));
                     PROTECT(v_F_efmi2_S2M2 = calcCapturabilite(v_F_efmi_S2M2 , effort));
                     PROTECT(v_F_efmi2_S2M3 = calcCapturabilite(v_F_efmi_S2M3 , effort));
                     PROTECT(v_F_efmi2_S2M4 = calcCapturabilite(v_F_efmi_S2M4 , effort));
                     PROTECT(v_F_efmi2_S3M1 = calcCapturabilite(v_F_efmi_S3M1 , effort));
                     PROTECT(v_F_efmi2_S3M2 = calcCapturabilite(v_F_efmi_S3M2 , effort));
                     PROTECT(v_F_efmi2_S3M3 = calcCapturabilite(v_F_efmi_S3M3 , effort));
                     PROTECT(v_F_efmi2_S3M4 = calcCapturabilite(v_F_efmi_S3M4 , effort));
                     PROTECT(v_F_efmi2_S4M1 = calcCapturabilite(v_F_efmi_S4M1 , effort));
                     PROTECT(v_F_efmi2_S4M2 = calcCapturabilite(v_F_efmi_S4M2 , effort));
                     PROTECT(v_F_efmi2_S4M3 = calcCapturabilite(v_F_efmi_S4M3 , effort));
                     PROTECT(v_F_efmi2_S4M4 = calcCapturabilite(v_F_efmi_S4M4 , effort));    //PROTECT + 15



                     PROTECT(v_FRWT_efmi2_S1M1 = calcCapturabilite(v_FRWT_efmi_S1M1 , effort));




//fichier2 << "STBBB1S1M1" << endl;
//
//for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
////Rprintf("ind_i %i ind_m %i v_F_efmi_S1M1 %f\n",ind_i,ind_m,REAL(v_F_efmi_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe]*10000000);
////Rprintf("ind_i %i ind_m %i effort %f\n",ind_i,ind_m,REAL(effort)[ind_f + ind_m*nbF]);
////Rprintf("ind_i %i ind_m %i v_F_efmi2_S1M1 %f\n",ind_i,ind_m,REAL(v_F_efmi2_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe]*10000000);
//
//std::stringstream tt4S1M1;
//tt4S1M1 << REAL(v_FRWT_efmi_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe];
//fichier2 << tt4S1M1.str() << endl;
//
//}
//
//fichier2 << "STBBB2S1M1" << endl;
//
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
//std::stringstream tt6S1M1;
//tt6S1M1 << REAL(effort)[ind_f + ind_m*nbF];
//fichier2 << tt6S1M1.str() << endl;
//
//}
//
//fichier2 << "STBBB3S1M1" << endl;
//
//for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
//for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//for (int ind_f = 0 ; ind_f < 1 ; ind_f++) {
//
//std::stringstream tt5S1M1;
//tt5S1M1 << REAL(v_FRWT_efmi2_S1M1)[ind_f + ind_m*nbF + ind_i*nbF*nbMe];
//fichier2 << tt5S1M1.str() << endl;
//
//}




                     PROTECT(v_FRWT_efmi2_S1M2 = calcCapturabilite(v_FRWT_efmi_S1M2 , effort));
                     PROTECT(v_FRWT_efmi2_S1M3 = calcCapturabilite(v_FRWT_efmi_S1M3 , effort));
                     PROTECT(v_FRWT_efmi2_S1M4 = calcCapturabilite(v_FRWT_efmi_S1M4 , effort));
                     PROTECT(v_FRWT_efmi2_S2M1 = calcCapturabilite(v_FRWT_efmi_S2M1 , effort));
                     PROTECT(v_FRWT_efmi2_S2M2 = calcCapturabilite(v_FRWT_efmi_S2M2 , effort));
                     PROTECT(v_FRWT_efmi2_S2M3 = calcCapturabilite(v_FRWT_efmi_S2M3 , effort));
                     PROTECT(v_FRWT_efmi2_S2M4 = calcCapturabilite(v_FRWT_efmi_S2M4 , effort));
                     PROTECT(v_FRWT_efmi2_S3M1 = calcCapturabilite(v_FRWT_efmi_S3M1 , effort));
                     PROTECT(v_FRWT_efmi2_S3M2 = calcCapturabilite(v_FRWT_efmi_S3M2 , effort));
                     PROTECT(v_FRWT_efmi2_S3M3 = calcCapturabilite(v_FRWT_efmi_S3M3 , effort));
                     PROTECT(v_FRWT_efmi2_S3M4 = calcCapturabilite(v_FRWT_efmi_S3M4 , effort));
                     PROTECT(v_FRWT_efmi2_S4M1 = calcCapturabilite(v_FRWT_efmi_S4M1 , effort));
                     PROTECT(v_FRWT_efmi2_S4M2 = calcCapturabilite(v_FRWT_efmi_S4M2 , effort));
                     PROTECT(v_FRWT_efmi2_S4M3 = calcCapturabilite(v_FRWT_efmi_S4M3 , effort));
                     PROTECT(v_FRWT_efmi2_S4M4 = calcCapturabilite(v_FRWT_efmi_S4M4 , effort));

                     PROTECT(v_FDWT_efmi2_S1M1 = calcCapturabilite(v_FDWT_efmi_S1M1 , effort));
                     PROTECT(v_FDWT_efmi2_S1M2 = calcCapturabilite(v_FDWT_efmi_S1M2 , effort));
                     PROTECT(v_FDWT_efmi2_S1M3 = calcCapturabilite(v_FDWT_efmi_S1M3 , effort));
                     PROTECT(v_FDWT_efmi2_S1M4 = calcCapturabilite(v_FDWT_efmi_S1M4 , effort));
                     PROTECT(v_FDWT_efmi2_S2M1 = calcCapturabilite(v_FDWT_efmi_S2M1 , effort));
                     PROTECT(v_FDWT_efmi2_S2M2 = calcCapturabilite(v_FDWT_efmi_S2M2 , effort));
                     PROTECT(v_FDWT_efmi2_S2M3 = calcCapturabilite(v_FDWT_efmi_S2M3 , effort));
                     PROTECT(v_FDWT_efmi2_S2M4 = calcCapturabilite(v_FDWT_efmi_S2M4 , effort));
                     PROTECT(v_FDWT_efmi2_S3M1 = calcCapturabilite(v_FDWT_efmi_S3M1 , effort));
                     PROTECT(v_FDWT_efmi2_S3M2 = calcCapturabilite(v_FDWT_efmi_S3M2 , effort));
                     PROTECT(v_FDWT_efmi2_S3M3 = calcCapturabilite(v_FDWT_efmi_S3M3 , effort));
                     PROTECT(v_FDWT_efmi2_S3M4 = calcCapturabilite(v_FDWT_efmi_S3M4 , effort));
                     PROTECT(v_FDWT_efmi2_S4M1 = calcCapturabilite(v_FDWT_efmi_S4M1 , effort));
                     PROTECT(v_FDWT_efmi2_S4M2 = calcCapturabilite(v_FDWT_efmi_S4M2 , effort));
                     PROTECT(v_FDWT_efmi2_S4M3 = calcCapturabilite(v_FDWT_efmi_S4M3 , effort));
                     PROTECT(v_FDWT_efmi2_S4M4 = calcCapturabilite(v_FDWT_efmi_S4M4 , effort));

                     PROTECT(v_iniF_efmi2_S1M1 = calcCapturabilite(v_iniF_efmi_S1M1 , effort));
                     PROTECT(v_iniF_efmi2_S1M2 = calcCapturabilite(v_iniF_efmi_S1M2 , effort));
                     PROTECT(v_iniF_efmi2_S1M3 = calcCapturabilite(v_iniF_efmi_S1M3 , effort));
                     PROTECT(v_iniF_efmi2_S1M4 = calcCapturabilite(v_iniF_efmi_S1M4 , effort));
                     PROTECT(v_iniF_efmi2_S2M1 = calcCapturabilite(v_iniF_efmi_S2M1 , effort));
                     PROTECT(v_iniF_efmi2_S2M2 = calcCapturabilite(v_iniF_efmi_S2M2 , effort));
                     PROTECT(v_iniF_efmi2_S2M3 = calcCapturabilite(v_iniF_efmi_S2M3 , effort));
                     PROTECT(v_iniF_efmi2_S2M4 = calcCapturabilite(v_iniF_efmi_S2M4 , effort));
                     PROTECT(v_iniF_efmi2_S3M1 = calcCapturabilite(v_iniF_efmi_S3M1 , effort));
                     PROTECT(v_iniF_efmi2_S3M2 = calcCapturabilite(v_iniF_efmi_S3M2 , effort));
                     PROTECT(v_iniF_efmi2_S3M3 = calcCapturabilite(v_iniF_efmi_S3M3 , effort));
                     PROTECT(v_iniF_efmi2_S3M4 = calcCapturabilite(v_iniF_efmi_S3M4 , effort));
                     PROTECT(v_iniF_efmi2_S4M1 = calcCapturabilite(v_iniF_efmi_S4M1 , effort));
                     PROTECT(v_iniF_efmi2_S4M2 = calcCapturabilite(v_iniF_efmi_S4M2 , effort));
                     PROTECT(v_iniF_efmi2_S4M3 = calcCapturabilite(v_iniF_efmi_S4M3 , effort));
                     PROTECT(v_iniF_efmi2_S4M4 = calcCapturabilite(v_iniF_efmi_S4M4 , effort));

                     PROTECT(v_iniFRWT_efmi2_S1M1 = calcCapturabilite(v_iniFRWT_efmi_S1M1 , effort));
                     PROTECT(v_iniFRWT_efmi2_S1M2 = calcCapturabilite(v_iniFRWT_efmi_S1M2 , effort));
                     PROTECT(v_iniFRWT_efmi2_S1M3 = calcCapturabilite(v_iniFRWT_efmi_S1M3 , effort));
                     PROTECT(v_iniFRWT_efmi2_S1M4 = calcCapturabilite(v_iniFRWT_efmi_S1M4 , effort));
                     PROTECT(v_iniFRWT_efmi2_S2M1 = calcCapturabilite(v_iniFRWT_efmi_S2M1 , effort));
                     PROTECT(v_iniFRWT_efmi2_S2M2 = calcCapturabilite(v_iniFRWT_efmi_S2M2 , effort));
                     PROTECT(v_iniFRWT_efmi2_S2M3 = calcCapturabilite(v_iniFRWT_efmi_S2M3 , effort));
                     PROTECT(v_iniFRWT_efmi2_S2M4 = calcCapturabilite(v_iniFRWT_efmi_S2M4 , effort));
                     PROTECT(v_iniFRWT_efmi2_S3M1 = calcCapturabilite(v_iniFRWT_efmi_S3M1 , effort));
                     PROTECT(v_iniFRWT_efmi2_S3M2 = calcCapturabilite(v_iniFRWT_efmi_S3M2 , effort));
                     PROTECT(v_iniFRWT_efmi2_S3M3 = calcCapturabilite(v_iniFRWT_efmi_S3M3 , effort));
                     PROTECT(v_iniFRWT_efmi2_S3M4 = calcCapturabilite(v_iniFRWT_efmi_S3M4 , effort));
                     PROTECT(v_iniFRWT_efmi2_S4M1 = calcCapturabilite(v_iniFRWT_efmi_S4M1 , effort));
                     PROTECT(v_iniFRWT_efmi2_S4M2 = calcCapturabilite(v_iniFRWT_efmi_S4M2 , effort));
                     PROTECT(v_iniFRWT_efmi2_S4M3 = calcCapturabilite(v_iniFRWT_efmi_S4M3 , effort));
                     PROTECT(v_iniFRWT_efmi2_S4M4 = calcCapturabilite(v_iniFRWT_efmi_S4M4 , effort));

                     PROTECT(v_iniFDWT_efmi2_S1M1 = calcCapturabilite(v_iniFDWT_efmi_S1M1 , effort));
                     PROTECT(v_iniFDWT_efmi2_S1M2 = calcCapturabilite(v_iniFDWT_efmi_S1M2 , effort));
                     PROTECT(v_iniFDWT_efmi2_S1M3 = calcCapturabilite(v_iniFDWT_efmi_S1M3 , effort));
                     PROTECT(v_iniFDWT_efmi2_S1M4 = calcCapturabilite(v_iniFDWT_efmi_S1M4 , effort));
                     PROTECT(v_iniFDWT_efmi2_S2M1 = calcCapturabilite(v_iniFDWT_efmi_S2M1 , effort));
                     PROTECT(v_iniFDWT_efmi2_S2M2 = calcCapturabilite(v_iniFDWT_efmi_S2M2 , effort));
                     PROTECT(v_iniFDWT_efmi2_S2M3 = calcCapturabilite(v_iniFDWT_efmi_S2M3 , effort));
                     PROTECT(v_iniFDWT_efmi2_S2M4 = calcCapturabilite(v_iniFDWT_efmi_S2M4 , effort));
                     PROTECT(v_iniFDWT_efmi2_S3M1 = calcCapturabilite(v_iniFDWT_efmi_S3M1 , effort));
                     PROTECT(v_iniFDWT_efmi2_S3M2 = calcCapturabilite(v_iniFDWT_efmi_S3M2 , effort));
                     PROTECT(v_iniFDWT_efmi2_S3M3 = calcCapturabilite(v_iniFDWT_efmi_S3M3 , effort));
                     PROTECT(v_iniFDWT_efmi2_S3M4 = calcCapturabilite(v_iniFDWT_efmi_S3M4 , effort));
                     PROTECT(v_iniFDWT_efmi2_S4M1 = calcCapturabilite(v_iniFDWT_efmi_S4M1 , effort));
                     PROTECT(v_iniFDWT_efmi2_S4M2 = calcCapturabilite(v_iniFDWT_efmi_S4M2 , effort));
                     PROTECT(v_iniFDWT_efmi2_S4M3 = calcCapturabilite(v_iniFDWT_efmi_S4M3 , effort));
                     PROTECT(v_iniFDWT_efmi2_S4M4 = calcCapturabilite(v_iniFDWT_efmi_S4M4 , effort));    //### +48

                    } else if((Qvec[e]==0) & (Svec[e]==0)){

                     PROTECT(v_F_efmi2 = calcCapturabilite(v_F_efmi , effort)); //PrintValue(v_F_efmi2);

                    } else if((Qvec[e]==0) & (Svec[e]==1)){

                     PROTECT(v_F_efmi2_G1 = calcCapturabilite(v_F_efmi_G1 , effort));
                     PROTECT(v_F_efmi2_G2 = calcCapturabilite(v_F_efmi_G2 , effort));

                    }

                    ////PrintValue(v_F_efmi2);
                            //et dans ce cas, l'effort � appliquer � la capturabilit� est...

                        dimE = INTEGER(dimEff);

                        if ((Qvec[e]==1) & (Svec[e]==0)) {
                            dimM = INTEGER(getAttrib(v_F_efmi2_S1M1, install("DimCst")));
                        } else if ((Qvec[e]==0) & (Svec[e]==0)){
                            dimM = INTEGER(getAttrib(v_F_efmi2, install("DimCst")));
                        } else if ((Qvec[e]==0) & (Svec[e]==1)){
                            dimM = INTEGER(getAttrib(v_F_efmi2_G1, install("DimCst")));
                        }

                        PROTECT(dimCstEff = allocVector(INTSXP,4));
                        dimEffort = INTEGER(dimCstEff);
                        for (int i = 0 ; i < 3 ; i++) dimEffort[i] = imin2( dimM[i] , dimE[i] );

                        //on conserve tout de m�me la dimension temporelle
                        dimEffort[3] = dimE[3];


                            PROTECT(formatEff = aggregObj(effort, dimCstEff));////PrintValue(formatEff);
                            rEff = REAL(formatEff);
//Rprintf("Mort5\n");
                        if ((Qvec[e]==1) & (Svec[e]==0)) {

//                         r_F_efmi_S1M1 = REAL(v_F_efmi2_S1M1);
//                         r_F_efmi_S1M2 = REAL(v_F_efmi2_S1M2);
//                         r_F_efmi_S1M3 = REAL(v_F_efmi2_S1M3);
//                         r_F_efmi_S1M4 = REAL(v_F_efmi2_S1M4);
//                         r_F_efmi_S2M1 = REAL(v_F_efmi2_S2M1);
//                         r_F_efmi_S2M2 = REAL(v_F_efmi2_S2M2);
//                         r_F_efmi_S2M3 = REAL(v_F_efmi2_S2M3);
//                         r_F_efmi_S2M4 = REAL(v_F_efmi2_S2M4);
//                         r_F_efmi_S3M1 = REAL(v_F_efmi2_S3M1);
//                         r_F_efmi_S3M2 = REAL(v_F_efmi2_S3M2);
//                         r_F_efmi_S3M3 = REAL(v_F_efmi2_S3M3);
//                         r_F_efmi_S3M4 = REAL(v_F_efmi2_S3M4);
//                         r_F_efmi_S4M1 = REAL(v_F_efmi2_S4M1);
//                         r_F_efmi_S4M2 = REAL(v_F_efmi2_S4M2);
//                         r_F_efmi_S4M3 = REAL(v_F_efmi2_S4M3);
//                         r_F_efmi_S4M4 = REAL(v_F_efmi2_S4M4);

//                         r_FRWT_efmi_S1M1 = REAL(v_FRWT_efmi2_S1M1);
//                         r_FRWT_efmi_S1M2 = REAL(v_FRWT_efmi2_S1M2);
//                         r_FRWT_efmi_S1M3 = REAL(v_FRWT_efmi2_S1M3);
//                         r_FRWT_efmi_S1M4 = REAL(v_FRWT_efmi2_S1M4);
//                         r_FRWT_efmi_S2M1 = REAL(v_FRWT_efmi2_S2M1);
//                         r_FRWT_efmi_S2M2 = REAL(v_FRWT_efmi2_S2M2);
//                         r_FRWT_efmi_S2M3 = REAL(v_FRWT_efmi2_S2M3);
//                         r_FRWT_efmi_S2M4 = REAL(v_FRWT_efmi2_S2M4);
//                         r_FRWT_efmi_S3M1 = REAL(v_FRWT_efmi2_S3M1);
//                         r_FRWT_efmi_S3M2 = REAL(v_FRWT_efmi2_S3M2);
//                         r_FRWT_efmi_S3M3 = REAL(v_FRWT_efmi2_S3M3);
//                         r_FRWT_efmi_S3M4 = REAL(v_FRWT_efmi2_S3M4);
//                         r_FRWT_efmi_S4M1 = REAL(v_FRWT_efmi2_S4M1);
//                         r_FRWT_efmi_S4M2 = REAL(v_FRWT_efmi2_S4M2);
//                         r_FRWT_efmi_S4M3 = REAL(v_FRWT_efmi2_S4M3);
//                         r_FRWT_efmi_S4M4 = REAL(v_FRWT_efmi2_S4M4);
//
//                         r_FDWT_efmi_S1M1 = REAL(v_FDWT_efmi2_S1M1);
//                         r_FDWT_efmi_S1M2 = REAL(v_FDWT_efmi2_S1M2);
//                         r_FDWT_efmi_S1M3 = REAL(v_FDWT_efmi2_S1M3);
//                         r_FDWT_efmi_S1M4 = REAL(v_FDWT_efmi2_S1M4);
//                         r_FDWT_efmi_S2M1 = REAL(v_FDWT_efmi2_S2M1);
//                         r_FDWT_efmi_S2M2 = REAL(v_FDWT_efmi2_S2M2);
//                         r_FDWT_efmi_S2M3 = REAL(v_FDWT_efmi2_S2M3);
//                         r_FDWT_efmi_S2M4 = REAL(v_FDWT_efmi2_S2M4);
//                         r_FDWT_efmi_S3M1 = REAL(v_FDWT_efmi2_S3M1);
//                         r_FDWT_efmi_S3M2 = REAL(v_FDWT_efmi2_S3M2);
//                         r_FDWT_efmi_S3M3 = REAL(v_FDWT_efmi2_S3M3);
//                         r_FDWT_efmi_S3M4 = REAL(v_FDWT_efmi2_S3M4);
//                         r_FDWT_efmi_S4M1 = REAL(v_FDWT_efmi2_S4M1);
//                         r_FDWT_efmi_S4M2 = REAL(v_FDWT_efmi2_S4M2);
//                         r_FDWT_efmi_S4M3 = REAL(v_FDWT_efmi2_S4M3);
//                         r_FDWT_efmi_S4M4 = REAL(v_FDWT_efmi2_S4M4);

                         r_iniF_efmi_S1M1 = REAL(v_iniF_efmi2_S1M1);
                         r_iniF_efmi_S1M2 = REAL(v_iniF_efmi2_S1M2);
                         r_iniF_efmi_S1M3 = REAL(v_iniF_efmi2_S1M3);
                         r_iniF_efmi_S1M4 = REAL(v_iniF_efmi2_S1M4);
                         r_iniF_efmi_S2M1 = REAL(v_iniF_efmi2_S2M1);
                         r_iniF_efmi_S2M2 = REAL(v_iniF_efmi2_S2M2);
                         r_iniF_efmi_S2M3 = REAL(v_iniF_efmi2_S2M3);
                         r_iniF_efmi_S2M4 = REAL(v_iniF_efmi2_S2M4);
                         r_iniF_efmi_S3M1 = REAL(v_iniF_efmi2_S3M1);
                         r_iniF_efmi_S3M2 = REAL(v_iniF_efmi2_S3M2);
                         r_iniF_efmi_S3M3 = REAL(v_iniF_efmi2_S3M3);
                         r_iniF_efmi_S3M4 = REAL(v_iniF_efmi2_S3M4);
                         r_iniF_efmi_S4M1 = REAL(v_iniF_efmi2_S4M1);
                         r_iniF_efmi_S4M2 = REAL(v_iniF_efmi2_S4M2);
                         r_iniF_efmi_S4M3 = REAL(v_iniF_efmi2_S4M3);
                         r_iniF_efmi_S4M4 = REAL(v_iniF_efmi2_S4M4);

                         r_iniFRWT_efmi_S1M1 = REAL(v_iniFRWT_efmi2_S1M1);
                         r_iniFRWT_efmi_S1M2 = REAL(v_iniFRWT_efmi2_S1M2);
                         r_iniFRWT_efmi_S1M3 = REAL(v_iniFRWT_efmi2_S1M3);
                         r_iniFRWT_efmi_S1M4 = REAL(v_iniFRWT_efmi2_S1M4);
                         r_iniFRWT_efmi_S2M1 = REAL(v_iniFRWT_efmi2_S2M1);
                         r_iniFRWT_efmi_S2M2 = REAL(v_iniFRWT_efmi2_S2M2);
                         r_iniFRWT_efmi_S2M3 = REAL(v_iniFRWT_efmi2_S2M3);
                         r_iniFRWT_efmi_S2M4 = REAL(v_iniFRWT_efmi2_S2M4);
                         r_iniFRWT_efmi_S3M1 = REAL(v_iniFRWT_efmi2_S3M1);
                         r_iniFRWT_efmi_S3M2 = REAL(v_iniFRWT_efmi2_S3M2);
                         r_iniFRWT_efmi_S3M3 = REAL(v_iniFRWT_efmi2_S3M3);
                         r_iniFRWT_efmi_S3M4 = REAL(v_iniFRWT_efmi2_S3M4);
                         r_iniFRWT_efmi_S4M1 = REAL(v_iniFRWT_efmi2_S4M1);
                         r_iniFRWT_efmi_S4M2 = REAL(v_iniFRWT_efmi2_S4M2);
                         r_iniFRWT_efmi_S4M3 = REAL(v_iniFRWT_efmi2_S4M3);
                         r_iniFRWT_efmi_S4M4 = REAL(v_iniFRWT_efmi2_S4M4);

                         r_iniFDWT_efmi_S1M1 = REAL(v_iniFDWT_efmi2_S1M1);
                         r_iniFDWT_efmi_S1M2 = REAL(v_iniFDWT_efmi2_S1M2);
                         r_iniFDWT_efmi_S1M3 = REAL(v_iniFDWT_efmi2_S1M3);
                         r_iniFDWT_efmi_S1M4 = REAL(v_iniFDWT_efmi2_S1M4);
                         r_iniFDWT_efmi_S2M1 = REAL(v_iniFDWT_efmi2_S2M1);
                         r_iniFDWT_efmi_S2M2 = REAL(v_iniFDWT_efmi2_S2M2);
                         r_iniFDWT_efmi_S2M3 = REAL(v_iniFDWT_efmi2_S2M3);
                         r_iniFDWT_efmi_S2M4 = REAL(v_iniFDWT_efmi2_S2M4);
                         r_iniFDWT_efmi_S3M1 = REAL(v_iniFDWT_efmi2_S3M1);
                         r_iniFDWT_efmi_S3M2 = REAL(v_iniFDWT_efmi2_S3M2);
                         r_iniFDWT_efmi_S3M3 = REAL(v_iniFDWT_efmi2_S3M3);
                         r_iniFDWT_efmi_S3M4 = REAL(v_iniFDWT_efmi2_S3M4);
                         r_iniFDWT_efmi_S4M1 = REAL(v_iniFDWT_efmi2_S4M1);
                         r_iniFDWT_efmi_S4M2 = REAL(v_iniFDWT_efmi2_S4M2);
                         r_iniFDWT_efmi_S4M3 = REAL(v_iniFDWT_efmi2_S4M3);
                         r_iniFDWT_efmi_S4M4 = REAL(v_iniFDWT_efmi2_S4M4);

                         PROTECT(dimCst_F_efmi = getAttrib(v_F_efmi2_S1M1, install("DimCst")));

                        } else if ((Qvec[e]==0) & (Svec[e]==0)){

                         r_F_efmi = REAL(v_F_efmi2);
                         PROTECT(dimCst_F_efmi = getAttrib(v_F_efmi2, install("DimCst")));

                        } else if ((Qvec[e]==0) & (Svec[e]==1)){

                         r_F_efmi_G1 = REAL(v_F_efmi2_G1);
                         r_F_efmi_G2 = REAL(v_F_efmi2_G2);
                         PROTECT(dimCst_F_efmi = getAttrib(v_F_efmi2_G1, install("DimCst")));

                        }

//Rprintf("Mort6\n");
                        //tests sur les dimensions
                        dim_Sr_e = INTEGER(dimCst_Sr_e);
                        if (((dim_Sr_e[0]!=0) & (dim_Sr_e[0]!=nbF)) | ((dim_Sr_e[1]!=0) & (dim_Sr_e[1]!=nbM)) |
                            ((dim_Sr_e[2]!=0) & (dim_Sr_e[2]!=nbI)) | ((dim_Sr_e[3]!=0) & (dim_Sr_e[3]!=nbT)))
                        {
                            error("Non_homogeneous dimensions in Sr_e element. Check .ini biological parameters files !!\n");
                        }

                        dim_d_efi = INTEGER(dimCst_d_efi);
                        if (((dim_d_efi[0]!=0) & (dim_d_efi[0]!=nbF)) | ((dim_d_efi[1]!=0) & (dim_d_efi[1]!=nbM)) |
                            ((dim_d_efi[2]!=0) & (dim_d_efi[2]!=nbI)) | ((dim_d_efi[3]!=0) & (dim_d_efi[3]!=nbT)))
                        {
                            error("Non_homogeneous dimensions in d_efi element. Check .ini biological parameters files !!\n");
                        }

                        dim_doth_ei = INTEGER(dimCst_doth_ei);
                        if (((dim_doth_ei[0]!=0) & (dim_doth_ei[0]!=nbF)) | ((dim_doth_ei[1]!=0) & (dim_doth_ei[1]!=nbM)) |
                            ((dim_doth_ei[2]!=0) & (dim_doth_ei[2]!=nbI)) | ((dim_doth_ei[3]!=0) & (dim_doth_ei[3]!=nbT)))
                        {
                            error("Non_homogeneous dimensions in doth_ei element. Check .ini biological parameters files !!\n");
                        }

                        dim_F_efmi = INTEGER(dimCst_F_efmi);
                        if (((dim_F_efmi[0]!=0) & (dim_F_efmi[0]!=nbF)) | ((dim_F_efmi[1]!=0) & (dim_F_efmi[1]!=nbM)) |
                            ((dim_F_efmi[2]!=0) & (dim_F_efmi[2]!=nbI)) | ((dim_F_efmi[3]!=0) & (dim_F_efmi[3]!=nbT)))
                        {
                            error("Non_homogeneous dimensions in F_efmi element. Check .ini biological parameters files !!\n");
                        }

                        //on d�termine l'attribut Dimension du tableau r�sultant -> dimCst (on en profite pour compter les dimensions r�elles + nombre de cellules)
                        PROTECT(dimCst = allocVector(INTSXP, 4));
                        dimC = INTEGER(dimCst);
                        int count = 0, prod = 1, count2 = 0, count3 = 0;

                        for (int k = 0 ; k < 4 ; k++) {

                            dimC[k] = imax2( imax2(dim_d_efi[k] , dim_F_efmi[k]) , dimEffort[k]);
                            if (k==3) dimC[3] = nbT; //on consid�re la donn�e temporellement
                            if (dimC[k]>0) {
                                count++;
                                prod = prod * dimC[k];
                            }

                        }
//Rprintf("Mort7\n");
                        PROTECT(Dim = allocVector(INTSXP, count));
                        int *dim = INTEGER(Dim);

                        for (int k = 0 ; k < 4 ; k++) {
                            if (dimC[k]>0) {
                                dim[count2] = dimC[k];
                                count2++;
                                }
                        }

//Rprintf("Mort8\n");
                    if (ind_t==0) {

                      if ((Qvec[e]==1) & (Svec[e]==0)) {

                        PROTECT(ans_11_S1M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S1M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S1M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S1M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S2M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S2M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S2M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S2M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S3M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S3M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S3M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S3M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S4M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S4M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S4M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_S4M4 = NEW_NUMERIC(prod));

                        PROTECT(ans_11l_S1M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S1M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S1M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S1M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S2M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S2M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S2M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S2M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S3M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S3M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S3M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S3M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S4M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S4M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S4M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_S4M4 = NEW_NUMERIC(prod));

                        PROTECT(ans_FRWT_S1M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S1M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S1M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S1M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S2M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S2M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S2M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S2M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S3M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S3M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S3M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S3M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S4M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S4M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S4M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FRWT_S4M4 = NEW_NUMERIC(prod));

                        PROTECT(ans_FDWT_S1M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S1M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S1M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S1M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S2M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S2M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S2M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S2M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S3M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S3M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S3M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S3M4 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S4M1 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S4M2 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S4M3 = NEW_NUMERIC(prod));
                        PROTECT(ans_FDWT_S4M4 = NEW_NUMERIC(prod));

                        setAttrib(ans_11_S1M1, R_DimSymbol, Dim);
                        setAttrib(ans_11_S1M2, R_DimSymbol, Dim);
                        setAttrib(ans_11_S1M3, R_DimSymbol, Dim);
                        setAttrib(ans_11_S1M4, R_DimSymbol, Dim);
                        setAttrib(ans_11_S2M1, R_DimSymbol, Dim);
                        setAttrib(ans_11_S2M2, R_DimSymbol, Dim);
                        setAttrib(ans_11_S2M3, R_DimSymbol, Dim);
                        setAttrib(ans_11_S2M4, R_DimSymbol, Dim);
                        setAttrib(ans_11_S3M1, R_DimSymbol, Dim);
                        setAttrib(ans_11_S3M2, R_DimSymbol, Dim);
                        setAttrib(ans_11_S3M3, R_DimSymbol, Dim);
                        setAttrib(ans_11_S3M4, R_DimSymbol, Dim);
                        setAttrib(ans_11_S4M1, R_DimSymbol, Dim);
                        setAttrib(ans_11_S4M2, R_DimSymbol, Dim);
                        setAttrib(ans_11_S4M3, R_DimSymbol, Dim);
                        setAttrib(ans_11_S4M4, R_DimSymbol, Dim);

                        setAttrib(ans_11l_S1M1, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S1M2, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S1M3, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S1M4, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S2M1, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S2M2, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S2M3, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S2M4, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S3M1, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S3M2, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S3M3, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S3M4, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S4M1, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S4M2, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S4M3, R_DimSymbol, Dim);
                        setAttrib(ans_11l_S4M4, R_DimSymbol, Dim);

                        setAttrib(ans_FRWT_S1M1, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S1M2, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S1M3, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S1M4, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S2M1, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S2M2, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S2M3, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S2M4, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S3M1, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S3M2, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S3M3, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S3M4, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S4M1, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S4M2, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S4M3, R_DimSymbol, Dim);
                        setAttrib(ans_FRWT_S4M4, R_DimSymbol, Dim);

                        setAttrib(ans_FDWT_S1M1, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S1M2, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S1M3, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S1M4, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S2M1, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S2M2, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S2M3, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S2M4, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S3M1, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S3M2, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S3M3, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S3M4, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S4M1, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S4M2, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S4M3, R_DimSymbol, Dim);
                        setAttrib(ans_FDWT_S4M4, R_DimSymbol, Dim);

                        rans_11_S1M1 = REAL(ans_11_S1M1);
                        rans_11_S1M2 = REAL(ans_11_S1M2);
                        rans_11_S1M3 = REAL(ans_11_S1M3);
                        rans_11_S1M4 = REAL(ans_11_S1M4);
                        rans_11_S2M1 = REAL(ans_11_S2M1);
                        rans_11_S2M2 = REAL(ans_11_S2M2);
                        rans_11_S2M3 = REAL(ans_11_S2M3);
                        rans_11_S2M4 = REAL(ans_11_S2M4);
                        rans_11_S3M1 = REAL(ans_11_S3M1);
                        rans_11_S3M2 = REAL(ans_11_S3M2);
                        rans_11_S3M3 = REAL(ans_11_S3M3);
                        rans_11_S3M4 = REAL(ans_11_S3M4);
                        rans_11_S4M1 = REAL(ans_11_S4M1);
                        rans_11_S4M2 = REAL(ans_11_S4M2);
                        rans_11_S4M3 = REAL(ans_11_S4M3);
                        rans_11_S4M4 = REAL(ans_11_S4M4);

                        rans_11l_S1M1 = REAL(ans_11l_S1M1);
                        rans_11l_S1M2 = REAL(ans_11l_S1M2);
                        rans_11l_S1M3 = REAL(ans_11l_S1M3);
                        rans_11l_S1M4 = REAL(ans_11l_S1M4);
                        rans_11l_S2M1 = REAL(ans_11l_S2M1);
                        rans_11l_S2M2 = REAL(ans_11l_S2M2);
                        rans_11l_S2M3 = REAL(ans_11l_S2M3);
                        rans_11l_S2M4 = REAL(ans_11l_S2M4);
                        rans_11l_S3M1 = REAL(ans_11l_S3M1);
                        rans_11l_S3M2 = REAL(ans_11l_S3M2);
                        rans_11l_S3M3 = REAL(ans_11l_S3M3);
                        rans_11l_S3M4 = REAL(ans_11l_S3M4);
                        rans_11l_S4M1 = REAL(ans_11l_S4M1);
                        rans_11l_S4M2 = REAL(ans_11l_S4M2);
                        rans_11l_S4M3 = REAL(ans_11l_S4M3);
                        rans_11l_S4M4 = REAL(ans_11l_S4M4);

                        rans_FDWT_S1M1 = REAL(ans_FDWT_S1M1);
                        rans_FDWT_S1M2 = REAL(ans_FDWT_S1M2);
                        rans_FDWT_S1M3 = REAL(ans_FDWT_S1M3);
                        rans_FDWT_S1M4 = REAL(ans_FDWT_S1M4);
                        rans_FDWT_S2M1 = REAL(ans_FDWT_S2M1);
                        rans_FDWT_S2M2 = REAL(ans_FDWT_S2M2);
                        rans_FDWT_S2M3 = REAL(ans_FDWT_S2M3);
                        rans_FDWT_S2M4 = REAL(ans_FDWT_S2M4);
                        rans_FDWT_S3M1 = REAL(ans_FDWT_S3M1);
                        rans_FDWT_S3M2 = REAL(ans_FDWT_S3M2);
                        rans_FDWT_S3M3 = REAL(ans_FDWT_S3M3);
                        rans_FDWT_S3M4 = REAL(ans_FDWT_S3M4);
                        rans_FDWT_S4M1 = REAL(ans_FDWT_S4M1);
                        rans_FDWT_S4M2 = REAL(ans_FDWT_S4M2);
                        rans_FDWT_S4M3 = REAL(ans_FDWT_S4M3);
                        rans_FDWT_S4M4 = REAL(ans_FDWT_S4M4);

                        rans_FRWT_S1M1 = REAL(ans_FRWT_S1M1);
                        rans_FRWT_S1M2 = REAL(ans_FRWT_S1M2);
                        rans_FRWT_S1M3 = REAL(ans_FRWT_S1M3);
                        rans_FRWT_S1M4 = REAL(ans_FRWT_S1M4);
                        rans_FRWT_S2M1 = REAL(ans_FRWT_S2M1);
                        rans_FRWT_S2M2 = REAL(ans_FRWT_S2M2);
                        rans_FRWT_S2M3 = REAL(ans_FRWT_S2M3);
                        rans_FRWT_S2M4 = REAL(ans_FRWT_S2M4);
                        rans_FRWT_S3M1 = REAL(ans_FRWT_S3M1);
                        rans_FRWT_S3M2 = REAL(ans_FRWT_S3M2);
                        rans_FRWT_S3M3 = REAL(ans_FRWT_S3M3);
                        rans_FRWT_S3M4 = REAL(ans_FRWT_S3M4);
                        rans_FRWT_S4M1 = REAL(ans_FRWT_S4M1);
                        rans_FRWT_S4M2 = REAL(ans_FRWT_S4M2);
                        rans_FRWT_S4M3 = REAL(ans_FRWT_S4M3);
                        rans_FRWT_S4M4 = REAL(ans_FRWT_S4M4);

                      } else if ((Qvec[e]==0) & (Svec[e]==0)) {
                        PROTECT(ans_11l = NEW_NUMERIC(prod));
                        PROTECT(ans_11 = NEW_NUMERIC(prod));


                        setAttrib(ans_11, R_DimSymbol, Dim);
                        setAttrib(ans_11l, R_DimSymbol, Dim);

                        rans_11 = REAL(ans_11);
                        rans_11l = REAL(ans_11l);

                      } else if ((Qvec[e]==0) & (Svec[e]==1)) {
                        PROTECT(ans_11l_G1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_G1 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l_G2 = NEW_NUMERIC(prod));
                        PROTECT(ans_11_G2 = NEW_NUMERIC(prod));

                        setAttrib(ans_11_G1, R_DimSymbol, Dim);
                        setAttrib(ans_11l_G1, R_DimSymbol, Dim);
                        setAttrib(ans_11_G2, R_DimSymbol, Dim);
                        setAttrib(ans_11l_G2, R_DimSymbol, Dim);

                        rans_11_G1 = REAL(ans_11_G1);
                        rans_11l_G1 = REAL(ans_11l_G1);
                        rans_11_G2 = REAL(ans_11_G2);
                        rans_11l_G2 = REAL(ans_11l_G2);
                      }




//Rprintf("Mort9\n");

                        PROTECT(dimnames = allocVector(VECSXP,count));
                        if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                        if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                        if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                        if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}


                    } else {

                      if ((Qvec[e]==1) & (Svec[e]==0)) {

                        rans_11_S1M1 = REAL(VECTOR_ELT(out_F_fmi_S1M1, e));
                        rans_11_S1M2 = REAL(VECTOR_ELT(out_F_fmi_S1M2, e));
                        rans_11_S1M3 = REAL(VECTOR_ELT(out_F_fmi_S1M3, e));
                        rans_11_S1M4 = REAL(VECTOR_ELT(out_F_fmi_S1M4, e));
                        rans_11_S2M1 = REAL(VECTOR_ELT(out_F_fmi_S2M1, e));
                        rans_11_S2M2 = REAL(VECTOR_ELT(out_F_fmi_S2M2, e));
                        rans_11_S2M3 = REAL(VECTOR_ELT(out_F_fmi_S2M3, e));
                        rans_11_S2M4 = REAL(VECTOR_ELT(out_F_fmi_S2M4, e));
                        rans_11_S3M1 = REAL(VECTOR_ELT(out_F_fmi_S3M1, e));
                        rans_11_S3M2 = REAL(VECTOR_ELT(out_F_fmi_S3M2, e));
                        rans_11_S3M3 = REAL(VECTOR_ELT(out_F_fmi_S3M3, e));
                        rans_11_S3M4 = REAL(VECTOR_ELT(out_F_fmi_S3M4, e));
                        rans_11_S4M1 = REAL(VECTOR_ELT(out_F_fmi_S4M1, e));
                        rans_11_S4M2 = REAL(VECTOR_ELT(out_F_fmi_S4M2, e));
                        rans_11_S4M3 = REAL(VECTOR_ELT(out_F_fmi_S4M3, e));
                        rans_11_S4M4 = REAL(VECTOR_ELT(out_F_fmi_S4M4, e));

                        rans_11l_S1M1 = REAL(VECTOR_ELT(out_Fr_fmi_S1M1, e));
                        rans_11l_S1M2 = REAL(VECTOR_ELT(out_Fr_fmi_S1M2, e));
                        rans_11l_S1M3 = REAL(VECTOR_ELT(out_Fr_fmi_S1M3, e));
                        rans_11l_S1M4 = REAL(VECTOR_ELT(out_Fr_fmi_S1M4, e));
                        rans_11l_S2M1 = REAL(VECTOR_ELT(out_Fr_fmi_S2M1, e));
                        rans_11l_S2M2 = REAL(VECTOR_ELT(out_Fr_fmi_S2M2, e));
                        rans_11l_S2M3 = REAL(VECTOR_ELT(out_Fr_fmi_S2M3, e));
                        rans_11l_S2M4 = REAL(VECTOR_ELT(out_Fr_fmi_S2M4, e));
                        rans_11l_S3M1 = REAL(VECTOR_ELT(out_Fr_fmi_S3M1, e));
                        rans_11l_S3M2 = REAL(VECTOR_ELT(out_Fr_fmi_S3M2, e));
                        rans_11l_S3M3 = REAL(VECTOR_ELT(out_Fr_fmi_S3M3, e));
                        rans_11l_S3M4 = REAL(VECTOR_ELT(out_Fr_fmi_S3M4, e));
                        rans_11l_S4M1 = REAL(VECTOR_ELT(out_Fr_fmi_S4M1, e));
                        rans_11l_S4M2 = REAL(VECTOR_ELT(out_Fr_fmi_S4M2, e));
                        rans_11l_S4M3 = REAL(VECTOR_ELT(out_Fr_fmi_S4M3, e));
                        rans_11l_S4M4 = REAL(VECTOR_ELT(out_Fr_fmi_S4M4, e));

                        rans_FRWT_S1M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M1, e));
                        rans_FRWT_S1M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M2, e));
                        rans_FRWT_S1M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M3, e));
                        rans_FRWT_S1M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M4, e));
                        rans_FRWT_S2M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M1, e));
                        rans_FRWT_S2M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M2, e));
                        rans_FRWT_S2M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M3, e));
                        rans_FRWT_S2M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M4, e));
                        rans_FRWT_S3M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M1, e));
                        rans_FRWT_S3M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M2, e));
                        rans_FRWT_S3M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M3, e));
                        rans_FRWT_S3M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M4, e));
                        rans_FRWT_S4M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M1, e));
                        rans_FRWT_S4M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M2, e));
                        rans_FRWT_S4M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M3, e));
                        rans_FRWT_S4M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M4, e));

                        rans_FDWT_S1M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M1, e));
                        rans_FDWT_S1M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M2, e));
                        rans_FDWT_S1M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M3, e));
                        rans_FDWT_S1M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M4, e));
                        rans_FDWT_S2M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M1, e));
                        rans_FDWT_S2M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M2, e));
                        rans_FDWT_S2M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M3, e));
                        rans_FDWT_S2M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M4, e));
                        rans_FDWT_S3M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M1, e));
                        rans_FDWT_S3M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M2, e));
                        rans_FDWT_S3M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M3, e));
                        rans_FDWT_S3M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M4, e));
                        rans_FDWT_S4M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M1, e));
                        rans_FDWT_S4M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M2, e));
                        rans_FDWT_S4M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M3, e));
                        rans_FDWT_S4M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M4, e));

                      } else if ((Qvec[e]==0) & (Svec[e]==0)) {

                        rans_11 = REAL(VECTOR_ELT(out_F_fmi, e));
                        rans_11l = REAL(VECTOR_ELT(out_Fr_fmi, e));

                      } else if ((Qvec[e]==0) & (Svec[e]==1)) {

                        rans_11_G1 = REAL(VECTOR_ELT(out_F_fmi_G1, e));
                        rans_11l_G1 = REAL(VECTOR_ELT(out_Fr_fmi_G1, e));
                        rans_11_G2 = REAL(VECTOR_ELT(out_F_fmi_G2, e));
                        rans_11l_G2 = REAL(VECTOR_ELT(out_Fr_fmi_G2, e));
                      }

                    }


                            //facteurs des indices pour gen�riciser le processus
//Rprintf("Mort10\n");
                            PROTECT(fFACT1 = iDim(dimC));
                            PROTECT(fFACT2 = iDim(dim_d_efi));
                            PROTECT(fFACT3 = iDim(dim_Sr_e));
                            PROTECT(fFACT4 = iDim(dim_F_efmi));
                            PROTECT(fFACT5 = iDim(dimEffort));
                            PROTECT(fFACTsup1 = iDim(dimNav));
                            PROTECT(fFACTsup2 = iDim(dimNbds));
                            PROTECT(fFACT6 = iDim(dim_doth_ei));

                            int *fFact1 = INTEGER(fFACT1);
                            int *fFact2 = INTEGER(fFACT2);
                            int *fFact3 = INTEGER(fFACT3);
                            int *fFact4 = INTEGER(fFACT4);
                            int *fFact5 = INTEGER(fFACT5);
                            int *fFact6 = INTEGER(fFACT6);

//Rprintf("Mort11\n");
                            //�quation

                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                        if ((Qvec[e]==1) & (Svec[e]==0)) {

                           rans_11_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_11_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_11_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                          // l'�valuation de rans_11 et rans_11l ne peut se faire que dans le module de dynamiques de pop car on a besoin de N

                           rans_11l_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] ;//*
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);
                           if (ind_i==0) rans_11l_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_11l_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                           rans_11l_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniF_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];// *
                            //(1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            //  r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);



                           rans_FRWT_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FRWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FRWT_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FRWT_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFRWT_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];


                           rans_FDWT_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];
                           if (ind_i==0) rans_FDWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                           rans_FDWT_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                           rans_FDWT_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_iniFDWT_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];


                        } else if ((Qvec[e]==0) & (Svec[e]==0)) {


                        rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                        rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] *
                            (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                              r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                        } else if ((Qvec[e]==0) & (Svec[e]==1)) {


                        rans_11_G1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi_G1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                        rans_11l_G1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi_G1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] *
                            (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                              r_d_efi_G1[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                        rans_11_G2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi_G2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                        rans_11l_G2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi_G2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] *
                            (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                              r_d_efi_G2[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                        }}

//Rprintf("Mort12\n");

                    if (ind_t==0) {

                      if ((Qvec[e]==1) & (Svec[e]==0)) {

                        setAttrib(ans_11_S1M1, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S1M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S1M1, install("DimCst"), dimCst); setAttrib(ans_11l_S1M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S1M1, e, ans_11_S1M1); SET_VECTOR_ELT(out_Fr_fmi_S1M1, e, ans_11l_S1M1);

                        setAttrib(ans_11_S1M2, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S1M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S1M2, install("DimCst"), dimCst); setAttrib(ans_11l_S1M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S1M2, e, ans_11_S1M2); SET_VECTOR_ELT(out_Fr_fmi_S1M2, e, ans_11l_S1M2);

                        setAttrib(ans_11_S1M3, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S1M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S1M3, install("DimCst"), dimCst); setAttrib(ans_11l_S1M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S1M3, e, ans_11_S1M3); SET_VECTOR_ELT(out_Fr_fmi_S1M3, e, ans_11l_S1M3);

                        setAttrib(ans_11_S1M4, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S1M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S1M4, install("DimCst"), dimCst); setAttrib(ans_11l_S1M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S1M4, e, ans_11_S1M4); SET_VECTOR_ELT(out_Fr_fmi_S1M4, e, ans_11l_S1M4);

                        setAttrib(ans_11_S2M1, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S2M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S2M1, install("DimCst"), dimCst); setAttrib(ans_11l_S2M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S2M1, e, ans_11_S2M1); SET_VECTOR_ELT(out_Fr_fmi_S2M1, e, ans_11l_S2M1);

                        setAttrib(ans_11_S2M2, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S2M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S2M2, install("DimCst"), dimCst); setAttrib(ans_11l_S2M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S2M2, e, ans_11_S2M2); SET_VECTOR_ELT(out_Fr_fmi_S2M2, e, ans_11l_S2M2);

                        setAttrib(ans_11_S2M3, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S2M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S2M3, install("DimCst"), dimCst); setAttrib(ans_11l_S2M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S2M3, e, ans_11_S2M3); SET_VECTOR_ELT(out_Fr_fmi_S2M3, e, ans_11l_S2M3);

                        setAttrib(ans_11_S2M4, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S2M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S2M4, install("DimCst"), dimCst); setAttrib(ans_11l_S2M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S2M4, e, ans_11_S2M4); SET_VECTOR_ELT(out_Fr_fmi_S2M4, e, ans_11l_S2M4);

                        setAttrib(ans_11_S3M1, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S3M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S3M1, install("DimCst"), dimCst); setAttrib(ans_11l_S3M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S3M1, e, ans_11_S3M1); SET_VECTOR_ELT(out_Fr_fmi_S3M1, e, ans_11l_S3M1);

                        setAttrib(ans_11_S3M2, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S3M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S3M2, install("DimCst"), dimCst); setAttrib(ans_11l_S3M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S3M2, e, ans_11_S3M2); SET_VECTOR_ELT(out_Fr_fmi_S3M2, e, ans_11l_S3M2);

                        setAttrib(ans_11_S3M3, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S3M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S3M3, install("DimCst"), dimCst); setAttrib(ans_11l_S3M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S3M3, e, ans_11_S3M3); SET_VECTOR_ELT(out_Fr_fmi_S3M3, e, ans_11l_S3M3);

                        setAttrib(ans_11_S3M4, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S3M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S3M4, install("DimCst"), dimCst); setAttrib(ans_11l_S3M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S3M4, e, ans_11_S3M4); SET_VECTOR_ELT(out_Fr_fmi_S3M4, e, ans_11l_S3M4);

                        setAttrib(ans_11_S4M1, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S4M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S4M1, install("DimCst"), dimCst); setAttrib(ans_11l_S4M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S4M1, e, ans_11_S4M1); SET_VECTOR_ELT(out_Fr_fmi_S4M1, e, ans_11l_S4M1);

                        setAttrib(ans_11_S4M2, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S4M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S4M2, install("DimCst"), dimCst); setAttrib(ans_11l_S4M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S4M2, e, ans_11_S4M2); SET_VECTOR_ELT(out_Fr_fmi_S4M2, e, ans_11l_S4M2);

                        setAttrib(ans_11_S4M3, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S4M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S4M3, install("DimCst"), dimCst); setAttrib(ans_11l_S4M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S4M3, e, ans_11_S4M3); SET_VECTOR_ELT(out_Fr_fmi_S4M3, e, ans_11l_S4M3);

                        setAttrib(ans_11_S4M4, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_S4M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_S4M4, install("DimCst"), dimCst); setAttrib(ans_11l_S4M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_S4M4, e, ans_11_S4M4); SET_VECTOR_ELT(out_Fr_fmi_S4M4, e, ans_11l_S4M4);



                        setAttrib(ans_FRWT_S1M1, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S1M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S1M1, install("DimCst"), dimCst); setAttrib(ans_FDWT_S1M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S1M1, e, ans_FRWT_S1M1); SET_VECTOR_ELT(out_FDWT_fmi_S1M1, e, ans_FDWT_S1M1);

                        setAttrib(ans_FRWT_S1M2, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S1M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S1M2, install("DimCst"), dimCst); setAttrib(ans_FDWT_S1M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S1M2, e, ans_FRWT_S1M2); SET_VECTOR_ELT(out_FDWT_fmi_S1M2, e, ans_FDWT_S1M2);

                        setAttrib(ans_FRWT_S1M3, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S1M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S1M3, install("DimCst"), dimCst); setAttrib(ans_FDWT_S1M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S1M3, e, ans_FRWT_S1M3); SET_VECTOR_ELT(out_FDWT_fmi_S1M3, e, ans_FDWT_S1M3);

                        setAttrib(ans_FRWT_S1M4, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S1M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S1M4, install("DimCst"), dimCst); setAttrib(ans_FDWT_S1M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S1M4, e, ans_FRWT_S1M4); SET_VECTOR_ELT(out_FDWT_fmi_S1M4, e, ans_FDWT_S1M4);

                        setAttrib(ans_FRWT_S2M1, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S2M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S2M1, install("DimCst"), dimCst); setAttrib(ans_FDWT_S2M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S2M1, e, ans_FRWT_S2M1); SET_VECTOR_ELT(out_FDWT_fmi_S2M1, e, ans_FDWT_S2M1);

                        setAttrib(ans_FRWT_S2M2, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S2M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S2M2, install("DimCst"), dimCst); setAttrib(ans_FDWT_S2M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S2M2, e, ans_FRWT_S2M2); SET_VECTOR_ELT(out_FDWT_fmi_S2M2, e, ans_FDWT_S2M2);

                        setAttrib(ans_FRWT_S2M3, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S2M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S2M3, install("DimCst"), dimCst); setAttrib(ans_FDWT_S2M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S2M3, e, ans_FRWT_S2M3); SET_VECTOR_ELT(out_FDWT_fmi_S2M3, e, ans_FDWT_S2M3);

                        setAttrib(ans_FRWT_S2M4, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S2M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S2M4, install("DimCst"), dimCst); setAttrib(ans_FDWT_S2M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S2M4, e, ans_FRWT_S2M4); SET_VECTOR_ELT(out_FDWT_fmi_S2M4, e, ans_FDWT_S2M4);

                        setAttrib(ans_FRWT_S3M1, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S3M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S3M1, install("DimCst"), dimCst); setAttrib(ans_FDWT_S3M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S3M1, e, ans_FRWT_S3M1); SET_VECTOR_ELT(out_FDWT_fmi_S3M1, e, ans_FDWT_S3M1);

                        setAttrib(ans_FRWT_S3M2, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S3M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S3M2, install("DimCst"), dimCst); setAttrib(ans_FDWT_S3M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S3M2, e, ans_FRWT_S3M2); SET_VECTOR_ELT(out_FDWT_fmi_S3M2, e, ans_FDWT_S3M2);

                        setAttrib(ans_FRWT_S3M3, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S3M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S3M3, install("DimCst"), dimCst); setAttrib(ans_FDWT_S3M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S3M3, e, ans_FRWT_S3M3); SET_VECTOR_ELT(out_FDWT_fmi_S3M3, e, ans_FDWT_S3M3);

                        setAttrib(ans_FRWT_S3M4, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S3M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S3M4, install("DimCst"), dimCst); setAttrib(ans_FDWT_S3M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S3M4, e, ans_FRWT_S3M4); SET_VECTOR_ELT(out_FDWT_fmi_S3M4, e, ans_FDWT_S3M4);

                        setAttrib(ans_FRWT_S4M1, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S4M1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S4M1, install("DimCst"), dimCst); setAttrib(ans_FDWT_S4M1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S4M1, e, ans_FRWT_S4M1); SET_VECTOR_ELT(out_FDWT_fmi_S4M1, e, ans_FDWT_S4M1);

                        setAttrib(ans_FRWT_S4M2, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S4M2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S4M2, install("DimCst"), dimCst); setAttrib(ans_FDWT_S4M2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S4M2, e, ans_FRWT_S4M2); SET_VECTOR_ELT(out_FDWT_fmi_S4M2, e, ans_FDWT_S4M2);

                        setAttrib(ans_FRWT_S4M3, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S4M3, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S4M3, install("DimCst"), dimCst); setAttrib(ans_FDWT_S4M3, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S4M3, e, ans_FRWT_S4M3); SET_VECTOR_ELT(out_FDWT_fmi_S4M3, e, ans_FDWT_S4M3);

                        setAttrib(ans_FRWT_S4M4, R_DimNamesSymbol, dimnames); setAttrib(ans_FDWT_S4M4, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_FRWT_S4M4, install("DimCst"), dimCst); setAttrib(ans_FDWT_S4M4, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_FRWT_fmi_S4M4, e, ans_FRWT_S4M4); SET_VECTOR_ELT(out_FDWT_fmi_S4M4, e, ans_FDWT_S4M4);

                      } else if ((Qvec[e]==0) & (Svec[e]==0)) {

                        setAttrib(ans_11, R_DimNamesSymbol, dimnames); setAttrib(ans_11l, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11, install("DimCst"), dimCst); setAttrib(ans_11l, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi, e, ans_11); SET_VECTOR_ELT(out_Fr_fmi, e, ans_11l);

                      } else if ((Qvec[e]==0) & (Svec[e]==1)) {

                        setAttrib(ans_11_G1, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_G1, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_G1, install("DimCst"), dimCst); setAttrib(ans_11l_G1, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_G1, e, ans_11_G1); SET_VECTOR_ELT(out_Fr_fmi_G1, e, ans_11l_G1);

                        setAttrib(ans_11_G2, R_DimNamesSymbol, dimnames); setAttrib(ans_11l_G2, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11_G2, install("DimCst"), dimCst); setAttrib(ans_11l_G2, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi_G2, e, ans_11_G2); SET_VECTOR_ELT(out_Fr_fmi_G2, e, ans_11l_G2);

                      }
                       SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
                        //UNPROTECT(3);
                        //if ((Qvec[e]==1) & (Svec[e]==0)) UNPROTECT(32);
                        //if ((Qvec[e]==1) & (Svec[e]==0)) UNPROTECT(32);

                    }


                    //il ne reste plus qu'� calculer Foth_i en soutrayant de Ftot_i la somme aux �ges de la mortalit� ventil�e non corrig�e, et Froth_i en lui appliquant doth_i


                        PROTECT(dimI = allocVector(INTSXP,4));
                        PROTECT(dimIT = allocVector(INTSXP,4));
                        PROTECT(DimIT = allocVector(INTSXP,2));
                        int *rdimI = INTEGER(dimI); rdimI[0] = 0; rdimI[1] = 0; rdimI[2] = nbI; rdimI[3] = dimF[3];
                        int *rdimIT = INTEGER(dimIT); rdimIT[0] = 0; rdimIT[1] = 0; rdimIT[2] = nbI; rdimIT[3] = nbT;
                        int *rDimIT = INTEGER(DimIT); rDimIT[0] = nbI; rDimIT[1] = nbT;

                        PROTECT(dimnamesIT = allocVector(VECSXP,2));
                        SET_VECTOR_ELT(dimnamesIT, 0, intAge);
                        SET_VECTOR_ELT(dimnamesIT, 1, times);

//Rprintf("Mort13\n");


                      if ((Qvec[e]==1) & (Svec[e]==0)) {

                        PROTECT(Foth_i_S1M1 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S1M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S1M1, R_DimSymbol, DimIT); setAttrib(Froth_i_S1M1, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S1M1, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S1M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S1M1, install("DimCst"), dimIT); setAttrib(Froth_i_S1M1, install("DimCst"), dimIT);
                        r_Foth_i_S1M1 = REAL(Foth_i_S1M1);
                        r_Froth_i_S1M1 = REAL(Froth_i_S1M1);

                        PROTECT(Foth_i_S1M2 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S1M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S1M2, R_DimSymbol, DimIT); setAttrib(Froth_i_S1M2, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S1M2, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S1M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S1M2, install("DimCst"), dimIT); setAttrib(Froth_i_S1M2, install("DimCst"), dimIT);
                        r_Foth_i_S1M2 = REAL(Foth_i_S1M2);
                        r_Froth_i_S1M2 = REAL(Froth_i_S1M2);

                        PROTECT(Foth_i_S1M3 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S1M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S1M3, R_DimSymbol, DimIT); setAttrib(Froth_i_S1M3, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S1M3, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S1M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S1M3, install("DimCst"), dimIT); setAttrib(Froth_i_S1M3, install("DimCst"), dimIT);
                        r_Foth_i_S1M3 = REAL(Foth_i_S1M3);
                        r_Froth_i_S1M3 = REAL(Froth_i_S1M3);

                        PROTECT(Foth_i_S1M4 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S1M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S1M4, R_DimSymbol, DimIT); setAttrib(Froth_i_S1M4, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S1M4, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S1M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S1M4, install("DimCst"), dimIT); setAttrib(Froth_i_S1M4, install("DimCst"), dimIT);
                        r_Foth_i_S1M4 = REAL(Foth_i_S1M4);
                        r_Froth_i_S1M4 = REAL(Froth_i_S1M4);

                        PROTECT(Foth_i_S2M1 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S2M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S2M1, R_DimSymbol, DimIT); setAttrib(Froth_i_S2M1, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S2M1, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S2M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S2M1, install("DimCst"), dimIT); setAttrib(Froth_i_S2M1, install("DimCst"), dimIT);
                        r_Foth_i_S2M1 = REAL(Foth_i_S2M1);
                        r_Froth_i_S2M1 = REAL(Froth_i_S2M1);

                        PROTECT(Foth_i_S2M2 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S2M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S2M2, R_DimSymbol, DimIT); setAttrib(Froth_i_S2M2, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S2M2, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S2M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S2M2, install("DimCst"), dimIT); setAttrib(Froth_i_S2M2, install("DimCst"), dimIT);
                        r_Foth_i_S2M2 = REAL(Foth_i_S2M2);
                        r_Froth_i_S2M2 = REAL(Froth_i_S2M2);

                        PROTECT(Foth_i_S2M3 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S2M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S2M3, R_DimSymbol, DimIT); setAttrib(Froth_i_S2M3, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S2M3, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S2M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S2M3, install("DimCst"), dimIT); setAttrib(Froth_i_S2M3, install("DimCst"), dimIT);
                        r_Foth_i_S2M3 = REAL(Foth_i_S2M3);
                        r_Froth_i_S2M3 = REAL(Froth_i_S2M3);

                        PROTECT(Foth_i_S2M4 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S2M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S2M4, R_DimSymbol, DimIT); setAttrib(Froth_i_S2M4, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S2M4, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S2M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S2M4, install("DimCst"), dimIT); setAttrib(Froth_i_S2M4, install("DimCst"), dimIT);
                        r_Foth_i_S2M4 = REAL(Foth_i_S2M4);
                        r_Froth_i_S2M4 = REAL(Froth_i_S2M4);

                        PROTECT(Foth_i_S3M1 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S3M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S3M1, R_DimSymbol, DimIT); setAttrib(Froth_i_S3M1, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S3M1, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S3M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S3M1, install("DimCst"), dimIT); setAttrib(Froth_i_S3M1, install("DimCst"), dimIT);
                        r_Foth_i_S3M1 = REAL(Foth_i_S3M1);
                        r_Froth_i_S3M1 = REAL(Froth_i_S3M1);

                        PROTECT(Foth_i_S3M2 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S3M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S3M2, R_DimSymbol, DimIT); setAttrib(Froth_i_S3M2, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S3M2, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S3M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S3M2, install("DimCst"), dimIT); setAttrib(Froth_i_S3M2, install("DimCst"), dimIT);
                        r_Foth_i_S3M2 = REAL(Foth_i_S3M2);
                        r_Froth_i_S3M2 = REAL(Froth_i_S3M2);

                        PROTECT(Foth_i_S3M3 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S3M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S3M3, R_DimSymbol, DimIT); setAttrib(Froth_i_S3M3, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S3M3, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S3M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S3M3, install("DimCst"), dimIT); setAttrib(Froth_i_S3M3, install("DimCst"), dimIT);
                        r_Foth_i_S3M3 = REAL(Foth_i_S3M3);
                        r_Froth_i_S3M3 = REAL(Froth_i_S3M3);

                        PROTECT(Foth_i_S3M4 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S3M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S3M4, R_DimSymbol, DimIT); setAttrib(Froth_i_S3M4, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S3M4, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S3M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S3M4, install("DimCst"), dimIT); setAttrib(Froth_i_S3M4, install("DimCst"), dimIT);
                        r_Foth_i_S3M4 = REAL(Foth_i_S3M4);
                        r_Froth_i_S3M4 = REAL(Froth_i_S3M4);

                        PROTECT(Foth_i_S4M1 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S4M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S4M1, R_DimSymbol, DimIT); setAttrib(Froth_i_S4M1, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S4M1, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S4M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S4M1, install("DimCst"), dimIT); setAttrib(Froth_i_S4M1, install("DimCst"), dimIT);
                        r_Foth_i_S4M1 = REAL(Foth_i_S4M1);
                        r_Froth_i_S4M1 = REAL(Froth_i_S4M1);

                        PROTECT(Foth_i_S4M2 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S4M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S4M2, R_DimSymbol, DimIT); setAttrib(Froth_i_S4M2, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S4M2, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S4M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S4M2, install("DimCst"), dimIT); setAttrib(Froth_i_S4M2, install("DimCst"), dimIT);
                        r_Foth_i_S4M2 = REAL(Foth_i_S4M2);
                        r_Froth_i_S4M2 = REAL(Froth_i_S4M2);

                        PROTECT(Foth_i_S4M3 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S4M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S4M3, R_DimSymbol, DimIT); setAttrib(Froth_i_S4M3, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S4M3, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S4M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S4M3, install("DimCst"), dimIT); setAttrib(Froth_i_S4M3, install("DimCst"), dimIT);
                        r_Foth_i_S4M3 = REAL(Foth_i_S4M3);
                        r_Froth_i_S4M3 = REAL(Froth_i_S4M3);

                        PROTECT(Foth_i_S4M4 = NEW_NUMERIC(nbI*nbT));
                        PROTECT(Froth_i_S4M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_S4M4, R_DimSymbol, DimIT); setAttrib(Froth_i_S4M4, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_S4M4, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_S4M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_S4M4, install("DimCst"), dimIT); setAttrib(Froth_i_S4M4, install("DimCst"), dimIT);
                        r_Foth_i_S4M4 = REAL(Foth_i_S4M4);
                        r_Froth_i_S4M4 = REAL(Froth_i_S4M4);


                        PROTECT(FRWToth_i_S1M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S1M1, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S1M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S1M1, install("DimCst"), dimIT);
                        r_FRWToth_i_S1M1 = REAL(FRWToth_i_S1M1);

                        PROTECT(FRWToth_i_S1M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S1M2, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S1M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S1M2, install("DimCst"), dimIT);
                        r_FRWToth_i_S1M2 = REAL(FRWToth_i_S1M2);

                        PROTECT(FRWToth_i_S1M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S1M3, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S1M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S1M3, install("DimCst"), dimIT);
                        r_FRWToth_i_S1M3 = REAL(FRWToth_i_S1M3);

                        PROTECT(FRWToth_i_S1M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S1M4, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S1M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S1M4, install("DimCst"), dimIT);
                        r_FRWToth_i_S1M4 = REAL(FRWToth_i_S1M4);

                        PROTECT(FRWToth_i_S2M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S2M1, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S2M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S2M1, install("DimCst"), dimIT);
                        r_FRWToth_i_S2M1 = REAL(FRWToth_i_S2M1);

                        PROTECT(FRWToth_i_S2M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S2M2, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S2M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S2M2, install("DimCst"), dimIT);
                        r_FRWToth_i_S2M2 = REAL(FRWToth_i_S2M2);

                        PROTECT(FRWToth_i_S2M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S2M3, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S2M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S2M3, install("DimCst"), dimIT);
                        r_FRWToth_i_S2M3 = REAL(FRWToth_i_S2M3);

                        PROTECT(FRWToth_i_S2M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S2M4, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S2M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S2M4, install("DimCst"), dimIT);
                        r_FRWToth_i_S2M4 = REAL(FRWToth_i_S2M4);

                        PROTECT(FRWToth_i_S3M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S3M1, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S3M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S3M1, install("DimCst"), dimIT);
                        r_FRWToth_i_S3M1 = REAL(FRWToth_i_S3M1);

                        PROTECT(FRWToth_i_S3M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S3M2, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S3M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S3M2, install("DimCst"), dimIT);
                        r_FRWToth_i_S3M2 = REAL(FRWToth_i_S3M2);

                        PROTECT(FRWToth_i_S3M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S3M3, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S3M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S3M3, install("DimCst"), dimIT);
                        r_FRWToth_i_S3M3 = REAL(FRWToth_i_S3M3);

                        PROTECT(FRWToth_i_S3M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S3M4, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S3M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S3M4, install("DimCst"), dimIT);
                        r_FRWToth_i_S3M4 = REAL(FRWToth_i_S3M4);

                        PROTECT(FRWToth_i_S4M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S4M1, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S4M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S4M1, install("DimCst"), dimIT);
                        r_FRWToth_i_S4M1 = REAL(FRWToth_i_S4M1);

                        PROTECT(FRWToth_i_S4M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S4M2, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S4M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S4M2, install("DimCst"), dimIT);
                        r_FRWToth_i_S4M2 = REAL(FRWToth_i_S4M2);

                        PROTECT(FRWToth_i_S4M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S4M3, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S4M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S4M3, install("DimCst"), dimIT);
                        r_FRWToth_i_S4M3 = REAL(FRWToth_i_S4M3);

                        PROTECT(FRWToth_i_S4M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FRWToth_i_S4M4, R_DimSymbol, DimIT);
                        setAttrib(FRWToth_i_S4M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FRWToth_i_S4M4, install("DimCst"), dimIT);
                        r_FRWToth_i_S4M4 = REAL(FRWToth_i_S4M4);


                        PROTECT(FDWToth_i_S1M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S1M1, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S1M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S1M1, install("DimCst"), dimIT);
                        r_FDWToth_i_S1M1 = REAL(FDWToth_i_S1M1);

                        PROTECT(FDWToth_i_S1M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S1M2, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S1M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S1M2, install("DimCst"), dimIT);
                        r_FDWToth_i_S1M2 = REAL(FDWToth_i_S1M2);

                        PROTECT(FDWToth_i_S1M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S1M3, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S1M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S1M3, install("DimCst"), dimIT);
                        r_FDWToth_i_S1M3 = REAL(FDWToth_i_S1M3);

                        PROTECT(FDWToth_i_S1M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S1M4, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S1M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S1M4, install("DimCst"), dimIT);
                        r_FDWToth_i_S1M4 = REAL(FDWToth_i_S1M4);

                        PROTECT(FDWToth_i_S2M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S2M1, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S2M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S2M1, install("DimCst"), dimIT);
                        r_FDWToth_i_S2M1 = REAL(FDWToth_i_S2M1);

                        PROTECT(FDWToth_i_S2M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S2M2, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S2M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S2M2, install("DimCst"), dimIT);
                        r_FDWToth_i_S2M2 = REAL(FDWToth_i_S2M2);

                        PROTECT(FDWToth_i_S2M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S2M3, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S2M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S2M3, install("DimCst"), dimIT);
                        r_FDWToth_i_S2M3 = REAL(FDWToth_i_S2M3);

                        PROTECT(FDWToth_i_S2M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S2M4, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S2M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S2M4, install("DimCst"), dimIT);
                        r_FDWToth_i_S2M4 = REAL(FDWToth_i_S2M4);

                        PROTECT(FDWToth_i_S3M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S3M1, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S3M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S3M1, install("DimCst"), dimIT);
                        r_FDWToth_i_S3M1 = REAL(FDWToth_i_S3M1);

                        PROTECT(FDWToth_i_S3M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S3M2, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S3M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S3M2, install("DimCst"), dimIT);
                        r_FDWToth_i_S3M2 = REAL(FDWToth_i_S3M2);

                        PROTECT(FDWToth_i_S3M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S3M3, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S3M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S3M3, install("DimCst"), dimIT);
                        r_FDWToth_i_S3M3 = REAL(FDWToth_i_S3M3);

                        PROTECT(FDWToth_i_S3M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S3M4, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S3M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S3M4, install("DimCst"), dimIT);
                        r_FDWToth_i_S3M4 = REAL(FDWToth_i_S3M4);

                        PROTECT(FDWToth_i_S4M1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S4M1, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S4M1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S4M1, install("DimCst"), dimIT);
                        r_FDWToth_i_S4M1 = REAL(FDWToth_i_S4M1);

                        PROTECT(FDWToth_i_S4M2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S4M2, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S4M2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S4M2, install("DimCst"), dimIT);
                        r_FDWToth_i_S4M2 = REAL(FDWToth_i_S4M2);

                        PROTECT(FDWToth_i_S4M3 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S4M3, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S4M3, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S4M3, install("DimCst"), dimIT);
                        r_FDWToth_i_S4M3 = REAL(FDWToth_i_S4M3);

                        PROTECT(FDWToth_i_S4M4 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(FDWToth_i_S4M4, R_DimSymbol, DimIT);
                        setAttrib(FDWToth_i_S4M4, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(FDWToth_i_S4M4, install("DimCst"), dimIT);
                        r_FDWToth_i_S4M4 = REAL(FDWToth_i_S4M4);


                      } else if ((Qvec[e]==0) & (Svec[e]==0)){
                        PROTECT(Foth_i = NEW_NUMERIC(nbI*nbT)); //attention, on consid�re la mortalit� initiale comme �tant d�finie sans dimension temporelle --> � revoir
                        PROTECT(Froth_i = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i, R_DimSymbol, DimIT); setAttrib(Froth_i, R_DimSymbol, DimIT);
                        setAttrib(Foth_i, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i, install("DimCst"), dimIT); setAttrib(Froth_i, install("DimCst"), dimIT);

                        r_Foth_i = REAL(Foth_i);
                        r_Froth_i = REAL(Froth_i);

                      } else if ((Qvec[e]==0) & (Svec[e]==1)){

                        PROTECT(Foth_i_G1 = NEW_NUMERIC(nbI*nbT)); //attention, on consid�re la mortalit� initiale comme �tant d�finie sans dimension temporelle --> � revoir
                        PROTECT(Froth_i_G1 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_G1, R_DimSymbol, DimIT); setAttrib(Froth_i_G1, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_G1, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_G1, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_G1, install("DimCst"), dimIT); setAttrib(Froth_i_G1, install("DimCst"), dimIT);
                        r_Foth_i_G1 = REAL(Foth_i_G1);
                        r_Froth_i_G1 = REAL(Froth_i_G1);

                        PROTECT(Foth_i_G2 = NEW_NUMERIC(nbI*nbT)); //attention, on consid�re la mortalit� initiale comme �tant d�finie sans dimension temporelle --> � revoir
                        PROTECT(Froth_i_G2 = NEW_NUMERIC(nbI*nbT));
                        setAttrib(Foth_i_G2, R_DimSymbol, DimIT); setAttrib(Froth_i_G2, R_DimSymbol, DimIT);
                        setAttrib(Foth_i_G2, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i_G2, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i_G2, install("DimCst"), dimIT); setAttrib(Froth_i_G2, install("DimCst"), dimIT);
                        r_Foth_i_G2 = REAL(Foth_i_G2);
                        r_Froth_i_G2 = REAL(Froth_i_G2);
                      }
//Rprintf("Mort14\n");

                      if ((Qvec[e]==1) & (Svec[e]==0)) {

                       double *fothi = REAL(getListElement(elmt, "iniFothi_S1M1")) ; double *fothi2 = REAL(getListElement(elmt, "Fothi_S1M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S1M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S1M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S1M1[ind_i+ind_t*nbI] = r_Foth_i_S1M1[ind_i+ind_t*nbI] ;//*

                       fothi = REAL(getListElement(elmt, "iniFothi_S1M2")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S1M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M2[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S1M2[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S1M2[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S1M2[ind_i+ind_t*nbI] = r_Foth_i_S1M2[ind_i+ind_t*nbI];// *
                        if (ind_i==0) r_Froth_i_S1M2[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S1M3")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S1M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S1M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S1M3[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S1M3[ind_i+ind_t*nbI] = r_Foth_i_S1M3[ind_i+ind_t*nbI];// *
                        if (ind_i==0) r_Froth_i_S1M3[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S1M4")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S1M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S1M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S1M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S1M4[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S1M4[ind_i+ind_t*nbI] = r_Foth_i_S1M4[ind_i+ind_t*nbI] ;//*
                        if (ind_i==0) r_Froth_i_S1M4[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S2M1")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S2M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S2M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S2M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S2M1[ind_i+ind_t*nbI] = r_Foth_i_S2M1[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S2M2")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S2M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S2M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S2M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S2M2[ind_i+ind_t*nbI] = r_Foth_i_S2M2[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S2M3")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S2M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S2M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S2M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S2M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S2M3[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S2M3[ind_i+ind_t*nbI] = r_Foth_i_S2M3[ind_i+ind_t*nbI];// *
                        if (ind_i==0) r_Froth_i_S2M3[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S2M4")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S2M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S2M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S2M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S2M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S2M4[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S2M4[ind_i+ind_t*nbI] = r_Foth_i_S2M4[ind_i+ind_t*nbI];// *
                        if (ind_i==0) r_Froth_i_S2M4[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S3M1")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S3M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S3M1[ind_i+ind_t*nbI] = r_Foth_i_S3M1[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S3M2")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S3M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S3M2[ind_i+ind_t*nbI] = r_Foth_i_S3M2[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S3M3")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S3M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S3M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S3M3[ind_i+ind_t*nbI] = r_Foth_i_S3M3[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S3M4")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S3M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S3M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_Foth_i_S3M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Foth_i_S3M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_Foth_i_S3M4[ind_i+(ind_t+1)*nbI] = 0.0;}
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_Froth_i_S3M4[ind_i+ind_t*nbI] = r_Foth_i_S3M4[ind_i+ind_t*nbI];// *
                        if (ind_i==0) r_Froth_i_S3M4[ind_i+ind_t*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFothi_S4M1")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S4M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S4M1[ind_i+ind_t*nbI] = r_Foth_i_S4M1[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S4M2")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S4M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S4M2[ind_i+ind_t*nbI] = r_Foth_i_S4M2[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S4M3")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S4M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S4M3[ind_i+ind_t*nbI] = r_Foth_i_S4M3[ind_i+ind_t*nbI];// *

                       fothi = REAL(getListElement(elmt, "iniFothi_S4M4")) ; fothi2 = REAL(getListElement(elmt, "Fothi_S4M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M4[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i_S4M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];
                       for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Froth_i_S4M4[ind_i+ind_t*nbI] = r_Foth_i_S4M4[ind_i+ind_t*nbI];// *



                       fothi = REAL(getListElement(elmt, "iniFLWothi_S1M1")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S1M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S1M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S1M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S1M2")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S1M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M2[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S1M2[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S1M2[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S1M3")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S1M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S1M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S1M3[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S1M4")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S1M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S1M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S1M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S1M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S2M1")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S2M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S2M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S2M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S2M2")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S2M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S2M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S2M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S2M3")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S2M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S2M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S2M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S2M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S2M3[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S2M4")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S2M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S2M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S2M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S2M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S2M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S3M1")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S3M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S3M2")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S3M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S3M3")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S3M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S3M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S3M4")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S3M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S3M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FRWToth_i_S3M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FRWToth_i_S3M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FRWToth_i_S3M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S4M1")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S4M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S4M2")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S4M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S4M3")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S4M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFLWothi_S4M4")) ; fothi2 = REAL(getListElement(elmt, "FLWothi_S4M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M4[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FRWToth_i_S4M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];



                       fothi = REAL(getListElement(elmt, "iniFDWothi_S1M1")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S1M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S1M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S1M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S1M2")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S1M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M2[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S1M2[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S1M2[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S1M3")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S1M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S1M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S1M3[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S1M4")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S1M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S1M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S1M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S1M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S2M1")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S2M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S2M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S2M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S2M2")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S2M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S2M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S2M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S2M3")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S2M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S2M3[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S2M3[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S2M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S2M3[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S2M4")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S2M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S2M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S2M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S2M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S2M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S3M1")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S3M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S3M2")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S3M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S3M3")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S3M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S3M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S3M4")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S3M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S3M4[ind_i+ind_t*nbI] = fothi[ind_i]; if (ind_i==0) r_FDWToth_i_S3M4[ind_i+ind_t*nbI] = 0.0;}
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {r_FDWToth_i_S3M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i]; if (ind_i==0) r_FDWToth_i_S3M4[ind_i+(ind_t+1)*nbI] = 0.0;}

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S4M1")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S4M1"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M1[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M1[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S4M2")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S4M2"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M2[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M2[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S4M3")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S4M3"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M3[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M3[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                       fothi = REAL(getListElement(elmt, "iniFDWothi_S4M4")) ; fothi2 = REAL(getListElement(elmt, "FDWothi_S4M4"));
                       if (ind_t==0) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M4[ind_i+ind_t*nbI] = fothi[ind_i];
                       if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_FDWToth_i_S4M4[ind_i+(ind_t+1)*nbI] = fothi2[ind_i];

                      } else if ((Qvec[e]==0) & (Svec[e]==0)){
                    //if (ind_t==0) {//PrintValue(v_F_efmi); //PrintValue(aggregObj(v_F_efmi, dimI)); }

                        double *sumFtot = REAL(getListElement(elmt, "F_i"));

                        rdimI[3] = dimC[3];

                        if (ind_t==0) { //on initialise
                            SEXP gg0 = R_NilValue;
                            PROTECT(gg0=aggregObj(ans_11, dimI));
                            double *sumFr = REAL(gg0);
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                r_Foth_i[ind_i+ind_t*nbI] = fmax2(0.0 , sumFtot[ind_i+ind_t*nbI*(dimF[3]>0)] - sumFr[ind_i+ind_t*nbI*(dimC[3]>0)]); //ON N'INTEGRE PAS DE MORTALITES NEGATIVES

                                if (FOTHoptim_use & (e==eTemp)) {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                } else {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = r_Foth_i[ind_i+ind_t*nbI];
                                }
                            }
                            ////PrintValue(v_F_efmi);//PrintValue(dimI);PrintValue(Foth_i);
                            //Rprintf("Initialisation l.3580, Foth_i = "); PrintValue(Foth_i);
                            UNPROTECT(1);
                        } else {
                           if (ind_t<(nbT-1)) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                if (FOTHoptim_use & (e==eTemp)) {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                } else {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = r_Foth_i[ind_i+ind_t*nbI];
                                }
                            }

                           }
                        }


                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                                r_Froth_i[ind_i+ind_t*nbI] = r_Foth_i[ind_i+ind_t*nbI] *
                                    (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                                    r_doth_ei[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);

                      } else if ((Qvec[e]==0) & (Svec[e]==1)){

                        double *sumFtot_G1 = REAL(getListElement(elmt, "F_i_G1"));
                        double *sumFtot_G2 = REAL(getListElement(elmt, "F_i_G2"));

                        rdimI[3] = dimC[3];
                        double *sumFr_G1 = REAL(aggregObj(ans_11_G1, dimI));
                        double *sumFr_G2 = REAL(aggregObj(ans_11_G2, dimI));

                        if (ind_t==0) { //on initialise
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                r_Foth_i_G1[ind_i+ind_t*nbI] = fmax2(0.0 , sumFtot_G1[ind_i+ind_t*nbI*(dimF[3]>0)] - sumFr_G1[ind_i+ind_t*nbI*(dimC[3]>0)]); //ON N'INTEGRE PAS DE MORTALITES NEGATIVES
                                r_Foth_i_G2[ind_i+ind_t*nbI] = fmax2(0.0 , sumFtot_G2[ind_i+ind_t*nbI*(dimF[3]>0)] - sumFr_G2[ind_i+ind_t*nbI*(dimC[3]>0)]);

                                if (FOTHoptim_use & (e==eTemp)) {
                                    r_Foth_i_G1[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                    r_Foth_i_G2[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                } else {
                                    r_Foth_i_G1[ind_i+(ind_t+1)*nbI] = r_Foth_i_G1[ind_i+ind_t*nbI];
                                    r_Foth_i_G2[ind_i+(ind_t+1)*nbI] = r_Foth_i_G2[ind_i+ind_t*nbI];
                                }
                            }
                            ////PrintValue(v_F_efmi);//PrintValue(dimI);//PrintValue(Foth_i);
                        } else {
                           if (ind_t<(nbT-1)) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                if (FOTHoptim_use & (e==eTemp)) {
                                    r_Foth_i_G1[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                    r_Foth_i_G2[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                } else {
                                    r_Foth_i_G1[ind_i+(ind_t+1)*nbI] = r_Foth_i_G1[ind_i+ind_t*nbI];
                                    r_Foth_i_G2[ind_i+(ind_t+1)*nbI] = r_Foth_i_G2[ind_i+ind_t*nbI];
                                }
                            }

                           }
                        }


                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                r_Froth_i_G1[ind_i+ind_t*nbI] = r_Foth_i_G1[ind_i+ind_t*nbI] *
                                    (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                                    r_doth_ei_G1[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);
                                r_Froth_i_G2[ind_i+ind_t*nbI] = r_Foth_i_G2[ind_i+ind_t*nbI] *
                                    (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                                    r_doth_ei_G2[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);
                            }

                      }


//Rprintf("Mort15\n");

                        //on n'oublie pas d'archiver dans eVar ce dont on aura besoin dans les it�rations suivantes
                        if ((Qvec[e]==0) & (Svec[e]==0)) SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 0, v_F_efmi2); //ESSENTIEL!!! : ne pas laisser d'ind�fini en premier �l�ment d'une liste ; il vaut mieux laisser la partie telle qu'initialis�e
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 1, formatEff);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 2, v_Sr_e);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 3, v_d_efi);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 4, fFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 5, fFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 6, fFACT3);//Rprintf("SOURCE\n\n\n");//PrintValue(fFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 7, fFACT4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 8, fFACT5);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 61, fFACT6);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 44, Foth_i); //Rprintf("Dans EVAR - Apres module mortalite, Foth_i = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 60, Froth_i);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 50, fFACTsup1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 51, fFACTsup2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 52, v_nbNav_f);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 53, v_nbds_f);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 100, v_F_efmi2_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 101, v_F_efmi2_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 102, v_F_efmi2_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 103, v_F_efmi2_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 104, v_F_efmi2_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 105, v_F_efmi2_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 106, v_F_efmi2_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 107, v_F_efmi2_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 108, v_F_efmi2_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 109, v_F_efmi2_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 110, v_F_efmi2_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 111, v_F_efmi2_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 112, v_F_efmi2_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 113, v_F_efmi2_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 114, v_F_efmi2_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 115, v_F_efmi2_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 116, Foth_i_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 117, Foth_i_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 118, Foth_i_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 119, Foth_i_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 120, Foth_i_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 121, Foth_i_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 122, Foth_i_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 123, Foth_i_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 124, Foth_i_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 125, Foth_i_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 126, Foth_i_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 127, Foth_i_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 128, Foth_i_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 129, Foth_i_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 130, Foth_i_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 131, Foth_i_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 132, Froth_i_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 133, Froth_i_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 134, Froth_i_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 135, Froth_i_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 136, Froth_i_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 137, Froth_i_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 138, Froth_i_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 139, Froth_i_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 140, Froth_i_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 141, Froth_i_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 142, Froth_i_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 143, Froth_i_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 144, Froth_i_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 145, Froth_i_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 146, Froth_i_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 147, Froth_i_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 160, v_FRWT_efmi2_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 161, v_FRWT_efmi2_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 162, v_FRWT_efmi2_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 163, v_FRWT_efmi2_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 164, v_FRWT_efmi2_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 165, v_FRWT_efmi2_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 166, v_FRWT_efmi2_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 167, v_FRWT_efmi2_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 168, v_FRWT_efmi2_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 169, v_FRWT_efmi2_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 170, v_FRWT_efmi2_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 171, v_FRWT_efmi2_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 172, v_FRWT_efmi2_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 173, v_FRWT_efmi2_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 174, v_FRWT_efmi2_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 175, v_FRWT_efmi2_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 176, FRWToth_i_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 177, FRWToth_i_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 178, FRWToth_i_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 179, FRWToth_i_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 180, FRWToth_i_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 181, FRWToth_i_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 182, FRWToth_i_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 183, FRWToth_i_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 184, FRWToth_i_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 185, FRWToth_i_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 186, FRWToth_i_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 187, FRWToth_i_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 188, FRWToth_i_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 189, FRWToth_i_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 190, FRWToth_i_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 191, FRWToth_i_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 192, v_FDWT_efmi2_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 193, v_FDWT_efmi2_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 194, v_FDWT_efmi2_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 195, v_FDWT_efmi2_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 196, v_FDWT_efmi2_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 197, v_FDWT_efmi2_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 198, v_FDWT_efmi2_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 199, v_FDWT_efmi2_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 200, v_FDWT_efmi2_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 201, v_FDWT_efmi2_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 202, v_FDWT_efmi2_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 203, v_FDWT_efmi2_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 204, v_FDWT_efmi2_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 205, v_FDWT_efmi2_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 206, v_FDWT_efmi2_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 207, v_FDWT_efmi2_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 208, FDWToth_i_S1M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 209, FDWToth_i_S1M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 210, FDWToth_i_S1M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 211, FDWToth_i_S1M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 212, FDWToth_i_S2M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 213, FDWToth_i_S2M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 214, FDWToth_i_S2M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 215, FDWToth_i_S2M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 216, FDWToth_i_S3M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 217, FDWToth_i_S3M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 218, FDWToth_i_S3M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 219, FDWToth_i_S3M4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 220, FDWToth_i_S4M1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 221, FDWToth_i_S4M2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 222, FDWToth_i_S4M3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 223, FDWToth_i_S4M4);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 224, Foth_i_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 225, Foth_i_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 226, Froth_i_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 227, Froth_i_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 228, v_F_efmi_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 229, v_F_efmi_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 230, v_F_efmi2_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 231, v_F_efmi2_G2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 232, v_d_efi_G1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 233, v_d_efi_G2);
//Rprintf("Mort16\n");

                    //if (indP==1) UNPROTECT(1);
                    //Rprintf("K3\n");
                    if ((Qvec[e]==1) & (Svec[e]==0)) {
                     UNPROTECT(103);
                     UNPROTECT(96);
                     UNPROTECT(64);
                    } else if ((Qvec[e]==0) & (Svec[e]==0)){
                     UNPROTECT(11);
                    }
                    else if ((Qvec[e]==0) & (Svec[e]==1)){
                     UNPROTECT(17);
                    }

                    if (ind_t==0){
                        if ((Qvec[e]==1) & (Svec[e]==0)) {
                                UNPROTECT(64+1);
                        } else if ((Qvec[e]==0) & (Svec[e]==0)) {
                                UNPROTECT(2+1);
                        } else if ((Qvec[e]==0) & (Svec[e]==1)) {
                                UNPROTECT(4+1);
                                }
                    }
                    UNPROTECT(23);
}

fUpdate = false;
//Rprintf("K4\n");
if (ind_t==0) UNPROTECT(1);
UNPROTECT(2);

} else {
    //fichier2 << "fUpdate = " << fUpdate << endl;
    //fichier2 << "ind_t = " << ind_t << endl;

for (int e = 0 ; e < nbE ; e++) {

//Rprintf("Mort20\n");fichier2 << "Mort20" << endl;

                    int nbI = length(VECTOR_ELT(namDC,e));

                    SEXP elmt;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                    double *r_nbv_f, *r_nbds_f, *r_nbds2_f;
                    r_nbv_f = REAL(getListElement(Flist, "nbv_f_m"));
                    r_nbds_f = REAL(getListElement(Flist, "effort1_f_m"));
                    r_nbds2_f = REAL(getListElement(Flist, "effort2_f_m"));
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 4));//Rprintf("Mort20.1\n");//Rprintf("%i %i lgth",length(VECTOR_ELT(EVAR, 0)),length(VECTOR_ELT(EVAR, 1)));
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 5));//Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 6));//Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 7));//Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 61));//Rprintf("Mort20.5\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 50));//Rprintf("Mort20.6\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 51));//Rprintf("Mort20.7\n");
//Rprintf("MortZ0\n");
                    int *fFact1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 4));//Rprintf("MortZ1\n");
                    int *fFact2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 5));//Rprintf("MortZ2\n");//PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 6));
                    int *fFact3 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 6));//Rprintf("MortZ3\n");
                    int *fFact4 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 7));//Rprintf("MortZ4\n");
                    int *fFact6 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 61));//Rprintf("MortZ5\n");
                    int *fFactSup1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 50));//Rprintf("MortZ6\n");
                    int *fFactSup2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 51));//Rprintf("MortZ7\n");
//Rprintf("Mort21\n"); fichier2 << "Mort21" << endl;
                    if ((Qvec[e]==1) & (Svec[e]==0)) {

                    double *rans_11_S1M1 = REAL(VECTOR_ELT(out_F_fmi_S1M1, e));
                    double *rans_11_S1M2 = REAL(VECTOR_ELT(out_F_fmi_S1M2, e));
                    double *rans_11_S1M3 = REAL(VECTOR_ELT(out_F_fmi_S1M3, e));
                    double *rans_11_S1M4 = REAL(VECTOR_ELT(out_F_fmi_S1M4, e));
                    double *rans_11_S2M1 = REAL(VECTOR_ELT(out_F_fmi_S2M1, e));
                    double *rans_11_S2M2 = REAL(VECTOR_ELT(out_F_fmi_S2M2, e));
                    double *rans_11_S2M3 = REAL(VECTOR_ELT(out_F_fmi_S2M3, e));
                    double *rans_11_S2M4 = REAL(VECTOR_ELT(out_F_fmi_S2M4, e));
                    double *rans_11_S3M1 = REAL(VECTOR_ELT(out_F_fmi_S3M1, e));
                    double *rans_11_S3M2 = REAL(VECTOR_ELT(out_F_fmi_S3M2, e));
                    double *rans_11_S3M3 = REAL(VECTOR_ELT(out_F_fmi_S3M3, e));
                    double *rans_11_S3M4 = REAL(VECTOR_ELT(out_F_fmi_S3M4, e));
                    double *rans_11_S4M1 = REAL(VECTOR_ELT(out_F_fmi_S4M1, e));
                    double *rans_11_S4M2 = REAL(VECTOR_ELT(out_F_fmi_S4M2, e));
                    double *rans_11_S4M3 = REAL(VECTOR_ELT(out_F_fmi_S4M3, e));
                    double *rans_11_S4M4 = REAL(VECTOR_ELT(out_F_fmi_S4M4, e));

                    double *rans_11l_S1M1 = REAL(VECTOR_ELT(out_Fr_fmi_S1M1, e));
                    double *rans_11l_S1M2 = REAL(VECTOR_ELT(out_Fr_fmi_S1M2, e));
                    double *rans_11l_S1M3 = REAL(VECTOR_ELT(out_Fr_fmi_S1M3, e));
                    double *rans_11l_S1M4 = REAL(VECTOR_ELT(out_Fr_fmi_S1M4, e));
                    double *rans_11l_S2M1 = REAL(VECTOR_ELT(out_Fr_fmi_S2M1, e));
                    double *rans_11l_S2M2 = REAL(VECTOR_ELT(out_Fr_fmi_S2M2, e));
                    double *rans_11l_S2M3 = REAL(VECTOR_ELT(out_Fr_fmi_S2M3, e));
                    double *rans_11l_S2M4 = REAL(VECTOR_ELT(out_Fr_fmi_S2M4, e));
                    double *rans_11l_S3M1 = REAL(VECTOR_ELT(out_Fr_fmi_S3M1, e));
                    double *rans_11l_S3M2 = REAL(VECTOR_ELT(out_Fr_fmi_S3M2, e));
                    double *rans_11l_S3M3 = REAL(VECTOR_ELT(out_Fr_fmi_S3M3, e));
                    double *rans_11l_S3M4 = REAL(VECTOR_ELT(out_Fr_fmi_S3M4, e));
                    double *rans_11l_S4M1 = REAL(VECTOR_ELT(out_Fr_fmi_S4M1, e));
                    double *rans_11l_S4M2 = REAL(VECTOR_ELT(out_Fr_fmi_S4M2, e));
                    double *rans_11l_S4M3 = REAL(VECTOR_ELT(out_Fr_fmi_S4M3, e));
                    double *rans_11l_S4M4 = REAL(VECTOR_ELT(out_Fr_fmi_S4M4, e));

                    double *rans_FRWT_S1M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M1, e));
                    double *rans_FRWT_S1M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M2, e));
                    double *rans_FRWT_S1M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M3, e));
                    double *rans_FRWT_S1M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S1M4, e));
                    double *rans_FRWT_S2M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M1, e));
                    double *rans_FRWT_S2M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M2, e));
                    double *rans_FRWT_S2M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M3, e));
                    double *rans_FRWT_S2M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S2M4, e));
                    double *rans_FRWT_S3M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M1, e));
                    double *rans_FRWT_S3M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M2, e));
                    double *rans_FRWT_S3M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M3, e));
                    double *rans_FRWT_S3M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S3M4, e));
                    double *rans_FRWT_S4M1 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M1, e));
                    double *rans_FRWT_S4M2 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M2, e));
                    double *rans_FRWT_S4M3 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M3, e));
                    double *rans_FRWT_S4M4 = REAL(VECTOR_ELT(out_FRWT_fmi_S4M4, e));

                    double *rans_FDWT_S1M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M1, e));
                    double *rans_FDWT_S1M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M2, e));
                    double *rans_FDWT_S1M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M3, e));
                    double *rans_FDWT_S1M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S1M4, e));
                    double *rans_FDWT_S2M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M1, e));
                    double *rans_FDWT_S2M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M2, e));
                    double *rans_FDWT_S2M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M3, e));
                    double *rans_FDWT_S2M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S2M4, e));
                    double *rans_FDWT_S3M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M1, e));
                    double *rans_FDWT_S3M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M2, e));
                    double *rans_FDWT_S3M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M3, e));
                    double *rans_FDWT_S3M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S3M4, e));
                    double *rans_FDWT_S4M1 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M1, e));
                    double *rans_FDWT_S4M2 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M2, e));
                    double *rans_FDWT_S4M3 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M3, e));
                    double *rans_FDWT_S4M4 = REAL(VECTOR_ELT(out_FDWT_fmi_S4M4, e));

                    double *r_F_efmi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 100));
                    double *r_F_efmi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 101));
                    double *r_F_efmi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 102));
                    double *r_F_efmi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 103));
                    double *r_F_efmi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 104));
                    double *r_F_efmi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 105));
                    double *r_F_efmi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 106));
                    double *r_F_efmi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 107));
                    double *r_F_efmi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 108));
                    double *r_F_efmi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 109));
                    double *r_F_efmi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 110));
                    double *r_F_efmi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 111));
                    double *r_F_efmi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 112));
                    double *r_F_efmi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 113));
                    double *r_F_efmi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 114));
                    double *r_F_efmi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 115));

                    double *r_Foth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 116));
                    double *r_Foth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 117));
                    double *r_Foth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 118));
                    double *r_Foth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 119));
                    double *r_Foth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 120));
                    double *r_Foth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 121));
                    double *r_Foth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 122));
                    double *r_Foth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 123));
                    double *r_Foth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 124));
                    double *r_Foth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 125));
                    double *r_Foth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 126));
                    double *r_Foth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 127));
                    double *r_Foth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 128));
                    double *r_Foth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 129));
                    double *r_Foth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 130));
                    double *r_Foth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 131));

                    double *r_Froth_it_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 132));
                    double *r_Froth_it_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 133));
                    double *r_Froth_it_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 134));
                    double *r_Froth_it_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 135));
                    double *r_Froth_it_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 136));
                    double *r_Froth_it_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 137));
                    double *r_Froth_it_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 138));
                    double *r_Froth_it_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 139));
                    double *r_Froth_it_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 140));
                    double *r_Froth_it_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 141));
                    double *r_Froth_it_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 142));
                    double *r_Froth_it_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 143));
                    double *r_Froth_it_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 144));
                    double *r_Froth_it_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 145));
                    double *r_Froth_it_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 146));
                    double *r_Froth_it_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 147));

                    double *r_FRWT_efmi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 160));
                    double *r_FRWT_efmi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 161));
                    double *r_FRWT_efmi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 162));
                    double *r_FRWT_efmi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 163));
                    double *r_FRWT_efmi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 164));
                    double *r_FRWT_efmi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 165));
                    double *r_FRWT_efmi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 166));
                    double *r_FRWT_efmi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 167));
                    double *r_FRWT_efmi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 168));
                    double *r_FRWT_efmi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 169));
                    double *r_FRWT_efmi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 170));
                    double *r_FRWT_efmi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 171));
                    double *r_FRWT_efmi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 172));
                    double *r_FRWT_efmi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 173));
                    double *r_FRWT_efmi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 174));
                    double *r_FRWT_efmi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 175));

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

                    double *r_FDWT_efmi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 192));
                    double *r_FDWT_efmi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 193));
                    double *r_FDWT_efmi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 194));
                    double *r_FDWT_efmi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 195));
                    double *r_FDWT_efmi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 196));
                    double *r_FDWT_efmi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 197));
                    double *r_FDWT_efmi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 198));
                    double *r_FDWT_efmi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 199));
                    double *r_FDWT_efmi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 200));
                    double *r_FDWT_efmi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 201));
                    double *r_FDWT_efmi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 202));
                    double *r_FDWT_efmi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 203));
                    double *r_FDWT_efmi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 204));
                    double *r_FDWT_efmi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 205));
                    double *r_FDWT_efmi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 206));
                    double *r_FDWT_efmi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 207));

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

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_11_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_11_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];


                    rans_11l_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];

                    rans_11l_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        rans_11_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];



                    rans_FRWT_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

//if ((ind_t==1) & (ind_f==0)) {
//
//std::stringstream tstt1,tstt2,tstt3,tstt4,tstt5;
//tstt1 << r_FRWT_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]];
//tstt2 << r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]];
//tstt3 << r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
//tstt4 << r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
//tstt5 << rans_FRWT_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]];
//
//fichier2 << "detail_T1" << tstt1.str() << " " << tstt2.str() << " " << tstt3.str() << " " << tstt4.str() << " " << tstt5.str() << endl;
//
//}


                    rans_FRWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FRWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FRWT_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FRWT_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FRWT_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];



                    rans_FDWT_S1M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S1M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S1M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S1M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S1M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S1M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S1M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S1M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S2M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S2M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S2M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S2M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S2M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S2M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S2M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S2M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S3M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S3M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S3M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S3M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S3M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S3M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S3M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];
                    if (ind_i==0) rans_FDWT_S3M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] = 0.0;

                    rans_FDWT_S4M1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S4M1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S4M2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S4M2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S4M3[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S4M3[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_FDWT_S4M4[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_FDWT_efmi_S4M4[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    }

                    if (ind_t<(nbT-1)) for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                            r_Foth_it_S1M1[ind_i+(ind_t+1)*nbI] = r_Foth_it_S1M1[ind_i+ind_t*nbI];
                            r_Foth_it_S1M2[ind_i+(ind_t+1)*nbI] = r_Foth_it_S1M2[ind_i+ind_t*nbI];
                            r_Foth_it_S1M3[ind_i+(ind_t+1)*nbI] = r_Foth_it_S1M3[ind_i+ind_t*nbI];
                            r_Foth_it_S1M4[ind_i+(ind_t+1)*nbI] = r_Foth_it_S1M4[ind_i+ind_t*nbI];
                            r_Foth_it_S2M1[ind_i+(ind_t+1)*nbI] = r_Foth_it_S2M1[ind_i+ind_t*nbI];
                            r_Foth_it_S2M2[ind_i+(ind_t+1)*nbI] = r_Foth_it_S2M2[ind_i+ind_t*nbI];
                            r_Foth_it_S2M3[ind_i+(ind_t+1)*nbI] = r_Foth_it_S2M3[ind_i+ind_t*nbI];
                            r_Foth_it_S2M4[ind_i+(ind_t+1)*nbI] = r_Foth_it_S2M4[ind_i+ind_t*nbI];
                            r_Foth_it_S3M1[ind_i+(ind_t+1)*nbI] = r_Foth_it_S3M1[ind_i+ind_t*nbI];
                            r_Foth_it_S3M2[ind_i+(ind_t+1)*nbI] = r_Foth_it_S3M2[ind_i+ind_t*nbI];
                            r_Foth_it_S3M3[ind_i+(ind_t+1)*nbI] = r_Foth_it_S3M3[ind_i+ind_t*nbI];
                            r_Foth_it_S3M4[ind_i+(ind_t+1)*nbI] = r_Foth_it_S3M4[ind_i+ind_t*nbI];
                            r_Foth_it_S4M1[ind_i+(ind_t+1)*nbI] = r_Foth_it_S4M1[ind_i+ind_t*nbI];
                            r_Foth_it_S4M2[ind_i+(ind_t+1)*nbI] = r_Foth_it_S4M2[ind_i+ind_t*nbI];
                            r_Foth_it_S4M3[ind_i+(ind_t+1)*nbI] = r_Foth_it_S4M3[ind_i+ind_t*nbI];
                            r_Foth_it_S4M4[ind_i+(ind_t+1)*nbI] = r_Foth_it_S4M4[ind_i+ind_t*nbI];

                            r_FRWToth_it_S1M1[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S1M1[ind_i+ind_t*nbI];
                            r_FRWToth_it_S1M2[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S1M2[ind_i+ind_t*nbI];
                            r_FRWToth_it_S1M3[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S1M3[ind_i+ind_t*nbI];
                            r_FRWToth_it_S1M4[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S1M4[ind_i+ind_t*nbI];
                            r_FRWToth_it_S2M1[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S2M1[ind_i+ind_t*nbI];
                            r_FRWToth_it_S2M2[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S2M2[ind_i+ind_t*nbI];
                            r_FRWToth_it_S2M3[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S2M3[ind_i+ind_t*nbI];
                            r_FRWToth_it_S2M4[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S2M4[ind_i+ind_t*nbI];
                            r_FRWToth_it_S3M1[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S3M1[ind_i+ind_t*nbI];
                            r_FRWToth_it_S3M2[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S3M2[ind_i+ind_t*nbI];
                            r_FRWToth_it_S3M3[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S3M3[ind_i+ind_t*nbI];
                            r_FRWToth_it_S3M4[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S3M4[ind_i+ind_t*nbI];
                            r_FRWToth_it_S4M1[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S4M1[ind_i+ind_t*nbI];
                            r_FRWToth_it_S4M2[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S4M2[ind_i+ind_t*nbI];
                            r_FRWToth_it_S4M3[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S4M3[ind_i+ind_t*nbI];
                            r_FRWToth_it_S4M4[ind_i+(ind_t+1)*nbI] = r_FRWToth_it_S4M4[ind_i+ind_t*nbI];

                            r_FDWToth_it_S1M1[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S1M1[ind_i+ind_t*nbI];
                            r_FDWToth_it_S1M2[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S1M2[ind_i+ind_t*nbI];
                            r_FDWToth_it_S1M3[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S1M3[ind_i+ind_t*nbI];
                            r_FDWToth_it_S1M4[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S1M4[ind_i+ind_t*nbI];
                            r_FDWToth_it_S2M1[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S2M1[ind_i+ind_t*nbI];
                            r_FDWToth_it_S2M2[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S2M2[ind_i+ind_t*nbI];
                            r_FDWToth_it_S2M3[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S2M3[ind_i+ind_t*nbI];
                            r_FDWToth_it_S2M4[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S2M4[ind_i+ind_t*nbI];
                            r_FDWToth_it_S3M1[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S3M1[ind_i+ind_t*nbI];
                            r_FDWToth_it_S3M2[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S3M2[ind_i+ind_t*nbI];
                            r_FDWToth_it_S3M3[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S3M3[ind_i+ind_t*nbI];
                            r_FDWToth_it_S3M4[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S3M4[ind_i+ind_t*nbI];
                            r_FDWToth_it_S4M1[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S4M1[ind_i+ind_t*nbI];
                            r_FDWToth_it_S4M2[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S4M2[ind_i+ind_t*nbI];
                            r_FDWToth_it_S4M3[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S4M3[ind_i+ind_t*nbI];
                            r_FDWToth_it_S4M4[ind_i+(ind_t+1)*nbI] = r_FDWToth_it_S4M4[ind_i+ind_t*nbI];

                    }

                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                        r_Froth_it_S1M1[ind_i+ind_t*nbI] = r_Foth_it_S1M1[ind_i+ind_t*nbI];
                        r_Froth_it_S1M2[ind_i+ind_t*nbI] = r_Foth_it_S1M2[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S1M2[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S1M3[ind_i+ind_t*nbI] = r_Foth_it_S1M3[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S1M3[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S1M4[ind_i+ind_t*nbI] = r_Foth_it_S1M4[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S1M4[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S2M1[ind_i+ind_t*nbI] = r_Foth_it_S2M1[ind_i+ind_t*nbI];
                        r_Froth_it_S2M2[ind_i+ind_t*nbI] = r_Foth_it_S2M2[ind_i+ind_t*nbI];
                        r_Froth_it_S2M3[ind_i+ind_t*nbI] = r_Foth_it_S2M3[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S2M3[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S2M4[ind_i+ind_t*nbI] = r_Foth_it_S2M4[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S2M4[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S3M1[ind_i+ind_t*nbI] = r_Foth_it_S3M1[ind_i+ind_t*nbI];
                        r_Froth_it_S3M2[ind_i+ind_t*nbI] = r_Foth_it_S3M2[ind_i+ind_t*nbI];
                        r_Froth_it_S3M3[ind_i+ind_t*nbI] = r_Foth_it_S3M3[ind_i+ind_t*nbI];
                        r_Froth_it_S3M4[ind_i+ind_t*nbI] = r_Foth_it_S3M4[ind_i+ind_t*nbI];
                        if (ind_i==0) r_Froth_it_S3M4[ind_i+ind_t*nbI] = 0.0;
                        r_Froth_it_S4M1[ind_i+ind_t*nbI] = r_Foth_it_S4M1[ind_i+ind_t*nbI];
                        r_Froth_it_S4M2[ind_i+ind_t*nbI] = r_Foth_it_S4M2[ind_i+ind_t*nbI];
                        r_Froth_it_S4M3[ind_i+ind_t*nbI] = r_Foth_it_S4M3[ind_i+ind_t*nbI];
                        r_Froth_it_S4M4[ind_i+ind_t*nbI] = r_Foth_it_S4M4[ind_i+ind_t*nbI];

                    }

                    } else if ((Qvec[e]==0) & (Svec[e]==0)) {
                    double *r_Sr_e = REAL(getListElement(elmt, "sr"));
                    double *r_d_efi = REAL(getListElement(elmt, "d_i"));
                    double *r_doth_ei = REAL(getListElement(elmt, "doth_i"));

                    double *rans_11 = REAL(VECTOR_ELT(out_F_fmi, e));
                    double *rans_11l = REAL(VECTOR_ELT(out_Fr_fmi, e));
                    double *r_F_efmi = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 0)); // attention capturabilite et non F
                    double *r_Foth_it = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44)); //Rprintf("Dans EVAR (l.4782), Fothi = "); PrintValue(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                    double *r_Froth_it = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];


                    rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                        r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                    }

                    if (ind_t<(nbT-1)) {
                        if (FOTHoptim_use & (e==eTemp)) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_it[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                        } else {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_it[ind_i+(ind_t+1)*nbI] = r_Foth_it[ind_i+ind_t*nbI];   //� modifier quand on consid�rera une mortalit� "autres" variable
                        }
                    }

                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                        r_Froth_it[ind_i+ind_t*nbI] = r_Foth_it[ind_i+ind_t*nbI] *
                            (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            r_doth_ei[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);

                    } else if ((Qvec[e]==0) & (Svec[e]==1)) {
                    double *r_Sr_e = REAL(getListElement(elmt, "sr"));
                    double *r_d_efi_G1 = REAL(getListElement(elmt, "d_i_G1"));
                    double *r_doth_ei_G1 = REAL(getListElement(elmt, "doth_i_G1"));
                    double *r_d_efi_G2 = REAL(getListElement(elmt, "d_i_G2"));
                    double *r_doth_ei_G2 = REAL(getListElement(elmt, "doth_i_G2"));

                    double *rans_11_G1 = REAL(VECTOR_ELT(out_F_fmi_G1, e));
                    double *rans_11l_G1 = REAL(VECTOR_ELT(out_Fr_fmi_G1, e));
                    double *r_F_efmi_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 230));// attention capturabilite et non F
                    double *r_Foth_it_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 224));
                    double *r_Froth_it_G1 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 226));

                    double *rans_11_G2 = REAL(VECTOR_ELT(out_F_fmi_G2, e));
                    double *rans_11l_G2 = REAL(VECTOR_ELT(out_Fr_fmi_G2, e));
                    double *r_F_efmi_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 231));// attention capturabilite et non F
                    double *r_Foth_it_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 225));
                    double *r_Froth_it_G2 = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 227));

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11_G1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_G1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];


                    rans_11l_G1[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_G1[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                        r_d_efi_G1[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);


                    rans_11_G2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_G2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]]*
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];


                    rans_11l_G2[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi_G2[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        r_nbds2_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                        r_d_efi_G2[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                    }

                    if (ind_t<(nbT-1)) {
                        if (FOTHoptim_use & (e==eTemp)) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                    r_Foth_it_G1[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                    r_Foth_it_G2[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                            }
                        } else {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                r_Foth_it_G1[ind_i+(ind_t+1)*nbI] = r_Foth_it_G1[ind_i+ind_t*nbI];   //� modifier quand on consid�rera une mortalit� "autres" variable
                                r_Foth_it_G2[ind_i+(ind_t+1)*nbI] = r_Foth_it_G2[ind_i+ind_t*nbI];
                                //fichier2 << "Age " << ind_i <<  " Foth dans Evar : " << r_Foth_it_G1[ind_i+ind_t*nbI] << endl;
                            }
                        }
                    }

                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                        r_Froth_it_G1[ind_i+ind_t*nbI] = r_Foth_it_G1[ind_i+ind_t*nbI] *
                            (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            r_doth_ei_G1[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);
                        r_Froth_it_G2[ind_i+ind_t*nbI] = r_Foth_it_G2[ind_i+ind_t*nbI] *
                            (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            r_doth_ei_G2[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);

                    }


                    }

//Rprintf("K5\n"); fichier2 << "K5" << endl;
                UNPROTECT(1);
//Rprintf("Mort\n");//PrintValue(out_Fbar_et);
}
}
UNPROTECT(1);
//Rprintf("End Mortalite\n");
//fichier2 << "End" << endl;

//fichier2.close();

}}

