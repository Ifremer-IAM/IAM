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

// #include "Modules.h" // Contient tout les modules
#include "BioEcoPar.h" // Class is defined in this file.

//using namespace Rcpp;
using namespace std;

//d�tection d'un caract�re donn� dans un objet SEXP de type SXPSTR

extern "C"
{
    bool isCharIn(SEXP names, const char *str)
    {
        int i;
        bool test = false;

        for (i = 0; i < length(names); i++)
            if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
            {
                test = true;
                break;
            }
        return test;
    }
}

extern "C" {
SEXP IAM(SEXP listInput, SEXP listSpec, SEXP listStochastic, SEXP listScen,
            SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, SEXP Bootstrp, SEXP nbBoot,
            SEXP GestInd, SEXP mOth, SEXP bounds, SEXP TAC, SEXP FBAR, SEXP othSpSup, SEXP effSup, SEXP GestParam, SEXP EcoDcf,
            SEXP EcoInd, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyF, SEXP Ftarg, SEXP W_Ftarg, SEXP MeanRec_Ftarg,
            SEXP parBHV, SEXP parQEX,
            SEXP tacCTRL, SEXP stochPrice, SEXP updateE, SEXP parOQD, SEXP bootVar = R_NilValue, SEXP verbose = 0)
{
int VERBOSE = INTEGER(verbose)[0];
//Rprintf("OO");
    if (INTEGER(Bootstrp)[0]==0) {
        
        if(VERBOSE){Rprintf("No bootstrap \n");}
        BioEcoPar *object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                            RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                            GestInd, mOth, bounds, TAC, FBAR, othSpSup, effSup, GestParam, EcoDcf,
                                            EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, Ftarg, W_Ftarg, MeanRec_Ftarg,
                                            parBHV, parQEX, tacCTRL, stochPrice, updateE, parOQD, VERBOSE);


        SEXP output, out_names, out_Foth, out_Foth_G1,out_Foth_G2;
        PROTECT(output = allocVector(VECSXP, 128)); //11/04/18 rajout de l'�l�ment reconcilSPP
        SET_VECTOR_ELT(output, 0, object->out_F_fmi);
        SET_VECTOR_ELT(output, 1, object->out_Z_eit);
        SET_VECTOR_ELT(output, 2, object->out_Fbar_et);
        SET_VECTOR_ELT(output, 3, object->out_N_eit);
        SET_VECTOR_ELT(output, 4, object->out_B_et);
        SET_VECTOR_ELT(output, 5, object->out_SSB_et);
        SET_VECTOR_ELT(output, 6, object->out_C_efmit);
        SET_VECTOR_ELT(output, 7, object->out_C_eit);
        SET_VECTOR_ELT(output, 8, object->out_Y_efmit);
        SET_VECTOR_ELT(output, 9, object->out_Y_eit);
        SET_VECTOR_ELT(output, 10, object->out_D_efmit);
        SET_VECTOR_ELT(output, 11, object->out_L_efmit);
        SET_VECTOR_ELT(output, 12, object->out_L_efmct);
        SET_VECTOR_ELT(output, 13, object->out_L_eit);
        SET_VECTOR_ELT(output, 14, object->out_P_t);
        //if (INTEGER(EcoDcf)[0]==0) SET_VECTOR_ELT(output, 15, object->out_Eco); else SET_VECTOR_ELT(output, 15, object->out_EcoDCF);
        SET_VECTOR_ELT(output, 15, object->out_EcoDCF);
        PROTECT(out_Foth = allocVector(VECSXP, object->nbE));
        setAttrib(out_Foth, R_NamesSymbol, object->sppList);
        if (object->nbE>0) {for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 44));}
        SET_VECTOR_ELT(output, 16, out_Foth);
        SET_VECTOR_ELT(output, 17, object->mu_nbds);
        SET_VECTOR_ELT(output, 18, object->mu_nbv);
        SET_VECTOR_ELT(output, 19, object->out_effort);
        SET_VECTOR_ELT(output, 20, object->out_Fr_fmi);
        SET_VECTOR_ELT(output, 21, VECTOR_ELT(object->fVar, 29));
        SET_VECTOR_ELT(output, 22, object->out_PQuot_et);
        SET_VECTOR_ELT(output, 23, object->out_typeGest);
        SET_VECTOR_ELT(output, 24, object->out_Ystat);
        SET_VECTOR_ELT(output, 25, object->out_Lstat);
        SET_VECTOR_ELT(output, 26, object->out_Dstat);
        SET_VECTOR_ELT(output, 27, object->out_Pstat);


        //peut-�tre � retirer ensuite --> export des indicateurs par morph/trimestre -------------------------------

        SET_VECTOR_ELT(output, 28, object->out_F_fmi_S1M1);
        SET_VECTOR_ELT(output, 29, object->out_F_fmi_S1M2);
        SET_VECTOR_ELT(output, 30, object->out_F_fmi_S1M3);
        SET_VECTOR_ELT(output, 31, object->out_F_fmi_S1M4);
        SET_VECTOR_ELT(output, 32, object->out_F_fmi_S2M1);
        SET_VECTOR_ELT(output, 33, object->out_F_fmi_S2M2);
        SET_VECTOR_ELT(output, 34, object->out_F_fmi_S2M3);
        SET_VECTOR_ELT(output, 35, object->out_F_fmi_S2M4);
        SET_VECTOR_ELT(output, 36, object->out_F_fmi_S3M1);
        SET_VECTOR_ELT(output, 37, object->out_F_fmi_S3M2);
        SET_VECTOR_ELT(output, 38, object->out_F_fmi_S3M3);
        SET_VECTOR_ELT(output, 39, object->out_F_fmi_S3M4);
        SET_VECTOR_ELT(output, 40, object->out_F_fmi_S4M1);
        SET_VECTOR_ELT(output, 41, object->out_F_fmi_S4M2);
        SET_VECTOR_ELT(output, 42, object->out_F_fmi_S4M3);
        SET_VECTOR_ELT(output, 43, object->out_F_fmi_S4M4);
        SET_VECTOR_ELT(output, 44, object->out_Fr_fmi_S1M1);
        SET_VECTOR_ELT(output, 45, object->out_Fr_fmi_S1M2);
        SET_VECTOR_ELT(output, 46, object->out_Fr_fmi_S1M3);
        SET_VECTOR_ELT(output, 47, object->out_Fr_fmi_S1M4);
        SET_VECTOR_ELT(output, 48, object->out_Fr_fmi_S2M1);
        SET_VECTOR_ELT(output, 49, object->out_Fr_fmi_S2M2);
        SET_VECTOR_ELT(output, 50, object->out_Fr_fmi_S2M3);
        SET_VECTOR_ELT(output, 51, object->out_Fr_fmi_S2M4);
        SET_VECTOR_ELT(output, 52, object->out_Fr_fmi_S3M1);
        SET_VECTOR_ELT(output, 53, object->out_Fr_fmi_S3M2);
        SET_VECTOR_ELT(output, 54, object->out_Fr_fmi_S3M3);
        SET_VECTOR_ELT(output, 55, object->out_Fr_fmi_S3M4);
        SET_VECTOR_ELT(output, 56, object->out_Fr_fmi_S4M1);
        SET_VECTOR_ELT(output, 57, object->out_Fr_fmi_S4M2);
        SET_VECTOR_ELT(output, 58, object->out_Fr_fmi_S4M3);
        SET_VECTOR_ELT(output, 59, object->out_Fr_fmi_S4M4);
        SET_VECTOR_ELT(output, 60, object->out_Z_eit_S1M1);
        SET_VECTOR_ELT(output, 61, object->out_Z_eit_S1M2);
        SET_VECTOR_ELT(output, 62, object->out_Z_eit_S1M3);
        SET_VECTOR_ELT(output, 63, object->out_Z_eit_S1M4);
        SET_VECTOR_ELT(output, 64, object->out_Z_eit_S2M1);
        SET_VECTOR_ELT(output, 65, object->out_Z_eit_S2M2);
        SET_VECTOR_ELT(output, 66, object->out_Z_eit_S2M3);
        SET_VECTOR_ELT(output, 67, object->out_Z_eit_S2M4);
        SET_VECTOR_ELT(output, 68, object->out_Z_eit_S3M1);
        SET_VECTOR_ELT(output, 69, object->out_Z_eit_S3M2);
        SET_VECTOR_ELT(output, 70, object->out_Z_eit_S3M3);
        SET_VECTOR_ELT(output, 71, object->out_Z_eit_S3M4);
        SET_VECTOR_ELT(output, 72, object->out_Z_eit_S4M1);
        SET_VECTOR_ELT(output, 73, object->out_Z_eit_S4M2);
        SET_VECTOR_ELT(output, 74, object->out_Z_eit_S4M3);
        SET_VECTOR_ELT(output, 75, object->out_Z_eit_S4M4);
        SET_VECTOR_ELT(output, 76, object->out_N_eit_S1M1);
        SET_VECTOR_ELT(output, 77, object->out_N_eit_S1M2);
        SET_VECTOR_ELT(output, 78, object->out_N_eit_S1M3);
        SET_VECTOR_ELT(output, 79, object->out_N_eit_S1M4);
        SET_VECTOR_ELT(output, 80, object->out_N_eit_S2M1);
        SET_VECTOR_ELT(output, 81, object->out_N_eit_S2M2);
        SET_VECTOR_ELT(output, 82, object->out_N_eit_S2M3);
        SET_VECTOR_ELT(output, 83, object->out_N_eit_S2M4);
        SET_VECTOR_ELT(output, 84, object->out_N_eit_S3M1);
        SET_VECTOR_ELT(output, 85, object->out_N_eit_S3M2);
        SET_VECTOR_ELT(output, 86, object->out_N_eit_S3M3);
        SET_VECTOR_ELT(output, 87, object->out_N_eit_S3M4);
        SET_VECTOR_ELT(output, 88, object->out_N_eit_S4M1);
        SET_VECTOR_ELT(output, 89, object->out_N_eit_S4M2);
        SET_VECTOR_ELT(output, 90, object->out_N_eit_S4M3);
        SET_VECTOR_ELT(output, 91, object->out_N_eit_S4M4);

        //01/04/2015 export des indicateurs suppl�mentaires : out_Ytot_fm, out_DD_efmi, out_DD_efmc, out_LD_efmi, out_LD_efmc,
        //out_statDD_efm, out_statLD_efm, out_statLDst_efm, out_statLDor_efm-------------------------------

        SET_VECTOR_ELT(output, 92, object->out_Ytot_fm);
        SET_VECTOR_ELT(output, 93, object->out_DD_efmi);
        SET_VECTOR_ELT(output, 94, object->out_DD_efmc);
        SET_VECTOR_ELT(output, 95, object->out_LD_efmi);
        SET_VECTOR_ELT(output, 96, object->out_LD_efmc);
        SET_VECTOR_ELT(output, 97, object->out_statDD_efm);
        SET_VECTOR_ELT(output, 98, object->out_statLD_efm);
        SET_VECTOR_ELT(output, 99, object->out_statLDst_efm);
        SET_VECTOR_ELT(output, 100, object->out_statLDor_efm);
        SET_VECTOR_ELT(output, 101, object->out_oqD_eft);
        SET_VECTOR_ELT(output, 102, object->out_oqD_et);
        SET_VECTOR_ELT(output, 103, object->out_oqDstat);
        SET_VECTOR_ELT(output, 104, object->reconcilSPP);
        SET_VECTOR_ELT(output, 105, object->TAC);
        SET_VECTOR_ELT(output, 106, object->TACbyF);

        //Ajout indicateurs sex-based
        PROTECT(out_Foth_G1 = allocVector(VECSXP, object->nbE));
        setAttrib(out_Foth_G1, R_NamesSymbol, object->sppList);
        if (object->nbE>0) {for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth_G1, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 224));}
        SET_VECTOR_ELT(output, 107, out_Foth_G1);
        PROTECT(out_Foth_G2 = allocVector(VECSXP, object->nbE));
        setAttrib(out_Foth_G2, R_NamesSymbol, object->sppList);
        if (object->nbE>0) {for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth_G2, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 225));}
        SET_VECTOR_ELT(output, 108, out_Foth_G2);
        SET_VECTOR_ELT(output, 109, object->out_F_fmi_G1);
        SET_VECTOR_ELT(output, 110, object->out_F_fmi_G2);
        SET_VECTOR_ELT(output, 111, object->out_Z_eit_G1);
        SET_VECTOR_ELT(output, 112, object->out_Z_eit_G2);
        SET_VECTOR_ELT(output, 113, object->out_N_eit_G1);
        SET_VECTOR_ELT(output, 114, object->out_N_eit_G2);
        SET_VECTOR_ELT(output, 115, object->out_Fr_fmi_G1);
        SET_VECTOR_ELT(output, 116, object->out_Fr_fmi_G2);
        SET_VECTOR_ELT(output, 117, object->out_C_efmit_G1);
        SET_VECTOR_ELT(output, 118, object->out_C_efmit_G2);
        SET_VECTOR_ELT(output, 119, object->out_C_eit_G1);
        SET_VECTOR_ELT(output, 120, object->out_C_eit_G2);
        SET_VECTOR_ELT(output, 121, object->out_L_et);
        SET_VECTOR_ELT(output, 122, object->out_L_pt);
        SET_VECTOR_ELT(output, 123, object->out_QuotaTrade_fe);
        SET_VECTOR_ELT(output, 124, object->out_allocEff_fm);
        SET_VECTOR_ELT(output, 125, object->out_PQuot_temp);
        SET_VECTOR_ELT(output, 126, object->out_diffLQ);
        SET_VECTOR_ELT(output, 127, object->intermGoFish);

//PrintValue(object->reconcilSPP);
        //----------------------------------------------------------------------------------------------------------

        //on nomme les �l�ments de output
        const char *namesOut[128] = {"F","Z","Fbar","N","B","SSB","C","Ctot","Y","Ytot","D","Li","Lc","Ltot","P","E","Fothi","mu_nbds","mu_nbv","Eff","Fr","GVLoths_f","PQuot","typeGest","Ystat","Lstat","Dstat","Pstat",//};
                                    "F_S1M1","F_S1M2","F_S1M3","F_S1M4","F_S2M1","F_S2M2","F_S2M3","F_S2M4","F_S3M1","F_S3M2","F_S3M3","F_S3M4","F_S4M1","F_S4M2","F_S4M3","F_S4M4",
                                    "Fr_S1M1","Fr_S1M2","Fr_S1M3","Fr_S1M4","Fr_S2M1","Fr_S2M2","Fr_S2M3","Fr_S2M4","Fr_S3M1","Fr_S3M2","Fr_S3M3","Fr_S3M4","Fr_S4M1","Fr_S4M2","Fr_S4M3","Fr_S4M4",
                                    "Z_S1M1","Z_S1M2","Z_S1M3","Z_S1M4","Z_S2M1","Z_S2M2","Z_S2M3","Z_S2M4","Z_S3M1","Z_S3M2","Z_S3M3","Z_S3M4","Z_S4M1","Z_S4M2","Z_S4M3","Z_S4M4",
                                    "N_S1M1","N_S1M2","N_S1M3","N_S1M4","N_S2M1","N_S2M2","N_S2M3","N_S2M4","N_S3M1","N_S3M2","N_S3M3","N_S3M4","N_S4M1","N_S4M2","N_S4M3","N_S4M4",
                                    "YTOT_fm", "DD_efmi", "DD_efmc", "LD_efmi", "LD_efmc", "statDD_efm", "statLD_efm", "statLDst_efm", "statLDor_efm",
                                    "oqD_ef","oqD_e","oqDstat_ef","reconcilSPP","TACtot","TACbyF","Fothi_G1","Fothi_G2","F_G1","F_G2","Z_G1","Z_G2","N_G1","N_G2","Fr_G1","Fr_G2","C_G1","C_G2","Ctot_G1","Ctot_G2","L_et","L_pt","TradedQ_f","allocEff_fm",
                                    "PQuot_conv","diffLQ_conv","GoFish"};
        PROTECT(out_names = allocVector(STRSXP, 128));

        for(int ct = 0; ct < 128; ct++) SET_STRING_ELT(out_names, ct, mkChar(namesOut[ct]));

        setAttrib(output, R_NamesSymbol, out_names);

        UNPROTECT(5);
        if(VERBOSE){Rprintf("---- Exit C++ ----\n");}
        return(output);
        delete object;


    } else {

        //on n'oublie pas d'activer les parties stochastiques pour que �a ait un sens

        //on commence par cr�er l'objet qui va accueillir la donn�e  (3 outputs pour l'instant : biomasse, SSB, captures --> � d�velopper selon les besoins)
        SEXP output, out_names, out_Foth, emptyObj;
        PROTECT(output = allocVector(VECSXP, 44)); //36
        SEXP eBoot;

        for (int ind = 0 ; ind < 44 ; ind++) { //36

            PROTECT(eBoot = allocVector(VECSXP, INTEGER(nbBoot)[0]));
            SET_VECTOR_ELT(output, ind, eBoot);

        }

        //on commence le bootstrap
        if(VERBOSE){Rprintf("Init bootstrap \n");}
        BioEcoPar *object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                    RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                    GestInd, mOth, bounds, TAC, FBAR, othSpSup, effSup, GestParam, EcoDcf,
                                    EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, Ftarg, W_Ftarg, MeanRec_Ftarg, parBHV, parQEX, tacCTRL, stochPrice, updateE, parOQD, VERBOSE);

        for (int it = 0 ; it < INTEGER(nbBoot)[0] ; it++) {
            if(VERBOSE && (it % 50) == 0){Rprintf("boot :");}
            if (it>0) object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                    RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                    GestInd, mOth, bounds, TAC, FBAR, othSpSup, effSup, GestParam, EcoDcf,
                                    EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, Ftarg, W_Ftarg, MeanRec_Ftarg, parBHV, parQEX, tacCTRL, stochPrice, updateE, parOQD, VERBOSE);

            //objet vide pour garder la structuration malgr� la non-s�lection de la variable en question
            PROTECT(emptyObj = allocVector(VECSXP, object->nbE));
            setAttrib(emptyObj, R_NamesSymbol, object->sppList);

            if (isCharIn(bootVar, "B")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 0), it, object->out_B_et);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 0), it, emptyObj);
            }

            if (isCharIn(bootVar, "SSB")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 1), it, object->out_SSB_et);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 1), it, emptyObj);
            }

            if (isCharIn(bootVar, "Ctot")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 2), it, object->out_C_eit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 2), it, emptyObj);
            }

            if (isCharIn(bootVar, "Ytot")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 3), it, object->out_Y_eit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 3), it, emptyObj);
            }

            if (isCharIn(bootVar, "Yfmi")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 4), it, object->out_Y_efmit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 4), it, emptyObj);
            }

            if (isCharIn(bootVar, "Ffmi")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 5), it, object->out_F_fmi);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 5), it, emptyObj);
            }

            if (isCharIn(bootVar, "Zeit")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 6), it, object->out_Z_eit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 6), it, emptyObj);
            }

            if (isCharIn(bootVar, "Fbar")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 7), it, object->out_Fbar_et);
                } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 7), it, emptyObj);
            }

            PROTECT(out_Foth = allocVector(VECSXP, object->nbE));
            setAttrib(out_Foth, R_NamesSymbol, object->sppList);
            if (object->nbE>0) {for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 44));}
            if (isCharIn(bootVar, "Foth")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 8), it, out_Foth);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 8), it, emptyObj);
            }

            if (isCharIn(bootVar, "mu_nbds")) SET_VECTOR_ELT(VECTOR_ELT(output, 9), it, object->mu_nbds);
            if (isCharIn(bootVar, "mu_nbv")) SET_VECTOR_ELT(VECTOR_ELT(output, 10), it, object->mu_nbv);

            if (isCharIn(bootVar, "N")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 11), it, object->out_N_eit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 11), it, emptyObj);
            }

            if (isCharIn(bootVar, "Eff")) SET_VECTOR_ELT(VECTOR_ELT(output, 12), it, object->out_effort);

            if (INTEGER(EcoDcf)[0]==0) {
                if (isCharIn(bootVar, "GVL_fme")) {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 13), it, VECTOR_ELT(object->out_Eco,1));
                } else {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 13), it, emptyObj);
                }

                if (isCharIn(bootVar, "GVLtot_fm")) SET_VECTOR_ELT(VECTOR_ELT(output, 14), it, VECTOR_ELT(object->out_Eco,2));
                if (isCharIn(bootVar, "GVLav_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 15), it, VECTOR_ELT(object->out_Eco,5));
                if (isCharIn(bootVar, "rtbs_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 16), it, VECTOR_ELT(object->out_Eco,11));
                if (isCharIn(bootVar, "gp_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 17), it, VECTOR_ELT(object->out_Eco,27));
                if (isCharIn(bootVar, "ps_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 18), it, VECTOR_ELT(object->out_Eco,29));
                if (isCharIn(bootVar, "gcf_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 19), it, VECTOR_ELT(object->out_Eco,25));
                if (isCharIn(bootVar, "gva_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 20), it, VECTOR_ELT(object->out_Eco,20));
                if (isCharIn(bootVar, "cs_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 21), it, VECTOR_ELT(object->out_Eco,18));
                if (isCharIn(bootVar, "sts_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 22), it, VECTOR_ELT(object->out_Eco,30));
                if (isCharIn(bootVar, "rtbsAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 23), it, VECTOR_ELT(object->out_Eco,54));
                if (isCharIn(bootVar, "csAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 24), it, VECTOR_ELT(object->out_Eco,55));
                if (isCharIn(bootVar, "gvaAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 25), it, VECTOR_ELT(object->out_Eco,56));
                if (isCharIn(bootVar, "gcfAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 26), it, VECTOR_ELT(object->out_Eco,57));
                if (isCharIn(bootVar, "psAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 27), it, VECTOR_ELT(object->out_Eco,58));
                if (isCharIn(bootVar, "stsAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 28), it, VECTOR_ELT(object->out_Eco,59));
                if (isCharIn(bootVar, "ccwCr_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 29), it, VECTOR_ELT(object->out_Eco,22));
                if (isCharIn(bootVar, "GVLtot_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 30), it, VECTOR_ELT(object->out_Eco,4));
                if (isCharIn(bootVar, "wagen_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 31), it, VECTOR_ELT(object->out_Eco,24));
                if (isCharIn(bootVar, "vcst_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 36), it, VECTOR_ELT(object->out_Eco,9));
                if (isCharIn(bootVar, "vcst_fm")) SET_VECTOR_ELT(VECTOR_ELT(output, 37), it, VECTOR_ELT(object->out_Eco,8));

                if (isCharIn(bootVar, "StatGVL_fme")) {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 43), it, VECTOR_ELT(object->out_Eco,60));
                } else {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 43), it, emptyObj);
                }



            } else {

                if (isCharIn(bootVar, "GVL_fme")) {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 13), it, VECTOR_ELT(object->out_EcoDCF,0));
                } else {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 13), it, emptyObj);
                }

                if (isCharIn(bootVar, "GVLtot_fm")) SET_VECTOR_ELT(VECTOR_ELT(output, 14), it, VECTOR_ELT(object->out_EcoDCF,1));
                if (isCharIn(bootVar, "GVLav_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 15), it, VECTOR_ELT(object->out_EcoDCF,4));
                if (isCharIn(bootVar, "rtbs_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 16), it, VECTOR_ELT(object->out_EcoDCF,7));
                if (isCharIn(bootVar, "gp_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 17), it, VECTOR_ELT(object->out_EcoDCF,22));
                if (isCharIn(bootVar, "ps_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 18), it, VECTOR_ELT(object->out_EcoDCF,23));
                if (isCharIn(bootVar, "gcf_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 19), it, VECTOR_ELT(object->out_EcoDCF,20));
                if (isCharIn(bootVar, "gva_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 20), it, VECTOR_ELT(object->out_EcoDCF,15));
                if (isCharIn(bootVar, "cs_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 21), it, VECTOR_ELT(object->out_EcoDCF,14));
                if (isCharIn(bootVar, "sts_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 22), it, VECTOR_ELT(object->out_EcoDCF,24));
                if (isCharIn(bootVar, "rtbsAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 23), it, VECTOR_ELT(object->out_EcoDCF,38));
                if (isCharIn(bootVar, "csAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 24), it, VECTOR_ELT(object->out_EcoDCF,39));
                if (isCharIn(bootVar, "gvaAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 25), it, VECTOR_ELT(object->out_EcoDCF,40));
                if (isCharIn(bootVar, "gcfAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 26), it, VECTOR_ELT(object->out_EcoDCF,41));
                if (isCharIn(bootVar, "psAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 27), it, VECTOR_ELT(object->out_EcoDCF,42));
                if (isCharIn(bootVar, "stsAct_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 28), it, VECTOR_ELT(object->out_EcoDCF,43));
                if (isCharIn(bootVar, "ccwCr_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 29), it, VECTOR_ELT(object->out_EcoDCF,17));
                if (isCharIn(bootVar, "GVLtot_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 30), it, VECTOR_ELT(object->out_EcoDCF,3));
                if (isCharIn(bootVar, "wagen_f")) SET_VECTOR_ELT(VECTOR_ELT(output, 31), it, VECTOR_ELT(object->out_EcoDCF,19));

                if (isCharIn(bootVar, "StatGVL_fme")) {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 43), it, VECTOR_ELT(object->out_EcoDCF,45));
                } else {
                    SET_VECTOR_ELT(VECTOR_ELT(output, 43), it, emptyObj);
                }


            }

            if (isCharIn(bootVar, "L_efmit")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 32), it, object->out_L_efmit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 32), it, emptyObj);
            }

            if (isCharIn(bootVar, "D_efmit")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 33), it, object->out_D_efmit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 33), it, emptyObj);
            }

            if (isCharIn(bootVar, "Fr_fmi")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 34), it, object->out_Fr_fmi);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 34), it, emptyObj);
            }

            if (isCharIn(bootVar, "C_efmit")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 35), it, object->out_C_efmit);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 35), it, emptyObj);
            }

            if (isCharIn(bootVar, "P")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 38), it, object->out_P_t);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 38), it, emptyObj);
            }

            if (isCharIn(bootVar, "Ystat")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 39), it, object->out_Ystat);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 39), it, emptyObj);
            }

           if (isCharIn(bootVar, "Lstat")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 40), it, object->out_Lstat);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 40), it, emptyObj);
            }

           if (isCharIn(bootVar, "Dstat")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 41), it, object->out_Dstat);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 41), it, emptyObj);
            }

           if (isCharIn(bootVar, "Pstat")) {
                SET_VECTOR_ELT(VECTOR_ELT(output, 42), it, object->out_Pstat);
            } else {
                SET_VECTOR_ELT(VECTOR_ELT(output, 42), it, emptyObj);
            }


            UNPROTECT(2);
        }

        //on nomme les �l�ments de output
        const char *namesOut[44] = {"B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth","mu_nbds","mu_nbv","N","Eff",
                                    "GVL_fme","GVLtot_fm","GVLav_f","rtbs_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
                                    "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f","GVLtot_f","wagen_f","L_efmit","D_efmit","Fr_fmi","C_efmit","vcst_f","vcst_fm","P",
                                    "Ystat","Lstat","Dstat","Pstat","StatGVL_fme"};

        PROTECT(out_names = allocVector(STRSXP, 44));

        for(int ct = 0; ct < 44; ct++) SET_STRING_ELT(out_names, ct, mkChar(namesOut[ct]));

        setAttrib(output, R_NamesSymbol, out_names);

        UNPROTECT(2+44);
        if(VERBOSE){Rprintf("---- Exit C++ ----\n");}
        return(output);
        delete object;

    }
}
}


int main()
{

    return 0;
}