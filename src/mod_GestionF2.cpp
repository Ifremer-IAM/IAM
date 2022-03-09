#include <Rdefines.h>
#include <Rmath.h>
//#include <Rcpp.h>

#include "BioEcoPar.h" // Class is defined in this file.
#include "array_fcts.h" // AggregObj function
// #include "Modules.h" // Contient tout les modules

//using namespace Rcpp;
// using namespace std;

//------------------------------------------
// Module 'Gestion F2'
//------------------------------------------

extern "C" {
  void BioEcoPar::abv_GestionF2(int ind_t, SEXP updateE, SEXP tacCTRL, SEXP FList, int VERBOSE)
  {
    if(VERBOSE){Rprintf("set t -1 - ");}
    int DELAY = INTEGER(updateE)[0];

    SPPstatOPT = INTEGER(getListElement(tacCTRL, "SPPstatOPT"));
    SPPspictOPT = INTEGER(getListElement(tacCTRL, "SPPspictOPT"));
    SPPdynOPT = INTEGER(getListElement(tacCTRL, "SPPdynOPT"));
    N_SPPstatOPT = length(getListElement(tacCTRL, "SPPstatOPT"));
    N_SPPspictOPT = length(getListElement(tacCTRL, "SPPspictOPT"));
    N_SPPdynOPT = length(getListElement(tacCTRL, "SPPdynOPT"));

    if ((delay<=ind_t) & (gestInd==1) & (DELAY>0)) { //DELAY = 1 -> on remet l'effort au niveau de l'instant initial
      if(VERBOSE){Rprintf("initial time ");}
      //on remet au niveau de l'instant pr�c�dent la mise en action du module Gestion

      double *nbdsFM3 = REAL(getListElement(FList, "effort1_f_m"));
      double *nbdsF3 = REAL(getListElement(FList, "effort1_f"));
      double *nbTripFM3 = REAL(getListElement(FList, "nbTrip_f_m"));
      double *nbTripF3 = REAL(getListElement(FList, "nbTrip_f"));
      double *nbvFM3 = REAL(getListElement(FList, "nbv_f_m"));
      double *nbvF3 = REAL(getListElement(FList, "nbv_f"));

      if (DELAY>delay) DELAY=delay; // TODO : should be triggered in R before launching module. why 2 var btw ?

      for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
        // assign pointeurs for output empty matrix
        if (var==1) {nbdsF3[ind_f] = REAL(NBDSF)[ind_f + nbF*(DELAY-1)];
          nbTripF3[ind_f] = REAL(NBDSF)[ind_f + nbF*(DELAY-1)];}
        if (var==2) nbvF3[ind_f] = REAL(NBVF)[ind_f + nbF*(DELAY-1)];

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
          if (var==1) {nbdsFM3[ind_f+nbF*ind_m] = REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];
            nbTripFM3[ind_f+nbF*ind_m] = REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];}
          if (var==2) nbvFM3[ind_f+nbF*ind_m] = REAL(NBVFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];
        }
      }

      if (N_SPPspictOPT>0) { //si esp�ce dynamique SPICT // TODO : if N == 0, loop will not trigger, useless.
        if(VERBOSE){Rprintf("spict sp ");}

        for (int i = 0; i < N_SPPspictOPT; i++){
          int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPspictOPT[i] - 1))), "modI")); //doit normalement etre egal a 1
          double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPspictOPT[i]), 44));
          for (int ag = 0; ag < nbi; ag++){
            Fothi2[ag + ind_t*nbi] = Fothi2[ag + (DELAY-1)*nbi];
          }
        }
      }

      if (N_SPPdynOPT>0) { //si esp�ce dynamique XSA ou SS3 // TODO : idem, 0 will not trigger loop
        if(VERBOSE){Rprintf("dyna sp ");}

        for (int i = 0; i < N_SPPdynOPT; i++){
          int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPdynOPT[i] - 1))), "modI"));
          int eTmp = SPPdynOPT[i] - 1;

          if (Qvec[eTmp]==0) {
            if(VERBOSE){Rprintf("xsa ");}
            double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 44));
            for (int ag = 0; ag < nbi; ag++){
              Fothi2[ag + ind_t*nbi] = Fothi2[ag + (DELAY-1)*nbi];
            }

          } else {
            if(VERBOSE){Rprintf("ss3 ");}
            double *Fothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 116));
            double *Fothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 117));
            double *Fothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 118));
            double *Fothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 119));
            double *Fothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 120));
            double *Fothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 121));
            double *Fothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 122));
            double *Fothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 123));
            double *Fothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 124));
            double *Fothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 125));
            double *Fothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 126));
            double *Fothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 127));
            double *Fothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 128));
            double *Fothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 129));
            double *Fothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 130));
            double *Fothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 131));

            for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + ind_t*nbi] = Fothi2_S1M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + ind_t*nbi] = Fothi2_S1M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + ind_t*nbi] = Fothi2_S1M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + ind_t*nbi] = Fothi2_S1M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + ind_t*nbi] = Fothi2_S2M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + ind_t*nbi] = Fothi2_S2M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + ind_t*nbi] = Fothi2_S2M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + ind_t*nbi] = Fothi2_S2M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + ind_t*nbi] = Fothi2_S3M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + ind_t*nbi] = Fothi2_S3M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + ind_t*nbi] = Fothi2_S3M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + ind_t*nbi] = Fothi2_S3M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + ind_t*nbi] = Fothi2_S4M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + ind_t*nbi] = Fothi2_S4M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + ind_t*nbi] = Fothi2_S4M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + ind_t*nbi] = Fothi2_S4M4[ag + (DELAY-1)*nbi];


            double *FRWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 176));
            double *FRWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 177));
            double *FRWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 178));
            double *FRWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 179));
            double *FRWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 180));
            double *FRWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 181));
            double *FRWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 182));
            double *FRWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 183));
            double *FRWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 184));
            double *FRWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 185));
            double *FRWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 186));
            double *FRWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 187));
            double *FRWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 188));
            double *FRWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 189));
            double *FRWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 190));
            double *FRWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 191));

            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + ind_t*nbi] = FRWTothi2_S1M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + ind_t*nbi] = FRWTothi2_S1M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + ind_t*nbi] = FRWTothi2_S1M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + ind_t*nbi] = FRWTothi2_S1M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + ind_t*nbi] = FRWTothi2_S2M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + ind_t*nbi] = FRWTothi2_S2M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + ind_t*nbi] = FRWTothi2_S2M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + ind_t*nbi] = FRWTothi2_S2M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + ind_t*nbi] = FRWTothi2_S3M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + ind_t*nbi] = FRWTothi2_S3M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + ind_t*nbi] = FRWTothi2_S3M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + ind_t*nbi] = FRWTothi2_S3M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + ind_t*nbi] = FRWTothi2_S4M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + ind_t*nbi] = FRWTothi2_S4M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + ind_t*nbi] = FRWTothi2_S4M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + ind_t*nbi] = FRWTothi2_S4M4[ag + (DELAY-1)*nbi];

            double *FDWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 208));
            double *FDWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 209));
            double *FDWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 210));
            double *FDWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 211));
            double *FDWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 212));
            double *FDWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 213));
            double *FDWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 214));
            double *FDWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 215));
            double *FDWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 216));
            double *FDWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 217));
            double *FDWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 218));
            double *FDWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 219));
            double *FDWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 220));
            double *FDWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 221));
            double *FDWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 222));
            double *FDWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 223));

            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + ind_t*nbi] = FDWTothi2_S1M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + ind_t*nbi] = FDWTothi2_S1M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + ind_t*nbi] = FDWTothi2_S1M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + ind_t*nbi] = FDWTothi2_S1M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + ind_t*nbi] = FDWTothi2_S2M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + ind_t*nbi] = FDWTothi2_S2M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + ind_t*nbi] = FDWTothi2_S2M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + ind_t*nbi] = FDWTothi2_S2M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + ind_t*nbi] = FDWTothi2_S3M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + ind_t*nbi] = FDWTothi2_S3M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + ind_t*nbi] = FDWTothi2_S3M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + ind_t*nbi] = FDWTothi2_S3M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + ind_t*nbi] = FDWTothi2_S4M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + ind_t*nbi] = FDWTothi2_S4M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + ind_t*nbi] = FDWTothi2_S4M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + ind_t*nbi] = FDWTothi2_S4M4[ag + (DELAY-1)*nbi];

          }
        }}

    }

    if(VERBOSE){Rprintf("launch module :");}
    GestionF2(ind_t, tacCTRL, VERBOSE);

  }
}

extern "C" {

  void BioEcoPar::GestionF2(int ind_t, SEXP tacCTRL, int VERBOSE)
  {
    if(VERBOSE){Rprintf("\n   Gest: ");}

    IND_T = ind_t; // TODO: WHY ??!!

    // int ITtot = maxIter; // TODO : why rename this variable ?
    int ITtot = INTEGER(getListElement(tacCTRL, "maxIter"))[0];
    double diffZmax = REAL(getListElement(tacCTRL, "diffZmax"))[0];
    double lambda = REAL(getListElement(tacCTRL, "lambda"))[0];
    bool goon = true;

    double *g_effSup = REAL(effSupMat);
    double *mpond_fm = REAL(out_allocEff_fm);

    SEXP listTempP, nDimFM, nDimF, nDim, copyEffort;

    // create DimCst elements
    PROTECT(nDimFM = allocVector(INTSXP,4)); // nbf * nbMe * nbT
    int *ndFM = INTEGER(nDimFM); ndFM[0] = nbF; ndFM[1] = nbMe; ndFM[2] = 0; ndFM[3] = nbT;
    PROTECT(nDimF = allocVector(INTSXP,4)); // nbf * nbT
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4)); // nbT
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;

    PROTECT(listTempP = duplicate(list));
    // PROTECT(eVarCopy = duplicate(eVar));
    //new
    double *g_effort1FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m"));
    double *g_effort1F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f"));
    double *g_nbTripFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f_m"));
    double *g_nbTripF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f"));
    double *g_nbvFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f_m"));
    double *g_nbvF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f"));
    double *g_effort2FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f_m"));
    double *g_effort2F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f"));


    for (int ind_f = 0 ; ind_f<nbF ; ind_f++) { //initialisation de l'effort de base
      g_effort1F[ind_f] = REAL(NBDSF)[ind_f];
      g_nbTripF[ind_f] = g_effort1F[ind_f];

      for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
        g_effort1FM[ind_f + nbF*ind_m] = (g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
          (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);
        g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m];
      }
    }

    PROTECT(copyEffort = duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f")));
    double *g_effort1F_copy = REAL(copyEffort); //duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f")));
    //new
    if(VERBOSE){Rprintf("Bio modules ");}
    Mortalite(listTempP, IND_T, eVar);
    DynamicPop(listTempP, IND_T, eVar,true);
    CatchDL(listTempP, IND_T, eVar, 0);
    if(VERBOSE){Rprintf("done\n");}

    //on initialise en remplissant Einterm_fm avec les valeurs maximales autorisees (fonction de effSupMat et la ponderation multiplicative U (mpond_fm)) : Eq A1
    //on rappelle que tout cela fonctionne uniquement si on est sur une base individuelle (dans le cas contraire, il faudrait descendre au niveau metier et ponderer par les nbv par metier

    for (int indF = 0 ; indF < nbF ; indF++) {
      if (ISNA(g_effSup[indF + nbF*ind_t])) {
        g_effSup[indF + nbF*ind_t] = 300; // Warning if is NA, max is 300 NBDS !
      }
      double denom_st1 = 0.0; double alpha_f_st1 = 0.0;
      for (int indM = 0 ; indM<nbMe ; indM++) {
        if(!ISNA(mpond_fm[indF + nbF*indM + nbF*nbMe*IND_T])){
          // TODO : denom = NA car sum(n,NA) == NA !!!
          denom_st1 = denom_st1 + (g_effort1F[indF] * g_effort2F[indF] * mpond_fm[indF + nbF*indM + nbF*nbMe*IND_T]); //c'est cense valoir eff_f normalement si la somme des facteurs de pond. vaut 1
        }
      }
      if (denom_st1>0) alpha_f_st1 = g_effSup[indF + nbF*ind_t] / denom_st1;
      for (int indM = 0 ; indM<nbMe ; indM++) {
        EffsupTMP_fm[1 + indF + nbF*indM] = alpha_f_st1 * g_effort1F[indF] * g_nbvF[indF] * g_effort2F[indF] * mpond_fm[indF + nbF*indM + nbF*nbMe*IND_T];
      }
    }

    for (int indF = 0 ; indF < nbF ; indF++){
      for (int indM = 0 ; indM<nbMe ; indM++){
        Einterm_fm[1 + indF + nbF*indM] = -1.0; //va servir a identifier a la fin les cellules non reconciliees, qui devront etre egales a l'effort courant!!!
      }
    }

    //� ce stade, on a initialise Einterm_fm
    //on commence l'ajustement avec les especes statiques indexees dans SPPstatOPT (s'il y en a)
    if (N_SPPstatOPT>0) { // TODO : loop is not run if 0
      if(VERBOSE){Rprintf("Static");}

      for (int ind = 0 ; ind < N_SPPstatOPT ; ind++) { //boucle sur ces especes statiques

        SEXP gg1=R_NilValue, gg2=R_NilValue;
        PROTECT(gg1=VECTOR_ELT(out_Lstat, SPPstatOPT[ind]-1));
        PROTECT(gg2=VECTOR_ELT(out_statLD_efm, SPPstatOPT[ind]-1));

        double *TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppListStat,SPPstatOPT[ind]-1))));
        double *totFM = REAL(aggregObj(gg1,nDimFM));
        double *totFM2 = REAL(aggregObj(gg2,nDimFM));

        for (int ind_f = 0 ; ind_f<nbF ; ind_f++) {

          double denom_st2 = 0.0; double alpha_f_st2 = 0.0;
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {

            if (g_effort1FM[ind_f + nbF*ind_m]>0)
            // TODO : denom = NA car sum(n,NA) == NA !!!
              denom_st2 = denom_st2 + mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * (totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
                (g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);

          }

          if ((denom_st2>0) & (g_effort1F[ind_f]>0))
            alpha_f_st2 = finite(TAC_byFleet[ind_f + nbF*IND_T] / (g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * denom_st2));

          //r�conciliation
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
            // TODO : denom = NA car sum(n,NA) == NA !!!
            double valTest = alpha_f_st2 * g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
            if ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0) {
              if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
                Einterm_fm[1 + ind_f + nbF*ind_m] = EffsupTMP_fm[1 + ind_f + nbF*ind_m]; //on �limine le -1
                SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
              }
              if (valTest < Einterm_fm[1 + ind_f + nbF*ind_m] ) {//r�conciliation seulement si la capture associ�e est non nulle
                Einterm_fm[1 + ind_f + nbF*ind_m] = valTest;
                SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppListStat,SPPstatOPT[ind]-1));
                SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppListStat,SPPstatOPT[ind]-1));
              }
            }
          }
          //Rprintf("T %i F %i \n",ind_t,ind_f);
        }
        UNPROTECT(2);
      }
      if(VERBOSE){Rprintf(". ");}
    }


    //on poursuit avec l'ajustement des especes dynamiques sous SPiCT indexees dans SPPspictOPT (s'il y en a)
    if (N_SPPspictOPT>0) { // TODO : loop is not run if 0
      if(VERBOSE){Rprintf("Spict");}

      for (int ind = 0 ; ind < N_SPPspictOPT ; ind++) { //boucle sur ces esp�ces Spict

        SEXP gg1=R_NilValue, gg2=R_NilValue /*, gg3=R_NilValue*/;
        PROTECT(gg1=VECTOR_ELT(out_L_efmit, SPPspictOPT[ind]-1));
        PROTECT(gg2=VECTOR_ELT(out_LD_efmi, SPPspictOPT[ind]-1));
        // PROTECT(gg3=VECTOR_ELT(out_L_eit, SPPspictOPT[ind]-1));

        TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))));
        // TAC_glob = REAL(getListElement(TAC, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))));
        double *totFM = REAL(aggregObj(gg1,nDimFM));
        double *totFM2 = REAL(aggregObj(gg2,nDimFM));
        // double *tot = REAL(aggregObj(gg3,nDim));
        // double *totMod = REAL(aggregObj(gg1,nDim));
        // double *totMod2 = REAL(aggregObj(gg2,nDim));

        for (int ind_f = 0 ; ind_f<nbF ; ind_f++) {

          double denom_st3 = 0.0; double alpha_f_st3 = 0.0;
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {

            if (g_effort1FM[ind_f + nbF*ind_m]>0)
              // TODO : denom = NA car sum(n,NA) == NA !!!
              denom_st3 = denom_st3 + mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * (totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
                (g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);
            //totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
            //totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
          }

          if ((denom_st3>0) & (g_effort1F[ind_f]>0)) {
            alpha_f_st3 = finite(TAC_byFleet[ind_f + nbF*IND_T] /(g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * denom_st3)); //Eq: A3
          }

          //r�conciliation
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
            // TODO : denom = NA car sum(n,NA) == NA !!!
            double valTest = alpha_f_st3 * g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];

            if ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0){
              if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
                Einterm_fm[1 + ind_f + nbF*ind_m] = EffsupTMP_fm[1 + ind_f + nbF*ind_m]; //on �limine le -1
                SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
              }
              if (valTest < Einterm_fm[1 + ind_f + nbF*ind_m] ) {//reconciliation seulement si la capture associee est non nulle
                Einterm_fm[1 + ind_f + nbF*ind_m] = valTest;
                SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPspictOPT[ind]-1));
                SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPspictOPT[ind]-1));
              }}
          }

        }
        //    double TACoth = TAC_glob[IND_T];
        //    for (int ind_ff = 0 ; ind_ff < nbF ; ind_ff++) TACoth = TACoth - TAC_byFleet[ind_ff + nbF*IND_T];
        //    //multF[ind_f+1] = finite(TACoth / (tot[IND_T] - totMod[IND_T] - totMod2[IND_T]));
        //    double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPspictOPT[ind]-1), 44));
        //    int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))), "modI"));
        //    for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*finite(TACoth / (tot[IND_T] - totMod[IND_T] - totMod2[IND_T])),0.0);
        UNPROTECT(3 - 1);
      }
      if(VERBOSE){Rprintf(". ");}
    }

    //on termine avec l'ajustement des esp�ces dynamiques sous XSA ou SS3 index�es dans SPPdynOPT (s'il y en a)
    if (N_SPPdynOPT>0) { // TODO : loop is not run if 0
      if(VERBOSE){Rprintf("Dyna");}

      //on finit d'initialiser Einterm_fm
      for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques restantes pour finaliser l'initialisation de Einterm_fm

        SEXP gg1=R_NilValue, gg2=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue;;
        PROTECT(Pgg1=VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1));
        PROTECT(Pgg2=VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1));
        PROTECT(gg1=aggregObj(Pgg1,nDimFM));
        PROTECT(gg2=aggregObj(Pgg2,nDimFM));

        double *totFM = REAL(gg1);
        double *totFM2 = REAL(gg2);

        for (int ind_f = 0 ; ind_f<nbF ; ind_f++){
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
            if ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0){
              if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
                Einterm_fm[1 + ind_f + nbF*ind_m] = EffsupTMP_fm[1 + ind_f + nbF*ind_m];//on �limine le -1
                SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
              }
            }
          }
        }
        UNPROTECT(4);
      }
      if(VERBOSE){Rprintf(". ");}
    }
    //fichier << "ST91" << endl;

    // Rprintf("\ncheck Einterm_fm : \n");
    // for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
    //   for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
    //     Rprintf("%4f ", Einterm_fm[1 + ind_f + nbF*ind_m]);
    //   }
    //   Rprintf("\n");
    // }

    // Rprintf("\ncheck Einterm_fm : \n");
    // for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
    //   for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
    //     Rprintf("%4f ", g_effort1FM[1 + ind_f + nbF*ind_m]);
    //   }
    //   Rprintf("\n");
    // }


    //les -1 restants sont remplac�s par les efforts initiaux --> les cellules fm non contraintes ne changent pas d'effort lors de l'ajustement
    for (int ind_f = 0 ; ind_f<nbF ; ind_f++){
      for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
        if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) { 
          // if(mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] == 0){// TODO : ici il ne faut pas mfm car 0 !
            // Einterm_fm[1 + ind_f + nbF*ind_m] = g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f];
          // } else {
            Einterm_fm[1 + ind_f + nbF*ind_m] = g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
          // }
          
          SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("---"));
          SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("---"));
        }
      }
    }
    if(VERBOSE){Rprintf(". ");}

    // Rprintf("effinitial Einterm_fm : \n");
    // for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
    //   for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
    //     Rprintf("%4f ", Einterm_fm[1 + ind_f + nbF*ind_m]);
    //   }
    //   Rprintf("\n");
    // }

    if (N_SPPdynOPT>0) {
      if(VERBOSE){Rprintf("Finding TAC");}

      int NBI = 0;

      ZoptSS3 = false;

      //SEXP reconcilSPP_copy;
      //on remet reconcilSPP au niveau qu'il avait avant l'engagement dans la boucle XSA/SS3

      //PROTECT(reconcilSPP_copy = duplicate(reconcilSPP));


      //on cr�e l'ensemble de matrices Einterm de longueur N_SPPdynOPT
      //double *EintermList_fm = NRvector(1,nbF*nbMe*N_SPPdynOPT);

      for (int indF = 0 ; indF < nbF ; indF++){  //initialisation de la copie de Einterm_fm
        for (int indM = 0 ; indM < nbMe ; indM++){
          Einterm_fm_copy[1 + indF + nbF*indM] = Einterm_fm[1 + indF + nbF*indM];
        }
      }

    // Rprintf("\neffort1F : \n");
    // for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
    //   // for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
    //     Rprintf("%4f ", g_effort1F[ind_f]);
    //   // }
    //   Rprintf("\n");
    // }

      //avant de passer aux ajustements suivants, on met a jour les efforts flottille-metier et flottille, puis on evalue Z associe
      for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        g_effort1F[ind_f] = 0.0;
        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
          g_effort1FM[ind_f + nbF*ind_m] = Einterm_fm[1 + ind_f + nbF*ind_m] / (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]) ;
          g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
          if (!ISNA(g_effort1FM[ind_f + nbF*ind_m])){
            g_effort1F[ind_f] = g_effort1F[ind_f] + g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m] ;
          }
        }
        g_effort1F[ind_f] = g_effort1F[ind_f] / (g_nbvF[ind_f] * g_effort2F[ind_f]) ;
        g_nbTripF[ind_f] = g_effort1F[ind_f] ;
      }

    // Rprintf("\neffort1F : \n");
    // for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
    //   // for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
    //     Rprintf("%4f ", g_effort1F[ind_f]);
    //   // }
    //   Rprintf("\n");
    // }

      Mortalite(listTempP, IND_T, eVar) ; //on genere les Z dans out_Z_eit (il faut s'assurer que l'application du Ztemp est bloquee -> utilisation de ZoptSS3)

      // on actualise les Ztemp ----
      for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) { //boucle sur les especes dynamiques XSA ou SS3 : initialisation des Z temporaires (ZtempList)

        NBI = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
        double *Zt = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;

        if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0) ) {
          for (int i = 0 ; i < NBI ; i++) {
            Zt[i] = REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] ; //XSA
          }

        } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1) ) {
          for (int i = 0 ; i < NBI ; i++){ // SEX
            Zt[i+(0*NBI)] = REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(1*NBI)] = REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
          }

        } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0) ){
          for (int i = 0 ; i < NBI ; i++) { //SS3

            Zt[i+(0*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(1*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(2*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(3*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(4*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(5*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(6*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(7*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(8*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(9*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(10*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(11*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(12*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(13*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(14*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            Zt[i+(15*NBI)] = REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;

          }
        }
      }

      if(VERBOSE){Rprintf("\n");}
      for (int IT = 0 ; (IT < ITtot) && (goon) ; IT++){
        if(VERBOSE){Rprintf("% ");}

        // on enclenche la boucle d'ajustement croise selon les phases par espece successives : fixation du Z, r�solution marginale, r�conciliation avec Einterm_fm, ajustement des Zfix pour chaque esp�ce
        for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) { //boucle sur les esp�ces dynamiques XSA ou SS3
          if (Qvec[SPPdynOPT[ind]-1]==1) ZoptSS3 = true; //SS3
        }

        for (int indF = 0 ; indF < nbF ; indF++){  //reinitialisation de Einterm_fm et de reconcilSPP
          for (int indM = 0 ; indM < nbMe ; indM++) {
            Einterm_fm[1 + indF + nbF*indM] = Einterm_fm_copy[1 + indF + nbF*indM] ;
            SET_STRING_ELT(reconcilSPP, indF + nbF*indM + nbF*nbMe*IND_T, STRING_ELT(reconcilSPP_copy,indF + nbF*indM + nbF*nbMe*IND_T));
          }
        }
        //on remet reconcilSPP_copy egal a reconcilSPP avant le lancement de l'ajustement XSA/SS3
        //reconcilSPP_copy = duplicate(reconcilSPP);

        // on reboote les efforts flottilles, on en deduit les efforts flottille-metier initiaux en fonction de U

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
          g_effort1F[ind_f] = g_effort1F_copy[ind_f] ;
          g_nbTripF[ind_f] = g_effort1F[ind_f] ;
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
            g_effort1FM[ind_f + nbF*ind_m] = g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m] ) ;
            g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
          }
        }


        Mortalite(listTempP, IND_T, eVar); //hors boucle espece a optimiser
        DynamicPop(listTempP, IND_T, eVar,true); //hors boucle espece a optimiser


        for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques XSA ou SS3

          int NBI2 = 0;
          NBI2 = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
          double *Zt = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;
          if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0)) {    //XSA
            for (int i = 0 ; i < NBI2 ; i++)
              REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i];

          } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1)) {    //sex-based
            for (int i = 0 ; i < NBI2 ; i++){
              REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(0*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(1*NBI2)];

            }

          } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0)) {                            //SS3 : sans doute inutile car boulot d�j� fait dans le module Dyn avec ZoptSS3 = true
            for (int i = 0 ; i < NBI2 ; i++)
            {
              REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(0*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(1*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(2*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(3*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(4*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(5*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(6*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(7*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(8*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(9*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(10*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(11*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(12*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(13*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(14*NBI2)];
              REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(15*NBI2)];
            }
          }
        }

        CatchDL(listTempP, IND_T, eVar, 0); //hors boucle esp�ce � optimiser

        for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les especes dynamiques XSA ou SS3 : determination de alpha, effort associe, reconciliation avec Einterm, maj Einterm
          SEXP gg1=R_NilValue, gg2=R_NilValue, /*gg3=R_NilValue,*/ gg4=R_NilValue, gg5=R_NilValue /*, gg6=R_NilValue, gg7=R_NilValue*/;

          PROTECT(gg1 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDimFM));
          PROTECT(gg2 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDimFM));
          // PROTECT(gg3 = aggregObj(VECTOR_ELT(out_L_eit, SPPdynOPT[ind]-1),nDim)); // TODO : remove
          PROTECT(gg4 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDimF));
          PROTECT(gg5 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDimF));
          // PROTECT(gg6 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDim)); // TODO : remove
          // PROTECT(gg7 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDim)); // TODO : remove

          TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))));
          // TAC_glob = REAL(getListElement(TAC, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))));
          double *totFM = REAL(gg1);
          double *totFM2 = REAL(gg2);
          double *totF = REAL(gg4);
          double *totF2 = REAL(gg5);
          // double *tot = REAL(gg3); // TODO : remove
          // double *totMod = REAL(gg6); // TODO : remove
          // double *totMod2 = REAL(gg7); // TODO : remove

          for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

            double alpha_f_st4 = 0.0 /*, alpha_f_st5 = 0.0*/;

            if (ind_f<nbF) {

              alpha_f_st4 = finite(TAC_byFleet[ind_f + nbF*IND_T] / (totF[ind_f + nbF*IND_T] + totF2[ind_f + nbF*IND_T])); //Eq: A4 ???????????????????

              //reconciliation
              for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
                double valTest = alpha_f_st4 * g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];

                if ( ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0) & (valTest < Einterm_fm[1 + ind_f + nbF*ind_m]) ) {//r�conciliation seulement si la capture associ�e est non nulle
                  Einterm_fm[1 + ind_f + nbF*ind_m] = valTest;
                  SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPdynOPT[ind]-1));
                }
              }

            } else {  //alpha_oth
              //fichier << "ST9.11" << endl;
              //                double TACoth = TAC_glob[IND_T];
              //                for (int ind_ff = 0 ; ind_ff < nbF ; ind_ff++){
              //                        TACoth = TACoth - TAC_byFleet[ind_ff + nbF*IND_T];
              //                }
              //                //fichier << "ST9.12" << endl;
              //                alpha_f_st5 = finite(TACoth / (tot[IND_T] - totMod[IND_T] - totMod2[IND_T]));
              //                //Rprintf("TACothIni %f TACoth %f tot %f totMod %f totMod2 %f\n",TAC_glob[IND_T],TACoth,tot[IND_T],totMod[IND_T],totMod2[IND_T]);
              //
              //                if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0)) {
              ////fichier << "ST9.13" << endl;
              //                    double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 44)); //PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 44)); Rprintf("alpha= %f \n",alpha_f_st5);
              //                    int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
              //                    for (int ag = 0; ag < ni; ag++) {
              //                            g_Fothi[ag + ni*IND_T] = g_Fothi[ag + ni*IND_T]*alpha_f_st5;
              //                            //Rprintf("Age: %f ; Fothi : %f \n",ag,g_Fothi[ag + ni*IND_T]);
              //                    }
              ////fichier << "ST9.14" << endl;
              //
              //                } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1)) {
              ////fichier << "ST9.13" << endl;
              //                    double *g_Fothi_G1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 224)); //PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 44)); Rprintf("alpha= %f \n",alpha_f_st5);
              //                    double *g_Fothi_G2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 225));
              //                    int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
              //                    for (int ag = 0; ag < ni; ag++) {
              //                            g_Fothi_G1[ag + ni*IND_T] = g_Fothi_G1[ag + ni*IND_T]*alpha_f_st5;
              //                            g_Fothi_G2[ag + ni*IND_T] = g_Fothi_G2[ag + ni*IND_T]*alpha_f_st5;
              //                            //Rprintf("Age: %f ; Fothi G1 : %f \n",ag,g_Fothi_G1[ag + ni*IND_T]);
              //                    }
              ////fichier << "ST9.14" << endl;
              //
              //                } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0)) {  //esp�ce SS3
              //
              ////fichier << "ST9.15" << endl;
              //                            double *Fothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 116));
              //                            double *Fothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 117));
              //                            double *Fothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 118));
              //                            double *Fothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 119));
              //                            double *Fothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 120));
              //                            double *Fothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 121));
              //                            double *Fothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 122));
              //                            double *Fothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 123));
              //                            double *Fothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 124));
              //                            double *Fothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 125));
              //                            double *Fothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 126));
              //                            double *Fothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 127));
              //                            double *Fothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 128));
              //                            double *Fothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 129));
              //                            double *Fothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 130));
              //                            double *Fothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 131));
              ////fichier << "ST9.16" << endl;
              //                            double *FRWTothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 176));
              //                            double *FRWTothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 177));
              //                            double *FRWTothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 178));
              //                            double *FRWTothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 179));
              //                            double *FRWTothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 180));
              //                            double *FRWTothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 181));
              //                            double *FRWTothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 182));
              //                            double *FRWTothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 183));
              //                            double *FRWTothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 184));
              //                            double *FRWTothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 185));
              //                            double *FRWTothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 186));
              //                            double *FRWTothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 187));
              //                            double *FRWTothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 188));
              //                            double *FRWTothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 189));
              //                            double *FRWTothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 190));
              //                            double *FRWTothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 191));
              ////fichier << "ST9.17" << endl;
              //                            double *FDWTothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 208));
              //                            double *FDWTothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 209));
              //                            double *FDWTothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 210));
              //                            double *FDWTothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 211));
              //                            double *FDWTothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 212));
              //                            double *FDWTothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 213));
              //                            double *FDWTothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 214));
              //                            double *FDWTothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 215));
              //                            double *FDWTothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 216));
              //                            double *FDWTothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 217));
              //                            double *FDWTothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 218));
              //                            double *FDWTothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 219));
              //                            double *FDWTothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 220));
              //                            double *FDWTothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 221));
              //                            double *FDWTothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 222));
              //                            double *FDWTothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPdynOPT[ind]-1), 223));
              ////fichier << "ST9.18" << endl;
              //                    int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
              //
              //                    for (int ag = 0; ag < ni; ag++) {
              //
              //                                Fothi_S1M1[ag + IND_T*ni] = Fothi_S1M1[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S1M2[ag + IND_T*ni] = Fothi_S1M2[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S1M3[ag + IND_T*ni] = Fothi_S1M3[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S1M4[ag + IND_T*ni] = Fothi_S1M4[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S2M1[ag + IND_T*ni] = Fothi_S2M1[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S2M2[ag + IND_T*ni] = Fothi_S2M2[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S2M3[ag + IND_T*ni] = Fothi_S2M3[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S2M4[ag + IND_T*ni] = Fothi_S2M4[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S3M1[ag + IND_T*ni] = Fothi_S3M1[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S3M2[ag + IND_T*ni] = Fothi_S3M2[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S3M3[ag + IND_T*ni] = Fothi_S3M3[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S3M4[ag + IND_T*ni] = Fothi_S3M4[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S4M1[ag + IND_T*ni] = Fothi_S4M1[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S4M2[ag + IND_T*ni] = Fothi_S4M2[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S4M3[ag + IND_T*ni] = Fothi_S4M3[ag + IND_T*ni]*alpha_f_st5;
              //                                Fothi_S4M4[ag + IND_T*ni] = Fothi_S4M4[ag + IND_T*ni]*alpha_f_st5;
              ////fichier << "ST9.19" << endl;
              //                                FRWTothi_S1M1[ag + IND_T*ni] = FRWTothi_S1M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S1M2[ag + IND_T*ni] = FRWTothi_S1M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S1M3[ag + IND_T*ni] = FRWTothi_S1M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S1M4[ag + IND_T*ni] = FRWTothi_S1M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S2M1[ag + IND_T*ni] = FRWTothi_S2M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S2M2[ag + IND_T*ni] = FRWTothi_S2M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S2M3[ag + IND_T*ni] = FRWTothi_S2M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S2M4[ag + IND_T*ni] = FRWTothi_S2M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S3M1[ag + IND_T*ni] = FRWTothi_S3M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S3M2[ag + IND_T*ni] = FRWTothi_S3M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S3M3[ag + IND_T*ni] = FRWTothi_S3M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S3M4[ag + IND_T*ni] = FRWTothi_S3M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S4M1[ag + IND_T*ni] = FRWTothi_S4M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S4M2[ag + IND_T*ni] = FRWTothi_S4M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S4M3[ag + IND_T*ni] = FRWTothi_S4M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FRWTothi_S4M4[ag + IND_T*ni] = FRWTothi_S4M4[ag + IND_T*ni]*alpha_f_st5;
              ////fichier << "ST9.21" << endl;
              //                                FDWTothi_S1M1[ag + IND_T*ni] = FDWTothi_S1M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S1M2[ag + IND_T*ni] = FDWTothi_S1M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S1M3[ag + IND_T*ni] = FDWTothi_S1M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S1M4[ag + IND_T*ni] = FDWTothi_S1M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S2M1[ag + IND_T*ni] = FDWTothi_S2M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S2M2[ag + IND_T*ni] = FDWTothi_S2M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S2M3[ag + IND_T*ni] = FDWTothi_S2M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S2M4[ag + IND_T*ni] = FDWTothi_S2M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S3M1[ag + IND_T*ni] = FDWTothi_S3M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S3M2[ag + IND_T*ni] = FDWTothi_S3M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S3M3[ag + IND_T*ni] = FDWTothi_S3M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S3M4[ag + IND_T*ni] = FDWTothi_S3M4[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S4M1[ag + IND_T*ni] = FDWTothi_S4M1[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S4M2[ag + IND_T*ni] = FDWTothi_S4M2[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S4M3[ag + IND_T*ni] = FDWTothi_S4M3[ag + IND_T*ni]*alpha_f_st5;
              //                                FDWTothi_S4M4[ag + IND_T*ni] = FDWTothi_S4M4[ag + IND_T*ni]*alpha_f_st5;
              //
              //                        }
              //
              //                    }
            }


          }
          UNPROTECT(7 - 3);

        }

        ZoptSS3 = false;
        //on evalue Z par espece associe a Einterm reconcilie
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
          g_effort1F[ind_f] = 0.0;
          for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
            g_effort1FM[ind_f + nbF*ind_m] = Einterm_fm[1 + ind_f + nbF*ind_m] / (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]) ;
            g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
            if (!ISNA(g_effort1FM[ind_f + nbF*ind_m]))
              g_effort1F[ind_f] = g_effort1F[ind_f] + g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m] ;
          }
          g_effort1F[ind_f] = g_effort1F[ind_f] / (g_nbvF[ind_f] * g_effort2F[ind_f]) ;
          g_nbTripF[ind_f] = g_effort1F[ind_f] ;
        }

        Mortalite(listTempP, IND_T, eVar) ; //on g�n�re les Z dans out_Z_eit (il faut s'assurer que l'application du Ztemp est bloqu�e -> utilisation de ZoptSS3)
        DynamicPop(listTempP, IND_T, eVar,true);

        //il faut maintenant en d�duire Z par esp�ce et modifier Zfix_e en cons�quence

        goon = false;

        for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les especes dynamiques XSA ou SS3 : determination de alpha, effort associe, reconciliation avec Einterm, maj Einterm
          int NBI = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
          double *Ztemp = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;

          if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0)) { // XSA
            for (int i = 0 ; i < NBI ; i++) {
              if (fabs(REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]) > diffZmax){
                //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                goon = true;
              }
              // TODO : replace by
              //goon = fabs(REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]) > diffZmax;
              Ztemp[i] = Ztemp[i] + lambda*(REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]);
            }

          } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1)) { // SEX
            for (int i = 0 ; i < NBI ; i++) {
              if (fabs(REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(0*NBI)] = Ztemp[i+(0*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(1*NBI)] = Ztemp[i+(1*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]);
            }

          } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0)){ //SS3
            for (int i = 0 ; i < NBI ; i++) {
              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(0*NBI)] = Ztemp[i+(0*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(1*NBI)] = Ztemp[i+(1*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(2*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(2*NBI)] = Ztemp[i+(2*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(2*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(3*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(3*NBI)] = Ztemp[i+(3*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(3*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(4*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(4*NBI)] = Ztemp[i+(4*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(4*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(5*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(5*NBI)] = Ztemp[i+(5*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(5*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(6*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(6*NBI)] = Ztemp[i+(6*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(6*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(7*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(7*NBI)] = Ztemp[i+(7*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(7*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(8*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(8*NBI)] = Ztemp[i+(8*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(8*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(9*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(9*NBI)] = Ztemp[i+(9*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(9*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(10*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(10*NBI)] = Ztemp[i+(10*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(10*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(11*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(11*NBI)] = Ztemp[i+(11*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(11*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(12*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(12*NBI)] = Ztemp[i+(12*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(12*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(13*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(13*NBI)] = Ztemp[i+(13*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(13*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(14*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(14*NBI)] = Ztemp[i+(14*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(14*NBI)]);

              if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(15*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
              Ztemp[i+(15*NBI)] = Ztemp[i+(15*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(15*NBI)]);
            }

          }
        }
      }
      //reconcilSPP = duplicate(reconcilSPP_copy);
      if(VERBOSE){Rprintf(".");}
    }



    if(VERBOSE){Rprintf("in Flist \n");}
    double *g_effort1FM_G = REAL(getListElement(FList, "effort1_f_m"));
    double *g_effort1F_G = REAL(getListElement(FList, "effort1_f"));
    double *g_nbTripF_G = REAL(getListElement(FList, "nbTrip_f"));
    double *g_nbvFM_G = REAL(getListElement(FList, "nbv_f_m"));
    double *g_nbvF_G = REAL(getListElement(FList, "nbv_f"));
    double *g_effort2FM_G = REAL(getListElement(FList, "effort2_f_m"));
    double *g_effort2F_G = REAL(getListElement(FList, "effort2_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
      g_effort1F_G[ind_f] = 0.0;
      for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
        g_effort1FM_G[ind_f + nbF*ind_m] = Einterm_fm[1 + ind_f + nbF*ind_m] / (g_effort2FM_G[ind_f + nbF*ind_m] * g_nbvFM_G[ind_f + nbF*ind_m]) ;
        g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
        if (!ISNA(g_effort1FM_G[ind_f + nbF*ind_m])){
          g_effort1F_G[ind_f] = g_effort1F_G[ind_f] + g_effort1FM_G[ind_f + nbF*ind_m] * g_effort2FM_G[ind_f + nbF*ind_m] * g_nbvFM_G[ind_f + nbF*ind_m] ;
        }
      }
      g_effort1F_G[ind_f] = g_effort1F_G[ind_f] / (g_nbvF_G[ind_f] * g_effort2F_G[ind_f]) ;
      g_nbTripF_G[ind_f] = g_effort1F_G[ind_f] ;
    }

    if(VERBOSE){Rprintf("End \n");}

    UNPROTECT(6 - 1);

    ZoptSS3 = false;
  }

}
