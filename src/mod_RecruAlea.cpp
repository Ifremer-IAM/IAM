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
// Module de recrutement al�atoire
//
//---------------------------------



extern "C" {

void BioEcoPar::RecAlea(SEXP list, SEXP listSto, int ind_t, int type, int *recTyp) //list : liste des param�tres d'entr�e ; listSto : liste des variables d'op�rations stochastiques ; type : 1 -> samples sur l'historique (temps variable), 2 -> samples sur l'historique (temps constant), 3 -> loi de distribution
{

if (type<3) {

       SEXP elmtIn, elmtMeanSto, elmtResSto, MeanSto, ResSto, Rec, dimRec;
    //on tire au sort pour chacune des esp�ces mod�lis�es un r�sidu et on l'ajoute � la moyenne g�om�trique pr�-calcul�e
       int index = 0;

    for (int e = 0 ; e < nbE ; e++) {

                    if (recTyp[e]==1) {

                        PROTECT(elmtIn = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                        if (type==1) {

                            PROTECT(elmtMeanSto = getListElement(listSto, "GeoMeanRec"));
                            PROTECT(elmtResSto = getListElement(listSto, "RecResiduals"));

                        } else {

                            PROTECT(elmtMeanSto = getListElement(listSto, "GeoMeanRecLink"));
                            PROTECT(elmtResSto = getListElement(listSto, "RecResidualsLink"));

                        }

                        PROTECT(MeanSto = getListElement(elmtMeanSto, CHAR(STRING_ELT(sppList,e))));
                        PROTECT(ResSto = getListElement(elmtResSto, CHAR(STRING_ELT(sppList,e))));

                        int ll = length(ResSto);

                    if (ll > 0) {

                        if (type==1) { //multiple tirage d'indice (1 par esp�ce)

                            index = ll;
                            while (index >= ll) index = (int)(rand() / (((double)RAND_MAX + 1)/ ll));

                        } else {        //unique tirage d'indice pour les esp�ces consid�r�es (historiques de m�me taille)

                            if (e==0) {

                                index = ll;
                                while (index >= ll) index = (int)(rand() / (((double)RAND_MAX + 1)/ ll));

                            }

                        }

                        PROTECT(Rec = getListElement(elmtIn, "N_i0t"));
                        PROTECT(dimRec = getAttrib(Rec, install("DimCst")));


                        int *dimrec;
                        double *rec = REAL(Rec), *geomean = REAL(MeanSto), *residuals = REAL(ResSto);


                        dimrec = INTEGER(dimRec);

                        if (dimrec[3]==0) rec[0] = geomean[0] + residuals[index]; else rec[ind_t] = geomean[0] + residuals[index];

                        UNPROTECT(2);

                    }

                    UNPROTECT(5);

                    }

    }



} else {



    SEXP elmtIn, elmtDist, elmtDistParOne, elmtDistParTwo, elmtDistParThree,
         elmtDistSp, elmtDistParOneSp, elmtDistParTwoSp, elmtDistParThreeSp, Rec, dimRec;
    //on g�n�re une variable al�atoire suivant une loi log-normale de param�tres sp�cifi�s

    for (int e = 0 ; e < nbE ; e++) {

                if (recTyp[e]==1) {

                PROTECT(elmtIn = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                PROTECT(elmtDist = getListElement(listSto, "RecDist"));
                PROTECT(elmtDistParOne = getListElement(listSto, "RecDistPar1"));
                PROTECT(elmtDistParTwo = getListElement(listSto, "RecDistPar2"));
                PROTECT(elmtDistParThree = getListElement(listSto, "RecDistPar3"));

                PROTECT(elmtDistSp = getListElement(elmtDist, CHAR(STRING_ELT(sppList,e))));
                PROTECT(elmtDistParOneSp = getListElement(elmtDistParOne, CHAR(STRING_ELT(sppList,e))));
                PROTECT(elmtDistParTwoSp = getListElement(elmtDistParTwo, CHAR(STRING_ELT(sppList,e))));
                PROTECT(elmtDistParThreeSp = getListElement(elmtDistParThree, CHAR(STRING_ELT(sppList,e))));

                double v_a = NA_REAL;

                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "beta") == 0) v_a = rbeta(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
        //        if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "nbeta") == 0) v_a = rnbeta(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0],REAL(elmtDistParThreeSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "binom") == 0) v_a = rbinom(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "cauchy") == 0) v_a = rcauchy(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "chisq") == 0) v_a = rchisq(REAL(elmtDistParOneSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "nchisq") == 0) v_a = rnchisq(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "exp") == 0) v_a = rexp(REAL(elmtDistParOneSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "f") == 0) v_a = rf(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
        //        if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "nf") == 0) v_a = rnf(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0],REAL(elmtDistParThreeSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "gamma") == 0) v_a = rgamma(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "geom") == 0) v_a = rgeom(REAL(elmtDistParOneSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "hyper") == 0) v_a = rhyper(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0],REAL(elmtDistParThreeSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "logis") == 0) v_a = rlogis(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "lnorm") == 0) v_a = rlnorm(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "nbinom") == 0) v_a = rnbinom(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "norm") == 0) v_a = rnorm(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "pois") == 0) v_a = rpois(REAL(elmtDistParOneSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "t") == 0) v_a = rt(REAL(elmtDistParOneSp)[0]);
         //       if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "nt") == 0) v_a = rnt(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
         //       if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "tukey") == 0) v_a = rtukey(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0],REAL(elmtDistParThreeSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "unif") == 0) v_a = runif(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "weibull") == 0) v_a = rweibull(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "wilcox") == 0) v_a = rwilcox(REAL(elmtDistParOneSp)[0],REAL(elmtDistParTwoSp)[0]);
                if (strcmp(CHAR(STRING_ELT(elmtDistSp, 0)), "signrank") == 0) v_a = rsignrank(REAL(elmtDistParOneSp)[0]);

                //if ... pour les autres lois --> � compl�ter

                PROTECT(Rec = getListElement(elmtIn, "N_i0t"));
                PROTECT(dimRec = getAttrib(Rec, install("DimCst")));

                int *dimrec;
                double *rec = REAL(Rec);

                dimrec = INTEGER(dimRec);

                if (!ISNA(v_a)) {
                    if (dimrec[3]==0) rec[0] = v_a; else rec[ind_t] = v_a;
                }

                UNPROTECT(11);

            }

    }
}

}}


