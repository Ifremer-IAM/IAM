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
// Module 'Report d'effort' selon une pond�ration des ratio profit par m�tier et effort par m�tier anticip�s
//------------------------------------------

extern "C" {

void BioEcoPar::FleetBehav(SEXP list, int ind_t, SEXP paramBehav) //ind_t>0
{

//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.FleetBehav.txt", ios::out | ios::trunc);

    SEXP Flist, nbds_f, nbds_f_m, fmtBhv, muBhv, alphaBhv, RTBS_f_m;

    PROTECT(Flist = getListElement(list, "Fleet"));
    PROTECT(nbds_f = getListElement(Flist, "nbds_f"));
    PROTECT(nbds_f_m = getListElement(Flist, "nbds_f_m"));
    PROTECT(fmtBhv = getListElement(paramBehav, "FMT"));
    PROTECT(muBhv = getListElement(paramBehav, "MU"));
    PROTECT(alphaBhv = getListElement(paramBehav, "ALPHA"));

    double *r_nbds_f = REAL(nbds_f), *r_nbds_f_m = REAL(nbds_f_m);
    int typeBhv = INTEGER(getListElement(paramBehav, "type"))[0];
    int posMuBhv = INTEGER(getListElement(paramBehav, "MUpos"))[0];
    bool isPos = (posMuBhv==1);

 //type n�1 : pas de report d'effort. Intervention sur l'effort au niveau flottille-m�tier via la matrice FMT
 // qui op�re additivement, avec redressement en cas d'effort r�sultant n�gatif ou sup�rieur � 365 somm� sur les m�tiers
 // L'effort au niveau flottille est ensuite r��valu� par agr�gation du niveau flottille-m�tier

    if ((typeBhv==1) & (fmtBhv != NULL)) {

       double *r_fmtBhv = REAL(fmtBhv);

       for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double sumEff_byF = 0.0;

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (!ISNA(r_nbds_f_m[ind_f + nbF*ind_m])) {

                r_nbds_f_m[ind_f + nbF*ind_m] = fmax2(r_nbds_f_m[ind_f + nbF*ind_m] + r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t],0.0); //NA+...=NA et max(NA,0)=NA
                sumEff_byF = sumEff_byF + r_nbds_f_m[ind_f + nbF*ind_m];

            }

        }

        //correction

        if (sumEff_byF>365.0) {

          for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

           r_nbds_f_m[ind_f + nbF*ind_m] = r_nbds_f_m[ind_f + nbF*ind_m]*365/sumEff_byF;

          }

        }

        // niveau flottille

        r_nbds_f[ind_f] = fmin2(sumEff_byF,365.0);

       }
    }


 //type n�2 : reports d'effort pilot�s. Intervention sur les m�tiers par flottille avec report conditionn� par une matrice FMT
 // de type :   | xx  xx   1 -0.5 -0.5   xx |
 //             | xx 0.7 0.3   xx -0.2 -0.8 |
 //             | ...                       |
 //
 // La quantit� brute de report par flottille-m�tier est ensuite �valu�e par multiplication de FMT par un vecteur MU de dimension nbF
 // MU est contraint pour que les reports soient coh�rents


    if ((typeBhv==2) & (fmtBhv != NULL) & (muBhv != NULL)) {

       double *r_fmtBhv = REAL(fmtBhv), *r_muBhv = REAL(muBhv);

       for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

           //d�termination de la validit� de MU_f et correction le cas �ch�ant

        double mu_limSup=-1.0, mu_limInf=0.0, finalMu=0.0;

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (!ISNA(r_nbds_f_m[ind_f + nbF*ind_m])) {

                if (mu_limSup<0) { //premi�re �valuation
                    if (r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]>0) {
                       mu_limSup = (365-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                       if (!isPos) mu_limInf = (0-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                    } else {
                       if (r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]<0) {
                        mu_limSup = (0-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                        if (!isPos) mu_limInf = (365-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t];
                       }
                    }
                } else {
                    if (r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]>0) {
                       mu_limSup = fmin2(mu_limSup , (365-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]);
                       if (!isPos) mu_limInf = fmax2(mu_limInf , (0-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]);
                    } else {
                       if (r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]<0) {
                        mu_limSup = fmin2(mu_limSup , (0-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]);
                        if (!isPos) mu_limInf = fmax2(mu_limInf , (365-r_nbds_f_m[ind_f + nbF*ind_m])/r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]);
                       }
                    }
                }
            }
        }

        mu_limSup = fmax2(mu_limSup,0.0);

        finalMu = fmax2(fmin2(r_muBhv[ind_f + nbF*ind_t],mu_limSup),mu_limInf);

        //calcul

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (!ISNA(r_nbds_f_m[ind_f + nbF*ind_m]))

                r_nbds_f_m[ind_f + nbF*ind_m] = r_nbds_f_m[ind_f + nbF*ind_m] + r_fmtBhv[ind_f + nbF*ind_m + nbF*nbMe*ind_t]*finalMu;

        }

        //normalement, pas besoin de r��valuer nbds_f car la conservation de l'effort est assur�e par la m�thodo
       }
    }




 //type n�3 : report d'effort orient� par pond�ration des ratio de profit et d'effort de l'ann�e pr�c�dente (cf P. Marchal).

    if ((typeBhv==3) & (ind_t>0) & (alphaBhv != NULL)) {

        if (ecodcf==0) {
            PROTECT(RTBS_f_m = VECTOR_ELT(out_Eco,10));
        } else {
            PROTECT(RTBS_f_m = VECTOR_ELT(out_EcoDCF,44));
        }

        double *r_RTBS_f_m = REAL(RTBS_f_m), *r_alphaBhv = REAL(alphaBhv);

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            r_alphaBhv[ind_f + nbF*ind_t] = fmax2(fmin2(r_alphaBhv[ind_f + nbF*ind_t],1.0),0.0);

            double totalRTBS_f = 0.0, totalEff_f = 0.0;

            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

                if (!ISNA(r_RTBS_f_m[ind_f + nbF*ind_m + nbF*nbMe*(ind_t-1)]))

                    totalRTBS_f = totalRTBS_f + fmax2(r_RTBS_f_m[ind_f + nbF*ind_m + nbF*nbMe*(ind_t-1)],0.0);

            }

            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

                if (!ISNA(r_nbds_f_m[ind_f + nbF*ind_m])) {

                    if (totalRTBS_f==0.0) {

                        r_nbds_f_m[ind_f + nbF*ind_m] = 0.0;

                    } else {

                        r_nbds_f_m[ind_f + nbF*ind_m] = r_nbds_f[ind_f] *
                                ((r_alphaBhv[ind_f + nbF*ind_t] * fmax2(r_RTBS_f_m[ind_f + nbF*ind_m + nbF*nbMe*(ind_t-1)],0.0) / totalRTBS_f) +
                                ((1 - r_alphaBhv[ind_f + nbF*ind_t]) * r_nbds_f_m[ind_f + nbF*ind_m] / r_nbds_f[ind_f]));

                        totalEff_f = totalEff_f + r_nbds_f_m[ind_f + nbF*ind_m];
                    }
                }

            }

            r_nbds_f[ind_f] = totalEff_f; //=0 si totalRTBS_f=0
        }

        UNPROTECT(1);

    }



    UNPROTECT(6);
}
}