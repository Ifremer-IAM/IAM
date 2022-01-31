#include <Rdefines.h>
#include <Rmath.h>

#include "BioEcoPar.h" // Class is defined in this file.



//---------------------------------
//
//   Module de gestion
//
//---------------------------------



extern "C" {

double BioEcoPar::fxTAC_glob(double mult) //par temps IND_T pour une esp�ce donn�e
{
    SEXP listTemp;
    Rprintf("Hello ");
    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));
////Rprintf("3");
    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f_m"));////Rprintf("4");
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f"));////Rprintf("5");
    double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));////Rprintf("6");
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));////Rprintf("7");
    double *g_eff2FM2 = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort2_f_m"));
    double *g_eff2F2 = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort2_f"));
    double *mpond_fm = REAL(m_fm);////Rprintf("9");
    double *mpond_oth = REAL(m_oth);////Rprintf("10");

    double result;

//Rprintf("mult %f \n",mult);
//Rprintf("var %i \n",var);
//Rprintf("gestyp %i \n",gestyp);

     for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double countEff = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m] + mult*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m]*(1 + mult*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (var==2) {
                if (gestyp==1) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m] + mult*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m]*(1 + mult*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (!ISNA(g_nbdsFM[ind_f+nbF*ind_m]) & !ISNA(g_nbvFM[ind_f+nbF*ind_m]) & !ISNA(g_eff2FM2[ind_f+nbF*ind_m])) countEff = countEff + g_nbdsFM[ind_f+nbF*ind_m]*g_nbvFM[ind_f+nbF*ind_m]*g_eff2FM2[ind_f+nbF*ind_m];
        }

        if (var==1) g_nbdsF[ind_f] = fmax2(countEff/(g_nbvF[ind_f]*g_eff2F2[ind_f]),0.0);
        if (var==2) g_nbvF[ind_f] = fmax2(countEff/(g_nbdsF[ind_f]*g_eff2F2[ind_f]),0.0);
     }

Rprintf("there ");
    //    for (int e = 0 ; e < nbE ; e++){  // e --> eTemp

                if (Qvec[eTemp]==0) {

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

                if (gestyp==1)
                for (int ag = 0; ag < ni; ag++)
                  g_Fothi[ag + IND_T*ni] = fmax2(g_Fothi[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                if (gestyp==2)
                for (int ag = 0; ag < ni; ag++)
                  g_Fothi[ag + IND_T*ni] = fmax2(g_Fothi[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);


                } else {  //esp�ce SS3


                        double *Fothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 116));
                        double *Fothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 117));
                        double *Fothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 118));
                        double *Fothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 119));
                        double *Fothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 120));
                        double *Fothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 121));
                        double *Fothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 122));
                        double *Fothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 123));
                        double *Fothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 124));
                        double *Fothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 125));
                        double *Fothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 126));
                        double *Fothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 127));
                        double *Fothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 128));
                        double *Fothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 129));
                        double *Fothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 130));
                        double *Fothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 131));

                        double *FRWTothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 176));
                        double *FRWTothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 177));
                        double *FRWTothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 178));
                        double *FRWTothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 179));
                        double *FRWTothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 180));
                        double *FRWTothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 181));
                        double *FRWTothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 182));
                        double *FRWTothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 183));
                        double *FRWTothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 184));
                        double *FRWTothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 185));
                        double *FRWTothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 186));
                        double *FRWTothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 187));
                        double *FRWTothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 188));
                        double *FRWTothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 189));
                        double *FRWTothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 190));
                        double *FRWTothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 191));

                        double *FDWTothi_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 208));
                        double *FDWTothi_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 209));
                        double *FDWTothi_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 210));
                        double *FDWTothi_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 211));
                        double *FDWTothi_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 212));
                        double *FDWTothi_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 213));
                        double *FDWTothi_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 214));
                        double *FDWTothi_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 215));
                        double *FDWTothi_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 216));
                        double *FDWTothi_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 217));
                        double *FDWTothi_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 218));
                        double *FDWTothi_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 219));
                        double *FDWTothi_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 220));
                        double *FDWTothi_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 221));
                        double *FDWTothi_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 222));
                        double *FDWTothi_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, eTemp), 223));

                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

                if (gestyp==1)
                    for (int ag = 0; ag < ni; ag++) {

                            Fothi_S1M1[ag + IND_T*ni] = fmax2(Fothi_S1M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S1M2[ag + IND_T*ni] = fmax2(Fothi_S1M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S1M3[ag + IND_T*ni] = fmax2(Fothi_S1M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S1M4[ag + IND_T*ni] = fmax2(Fothi_S1M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S2M1[ag + IND_T*ni] = fmax2(Fothi_S2M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S2M2[ag + IND_T*ni] = fmax2(Fothi_S2M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S2M3[ag + IND_T*ni] = fmax2(Fothi_S2M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S2M4[ag + IND_T*ni] = fmax2(Fothi_S2M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S3M1[ag + IND_T*ni] = fmax2(Fothi_S3M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S3M2[ag + IND_T*ni] = fmax2(Fothi_S3M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S3M3[ag + IND_T*ni] = fmax2(Fothi_S3M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S3M4[ag + IND_T*ni] = fmax2(Fothi_S3M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S4M1[ag + IND_T*ni] = fmax2(Fothi_S4M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S4M2[ag + IND_T*ni] = fmax2(Fothi_S4M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S4M3[ag + IND_T*ni] = fmax2(Fothi_S4M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            Fothi_S4M4[ag + IND_T*ni] = fmax2(Fothi_S4M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);

                            FRWTothi_S1M1[ag + IND_T*ni] = fmax2(FRWTothi_S1M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S1M2[ag + IND_T*ni] = fmax2(FRWTothi_S1M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S1M3[ag + IND_T*ni] = fmax2(FRWTothi_S1M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S1M4[ag + IND_T*ni] = fmax2(FRWTothi_S1M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S2M1[ag + IND_T*ni] = fmax2(FRWTothi_S2M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S2M2[ag + IND_T*ni] = fmax2(FRWTothi_S2M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S2M3[ag + IND_T*ni] = fmax2(FRWTothi_S2M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S2M4[ag + IND_T*ni] = fmax2(FRWTothi_S2M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S3M1[ag + IND_T*ni] = fmax2(FRWTothi_S3M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S3M2[ag + IND_T*ni] = fmax2(FRWTothi_S3M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S3M3[ag + IND_T*ni] = fmax2(FRWTothi_S3M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S3M4[ag + IND_T*ni] = fmax2(FRWTothi_S3M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S4M1[ag + IND_T*ni] = fmax2(FRWTothi_S4M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S4M2[ag + IND_T*ni] = fmax2(FRWTothi_S4M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S4M3[ag + IND_T*ni] = fmax2(FRWTothi_S4M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FRWTothi_S4M4[ag + IND_T*ni] = fmax2(FRWTothi_S4M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);

                            FDWTothi_S1M1[ag + IND_T*ni] = fmax2(FDWTothi_S1M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S1M2[ag + IND_T*ni] = fmax2(FDWTothi_S1M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S1M3[ag + IND_T*ni] = fmax2(FDWTothi_S1M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S1M4[ag + IND_T*ni] = fmax2(FDWTothi_S1M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S2M1[ag + IND_T*ni] = fmax2(FDWTothi_S2M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S2M2[ag + IND_T*ni] = fmax2(FDWTothi_S2M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S2M3[ag + IND_T*ni] = fmax2(FDWTothi_S2M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S2M4[ag + IND_T*ni] = fmax2(FDWTothi_S2M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S3M1[ag + IND_T*ni] = fmax2(FDWTothi_S3M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S3M2[ag + IND_T*ni] = fmax2(FDWTothi_S3M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S3M3[ag + IND_T*ni] = fmax2(FDWTothi_S3M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S3M4[ag + IND_T*ni] = fmax2(FDWTothi_S3M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S4M1[ag + IND_T*ni] = fmax2(FDWTothi_S4M1[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S4M2[ag + IND_T*ni] = fmax2(FDWTothi_S4M2[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S4M3[ag + IND_T*ni] = fmax2(FDWTothi_S4M3[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);
                            FDWTothi_S4M4[ag + IND_T*ni] = fmax2(FDWTothi_S4M4[ag + IND_T*ni] + mult*mpond_oth[eTemp],0.0);

                    }

                if (gestyp==2)
                    for (int ag = 0; ag < ni; ag++) {

                            Fothi_S1M1[ag + IND_T*ni] = fmax2(Fothi_S1M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S1M2[ag + IND_T*ni] = fmax2(Fothi_S1M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S1M3[ag + IND_T*ni] = fmax2(Fothi_S1M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S1M4[ag + IND_T*ni] = fmax2(Fothi_S1M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S2M1[ag + IND_T*ni] = fmax2(Fothi_S2M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S2M2[ag + IND_T*ni] = fmax2(Fothi_S2M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S2M3[ag + IND_T*ni] = fmax2(Fothi_S2M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S2M4[ag + IND_T*ni] = fmax2(Fothi_S2M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S3M1[ag + IND_T*ni] = fmax2(Fothi_S3M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S3M2[ag + IND_T*ni] = fmax2(Fothi_S3M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S3M3[ag + IND_T*ni] = fmax2(Fothi_S3M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S3M4[ag + IND_T*ni] = fmax2(Fothi_S3M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S4M1[ag + IND_T*ni] = fmax2(Fothi_S4M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S4M2[ag + IND_T*ni] = fmax2(Fothi_S4M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S4M3[ag + IND_T*ni] = fmax2(Fothi_S4M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            Fothi_S4M4[ag + IND_T*ni] = fmax2(Fothi_S4M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);

                            FRWTothi_S1M1[ag + IND_T*ni] = fmax2(FRWTothi_S1M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S1M2[ag + IND_T*ni] = fmax2(FRWTothi_S1M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S1M3[ag + IND_T*ni] = fmax2(FRWTothi_S1M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S1M4[ag + IND_T*ni] = fmax2(FRWTothi_S1M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S2M1[ag + IND_T*ni] = fmax2(FRWTothi_S2M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S2M2[ag + IND_T*ni] = fmax2(FRWTothi_S2M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S2M3[ag + IND_T*ni] = fmax2(FRWTothi_S2M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S2M4[ag + IND_T*ni] = fmax2(FRWTothi_S2M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S3M1[ag + IND_T*ni] = fmax2(FRWTothi_S3M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S3M2[ag + IND_T*ni] = fmax2(FRWTothi_S3M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S3M3[ag + IND_T*ni] = fmax2(FRWTothi_S3M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S3M4[ag + IND_T*ni] = fmax2(FRWTothi_S3M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S4M1[ag + IND_T*ni] = fmax2(FRWTothi_S4M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S4M2[ag + IND_T*ni] = fmax2(FRWTothi_S4M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S4M3[ag + IND_T*ni] = fmax2(FRWTothi_S4M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FRWTothi_S4M4[ag + IND_T*ni] = fmax2(FRWTothi_S4M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);

                            FDWTothi_S1M1[ag + IND_T*ni] = fmax2(FDWTothi_S1M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S1M2[ag + IND_T*ni] = fmax2(FDWTothi_S1M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S1M3[ag + IND_T*ni] = fmax2(FDWTothi_S1M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S1M4[ag + IND_T*ni] = fmax2(FDWTothi_S1M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S2M1[ag + IND_T*ni] = fmax2(FDWTothi_S2M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S2M2[ag + IND_T*ni] = fmax2(FDWTothi_S2M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S2M3[ag + IND_T*ni] = fmax2(FDWTothi_S2M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S2M4[ag + IND_T*ni] = fmax2(FDWTothi_S2M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S3M1[ag + IND_T*ni] = fmax2(FDWTothi_S3M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S3M2[ag + IND_T*ni] = fmax2(FDWTothi_S3M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S3M3[ag + IND_T*ni] = fmax2(FDWTothi_S3M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S3M4[ag + IND_T*ni] = fmax2(FDWTothi_S3M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S4M1[ag + IND_T*ni] = fmax2(FDWTothi_S4M1[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S4M2[ag + IND_T*ni] = fmax2(FDWTothi_S4M2[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S4M3[ag + IND_T*ni] = fmax2(FDWTothi_S4M3[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                            FDWTothi_S4M4[ag + IND_T*ni] = fmax2(FDWTothi_S4M4[ag + IND_T*ni]*(1+mult*mpond_oth[eTemp]),0.0);
                    }

                }


        // }
Rprintf("general ");

if ((trgt==1) | (trgt==3) | (trgt==999)) {//on vise un TAC ou une biomasse

    Mortalite(listTemp, IND_T, eVarCopy);
    DynamicPop(listTemp, IND_T, eVarCopy,true);
    CatchDL(listTemp, IND_T, eVarCopy);
////PrintValue(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
////PrintValue(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));
//Rprintf("16");
Rprintf("le survivant ");
    double *tot ;
    if (trgt==999) {
     Mortalite(listTemp, IND_T+1, eVarCopy);
     DynamicPop(listTemp, IND_T+1, eVarCopy,true);  //� revoir !!!!!!!!!!!!!!!!
     tot = REAL(VECTOR_ELT(out_B_et, eTemp));
     result = TAC_glob[IND_T+1]-tot[IND_T+1]; //Rprintf("%f %f %f %f\n",mult,TAC_glob[IND_T+1],tot[IND_T+1],result);
    } else {
Rprintf("de tout les temps ");
     SEXP nDim = allocVector(INTSXP,4);
Rprintf(". ");
     int *nd = INTEGER(nDim); for (int i = 0; i<3; i++) nd[i] = 0; nd[3] = nbT;
     Rprintf(". ");
     tot = REAL(aggregObj(VECTOR_ELT(out_L_eit, eTemp),nDim));
     Rprintf(". ");
     result = TAC_glob[IND_T]-tot[IND_T]; //Rprintf("%f %f %f %f\n",mult,TAC_glob[IND_T],tot[IND_T],result);
     Rprintf("ouf ");
    }
//Rprintf("fxtac : TAC_glob %f TOT %f\n",TAC_glob[IND_T],tot[IND_T]);

} else { //on vise un Fbar
    Mortalite(listTemp, IND_T, eVarCopy);
    DynamicPop(listTemp, IND_T, eVarCopy,true);

    double *tot = REAL(VECTOR_ELT(out_Fbar_et, eTemp));
    result = Fbar_trgt[IND_T]-tot[IND_T]; //� interchanger avec la ligne du dessous pour une limitation plus restrictive selon SSB
    if ((Blim_trigger!=0) & !ISNA(Blim_CPP)) {
        result = Fbar_trgt[IND_T]*fmin2(REAL(VECTOR_ELT(out_SSB_et, eTemp))[IND_T]/Blim_CPP , 1.0) - tot[IND_T];
    } else {
        result = Fbar_trgt[IND_T] - tot[IND_T];
    }
    //Rprintf("%f %f %f %f\n",mult,Fbar_trgt[IND_T],tot[IND_T],result);
    //Rprintf("fxtac : FBARtarget %f TOT %f\n",Fbar_trgt[IND_T],tot[IND_T]);

}
Rprintf("kenobi ");
    UNPROTECT(2);

    return result;

}
}



//------------------------------------------
// Module de gestion : ajustement des variables d'effort (nbds (param�tre "var" = 1) ou nbv (param�tre "var" = 2))
// avec objectif d'atteinte du TAC (param�tre "trgt = 1") OU du Fbar (param�tre "trgt = 2")
// Un 3�me param�tre "delay" sp�cifie le d�lai de premi�re applicaton de l'ajustement (valeur par d�faut et minimale = 1).
// Enfin, un 4�me param�tre "upd" (update) permet de sp�cifier si le multiplicateur s'applique � la donn�e initiale � chaque pas de temps ("upd" = 1),
// ou si elle s'applique � la donn�e � l'instant pr�c�dent ("upd" = 2).
//------------------------------------------



extern "C" {

void BioEcoPar::Gestion(SEXP list, int ind_t, int VERBOSE) //param�tres en entr�e pas forc�ment utiles dans la mesure o� ils doivent rester constant tout au long de la simulation
{                                               //ajout de trgt = 22 pour consid�rer un ajustement inf�rieur � Fmsy, ie Fmsy*SSB/MSYBtrigger
                                                //ajout de trgt = 4 la biomasse limite sup�rieure (~Bmax) ?? : inactif pour le moment
//on teste la validit� des param�tres d'entr�e

    if ((var!=1) & (var!=2)) error("Wrong 'var' parameter in 'Gestion' module!!\n");
    if ((trgt!=1) & (trgt!=2) & (trgt!=3) & (trgt!=22) & (trgt!=999)) error("Wrong 'trgt' parameter in 'Gestion' module!!\n");
    if ((delay<1) | (delay>nbT)) error("Wrong 'delay' parameter in 'Gestion' module!!\n");
    if ((upd!=1) & (upd!=2)) error("Wrong 'upd' parameter in 'Gestion' module!!\n");

    if ((ind_t<delay) | ((trgt==999) & (ind_t==nbT-1))) {

    } else {

    IND_T = ind_t;//Rprintf("t %i",IND_T);
    double *mu_;//Rprintf("1");
    if (var==1){
        mu_ = REAL(mu_nbds);
    } else{
        mu_ = REAL(mu_nbv);
    }

    int NBMAX=1;
	int nb=NBMAX;
	float tol;
	int NbInter = 10;

    BEfn1 p = &BioEcoPar::fxTAC_glob;

    double *xb1,*xb2;
    xb1 = new double[NBMAX+1];
    xb2 = new double[NBMAX+1];
//Rprintf("2");
if(VERBOSE){Rprintf(" zbrak ");}
    zbrak(p,X1,X2,NbInter,xb1,xb2,&nb);
    if(VERBOSE){Rprintf("loop zbrent ");}
    for (int i=1;i<=nb;i++) {
        if(VERBOSE){Rprintf(".");}
        tol=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        double result=zbrent(p,xb1[i],xb2[i],tol);
        //Rprintf("result : %f",result);//Rprintf("t %i",IND_T);
        mu_[IND_T] = result;
    }
    delete xb1;
    delete xb2;

   }

}}




// Numerical Recipes //----------------------------------------------------------------------------------------

// --------  d�termination racine (unidimensionnel)

void BioEcoPar::zbrak(BEfn1 fx, double x1, double x2, int n, double xb1[],
	double xb2[], int *nb)
{
	int nbb,i;
	double x,fp,fc,dx;

	nbb=0;
	dx=(x2-x1)/n;
	fp=(this->*fx)(x=x1);
	for (i=1;i<=n;i++) {
		fc=(this->*fx)(x += dx);
		if (fc*fp <= 0.0) {
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb == nbb) return;

		}
		fp=fc;
	}
	*nb = nbb;
}
/* (C) Copr. 1986-92 Numerical Recipes Software *pA24. */


#define ITMAX 100
#define EPS 3.0e-6
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double  BioEcoPar::zbrent(BEfn1 fx, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
	double fa=(this->*fx)(a),fb=(this->*fx)(b),fc,p,q,r,s,tol1,xm;

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if (((fb > 0.0) && (fc > 0.0)) || ((fb < 0.0) && (fc < 0.0))) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if ((fabs(xm) <= tol1) || (fb == 0.0)) return b;
		if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb))) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(this->*fx)(b);
	}
//Rprintf("Maximum number of iterations exceeded in zbrent \n");
	return b;//0.0;  //modif 17/05/2013
}

#undef ITMAX
#undef EPS
#undef SIGN

