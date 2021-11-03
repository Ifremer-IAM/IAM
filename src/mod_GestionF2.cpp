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
// Module 'Gestion F2'
//------------------------------------------


// TODO : change from to int to void
extern "C" {

int BioEcoPar::GestionF2(int ind_t)
{

//string str1, str2, str3;
//str1 = "testGestion";//"\\home1\\datahome\\fbriton\\AMURE\\Sc_bug_hke\\debugHKE_V";
//str3 = "_V";
//str2 = ".txt";
//
//std::stringstream ss, mp;
//mp << ind_t;
//ss << EcoIndCopy[0];
//str1 = str1 + mp.str()+ str3 + ss.str() + str2;
//
//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test_GestionF2.txt", ios::out | ios::trunc);

//
//fichier << "D�but" << endl;

    IND_T = ind_t;
//spQ = spp;

	int ITtot = maxIter;
//if (spp>=nbE | length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,spQ))), "modI"))==1) ITtot = 1; //si esp�ce statique ou SPiCT, pas besoin d'it�rations successives
	// si qqs i, |Z_i - Ztemp_i|<diffZmax, on arr�te, sinon on continue sous r�serve que IT<ITtot
	bool goon = true;

	//double ftol = 0.00000001;

    double *g_effSup = REAL(effSupMat);
    double *mpond_fm = REAL(out_allocEff_fm);

    //double *totFM, *totFM2, *totF, *totF2, *totFF, *totFF2, *tot, *totMod, *totMod2;

    SEXP listTempP, nDimFM, nDimF, nDim, copyEffort;

    PROTECT(nDimFM = allocVector(INTSXP,4));
    int *ndFM = INTEGER(nDimFM); ndFM[0] = nbF; ndFM[1] = nbMe; ndFM[2] = 0; ndFM[3] = nbT;
    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;

    PROTECT(listTempP = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));
//new
    double *g_effort1FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m"));
    double *g_effort1F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f"));
    double *g_nbTripFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f_m"));
    double *g_nbTripF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f"));
    double *g_nbvFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f_m"));
    double *g_nbvF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f"));
    double *g_effort2FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f_m"));
    double *g_effort2F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f"));


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



    for (int ind_f = 0 ; ind_f<nbF ; ind_f++) {//initialisation de l'effort de base

            g_effort1F[ind_f] = REAL(NBDSF)[ind_f + nbF*0];
            g_nbTripF[ind_f] = g_effort1F[ind_f];

        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {

            g_effort1FM[ind_f + nbF*ind_m] = (g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
                                                      (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);
            g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m];

        }

    }

//for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
// if ((ind_t==1) & (ind_f==0)) {

//std::stringstream ggg2;
 //       ggg2 << g_effort1FM[ind_f + nbF*ind_m];

 //       fichier << "effort_step2T1" << ggg2.str() << endl;
//
//    }
//}

    PROTECT(copyEffort = duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f")));
    //double *g_effort1FM_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m")));
    double *g_effort1F_copy = REAL(copyEffort); //duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f")));
    //double *g_effort2FM_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort2_f_m")));
    //double *g_effort2F_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort2_f")));
//new
    //Rprintf("call.Mortalite.1\n");fichier << "call.Mortalite.1" << endl;
    Mortalite(listTempP, IND_T, eVar);
    //Rprintf("end.Mortalite.1\n");fichier << "end.Mortalite.1" << endl;

    //Rprintf("call.DynamicPop.1\n");fichier << "call.DynamicPop.1" << endl;
    DynamicPop(listTempP, IND_T, eVar,true);
    //Rprintf("end.DynamicPop.1\n");fichier << "end.DynamicPop.1" << endl;

    //Rprintf("call.CatchDL.1\n");fichier << "call.CatchDL.1" << endl;
    CatchDL(listTempP, IND_T, eVar, 0);
    //Rprintf("end.CatchDL.1\n");fichier << "end.CatchDL.1" << endl;


//    double *g_effort1FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m"));
//    double *g_effort1FM_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m")));
//    double *g_effort1F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f"));
//    double *g_effort1F_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort1_f")));
//    double *g_nbTripFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f_m"));
//    double *g_nbTripF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f"));
//    double *g_nbvFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f_m"));
//    double *g_nbvF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbv_f"));
//    double *g_effort2FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f_m"));
//    double *g_effort2FM_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort2_f_m")));
//    double *g_effort2F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort2_f"));
//    double *g_effort2F_copy = REAL(duplicate(getListElement(getListElement(listTempP, "Fleet"), "effort2_f")));

//fichier << "ST1" << endl;
//on initialise en remplissant Einterm_fm avec les valeurs maximales autoris�es (fonction de effSupMat et la pond�ration multiplicative U (mpond_fm)) : Eq A1
//on rappelle que tout cela fonctionne uniquement si on est sur une base individuelle (dans le cas contraire, il faudrait descendre au niveau m�tier et pond�rer par les nbv par m�tier

    for (int indF = 0 ; indF < nbF ; indF++) {
        if (ISNA(g_effSup[indF + nbF*ind_t])) g_effSup[indF + nbF*ind_t] = 300; //ATTENTION : important � noter, si un navire n'est contraint par rien, c'est direct 300 jours de mer
        double denom_st1 = 0.0; double alpha_f_st1 = 0.0;
        for (int indM = 0 ; indM<nbMe ; indM++) {
                denom_st1 = denom_st1 + (g_effort1F[indF] * g_effort2F[indF] * mpond_fm[indF + nbF*indM + nbF*nbMe*IND_T]); //c'est cens� valoir eff_f normalement si la somme des facteurs de pond. vaut 1
        }
        if (denom_st1>0) alpha_f_st1 = g_effSup[indF + nbF*ind_t] / denom_st1;
        for (int indM = 0 ; indM<nbMe ; indM++) {
                EffsupTMP_fm[1 + indF + nbF*indM] = alpha_f_st1 * g_effort1F[indF] * g_nbvF[indF] * g_effort2F[indF] * mpond_fm[indF + nbF*indM + nbF*nbMe*IND_T];
    }
    }
//fichier << "ST2" << endl;
    for (int indF = 0 ; indF < nbF ; indF++)
    for (int indM = 0 ; indM<nbMe ; indM++)
      Einterm_fm[1 + indF + nbF*indM] = -1.0; //va servir � identifier � la fin les cellules non r�concili�es, qui devront �tre �gales � l'effort courant!!!

//� ce stade, on a initialis� Einterm_fm
//fichier << "ST3" << endl;
//on commence l'ajustement avec les esp�ces statiques index�es dans SPPstatOPT (s'il y en a)
if (N_SPPstatOPT>0) {

    for (int ind = 0 ; ind < N_SPPstatOPT ; ind++) { //boucle sur ces esp�ces statiques

        SEXP gg1=R_NilValue, gg2=R_NilValue;
        PROTECT(gg1=VECTOR_ELT(out_Lstat, SPPstatOPT[ind]-1));
        PROTECT(gg2=VECTOR_ELT(out_statLD_efm, SPPstatOPT[ind]-1));

        double *TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppListStat,SPPstatOPT[ind]-1))));
        double *totFM = REAL(aggregObj(gg1,nDimFM)); //PrintValue(aggregObj(VECTOR_ELT(out_L_efmit, eTemp),nDimFM));
        double *totFM2 = REAL(aggregObj(gg2,nDimFM));
        //double *totF = REAL(aggregObj(VECTOR_ELT(out_Lstat, SPPstatOPT[ind]-1),nDimF));
        //double *totF2 = REAL(aggregObj(VECTOR_ELT(out_statLD_efm, SPPstatOPT[ind]-1),nDimF));
//fichier << "ST4" << endl;
        for (int ind_f = 0 ; ind_f<nbF ; ind_f++) {

            double denom_st2 = 0.0; double alpha_f_st2 = 0.0;
            for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {

              if (g_effort1FM[ind_f + nbF*ind_m]>0)
                 denom_st2 = denom_st2 + mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * (totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
                                  (g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);
                    //Rprintf("%f\n", totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T]);
                //totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
                    //Rprintf("%f\n", totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T]);
                //totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
            }

            if ((denom_st2>0) & (g_effort1F[ind_f]>0))
                    alpha_f_st2 = finite(TAC_byFleet[ind_f + nbF*IND_T] / (g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * denom_st2));

            //r�conciliation
            for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
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
}


//on poursuit avec l'ajustement des esp�ces dynamiques sous SPiCT index�es dans SPPspictOPT (s'il y en a)
if (N_SPPspictOPT>0) {

    for (int ind = 0 ; ind < N_SPPspictOPT ; ind++) { //boucle sur ces esp�ces Spict
//fichier << "ST5" << endl;
        SEXP gg1=R_NilValue, gg2=R_NilValue, gg3=R_NilValue;
        PROTECT(gg1=VECTOR_ELT(out_L_efmit, SPPspictOPT[ind]-1));
        PROTECT(gg2=VECTOR_ELT(out_LD_efmi, SPPspictOPT[ind]-1));
        PROTECT(gg3=VECTOR_ELT(out_L_eit, SPPspictOPT[ind]-1));

        TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))));//fichier << "ST51" << endl;
        TAC_glob = REAL(getListElement(TAC, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))));//fichier << "ST52" << endl;
        double *totFM = REAL(aggregObj(gg1,nDimFM));//fichier << "ST53" << endl; //PrintValue(VECTOR_ELT(out_L_efmit, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_L_efmit, eTemp),nDimFM));
        double *totFM2 = REAL(aggregObj(gg2,nDimFM));//fichier << "ST54" << endl; //PrintValue(VECTOR_ELT(out_LD_efmi, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_LD_efmi, eTemp),nDimFM));
        //double *totF = REAL(aggregObj(VECTOR_ELT(out_L_efmit, SPPspictOPT[ind]-1),nDimF));//fichier << "ST55" << endl;
        //double *totF2 = REAL(aggregObj(VECTOR_ELT(out_LD_efmi, SPPspictOPT[ind]-1),nDimF));//fichier << "ST56" << endl;
        double *tot = REAL(aggregObj(gg3,nDim));//fichier << "ST57" << endl;
        double *totMod = REAL(aggregObj(gg1,nDim));//fichier << "ST58" << endl;
        double *totMod2 = REAL(aggregObj(gg2,nDim));//fichier << "ST59" << endl;
        //fichier << "ST6" << endl;
        for (int ind_f = 0 ; ind_f<nbF ; ind_f++) {

            double denom_st3 = 0.0; double alpha_f_st3 = 0.0;
            for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {

              if (g_effort1FM[ind_f + nbF*ind_m]>0)
                 //fichier << "ST6.1" << endl;
                 denom_st3 = denom_st3 + mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * (totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T]) /
                                 (g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]);
                    //Rprintf("%f\n", totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T]);
                //totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
                    //Rprintf("%f\n", totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T]);
                //totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] = mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] * totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / g_effort1FM[ind_f + nbF*ind_m];
            }

            if ((denom_st3>0) & (g_effort1F[ind_f]>0)) {
                //fichier << "ST6.2" << endl;
                alpha_f_st3 = finite(TAC_byFleet[ind_f + nbF*IND_T] /(g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * denom_st3)); //Eq: A3
            }

            //r�conciliation
            for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
                double valTest = alpha_f_st3 * g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
                if ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0){
                    //fichier << "ST6.3" << endl;
                    if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
                    //fichier << "ST6.4" << endl;
                    Einterm_fm[1 + ind_f + nbF*ind_m] = EffsupTMP_fm[1 + ind_f + nbF*ind_m]; //on �limine le -1
                    SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                    SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                    }
                    if (valTest < Einterm_fm[1 + ind_f + nbF*ind_m] ) {//r�conciliation seulement si la capture associ�e est non nulle
                    //fichier << "ST6.5" << endl;
                    Einterm_fm[1 + ind_f + nbF*ind_m] = valTest;
                    SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPspictOPT[ind]-1));
                    SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPspictOPT[ind]-1));
                }}
            }
            //Rprintf("T %i F %i \n",ind_t,ind_f);

        }
//fichier << "ST7" << endl;
//    double TACoth = TAC_glob[IND_T];
//    for (int ind_ff = 0 ; ind_ff < nbF ; ind_ff++) TACoth = TACoth - TAC_byFleet[ind_ff + nbF*IND_T];
//    //multF[ind_f+1] = finite(TACoth / (tot[IND_T] - totMod[IND_T] - totMod2[IND_T]));
//    double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPspictOPT[ind]-1), 44));
//    int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPspictOPT[ind]-1))), "modI"));
//    for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*finite(TACoth / (tot[IND_T] - totMod[IND_T] - totMod2[IND_T])),0.0);
//fichier << "ST8" << endl;
    UNPROTECT(3);
    }
}



//� ce stade, on garde la trace de Einterm_fm et listTempP original
//SEXP EintermTMP, listTempPP;
//PROTECT(EintermTMP = duplicate(Einterm_fm));
//PROTECT(listTempPP = duplicate(listTempP));

//on termine avec l'ajustement des esp�ces dynamiques sous XSA ou SS3 index�es dans SPPdynOPT (s'il y en a)

if (N_SPPdynOPT>0) {
//fichier << "ST7" << endl;
    //on finit d'initialiser Einterm_fm
     for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques restantes pour finaliser l'initialisation de Einterm_fm
//fichier << "ST9" << endl;
            SEXP gg1=R_NilValue, gg2=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue;;
            PROTECT(Pgg1=VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1));
            PROTECT(Pgg2=VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1));
            PROTECT(gg1=aggregObj(Pgg1,nDimFM));
            PROTECT(gg2=aggregObj(Pgg2,nDimFM));

            double *totFM = REAL(gg1); //fichier << "ST9a" << endl;//PrintValue(VECTOR_ELT(out_L_efmit, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_L_efmit, eTemp),nDimFM));
            double *totFM2 = REAL(gg2); //fichier << "ST9b" << endl;//PrintValue(VECTOR_ELT(out_LD_efmi, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_LD_efmi, eTemp),nDimFM));

            for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
            for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
               if ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0){
                        //fichier << "ST7.1" << endl;
                        if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
                          //fichier << "ST7.2" << endl;
                          Einterm_fm[1 + ind_f + nbF*ind_m] = EffsupTMP_fm[1 + ind_f + nbF*ind_m]; //fichier << "ST9c" << endl;//on �limine le -1
                          //if (ind_f==0 & ind_m==0) PrintValue(reconcilSPP);
                          SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                          SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("MAX"));
                        }
                }
            }
             UNPROTECT(4);
     }
}
//fichier << "ST91" << endl;
 //les -1 restants sont remplac�s par les efforts initiaux --> les cellules fm non contraintes ne changent pas d'effort lors de l'ajustement
for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
      if (Einterm_fm[1 + ind_f + nbF*ind_m]<-0.5) {
      //fichier << "ST7.3" << endl;
      Einterm_fm[1 + ind_f + nbF*ind_m] = g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
      SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("---"));
      SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*IND_T, mkChar("---"));
      }
}


//fichier << "ST10" << endl;

if (N_SPPdynOPT>0) {

    int NBI = 0;

    ZoptSS3 = false;

    //SEXP reconcilSPP_copy;
    //on remet reconcilSPP au niveau qu'il avait avant l'engagement dans la boucle XSA/SS3

    //PROTECT(reconcilSPP_copy = duplicate(reconcilSPP));


    //on cr�e l'ensemble de matrices Einterm de longueur N_SPPdynOPT
//double *EintermList_fm = NRvector(1,nbF*nbMe*N_SPPdynOPT);

//fichier << "ST7.31" << endl;
for (int indF = 0 ; indF < nbF ; indF++)  //initialisation de la copie de Einterm_fm
for (int indM = 0 ; indM < nbMe ; indM++) Einterm_fm_copy[1 + indF + nbF*indM] = Einterm_fm[1 + indF + nbF*indM] ;

//for (int indF = 0 ; indF < nbF ; indF++)  //initialisation des copies de Einterm_fm qui serviront de base de r�conciliation sur le'ensemble des esp�ces dynamiques (hors Spict)
//for (int indM = 0 ; indM < nbMe ; indM++)
//for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) EintermList_fm[1 + indF + nbF*indM+ nbF*nbMe*ind] = Einterm_fm[1 + indF + nbF*indM] ;

//fichier << "ST11" << endl;
    //avant de passer aux ajustements suivants, on met � jour les efforts flottille-m�tier et flottille, puis on �value Z associ�
//fichier << "ST7.32" << endl;
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        //fichier << "ST7.34" << endl;
        g_effort1F[ind_f] = 0.0;
        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
                //fichier << "ST7.35" << endl;
                g_effort1FM[ind_f + nbF*ind_m] = Einterm_fm[1 + ind_f + nbF*ind_m] / (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m]) ;
                g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
                if (!ISNA(g_effort1FM[ind_f + nbF*ind_m])){
                     g_effort1F[ind_f] = g_effort1F[ind_f] + g_effort1FM[ind_f + nbF*ind_m] * g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m] ;
                     //fichier << "ST7.4" << endl;
                }
        }
        //fichier << "ST7.41" << endl;
        g_effort1F[ind_f] = g_effort1F[ind_f] / (g_nbvF[ind_f] * g_effort2F[ind_f]) ;
        g_nbTripF[ind_f] = g_effort1F[ind_f] ;
        }


//        for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg3;
//                ggg3 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step3T1 " << ggg3.str() << endl;
//
//            }
//        }


//fichier << "ST7.42" << endl;
//fichier << "ST12" << endl;
    //Rprintf("call.Mortalite.SPPdynOPT.1\n");fichier << "call.Mortalite.SPPdynOPT.1" << endl;
    Mortalite(listTempP, IND_T, eVar) ; //on g�n�re les Z dans out_Z_eit (il faut s'assurer que l'application du Ztemp est bloqu�e -> utilisation de ZoptSS3)
    //Rprintf("end.Mortalite.SPPdynOPT.1\n");fichier << "end.Mortalite.SPPdynOPT.1" << endl;
//fichier << "ST7.43" << endl;
// on actualise les Ztemp ----
    for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) { //boucle sur les esp�ces dynamiques XSA ou SS3 : initialisation des Z temporaires (ZtempList)

        NBI = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
//fichier << "ST13" << endl;
        double *Zt = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;
//fichier << "ST7.44" << endl;
        if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0) ) {
            //fichier << "ST8.1" << endl;
            for (int i = 0 ; i < NBI ; i++) Zt[i] = REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] ; //XSA

        } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1) ) {
            //fichier << "ST8.1" << endl;
            for (int i = 0 ; i < NBI ; i++){
                    Zt[i+(0*NBI)] = REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
                    Zt[i+(1*NBI)] = REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] ;
            }

        } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0) ){
            //fichier << "ST8.2" << endl;
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
        //fichier << "ST14" << endl;
    }



    for (int IT = 0 ; IT < ITtot ; IT++){
//fichier << "ST8.20" << endl;
       if (goon) {

        // on enclenche la boucle d'ajustement crois� selon les phases par esp�ce successives : fixation du Z, r�solution marginale, r�conciliation avec Einterm_fm, ajustement des Zfix pour chaque esp�ce
//fichier << "ST8.21" << endl;
        for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) { //boucle sur les esp�ces dynamiques XSA ou SS3
            if (Qvec[SPPdynOPT[ind]-1]==1) ZoptSS3 = true; //SS3
        }

//
//        for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg00;
//                ggg00 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step00T1 " << ggg00.str() << endl;
//
//            }
//        }

//fichier << "ST8.22" << endl;
for (int indF = 0 ; indF < nbF ; indF++)  //r�initialisation de Einterm_fm et de reconcilSPP
for (int indM = 0 ; indM < nbMe ; indM++) {
        Einterm_fm[1 + indF + nbF*indM] = Einterm_fm_copy[1 + indF + nbF*indM] ;
        SET_STRING_ELT(reconcilSPP, indF + nbF*indM + nbF*nbMe*IND_T, STRING_ELT(reconcilSPP_copy,indF + nbF*indM + nbF*nbMe*IND_T));
}
//fichier << "ST8.23" << endl;
//on remet reconcilSPP_copy �gal � reconcilSPP avant le lancement de l'ajustement XSA/SS3
//reconcilSPP_copy = duplicate(reconcilSPP);

    // on reboote les efforts flottilles, on en d�duit les efforts flottille-m�tier initiaux en fonction de U


//        for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg0, ggg01, ggg02, ggg03, ggg04, ggg05, ggg06, ggg07;
//                ggg0 << g_effort1FM[ind_f + nbF*ind_m];
//                ggg01 << g_effort1F[ind_f];
//                ggg02 << g_effort2F[ind_f];
//                ggg03 << g_nbvF[ind_f];
//                ggg04 << mpond_fm[ind_f + nbF*ind_m];
//                ggg05 << g_effort2FM[ind_f + nbF*ind_m];
//                ggg06 << g_nbvFM[ind_f + nbF*ind_m];
//                ggg07 << g_effort1F_copy[ind_f];
//
//                fichier << "effort_step0T1 " << ggg0.str() << endl;
//                fichier << " ef1F " << ggg01.str() << " ef2F " << ggg02.str() << " nbvF " << ggg03.str() <<" pond " << ggg04.str() << " ef2FM " << ggg05.str() << " nbvFM " << ggg06.str() << " eff1F_copy " << ggg07.str() << endl;
//
//            }
//        }


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        g_effort1F[ind_f] = g_effort1F_copy[ind_f] ;
        g_nbTripF[ind_f] = g_effort1F[ind_f] ;
        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
                    g_effort1FM[ind_f + nbF*ind_m] = g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T] / (g_effort2FM[ind_f + nbF*ind_m] * g_nbvFM[ind_f + nbF*ind_m] ) ;
                    g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
        }
    }


//            for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg4;
//                ggg4 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step4T1 " << ggg4.str() << endl;
//
//            }
//        }


//fichier << "ST8.24" << endl;
//fichier << "ST15" << endl;
        //Rprintf("call.Mortalite.SPPdynOPT.2\n");fichier << "call.Mortalite.SPPdynOPT.2" << endl;
        Mortalite(listTempP, IND_T, eVar); //hors boucle esp�ce � optimiser
        //Rprintf("end.Mortalite.SPPdynOPT.2\n");fichier << "end.Mortalite.SPPdynOPT.2" << endl;

        //Rprintf("call.DynamicPop.SPPdynOPT.2\n");fichier << "call.DynamicPop.SPPdynOPT.2" << endl;
        DynamicPop(listTempP, IND_T, eVar,true); //hors boucle esp�ce � optimiser
        //Rprintf("end.DynamicPop.SPPdynOPT.2\n");fichier << "end.DynamicPop.SPPdynOPT.2" << endl;
//fichier << "ST8.25" << endl;
    for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques XSA ou SS3

        int NBI2 = 0;
//fichier << "ST8.26" << endl;
        NBI2 = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
        double *Zt = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;
//fichier << "ST8.27" << endl;
//fichier << "ST16" << endl;
        if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0)) {    //XSA
          //fichier << "ST8.28" << endl;
          for (int i = 0 ; i < NBI2 ; i++)
            REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i];

        } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1)) {    //sex-based
          //fichier << "ST8.28" << endl;
          for (int i = 0 ; i < NBI2 ; i++){
            REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(0*NBI2)];
            REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI2*IND_T] = Zt[i+(1*NBI2)];

          }

        } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0)) {                            //SS3 : sans doute inutile car boulot d�j� fait dans le module Dyn avec ZoptSS3 = true
          //fichier << "ST8.29" << endl;
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
//fichier << "ST17" << endl;
    }

        //Rprintf("call.CatchDL.SPPdynOPT.2\n");fichier << "call.CatchDL.SPPdynOPT.2" << endl;
        CatchDL(listTempP, IND_T, eVar, 0); //hors boucle esp�ce � optimiser
        //Rprintf("end.CatchDL.SPPdynOPT.2\n");fichier << "end.CatchDL.SPPdynOPT.2" << endl;

    for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques XSA ou SS3 : d�termination de alpha, effort associ�, r�conciliation avec Einterm, m�j Einterm
//fichier << "ST18" << endl;
//fichier << "ST8.29" << endl;
       SEXP gg1=R_NilValue, gg2=R_NilValue, gg3=R_NilValue, gg4=R_NilValue, gg5=R_NilValue, gg6=R_NilValue, gg7=R_NilValue;

        PROTECT(gg1 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDimFM));
        PROTECT(gg2 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDimFM));
        PROTECT(gg3 = aggregObj(VECTOR_ELT(out_L_eit, SPPdynOPT[ind]-1),nDim));
        PROTECT(gg4 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDimF));
        PROTECT(gg5 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDimF));
        PROTECT(gg6 = aggregObj(VECTOR_ELT(out_L_efmit, SPPdynOPT[ind]-1),nDim));
        PROTECT(gg7 = aggregObj(VECTOR_ELT(out_LD_efmi, SPPdynOPT[ind]-1),nDim));

        TAC_byFleet = REAL(getListElement(TACbyF, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))));
        TAC_glob = REAL(getListElement(TAC, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))));
        double *totFM = REAL(gg1); //PrintValue(VECTOR_ELT(out_L_efmit, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_L_efmit, eTemp),nDimFM));
        double *totFM2 = REAL(gg2); //PrintValue(VECTOR_ELT(out_LD_efmi, eTemp)) ; PrintValue(aggregObj(VECTOR_ELT(out_LD_efmi, eTemp),nDimFM));
        double *totF = REAL(gg4);
        double *totF2 = REAL(gg5);
        double *tot = REAL(gg3);
        double *totMod = REAL(gg6);
        double *totMod2 = REAL(gg7);
//fichier << "ST8.291" << endl;

//fichier << "ST19" << endl;
        for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

                    ////Rprintf("TACf %f totF %f totF2 %f denom %f mult %f\n",TAC_byFleet[ind_f + nbF*IND_T],totF[ind_f + nbF*IND_T],totF2[ind_f + nbF*IND_T],denom,multF[ind_f+1]);
                    ////Rprintf("ind_f %i multF %f \n",ind_f,multF[ind_f+1]);
                    //multF[ind_f+1] = 1; }
//fichier << "ST8.292" << endl;
            double alpha_f_st4 = 0.0, alpha_f_st5 = 0.0;

            if (ind_f<nbF) {

                alpha_f_st4 = finite(TAC_byFleet[ind_f + nbF*IND_T] / (totF[ind_f + nbF*IND_T] + totF2[ind_f + nbF*IND_T])); //Eq: A4 ???????????????????
//fichier << "ST8.293" << endl;
//std::stringstream ff;
//ff << alpha_f_st4;
//fichier << ff.str() << endl;

                //r�conciliation
                for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//                        fichier << "ST8.294" << endl;
                    double valTest = alpha_f_st4 * g_effort1F[ind_f] * g_effort2F[ind_f] * g_nbvF[ind_f] * mpond_fm[ind_f + nbF*ind_m + nbF*nbMe*IND_T];

//                    std::stringstream gg;
//                    gg << totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
//                    std::stringstream hh;
//                    hh << totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T];
//                    fichier << gg.str() << endl;
//                    fichier << hh.str() << endl;

                    if ( ((totFM[ind_f + nbF*ind_m + nbF*nbMe*IND_T] + totFM2[ind_f + nbF*ind_m + nbF*nbMe*IND_T])>0) & (valTest < Einterm_fm[1 + ind_f + nbF*ind_m]) ) {//r�conciliation seulement si la capture associ�e est non nulle
//                        fichier << "ST9.1" << endl;
                        Einterm_fm[1 + ind_f + nbF*ind_m] = valTest;
                        SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*IND_T, STRING_ELT(sppList,SPPdynOPT[ind]-1));
                }
                }
                //Rprintf("T %i F %i \n",ind_t,ind_f);

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
            UNPROTECT(7);

           }

    ZoptSS3 = false;
//fichier << "ST20" << endl;
    //on �value Z par esp�ce associ� � Einterm r�concili�
//fichier << "ST9.22" << endl;
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

//
//        for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg5;
//                ggg5 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step5T1 " << ggg5.str() << endl;
//
//            }
//        }

//fichier << "ST9.23" << endl;
//fichier << "ST21" << endl;
    //Rprintf("call.Mortalite.SPPdynOPT.3\n");fichier << "call.Mortalite.SPPdynOPT.3" << endl;
    Mortalite(listTempP, IND_T, eVar) ; //on g�n�re les Z dans out_Z_eit (il faut s'assurer que l'application du Ztemp est bloqu�e -> utilisation de ZoptSS3)
    //Rprintf("end.Mortalite.SPPdynOPT.3\n");fichier << "end.Mortalite.SPPdynOPT.3" << endl;

//       for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg7;
//                ggg7 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step7T1 " << ggg7.str() << endl;
//
//            }
//        }
    //Rprintf("call.DynamicPop.SPPdynOPT.3\n");fichier << "call.DynamicPop.SPPdynOPT.3" << endl;
    DynamicPop(listTempP, IND_T, eVar,true);
    //Rprintf("end.DynamicPop.SPPdynOPT.3\n");fichier << "end.DynamicPop.SPPdynOPT.3" << endl;




//      for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg8;
//                ggg8 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step8T1 " << ggg8.str() << endl;
//
//            }
//        }


//fichier << "ST9.24" << endl;
//il faut maintenant en d�duire Z par esp�ce et modifier Zfix_e en cons�quence

	goon = false;

	for (int ind = 0 ; ind < N_SPPdynOPT ; ind++) {  //boucle sur les esp�ces dynamiques XSA ou SS3 : d�termination de alpha, effort associ�, r�conciliation avec Einterm, m�j Einterm
//fichier << "ST10.1" << endl;
      int NBI = length(getListElement(getListElement(listTempP, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1))), "modI"));
      double *Ztemp = REAL(getListElement(ZtempList, CHAR(STRING_ELT(sppList,SPPdynOPT[ind]-1)))) ;
//fichier << "ST22" << endl;
      if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==0)) {
//fichier << "ST221" << endl;
//fichier << "ST10.2" << endl;
            for (int i = 0 ; i < NBI ; i++) {
//fichier << "ST222" << endl;
                //Rprintf("IT %i time %i indiv %i\n",IT,IND_T,ind);
                //Rprintf("Z %f Ztmp %f diff ZZ%f\n", REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T], Ztemp[i], REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]);fichier << "ST223" << endl;
                //fichier << "ST10.3" << endl;
                if (fabs(REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i] = Ztemp[i] + lambda*(REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]);
//fichier << "ST10.4" << endl;
            }

       } else if ((Qvec[SPPdynOPT[ind]-1]==0) & (Svec[SPPdynOPT[ind]-1]==1)) {
//fichier << "ST221" << endl;
//fichier << "ST10.2" << endl;
            for (int i = 0 ; i < NBI ; i++) {
//fichier << "ST222" << endl;
                //Rprintf("IT %i time %i indiv %i\n",IT,IND_T,ind);
                //Rprintf("Z %f Ztmp %f diff ZZ%f\n", REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T], Ztemp[i], REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]);fichier << "ST223" << endl;
                //fichier << "ST10.3" << endl;
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(0*NBI)] = Ztemp[i+(0*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_G1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]);

                //Rprintf("IT %i time %i indiv %i\n",IT,IND_T,ind);
                //Rprintf("Z %f Ztmp %f diff ZZ%f\n", REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T], Ztemp[i], REAL(VECTOR_ELT(out_Z_eit,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i]);fichier << "ST223" << endl;
                //fichier << "ST10.3" << endl;
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(1*NBI)] = Ztemp[i+(1*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_G2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]);
//fichier << "ST10.4" << endl;
            }

       } else if ((Qvec[SPPdynOPT[ind]-1]==1) & (Svec[SPPdynOPT[ind]-1]==0)){ //SS3
//fichier << "ST224" << endl;
//fichier << "ST10.5" << endl;
            for (int i = 0 ; i < NBI ; i++) {
                    //fichier << "ST10.6" << endl;
//fichier << "ST225" << endl;
                //Rprintf("IT %i time %i indiv %i\n",IT,IND_T,ind);
                //Rprintf("diffZZ S1M1 %f \n", REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]);fichier << "ST226" << endl;
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(0*NBI)] = Ztemp[i+(0*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(0*NBI)]);

                //Rprintf("diffZZ S1M2 %f \n", REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(1*NBI)] = Ztemp[i+(1*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(1*NBI)]);

                //Rprintf("diffZZ S1M3 %f \n", REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(2*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(2*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(2*NBI)] = Ztemp[i+(2*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(2*NBI)]);

                //Rprintf("diffZZ S1M4 %f \n", REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(3*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(3*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(3*NBI)] = Ztemp[i+(3*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S1M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(3*NBI)]);

                //Rprintf("diffZZ S2M1 %f \n", REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(4*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(4*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(4*NBI)] = Ztemp[i+(4*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(4*NBI)]);

                //Rprintf("diffZZ S2M2 %f \n", REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(5*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(5*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(5*NBI)] = Ztemp[i+(5*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(5*NBI)]);

                //Rprintf("diffZZ S2M3 %f \n", REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(6*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(6*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(6*NBI)] = Ztemp[i+(6*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(6*NBI)]);

                //Rprintf("diffZZ S2M4 %f \n", REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(7*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(7*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(7*NBI)] = Ztemp[i+(7*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S2M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(7*NBI)]);

                //Rprintf("diffZZ S3M1 %f \n", REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(8*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(8*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(8*NBI)] = Ztemp[i+(8*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(8*NBI)]);

                //Rprintf("diffZZ S3M2 %f \n", REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(9*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(9*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(9*NBI)] = Ztemp[i+(9*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(9*NBI)]);

                //Rprintf("diffZZ S3M3 %f \n", REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(10*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(10*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(10*NBI)] = Ztemp[i+(10*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(10*NBI)]);

                //Rprintf("diffZZ S3M4 %f \n", REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(11*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(11*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(11*NBI)] = Ztemp[i+(11*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S3M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(11*NBI)]);

                //Rprintf("diffZZ S4M1 %f \n", REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(12*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(12*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(12*NBI)] = Ztemp[i+(12*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M1,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(12*NBI)]);

                //Rprintf("diffZZ S4M2 %f \n", REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(13*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(13*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(13*NBI)] = Ztemp[i+(13*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M2,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(13*NBI)]);

                //Rprintf("diffZZ S4M3 %f \n", REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(14*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(14*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(14*NBI)] = Ztemp[i+(14*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M3,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(14*NBI)]);

                //Rprintf("diffZZ S4M4 %f \n", REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(15*NBI)]);
                if (fabs(REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(15*NBI)]) > diffZmax) goon = true; //on continue si l'une des valeurs est sup�rieurs � diffZiMax
                Ztemp[i+(15*NBI)] = Ztemp[i+(15*NBI)] + lambda*(REAL(VECTOR_ELT(out_Z_eit_S4M4,SPPdynOPT[ind]-1))[i+NBI*IND_T] - Ztemp[i+(15*NBI)]);
//fichier << "ST10.7" << endl;
            }

      }
	}


//	     for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg9;
//                ggg9 << g_effort1FM[ind_f + nbF*ind_m];
//
//                fichier << "effort_step9T1 " << ggg9.str() << endl;
//
//            }
//        }


    }

}
//reconcilSPP = duplicate(reconcilSPP_copy);
//UNPROTECT(1);

//free_vector(EintermList_fm,1,nbF*nbMe*N_SPPdynOPT);

}

//fichier << "ST23" << endl;

//fichier << "ST10.8" << endl;
        double *g_effort1FM_G = REAL(getListElement(FList, "effort1_f_m"));
//        double *g_effort1FM_Gcopy = REAL(duplicate(getListElement(FList, "effort1_f_m")));
        double *g_effort1F_G = REAL(getListElement(FList, "effort1_f"));
//        double *g_effort1F_Gcopy = REAL(duplicate(getListElement(FList, "effort1_f")));
        //double *g_nbTripFM_G = REAL(getListElement(FList, "nbTrip_f_m"));
        double *g_nbTripF_G = REAL(getListElement(FList, "nbTrip_f"));
        double *g_nbvFM_G = REAL(getListElement(FList, "nbv_f_m"));
        double *g_nbvF_G = REAL(getListElement(FList, "nbv_f"));
        double *g_effort2FM_G = REAL(getListElement(FList, "effort2_f_m"));
//        double *g_effort1FM_Gcopy = REAL(duplicate(getListElement(FList, "effort1_f_m")));
        double *g_effort2F_G = REAL(getListElement(FList, "effort2_f"));
//fichier << "ST10.9" << endl;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
                //fichier << "ST11.1" << endl;
        g_effort1F_G[ind_f] = 0.0;
        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
                g_effort1FM_G[ind_f + nbF*ind_m] = Einterm_fm[1 + ind_f + nbF*ind_m] / (g_effort2FM_G[ind_f + nbF*ind_m] * g_nbvFM_G[ind_f + nbF*ind_m]) ;
                g_nbTripFM[ind_f + nbF*ind_m] = g_effort1FM[ind_f + nbF*ind_m] ;
                if (!ISNA(g_effort1FM_G[ind_f + nbF*ind_m])){
                     g_effort1F_G[ind_f] = g_effort1F_G[ind_f] + g_effort1FM_G[ind_f + nbF*ind_m] * g_effort2FM_G[ind_f + nbF*ind_m] * g_nbvFM_G[ind_f + nbF*ind_m] ;
                    //fichier << "ST11.2" << endl;
                }
        }
        g_effort1F_G[ind_f] = g_effort1F_G[ind_f] / (g_nbvF_G[ind_f] * g_effort2F_G[ind_f]) ;
        g_nbTripF_G[ind_f] = g_effort1F_G[ind_f] ;
        }

//        for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//            if ((ind_t==1) & (ind_f==0)) {
//
//                std::stringstream ggg6;
//                ggg6 << g_effort1FM_G[ind_f + nbF*ind_m];
//
//                fichier << "effort_step6T1 " << ggg6.str() << endl;
//
//            }
//        }




//  free_matrix(q,1,2,1,1);
//	free_vector(z,1,2);
//	free_vector(x,1,1);
//  free_vector(multF,1,nbF+1);

 UNPROTECT(6);


    ZoptSS3 = false;
	return 0;

//fichier.close();
}

}
