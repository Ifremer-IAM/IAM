// #include <stdlib.h> // loaded with R.h
// #include <stdio.h> // loaded with R.h
// #include <time.h>
// #include <vector> // loaded with R.h
// #include <math.h> // loaded with R.h
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

//------------------------------------------------------------------------------------
//constructeur de la classe Param (voir 'param.h' pour les descriptions de variables)
//------------------------------------------------------------------------------------

BioEcoPar::BioEcoPar(SEXP listInput /* object@input */, SEXP listSpec /* object@specific */, SEXP listStochastic /* object@stochastic */,
                     SEXP listScen /* object@scenario */, SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, /*SEXP Bootstrp, SEXP nbBootstrp, */ // TODO : remove unused arg
                     SEXP GestInd, SEXP mOth, SEXP bounds, SEXP TACL, SEXP TACtot, SEXP FBAR, /*SEXP othSpSup,*/ SEXP effSup, SEXP GestParam, /*SEXP EcoDcf,*/
                     SEXP persCalc, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyFL, SEXP Ftarg, SEXP W_Ftarg, SEXP MeanRec_Ftarg,
                     SEXP parBHV, SEXP parQEX,
                     SEXP tacCTRL, SEXP stochPrice, SEXP updateE, SEXP parOQD, int VERBOSE, int force_T)
{

//ofstream fichier("C:\\Users\\fbriton\\Dropbox\\These\\IAM_Dvt\\test.txt", ios::out | ios::trunc);
//fichier << "D�but " << endl;


PROTECT(listSpec);
PROTECT(parOQD); //+1
effortIni = REAL(getListElement(getListElement(listInput, "Fleet"), "nbds_f"));
effort1Ini = REAL(getListElement(getListElement(listInput, "Fleet"), "effort1_f"));
PROTECT_INDEX ipx_list;
PROTECT_INDEX ipx_FList;
PROTECT_INDEX ipx_eVar;
PROTECT_INDEX ipx_eStatVar;
PROTECT_INDEX ipx_fVar;
PROTECT_INDEX ipx_list_copy;
PROTECT_INDEX ipx_FList_copy;
PROTECT_INDEX ipx_eVar_copy;
PROTECT_INDEX ipx_eStatVar_copy;
PROTECT_INDEX ipx_fVar_copy;

if(VERBOSE){Rprintf("Step 0 | ");}

PROTECT_WITH_INDEX(list = duplicate(listInput),&ipx_list);
PROTECT_WITH_INDEX(list_copy = duplicate(list),&ipx_list_copy);
PROTECT_WITH_INDEX(FList = getListElement(list, "Fleet"),&ipx_FList);
PROTECT_WITH_INDEX(FList_copy = getListElement(list_copy, "Fleet"),&ipx_FList_copy);


PROTECT(listQR = getListElement(parOQD, "listQR")); //+1
PROTECT(listQR_f = getListElement(parOQD, "listQR_f")); //+1
activeQR = INTEGER(getListElement(parOQD, "activeQR"))[0];

PROTECT(sppList = getListElement(listSpec, "Species"));
PROTECT(sppListStat = getListElement(listSpec, "StaticSpp"));
PROTECT(fleetList = getListElement(listSpec, "Fleet"));
PROTECT(metierList = getListElement(listSpec, "MetierEco"));
PROTECT(metierListEco = getListElement(listSpec, "MetierEco"));
PROTECT(namDC = getListElement(listSpec, "Ages"));
PROTECT(t_init = getListElement(listSpec, "t_init"));
PROTECT(times = getListElement(listSpec, "times"));
PROTECT(Q = getListElement(listSpec, "Q"));
PROTECT(S = getListElement(listSpec, "S")); //25

Qvec = INTEGER(Q);
Svec = INTEGER(S);

Zoptim = effortIni;
FOTHoptim = effortIni;
Zoptim_use = false;
FOTHoptim_use = false;


if(VERBOSE){Rprintf("Step 0.1 | ");}
//Rprintf("step0\n");fichier << "step0 " << endl;

nbT = INTEGER(getListElement(listSpec, "NbSteps"))[0];
nbF = length(fleetList);
nbM = length(metierList);
nbMe = length(metierListEco);
nbE = length(sppList);
nbEstat = length(sppListStat);
nbP = 0;

SEXP dimCstF, dimCstFM;
int *dCF, *dCFM;
Rf_protect(dimCstF = Rf_allocVector(INTSXP, 4));
Rf_protect(dimCstFM = Rf_allocVector(INTSXP, 4));

dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;

drCopy = REAL(dr)[0];

int conform = 0;
typeGest = 0;
SEXP FList_copy, list_copy, eVar_copy, eStatVar_copy, fVar_copy;

recType1 = INTEGER(RecType1);//vecteur d'entiers de longueur nbE
recType2 = INTEGER(RecType2);//vecteur d'entiers de longueur nbE
recType3 = INTEGER(RecType3);//vecteur d'entiers de longueur nbE
// boot = INTEGER(Bootstrp)[0]; // TODO : remove unused arg
// nbBoot = INTEGER(nbBootstrp)[0];
boolQ = true;// param�tre vou� � rester fixe -> on calculera toujours la capturabilit� afin de moduler la mortalit� en fonction de l'effort de p�che
constMM = true; //on calcule la capturabilit� via l'effort par flottille (incompatibilit� des niveaux m�tiers entre bio et �co)
fUpdate = true;    // � t=0, on remet � jour
dUpdate = true;    //
cUpdate = true;    //
pflex = !isNull(getListElement(list, "Market"));  //Rprintf("pflex = %d \n",pflex) ; //
eUpdate = true;    //

if (pflex) {
    if(VERBOSE){Rprintf("pflex | ");}
    PROTECT(sppListAll = getListElement(getListElement(list, "Market"),"modE")); //PrintValue(sppListAll);
    PROTECT(pList = getListElement(getListElement(list, "Market"),"modP"));
    nbP = length(pList);
} else {
    PROTECT(sppListAll = getListElement(listSpec, "AllSpp")); // TODO : this should not require Market sheet to be input
}
nbEall = length(sppListAll);

scen = INTEGER(Scenarii)[0];
bhv_active = INTEGER(getListElement(parBHV, "active"))[0];

tolVarTACinf_CPP = REAL(getListElement(tacCTRL, "tolVarTACinf"))[0];
tolVarTACsup_CPP = REAL(getListElement(tacCTRL, "tolVarTACsup"))[0];
corVarTACval_CPP = REAL(getListElement(tacCTRL, "corVarTACval"))[0];
corVarTACnby_CPP = INTEGER(getListElement(tacCTRL, "corVarTACnby"))[0];
Blim_CPP = REAL(getListElement(tacCTRL, "Blim"))[0];
Bmax_CPP = REAL(getListElement(tacCTRL, "Bmax"))[0];
Blim_trigger = INTEGER(getListElement(tacCTRL, "BlimTrigger"))[0];

//gestyp = INTEGER(getListElement(tacCTRL, "typeMng"))[0]; //gestyp est pris prioritairement dans tacCTRL, sinon par l'interface
//if (ISNA(gestyp))
gestyp = INTEGER(GestParam)[5] + 1; //gestyp=1(+) ou 2(x)

// PROTECT(othSpSupList = othSpSup);
PROTECT(effSupMat = effSup);

gestInd = INTEGER(GestInd)[0];
PROTECT(m_fm = duplicate(mFM)); if (length(m_fm)!=nbF*nbM) error("Check dimension of array 'mFleetMetier'!!\n");
PROTECT(m_oth = duplicate(mOth)); if (length(m_oth)!=nbE) error("Check dimension of array 'mOth'!!\n");
X1 = REAL(bounds)[0];
X2 = REAL(bounds)[1];

// TAC_glob = REAL(TAC);  //� corriger
Fbar_trgt = REAL(FBAR);  //� corriger
PROTECT(TACbyF = TACbyFL);  //� corriger
PROTECT(TAC = TACL); // TODO :pourquoi juste un changement de nom la en fait ?
TAC_glob = REAL(TACtot);  // TODO : TAC est une list alors qu'avant il etait unique.

PROTECT(inpFtarg = Ftarg);
PROTECT(inpW_Ftarg = W_Ftarg);
PROTECT(inpMeanRec_Ftarg = MeanRec_Ftarg);

nbEQuota = length(getAttrib(TAC, R_NamesSymbol));
PROTECT(sppListQ = getAttrib(TAC, R_NamesSymbol));

nbEQuotaMarket = length(getAttrib(VECTOR_ELT(parQEX,1), R_NamesSymbol));
PROTECT(sppListQM = getAttrib(VECTOR_ELT(parQEX,1), R_NamesSymbol));
nbEQuotaMarket_dyn = nbEQuota;
PROTECT(sppListQM_dyn = sppListQ);

PROTECT(Qholdings = VECTOR_ELT(parQEX,8));

//tac_ctrl = tacCTRL;
recList = getListElement(tacCTRL, "recList");//PrintValue(recList);
recParamList = getListElement(tacCTRL, "recParamList");
ParamSPMlist = getListElement(tacCTRL, "ParamSPMList");
maxIter = INTEGER(getListElement(tacCTRL, "maxIter"))[0]; // TODO : only used in GestionF2
diffZmax = REAL(getListElement(tacCTRL, "diffZmax"))[0];
lambda = REAL(getListElement(tacCTRL, "lambda"))[0];
int t_stop = INTEGER(getListElement(tacCTRL, "t_stop"))[0];
ZoptSS3 = false;

eTemp = INTEGER(GestParam)[0];//2;
var = INTEGER(GestParam)[1];//1;
trgt = INTEGER(GestParam)[2];//1;  si NA, objectif=biomasse de l'ann�e suivante --> trgt = 999
if (ISNA(trgt)) trgt = 999;
delay = INTEGER(GestParam)[3];//2;
upd = INTEGER(GestParam)[4];//2;

//if (eTemp<nbE) {
//  if (Qvec[eTemp]==0) {
//    Ztemp = NRvector(1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI")));
//  } else {
//    Ztemp = NRvector(1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"))*16); //les Z aux �ges par morph et saison
//  }
//} else {
//  Ztemp = NRvector(1,2);
//}
Etemp = NRvector(1,nbF); //effort moyen par flottille en nombre de mar�es
expEff = 1.2; //1;  //facteur d'expansion de l'effort maximal par flottille autoris� dans le cadre de l'optimisation GestionF2 et QuotaExch

Einterm_fm = NRvector(1,nbF*nbMe);
Einterm_fm_copy = NRvector(1,nbF*nbMe);
EffsupTMP_fm = NRvector(1,nbF*nbMe);
multFOTHinterm_e = NRvector(1,nbE);
PROTECT(reconcilSPP = alloc3DArray(STRSXP,nbF,nbMe,nbT));
PROTECT(reconcilSPP_copy = alloc3DArray(STRSXP,nbF,nbMe,nbT));
//on l'initialise
for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
for (int ind_m = 0 ; ind_m < nbMe ; ind_m++)
for (int ind_t = 0 ; ind_t < nbT ; ind_t++) {
SET_STRING_ELT(reconcilSPP, ind_f + nbF*ind_m + nbF*nbMe*ind_t, NA_STRING); //SET_STRING_ELT(reconcilSPP, ind_f + nbF * ind_m, mkChar(STRING_ELT(sppList,eTemp)));
SET_STRING_ELT(reconcilSPP_copy, ind_f + nbF*ind_m + nbF*nbMe*ind_t, NA_STRING);
}
PROTECT(ZtempList = getListElement(tacCTRL, "Ztemp"));


SRInd = INTEGER(SRind);

bool door = true;

//il faut initialiser 'eVar' dans lequel on integrera toutes les variables intermediaires a decliner par espece
SEXP eltE;
PROTECT_WITH_INDEX(eVar = allocVector(VECSXP, nbE),&ipx_eVar);
if (nbE>0) {
    setAttrib(eVar, R_NamesSymbol, sppList);
    for (int e = 0 ; e < nbE ; e++) {
      PROTECT(eltE = allocVector(VECSXP,250)); //ex 62
      SET_VECTOR_ELT(eVar, e, eltE);
      UNPROTECT(1);
    }
}
PROTECT_WITH_INDEX(eVar_copy = duplicate(eVar),&ipx_eVar_copy); //19

PROTECT_WITH_INDEX(eStatVar = allocVector(VECSXP, nbEstat),&ipx_eStatVar);
if (nbEstat>0) {
    setAttrib(eStatVar, R_NamesSymbol, sppListStat);
    for (int e = 0 ; e < nbEstat ; e++) {
        PROTECT(eltE = allocVector(VECSXP,10));
        SET_VECTOR_ELT(eStatVar, e, eltE);
        UNPROTECT(1);
    }
}
PROTECT_WITH_INDEX(eStatVar_copy = duplicate(eStatVar),&ipx_eStatVar_copy); //21


if(VERBOSE){Rprintf("Protect out | ");}
{
PROTECT_WITH_INDEX(fVar = allocVector(VECSXP, 34),&ipx_fVar); //32= rtbsIni_f & 33=rtbsIni_f_m & 33=ETini_f_m
PROTECT_WITH_INDEX(fVar_copy = duplicate(fVar),&ipx_fVar_copy);

PROTECT(out_F_fmi = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_G1 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_G2 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_G1 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_G2 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_G1 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_G2 = allocVector(VECSXP, nbE));
PROTECT(out_SRmod = allocVector(VECSXP, nbE));
PROTECT(out_Ystat = allocVector(VECSXP, nbEstat));//PROTECT(out_N_eitQ = allocVector(VECSXP, nbE));
PROTECT(out_Lstat = allocVector(VECSXP, nbEstat));//PROTECT(out_F_itQ = allocVector(VECSXP, nbE));
PROTECT(out_Dstat = allocVector(VECSXP, nbEstat));//PROTECT(out_SSB_etQ = allocVector(VECSXP, nbE));
PROTECT(out_PQuot_et = allocVector(VECSXP, nbEQuotaMarket));
PROTECT(out_QuotaTrade_fe = allocVector(VECSXP, nbEQuotaMarket_dyn));
PROTECT(out_diffLQ = allocVector(VECSXP, nbEQuotaMarket_dyn));
PROTECT(out_PQuot_temp = allocVector(VECSXP, nbEQuotaMarket_dyn));

PROTECT(out_Fbar_et = allocVector(VECSXP, nbE));
PROTECT(out_N_eit = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_G1 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_G2 = allocVector(VECSXP, nbE));
PROTECT(out_B_et = allocVector(VECSXP, nbE));
PROTECT(out_SSB_et = allocVector(VECSXP, nbE));
PROTECT(out_C_efmit = allocVector(VECSXP, nbE));
PROTECT(out_C_eit = allocVector(VECSXP, nbE));
PROTECT(out_C_efmit_G1 = allocVector(VECSXP, nbE));
PROTECT(out_C_eit_G1 = allocVector(VECSXP, nbE));
PROTECT(out_C_efmit_G2 = allocVector(VECSXP, nbE));
PROTECT(out_C_eit_G2 = allocVector(VECSXP, nbE));
PROTECT(out_Y_efmit = allocVector(VECSXP, nbE));
PROTECT(out_Y_eit = allocVector(VECSXP, nbE));
PROTECT(out_D_efmit = allocVector(VECSXP, nbE));
PROTECT(out_L_efmit = allocVector(VECSXP, nbE));
PROTECT(out_L_efmct = allocVector(VECSXP, nbE));
PROTECT(out_L_eit = allocVector(VECSXP, nbE));
PROTECT(out_L_et = allocVector(VECSXP, nbEall));
PROTECT(out_L_pt = allocVector(VECSXP, nbP));//43

PROTECT(out_oqDstat = allocVector(VECSXP, nbEstat));
PROTECT(out_oqD_eft = allocVector(VECSXP, nbE));
PROTECT(out_oqD_et = allocVector(VECSXP, nbE));

PROTECT(intermBIOMspict = allocVector(VECSXP, nbE)); //+1 ajout pour ins�rer les 16 valeurs de biomasses interm�diaires de chacune des esp�ces SPiCT lors de l'�valuation des Bt+1 et des captures

PROTECT(out_F_fmi_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_F_fmi_S4M4 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi_S4M4 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_FRWT_fmi_S4M4 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_FDWT_fmi_S4M4 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit_S4M4 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S1M1 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S1M2 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S1M3 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S1M4 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S2M1 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S2M2 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S2M3 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S2M4 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S3M1 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S3M2 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S3M3 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S3M4 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S4M1 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S4M2 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S4M3 = allocVector(VECSXP, nbE));
PROTECT(out_N_eit_S4M4 = allocVector(VECSXP, nbE)); //+64 = 107

PROTECT(out_P_t = allocVector(VECSXP, nbE));
PROTECT(out_Pstat = allocVector(VECSXP, nbEstat));
// PROTECT(out_Eco = allocVector(VECSXP, 69)); // TODO useless ?
PROTECT(out_EcoDCF = allocVector(VECSXP, 60));
PROTECT(out_effort = allocVector(VECSXP, 6)); //nbv_f, effort1_f, effort2_f, nbv_f_m, effort1_f_m, effort2_f_m

PROTECT(mu_nbds = allocVector(REALSXP, nbT)); //il reste la mise en forme � op�rer
PROTECT(mu_nbv = allocVector(REALSXP, nbT));
PROTECT(out_typeGest = allocVector(INTSXP, nbT)); // TODO : Florence remove this, is it useless ?

PROTECT(out_Ytot_fm = NEW_NUMERIC(nbF*nbMe*nbT));
PROTECT(out_DD_efmi = allocVector(VECSXP, nbE));
PROTECT(out_DD_efmc = allocVector(VECSXP, nbE));
PROTECT(out_LD_efmi = allocVector(VECSXP, nbE));
PROTECT(out_LD_efmc = allocVector(VECSXP, nbE));
PROTECT(out_statDD_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLD_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLDst_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLDor_efm = allocVector(VECSXP, nbEstat));
}

if(VERBOSE){Rprintf("Step 0.3 | ");}

int *typegest = INTEGER(out_typeGest);
double *mu_nbds_t = REAL(mu_nbds); for (int i=0; i<nbT; i++) mu_nbds_t[i] = 0.0; //initialisation
double *mu_nbv_t = REAL(mu_nbv); for (int i=0; i<nbT; i++) mu_nbv_t[i] = 0.0;    //

double *mpond_fm = REAL(m_fm);
double *mpond_oth = REAL(m_oth);

//on n'oublie pas de composer l'objet de sortie d�crivant les variables 'nbv' et 'nbds'
//SEXP NBVF, NBVFM, NBDSF, NBDSFM, dnmsF, dnmsFM, nmsEF;
PROTECT(NBVF = allocMatrix(REALSXP,nbF,nbT));
PROTECT(NBVFM = alloc3DArray(REALSXP,nbF,nbMe,nbT));
PROTECT(NBDSF = allocMatrix(REALSXP,nbF,nbT));
PROTECT(NBDSFM = alloc3DArray(REALSXP,nbF,nbMe,nbT));
PROTECT(EFF2F = allocMatrix(REALSXP,nbF,nbT));
PROTECT(EFF2FM = alloc3DArray(REALSXP,nbF,nbMe,nbT));
PROTECT(dnmsF = allocVector(VECSXP,2));
PROTECT(dnmsFM = allocVector(VECSXP,3));


PROTECT(dnmsIter = allocVector(VECSXP,2));
int itmaxQ = INTEGER(VECTOR_ELT(parQEX,7))[0];
PROTECT(itListQ = NEW_INTEGER(itmaxQ));
for (int i=0; i<itmaxQ; i++) INTEGER(itListQ)[i] = i;


SET_VECTOR_ELT(dnmsF, 0, fleetList); SET_VECTOR_ELT(dnmsF, 1, times);
SET_VECTOR_ELT(dnmsFM, 0, fleetList); SET_VECTOR_ELT(dnmsFM, 1, metierListEco); SET_VECTOR_ELT(dnmsFM, 2, times);
SET_VECTOR_ELT(dnmsIter, 0, itListQ); SET_VECTOR_ELT(dnmsIter, 1, times);
setAttrib(NBVF, R_DimNamesSymbol, dnmsF); setAttrib(NBVFM, R_DimNamesSymbol, dnmsFM);
setAttrib(NBVF, install("DimCst"), dimCstF); setAttrib(NBVFM, install("DimCst"), dimCstFM);
setAttrib(NBDSF, R_DimNamesSymbol, dnmsF); setAttrib(NBDSFM, R_DimNamesSymbol, dnmsFM);
setAttrib(NBDSF, install("DimCst"), dimCstF); setAttrib(NBDSFM, install("DimCst"), dimCstFM);
setAttrib(EFF2F, R_DimNamesSymbol, dnmsF); setAttrib(EFF2FM, R_DimNamesSymbol, dnmsFM);
setAttrib(EFF2F, install("DimCst"), dimCstF); setAttrib(EFF2FM, install("DimCst"), dimCstFM);

PROTECT(out_allocEff_fm = alloc3DArray(REALSXP,nbF,nbMe,nbT));
setAttrib(out_allocEff_fm, R_DimNamesSymbol, dnmsFM);
setAttrib(out_allocEff_fm, install("DimCst"), dimCstFM);
// Intitialisation avec mfm mais modifie dans le module QuotaExchange si appele
for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
for (int ind_m = 0 ; ind_m < nbMe ; ind_m++)
for (int ind_t = 0 ; ind_t < nbT ; ind_t++) {
        REAL(out_allocEff_fm)[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = REAL(m_fm)[ind_f + nbF*ind_m];}
//PrintValue(out_allocEff_fm);

setAttrib(reconcilSPP, R_DimNamesSymbol, dnmsFM); setAttrib(reconcilSPP, install("DimCst"), dimCstFM);
setAttrib(reconcilSPP_copy, R_DimNamesSymbol, dnmsFM); setAttrib(reconcilSPP_copy, install("DimCst"), dimCstFM);

double *NBVf = REAL(NBVF);
double *NBVfm = REAL(NBVFM);
double *NBDSf = REAL(NBDSF);
double *NBDSfm = REAL(NBDSFM);
double *EFF2f = REAL(EFF2F);
double *EFF2fm = REAL(EFF2FM);

SEXP ans_PQuot_et, pQuotaIni;
if (nbEQuotaMarket>0) {
 if(VERBOSE){Rprintf("Step 0.3.1");}
 setAttrib(out_PQuot_et, R_NamesSymbol, sppListQM);
 for (int e = 0 ; e < nbEQuotaMarket ; e++) {
    PROTECT(ans_PQuot_et = NEW_NUMERIC(nbT));
    setAttrib(ans_PQuot_et, R_NamesSymbol, times);
    PROTECT(pQuotaIni = VECTOR_ELT(VECTOR_ELT(parQEX,1),e));
    for (int ind_t=0; ind_t < nbT ; ind_t++) REAL(ans_PQuot_et)[ind_t] = REAL(pQuotaIni)[ind_t];
    SET_VECTOR_ELT(out_PQuot_et, e, ans_PQuot_et);
    UNPROTECT(2);
 }
 if(VERBOSE){Rprintf(" | ");}
}

SEXP ans_QuotaTrade_fe;
if (nbEQuotaMarket_dyn>0) {
 if(VERBOSE){Rprintf("Step 0.3.2");}
 setAttrib(out_QuotaTrade_fe, R_NamesSymbol, sppListQM_dyn);
 for (int e = 0 ; e < nbEQuotaMarket_dyn ; e++) {
    PROTECT(ans_QuotaTrade_fe = allocMatrix(REALSXP,nbF,nbT));
    setAttrib(ans_QuotaTrade_fe, R_DimNamesSymbol, dnmsF);
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
        for (int ind_t = 0 ; ind_t < nbT ; ind_t++) REAL(ans_QuotaTrade_fe)[ind_f + nbF*ind_t] = 0.0;}
    SET_VECTOR_ELT(out_QuotaTrade_fe, e, ans_QuotaTrade_fe);
    UNPROTECT(1);
 }
 if(VERBOSE){Rprintf(" | ");}
}

SEXP ans_diffLQ;
if (nbEQuotaMarket_dyn>0) {
 if(VERBOSE){Rprintf("Step 0.3.3");}
 setAttrib(out_diffLQ, R_NamesSymbol, sppListQM_dyn);
 for (int e = 0 ; e < nbEQuotaMarket_dyn ; e++) {
    PROTECT(ans_diffLQ = allocMatrix(REALSXP,itmaxQ,nbT));
    setAttrib(ans_diffLQ, R_DimNamesSymbol, dnmsIter);
    for (int i = 0 ; i < itmaxQ ; i++){
        for (int ind_t = 0 ; ind_t < nbT ; ind_t++) REAL(ans_diffLQ)[i + itmaxQ*ind_t] = NA_REAL;}
    SET_VECTOR_ELT(out_diffLQ, e, ans_diffLQ);
    UNPROTECT(1);
 }
 if(VERBOSE){Rprintf(" | ");}
}

SEXP ans_PQuot_temp;
if (nbEQuotaMarket_dyn>0) {
 if(VERBOSE){Rprintf("Step 0.3.4");}
 setAttrib(out_PQuot_temp, R_NamesSymbol, sppListQM_dyn);
 for (int e = 0 ; e < nbEQuotaMarket_dyn ; e++) {
    PROTECT(ans_PQuot_temp = allocMatrix(REALSXP,itmaxQ,nbT));
    setAttrib(ans_PQuot_temp, R_DimNamesSymbol, dnmsIter);
    for (int i = 0 ; i < itmaxQ ; i++){
        for (int ind_t = 0 ; ind_t < nbT ; ind_t++) REAL(ans_PQuot_temp)[i + itmaxQ*ind_t] = NA_REAL;}
    SET_VECTOR_ELT(out_PQuot_temp, e, ans_PQuot_temp);
    UNPROTECT(1);
 }
 if(VERBOSE){Rprintf(" | ");}
}

if(VERBOSE){Rprintf("Step 0.3.5");}
PROTECT(intermGoFish = allocMatrix(REALSXP,nbF,nbT));
setAttrib(intermGoFish, R_DimNamesSymbol, dnmsF);
setAttrib(intermGoFish, install("DimCst"), dimCstF);
for (int ind_t = 0 ; ind_t < nbT ; ind_t++){
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) REAL(intermGoFish)[ind_t*nbF + ind_f] = 1.0; //initialisation
}
if(VERBOSE){Rprintf(" | ");}

PROTECT(multPrice = allocVector(VECSXP, nbE+nbEstat));
SEXP ans_multPrice;
if ((nbE+nbEstat)>0) {
 if(VERBOSE){Rprintf("Step 0.3.6");}
 setAttrib(multPrice, R_NamesSymbol, sppListAll);
 for (int e = 0 ; e < (nbE+nbEstat) ; e++) {
    PROTECT(ans_multPrice = NEW_NUMERIC(nbT));
    setAttrib(ans_multPrice, R_NamesSymbol, times);
    for (int ind_t=0; ind_t < nbT ; ind_t++) REAL(ans_multPrice)[ind_t] = 1.0;
    SET_VECTOR_ELT(multPrice, e, ans_multPrice);
    UNPROTECT(1);
 }
 if(VERBOSE){Rprintf(" | ");}
}

if(VERBOSE){Rprintf("\nLoop :");}

for (int it = 0; (it < nbT) & (it < force_T) ; it++) {
    if(VERBOSE){Rprintf("\n========================== T = %d => \n", it);}

//Rprintf("ini1");fichier << "ini1" << endl;

//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 5));////Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));////Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 7));////Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 61));////Rprintf("Mort20.5\n");


//    if (it==0) 	{
////        //float Y=;
////        //float **P=&NRmatrix(1,4,1,3);
//        MinimizeF(NRmatrix(1,4,1,3), NRvector(1,4), 3, 0.0001);
//    }



if (nbE>0) {

    if (it>=1) {
        RecAlea(list, listStochastic, it, 1, recType1); //IMPORTANT : a effectuer AVANT la procedure d'optimisation
        RecAlea(list, listStochastic, it, 2, recType2);
        RecAlea(list, listStochastic, it, 3, recType3);
    }

        SRmod(list, listSR, it, TypeSR, SRInd); //important : a envoyer avant 'DynamicPop'

}
//Rprintf("C");fichier << "C" << endl;
if (scen & (it>=1)) Scenario(list, listScen, it); //modif MM 16/01/2012

//if (bhv_active & (it>=1)) FleetBehav(list, it, parBHV);


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//------------------------ Ajustement initial ---------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Rprintf("intro0\n");fichier << "intro0" << endl;
 //copie de sauvegarde car il faudra revenir aux valeurs initiales s'il y a une correction des ajustements

REPROTECT(list_copy = duplicate(list),ipx_list_copy);//Rprintf("intro0.1\n");
REPROTECT(FList_copy = getListElement(list_copy, "Fleet"),ipx_FList_copy); //Rprintf("intro0.2\n");




////Rprintf("ini2");fichier << "ini2" << endl;

//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 5));//Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));//Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 7));//Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 61));//Rprintf("Mort20.5\n");

////Rprintf("ini2");

//PrintValue(VECTOR_ELT(VECTOR_ELT(duplicate(eVar), 0), 5));//Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(duplicate(eVar), 0), 6));//Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(duplicate(eVar), 0), 7));//Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(duplicate(eVar), 0), 61));//Rprintf("Mort20.5\n");


//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar,0),0));//Rprintf("a0");
//PrintValue(VECTOR_ELT(eVar,0));//Rprintf("a1");//if (it>0) error("Unexpectedz scondition occurred");

////Rprintf("AAAAA3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));
//SEXP eVarDup;
//PROTECT(eVarDup = duplicate(eVar));
////Rprintf("AAAAA31\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVarDup, 0), 6));
//UNPROTECT(1);
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));
////Rprintf("AAAAA32\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVarDup, 0), 6));


REPROTECT(eVar_copy = duplicate(eVar),ipx_eVar_copy);//Rprintf("intro0.3\n");fichier << "intro0.3" << endl;
//SEXP jojo = duplicate(eVar);//Rprintf("intro0.2.1\n");
//REPROTECT(eVar_copy = eVar,ipx_eVar_copy);//Rprintf("intro0.3\n");

////Rprintf("ini3");

//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 5));//Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));//Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 7));//Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 61));//Rprintf("Mort20.5\n");


REPROTECT(eStatVar_copy = duplicate(eStatVar),ipx_eStatVar_copy);//Rprintf("intro0.4\n");fichier << "intro.4" << endl;
REPROTECT(fVar_copy = duplicate(fVar),ipx_fVar_copy);//Rprintf("intro0.5\n");fichier << "intro0.5" << endl;
//Rprintf("intro1\n");fichier << "intro1" << endl;


if ((delay<=it) & (gestInd==1) & (it>=1) & isNull(TACbyF) & (eTemp<nbE)) { //seulement si esp�ce dynamique

    if(VERBOSE){Rprintf(" | Gestion");}

//------------------------------------------------------------------------------------------ ajout updateE : d�but
int DELAY = INTEGER(updateE)[0];

if ((delay<=it) & (gestInd==1) & (DELAY>0)) { //DELAY = 1 -> on remet l'effort au niveau de l'instant initial

//on remet au niveau de l'instant pr�c�dent la mise en action du module Gestion

    double *nbdsFM3 = REAL(getListElement(FList, "effort1_f_m"));
    double *nbdsF3 = REAL(getListElement(FList, "effort1_f"));
    double *nbTripFM3 = REAL(getListElement(FList, "nbTrip_f_m"));
    double *nbTripF3 = REAL(getListElement(FList, "nbTrip_f"));
    double *nbvFM3 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF3 = REAL(getListElement(FList, "nbv_f"));

    if (DELAY>delay) DELAY=delay;

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){


        if (var==1) {nbdsF3[ind_f] = REAL(NBDSF)[ind_f + nbF*(DELAY-1)];
                     nbTripF3[ind_f] = REAL(NBDSF)[ind_f + nbF*(DELAY-1)];}
        if (var==2) nbvF3[ind_f] = REAL(NBVF)[ind_f + nbF*(DELAY-1)];

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) {nbdsFM3[ind_f+nbF*ind_m] = REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];
                         nbTripFM3[ind_f+nbF*ind_m] = REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];}
            if (var==2) nbvFM3[ind_f+nbF*ind_m] = REAL(NBVFM)[ind_f + nbF*ind_m + nbF*nbMe*(DELAY-1)];

        }
    }
    //Rprintf("AA");fichier << "AA" << endl;
    //PrintValue(getListElement(FList, "effort1_f_m"));
            int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    if (Qvec[eTemp]==0) {

            double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 44));
            for (int ag = 0; ag < nbi; ag++) Fothi2[ag + it*nbi] = Fothi2[ag + (DELAY-1)*nbi];

    } else {

            double *Fothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 116)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + it*nbi] = Fothi2_S1M1[ag + (DELAY-1)*nbi];
            double *Fothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 117)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + it*nbi] = Fothi2_S1M2[ag + (DELAY-1)*nbi];
            double *Fothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 118)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + it*nbi] = Fothi2_S1M3[ag + (DELAY-1)*nbi];
            double *Fothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 119)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + it*nbi] = Fothi2_S1M4[ag + (DELAY-1)*nbi];
            double *Fothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 120)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + it*nbi] = Fothi2_S2M1[ag + (DELAY-1)*nbi];
            double *Fothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 121)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + it*nbi] = Fothi2_S2M2[ag + (DELAY-1)*nbi];
            double *Fothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 122)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + it*nbi] = Fothi2_S2M3[ag + (DELAY-1)*nbi];
            double *Fothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 123)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + it*nbi] = Fothi2_S2M4[ag + (DELAY-1)*nbi];
            double *Fothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 124)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + it*nbi] = Fothi2_S3M1[ag + (DELAY-1)*nbi];
            double *Fothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 125)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + it*nbi] = Fothi2_S3M2[ag + (DELAY-1)*nbi];
            double *Fothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 126)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + it*nbi] = Fothi2_S3M3[ag + (DELAY-1)*nbi];
            double *Fothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 127)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + it*nbi] = Fothi2_S3M4[ag + (DELAY-1)*nbi];
            double *Fothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 128)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + it*nbi] = Fothi2_S4M1[ag + (DELAY-1)*nbi];
            double *Fothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 129)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + it*nbi] = Fothi2_S4M2[ag + (DELAY-1)*nbi];
            double *Fothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 130)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + it*nbi] = Fothi2_S4M3[ag + (DELAY-1)*nbi];
            double *Fothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 131)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + it*nbi] = Fothi2_S4M4[ag + (DELAY-1)*nbi];

            double *FRWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 176));
            double *FRWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 177));
            double *FRWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 178));
            double *FRWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 179));
            double *FRWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 180));
            double *FRWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 181));
            double *FRWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 182));
            double *FRWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 183));
            double *FRWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 184));
            double *FRWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 185));
            double *FRWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 186));
            double *FRWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 187));
            double *FRWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 188));
            double *FRWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 189));
            double *FRWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 190));
            double *FRWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 191));

            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + it*nbi] = FRWTothi2_S1M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + it*nbi] = FRWTothi2_S1M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + it*nbi] = FRWTothi2_S1M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + it*nbi] = FRWTothi2_S1M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + it*nbi] = FRWTothi2_S2M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + it*nbi] = FRWTothi2_S2M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + it*nbi] = FRWTothi2_S2M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + it*nbi] = FRWTothi2_S2M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + it*nbi] = FRWTothi2_S3M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + it*nbi] = FRWTothi2_S3M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + it*nbi] = FRWTothi2_S3M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + it*nbi] = FRWTothi2_S3M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + it*nbi] = FRWTothi2_S4M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + it*nbi] = FRWTothi2_S4M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + it*nbi] = FRWTothi2_S4M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + it*nbi] = FRWTothi2_S4M4[ag + (DELAY-1)*nbi];

            double *FDWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 208));
            double *FDWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 209));
            double *FDWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 210));
            double *FDWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 211));
            double *FDWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 212));
            double *FDWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 213));
            double *FDWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 214));
            double *FDWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 215));
            double *FDWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 216));
            double *FDWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 217));
            double *FDWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 218));
            double *FDWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 219));
            double *FDWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 220));
            double *FDWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 221));
            double *FDWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 222));
            double *FDWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 223));

            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + it*nbi] = FDWTothi2_S1M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + it*nbi] = FDWTothi2_S1M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + it*nbi] = FDWTothi2_S1M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + it*nbi] = FDWTothi2_S1M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + it*nbi] = FDWTothi2_S2M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + it*nbi] = FDWTothi2_S2M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + it*nbi] = FDWTothi2_S2M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + it*nbi] = FDWTothi2_S2M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + it*nbi] = FDWTothi2_S3M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + it*nbi] = FDWTothi2_S3M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + it*nbi] = FDWTothi2_S3M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + it*nbi] = FDWTothi2_S3M4[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + it*nbi] = FDWTothi2_S4M1[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + it*nbi] = FDWTothi2_S4M2[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + it*nbi] = FDWTothi2_S4M3[ag + (DELAY-1)*nbi];
            for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + it*nbi] = FDWTothi2_S4M4[ag + (DELAY-1)*nbi];

    }

}
//------------------------------------------------------------------------------------------ ajout updateE : fin
if(VERBOSE){Rprintf("\n Before module");}

////Rprintf("AAAAAAA\n");
Gestion(list, it, VERBOSE);

if(VERBOSE){Rprintf("after Module");}


////Rprintf("bb");
//PrintValue(getListElement(getListElement(list, "Fleet"), "nbds_f"));
//PrintValue(getListElement(FList, "nbds_f"));
////PrintValue(mu_nbds);
////PrintValue(mu_nbv);

////Rprintf("D");
//on met � jour les variables sur lesquelles op�re le multiplicateur
////Rprintf("E");
////PrintValue(mu_nbds);
////PrintValue(mu_nbv);
    double *mu_nbds_t2 = REAL(mu_nbds);
    double *mu_nbv_t2 = REAL(mu_nbv);
    double *nbdsFM2 = REAL(getListElement(FList, "effort1_f_m"));
    double *nbdsF2 = REAL(getListElement(FList, "effort1_f"));
    double *nbvFM2 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF2 = REAL(getListElement(FList, "nbv_f"));
    double *eff2FM2 = REAL(getListElement(FList, "effort2_f_m"));
    double *eff2F2 = REAL(getListElement(FList, "effort2_f"));

// ATTENTION : dor�navant, on doit avoir POUR CHAQUE FLOTTILLE des niveaux m�tiers exhaustifs ie sum_m ind_fm = ind_f --> il faudra des proc�dures de v�rifications dans les routines d'importation


     for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double countEff = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m] + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m]*(1 + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (var==2) {
                if (gestyp==1) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m] + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m]*(1 + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (!ISNA(nbdsFM2[ind_f+nbF*ind_m]) & !ISNA(nbvFM2[ind_f+nbF*ind_m]) & !ISNA(eff2FM2[ind_f+nbF*ind_m])) countEff = countEff + nbdsFM2[ind_f+nbF*ind_m]*nbvFM2[ind_f+nbF*ind_m]*eff2FM2[ind_f+nbF*ind_m];
        }

        if (var==1) nbdsF2[ind_f] = fmax2(countEff/(nbvF2[ind_f]*eff2F2[ind_f]),0.0); //NBDSf
        if (var==2) nbvF2[ind_f] = fmax2(countEff/(nbdsF2[ind_f]*eff2F2[ind_f]),0.0);
     }

      //  for (int e = 0 ; e < nbE ; e++){

                if (Qvec[eTemp]==0) {

                double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

                if (var==1) {
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbds_t2[it]*mpond_oth[eTemp],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbds_t2[it]*mpond_oth[eTemp]),0.0);
                }

                if (var==2){
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbv_t2[it]*mpond_oth[eTemp],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbv_t2[it]*mpond_oth[eTemp]),0.0);
                }


                } else {  //esp�ce SS3


                        double *Fothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 116));
                        double *Fothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 117));
                        double *Fothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 118));
                        double *Fothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 119));
                        double *Fothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 120));
                        double *Fothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 121));
                        double *Fothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 122));
                        double *Fothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 123));
                        double *Fothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 124));
                        double *Fothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 125));
                        double *Fothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 126));
                        double *Fothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 127));
                        double *Fothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 128));
                        double *Fothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 129));
                        double *Fothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 130));
                        double *Fothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 131));

                        double *FRWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 176));
                        double *FRWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 177));
                        double *FRWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 178));
                        double *FRWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 179));
                        double *FRWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 180));
                        double *FRWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 181));
                        double *FRWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 182));
                        double *FRWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 183));
                        double *FRWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 184));
                        double *FRWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 185));
                        double *FRWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 186));
                        double *FRWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 187));
                        double *FRWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 188));
                        double *FRWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 189));
                        double *FRWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 190));
                        double *FRWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 191));

                        double *FDWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 208));
                        double *FDWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 209));
                        double *FDWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 210));
                        double *FDWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 211));
                        double *FDWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 212));
                        double *FDWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 213));
                        double *FDWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 214));
                        double *FDWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 215));
                        double *FDWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 216));
                        double *FDWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 217));
                        double *FDWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 218));
                        double *FDWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 219));
                        double *FDWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 220));
                        double *FDWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 221));
                        double *FDWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 222));
                        double *FDWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 223));

                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));
                double multip = 0.0;
                if (var==1) multip = mu_nbds_t2[it];
                if (var==2) multip = mu_nbv_t2[it];

                if (gestyp==1)
                    for (int ag = 0; ag < nbi; ag++) {

                            Fothi2_S1M1[ag + it*nbi] = fmax2(Fothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M2[ag + it*nbi] = fmax2(Fothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M3[ag + it*nbi] = fmax2(Fothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M4[ag + it*nbi] = fmax2(Fothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M1[ag + it*nbi] = fmax2(Fothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M2[ag + it*nbi] = fmax2(Fothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M3[ag + it*nbi] = fmax2(Fothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M4[ag + it*nbi] = fmax2(Fothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M1[ag + it*nbi] = fmax2(Fothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M2[ag + it*nbi] = fmax2(Fothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M3[ag + it*nbi] = fmax2(Fothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M4[ag + it*nbi] = fmax2(Fothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M1[ag + it*nbi] = fmax2(Fothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M2[ag + it*nbi] = fmax2(Fothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M3[ag + it*nbi] = fmax2(Fothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M4[ag + it*nbi] = fmax2(Fothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                            FRWTothi2_S1M1[ag + it*nbi] = fmax2(FRWTothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M2[ag + it*nbi] = fmax2(FRWTothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M3[ag + it*nbi] = fmax2(FRWTothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M4[ag + it*nbi] = fmax2(FRWTothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M1[ag + it*nbi] = fmax2(FRWTothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M2[ag + it*nbi] = fmax2(FRWTothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M3[ag + it*nbi] = fmax2(FRWTothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M4[ag + it*nbi] = fmax2(FRWTothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M1[ag + it*nbi] = fmax2(FRWTothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M2[ag + it*nbi] = fmax2(FRWTothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M3[ag + it*nbi] = fmax2(FRWTothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M4[ag + it*nbi] = fmax2(FRWTothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M1[ag + it*nbi] = fmax2(FRWTothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M2[ag + it*nbi] = fmax2(FRWTothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M3[ag + it*nbi] = fmax2(FRWTothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M4[ag + it*nbi] = fmax2(FRWTothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                            FDWTothi2_S1M1[ag + it*nbi] = fmax2(FDWTothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M2[ag + it*nbi] = fmax2(FDWTothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M3[ag + it*nbi] = fmax2(FDWTothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M4[ag + it*nbi] = fmax2(FDWTothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M1[ag + it*nbi] = fmax2(FDWTothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M2[ag + it*nbi] = fmax2(FDWTothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M3[ag + it*nbi] = fmax2(FDWTothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M4[ag + it*nbi] = fmax2(FDWTothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M1[ag + it*nbi] = fmax2(FDWTothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M2[ag + it*nbi] = fmax2(FDWTothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M3[ag + it*nbi] = fmax2(FDWTothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M4[ag + it*nbi] = fmax2(FDWTothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M1[ag + it*nbi] = fmax2(FDWTothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M2[ag + it*nbi] = fmax2(FDWTothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M3[ag + it*nbi] = fmax2(FDWTothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M4[ag + it*nbi] = fmax2(FDWTothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                    }

                if (gestyp==2)
                    for (int ag = 0; ag < nbi; ag++) {

                            Fothi2_S1M1[ag + it*nbi] = fmax2(Fothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M2[ag + it*nbi] = fmax2(Fothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M3[ag + it*nbi] = fmax2(Fothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M4[ag + it*nbi] = fmax2(Fothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M1[ag + it*nbi] = fmax2(Fothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M2[ag + it*nbi] = fmax2(Fothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M3[ag + it*nbi] = fmax2(Fothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M4[ag + it*nbi] = fmax2(Fothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M1[ag + it*nbi] = fmax2(Fothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M2[ag + it*nbi] = fmax2(Fothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M3[ag + it*nbi] = fmax2(Fothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M4[ag + it*nbi] = fmax2(Fothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M1[ag + it*nbi] = fmax2(Fothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M2[ag + it*nbi] = fmax2(Fothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M3[ag + it*nbi] = fmax2(Fothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M4[ag + it*nbi] = fmax2(Fothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);

                            FRWTothi2_S1M1[ag + it*nbi] = fmax2(FRWTothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M2[ag + it*nbi] = fmax2(FRWTothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M3[ag + it*nbi] = fmax2(FRWTothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M4[ag + it*nbi] = fmax2(FRWTothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M1[ag + it*nbi] = fmax2(FRWTothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M2[ag + it*nbi] = fmax2(FRWTothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M3[ag + it*nbi] = fmax2(FRWTothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M4[ag + it*nbi] = fmax2(FRWTothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M1[ag + it*nbi] = fmax2(FRWTothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M2[ag + it*nbi] = fmax2(FRWTothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M3[ag + it*nbi] = fmax2(FRWTothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M4[ag + it*nbi] = fmax2(FRWTothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M1[ag + it*nbi] = fmax2(FRWTothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M2[ag + it*nbi] = fmax2(FRWTothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M3[ag + it*nbi] = fmax2(FRWTothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M4[ag + it*nbi] = fmax2(FRWTothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);

                            FDWTothi2_S1M1[ag + it*nbi] = fmax2(FDWTothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M2[ag + it*nbi] = fmax2(FDWTothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M3[ag + it*nbi] = fmax2(FDWTothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M4[ag + it*nbi] = fmax2(FDWTothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M1[ag + it*nbi] = fmax2(FDWTothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M2[ag + it*nbi] = fmax2(FDWTothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M3[ag + it*nbi] = fmax2(FDWTothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M4[ag + it*nbi] = fmax2(FDWTothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M1[ag + it*nbi] = fmax2(FDWTothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M2[ag + it*nbi] = fmax2(FDWTothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M3[ag + it*nbi] = fmax2(FDWTothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M4[ag + it*nbi] = fmax2(FDWTothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M1[ag + it*nbi] = fmax2(FDWTothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M2[ag + it*nbi] = fmax2(FDWTothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M3[ag + it*nbi] = fmax2(FDWTothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M4[ag + it*nbi] = fmax2(FDWTothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                    }


                }


         //}
}


//Rprintf("intro2\n");fichier << "intro2" << endl;
// this function is only present in IAM20 and florence work...where does it comes from ?
if ( (delay<=it) & !isNull(Ftarg) & !isNull(W_Ftarg) & (it>=1) & ((t_stop==0) | (t_stop>it)) ) {
    if(VERBOSE){Rprintf(" | EstimationTACfromF");}
    int oooo = EstimationTACfromF(it) ;
    oooo = oooo * 2;
    if(VERBOSE){Rprintf(" | ");}
}

// Florence QuotaMarket avec it > 1
if ((INTEGER(VECTOR_ELT(parQEX,0))[0]==1) & (delay<=it) & !isNull(TACbyF) & !isNull(TAC) & (it>=1) & ((t_stop==0) | (t_stop>it))) {
   if(VERBOSE){Rprintf("QuotaMarket");}
   QuotaMarket(list, VECTOR_ELT(parQEX,1), VECTOR_ELT(parQEX,2), VECTOR_ELT(parQEX,3), REAL(VECTOR_ELT(parQEX,4))[0],REAL(VECTOR_ELT(parQEX,5))[0], REAL(VECTOR_ELT(parQEX,6))[0],INTEGER(VECTOR_ELT(parQEX,7))[0], parBHV, it, INTEGER(persCalc)[0]);
   if(VERBOSE){Rprintf(" | ");}
}


//if ((INTEGER(VECTOR_ELT(parQEX,0))[0]==0) & (delay<=it) & !all_is_na(TACbyF) & !all_is_na(TAC) & (it>=1) & (gestInd==1) & (t_stop==0 | t_stop>it)) {  //optimisation TAC par flottille activ�e si au moins un �l�ment de TACbyF est renseign�
if ( (delay<=it) & !isNull(TACbyF) & !isNull(TAC) & (it>=1) & ((t_stop==0) | (t_stop>it))) {
    if(VERBOSE){Rprintf("GestionF2");}
    abv_GestionF2(it, updateE, tacCTRL, FList, VERBOSE);
    if(VERBOSE){Rprintf(" | ");}
}


//if ((INTEGER(Bootstrp)[0]==0) & (INTEGER(VECTOR_ELT(parQEX,0))[0]==1) & (gestInd==0) & (delay<=it)) {
//
//
//    QuotaExchV2(REAL(VECTOR_ELT(parQEX,1))[0],REAL(VECTOR_ELT(parQEX,2))[0],REAL(VECTOR_ELT(parQEX,3))[0],
//                        REAL(VECTOR_ELT(parQEX,4))[0], eTemp, REAL(VECTOR_ELT(parQEX,5))[0], it);
//
//
//}

//Rprintf("intro4\n");fichier << "intro4" << endl;

//if (scen & it>=1) Scenario(list, listScen, it);
//Rprintf("E\n");
//on remplit l'objet de sortie d�crivant les variables 'nbv' et 'nbds'
double *nbdsFM4 = REAL(getListElement(FList, "effort1_f_m"));
double *nbdsF4 = REAL(getListElement(FList, "effort1_f"));
double *nbvFM4 = REAL(getListElement(FList, "nbv_f_m"));
double *nbvF4 = REAL(getListElement(FList, "nbv_f"));
double *eff2FM4 = REAL(getListElement(FList, "effort2_f_m"));
double *eff2F4 = REAL(getListElement(FList, "effort2_f"));

for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

    NBVf[ind_f + nbF*it] = nbvF4[ind_f];
    NBDSf[ind_f + nbF*it] = nbdsF4[ind_f];
    EFF2f[ind_f + nbF*it] = eff2F4[ind_f];

    for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

        NBVfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbvFM4[ind_f + nbF*ind_m];
        NBDSfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbdsFM4[ind_f + nbF*ind_m];
        EFF2fm[ind_f + nbF*ind_m + nbF*nbMe*it] = eff2FM4[ind_f + nbF*ind_m];
    }
}


//Rprintf("F\n");fichier << "F" << endl;

//3 modules avec pas de temps diff�renci� au niveau trimestre

if (nbE>0) {
if(VERBOSE){Rprintf("\n  Mortalite");}
 Mortalite(list, it, eVar, VERBOSE);//Rprintf("\nG");fichier << "G" << endl;//if (it>4) error("BBBhh");////PrintValue(out_Fr_fmi);//PrintValue(VECTOR_ELT(eVar,60));
 if(VERBOSE){Rprintf(" | ");}

if(VERBOSE){Rprintf("\n  DynamicPop");}
 DynamicPop(list, it, eVar, true, 0);//Rprintf("\nH");fichier << "H" << endl;////PrintValue(out_Z_eit);//PrintValue(out_N_eitQ);//PrintValue(out_N_eit);
 if(VERBOSE){Rprintf(" | ");}
}

if(VERBOSE){Rprintf("\n  CatchDL");}
CatchDL(list, it, eVar, 0);//Rprintf("\nI");fichier << "I" << endl;////PrintValue(out_Y_eit);
if(VERBOSE){Rprintf(" | ");}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//------------------------ Correction de l'ajustement initial si vise une esp�ce dynamique -------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

conform = 0; //0=OK

if (eTemp<nbE) { //esp�ce dynamique
    if(VERBOSE){Rprintf("\n  Correct");}

    SEXP NDIM = allocVector(INTSXP,4);
    int *ndi = INTEGER(NDIM); for (int q=0;q<3;q++) ndi[q]=0; ndi[3]=nbT;
    double YtotN = REAL(aggregObj(VECTOR_ELT(out_Y_eit,eTemp),NDIM))[it];
    double YtotNm1 = REAL(aggregObj(VECTOR_ELT(out_Y_eit,eTemp),NDIM))[it-1];
    double FbarNm1=NA_REAL, FbarNm2=NA_REAL, FbarNm3=NA_REAL;
    if (it>=3){
        //FbarN = REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it];
        FbarNm1 = REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it-1];
        FbarNm2 = REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it-2];
        FbarNm3 = REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it-3];
    }
    double SsbN = REAL(VECTOR_ELT(out_SSB_et, eTemp))[it];

    ////Rprintf("time %i\n",it);

    typeGest = 0;
    if ((trgt==2) & (gestInd==1)) typeGest = 3;
    if (((trgt==1) | (trgt==3)) & (gestInd==1)) typeGest = 2;

    //point n�6 (point dominant sur les deux autres)
    if ((delay<=it) & (gestInd==1) & ((trgt==1) | (trgt==3)) & !ISNA(Blim_CPP) & (it>=1)) {
        if (SsbN<Blim_CPP) {
            trgt=2;////Rprintf("conform 6\n");
            conform = 6;
            typeGest = 6;
            door=false;
        }
    }

    //point n�4
    if ((delay<=it) & (gestInd==1) & (trgt==2) & (SsbN>=Blim_CPP) & door & !ISNA(tolVarTACinf_CPP) & !ISNA(tolVarTACsup_CPP) & (it>=1)) {
        if ((YtotN<YtotNm1*tolVarTACinf_CPP) | (YtotN>YtotNm1*tolVarTACsup_CPP)) {
          trgt=1;
          door=false;//pour hi�rarchiser les actions
          if (YtotN<YtotNm1*tolVarTACinf_CPP) TAC_glob[it] = YtotNm1*tolVarTACinf_CPP; ////Rprintf("Ytot %f\n",TAC_glob[it]);
          if (YtotN>YtotNm1*tolVarTACsup_CPP) TAC_glob[it] = YtotNm1*tolVarTACsup_CPP;
          conform = 4; //point de blocage
          typeGest = 4;
        }
    }

    //point n�5 (corVarTACbby n'est pas utilis�)
    if ((delay<=it) & (gestInd==1) & ((trgt==1) | (trgt==3)) & door & !ISNA(corVarTACval_CPP) & (it>=3)) { //HYP: TAC suivi fixe au cours du temps
        if ((FbarNm2<FbarNm1) & (FbarNm3<FbarNm2)) {
          for (int yr=it;yr<nbT;yr++) TAC_glob[yr] = TAC_glob[it-1]*corVarTACval_CPP; ////Rprintf("Ytot %f",TAC_glob[yr]);}
          ////Rprintf("conform 5\n");
          conform = 5; //point de blocage
          typeGest = 5;
        }
    }

    typegest[it] = typeGest;
    door = true;

}

////Rprintf("trgt %i\n",trgt);
//selon les valeurs de 'conform', on engage ou pas le processus de correction

if (conform>0) {
    if(VERBOSE){Rprintf("yes");}

//if (it==3) {//Rprintf("AAAAAAAAAAAA\n"); //PrintValue(FList_copy);}
//on r�initialise Flist, list, eVar (incluant Fothi) et fVar (pas s�r que ce dernier soit utile, mais par pr�caution...)
REPROTECT(list = duplicate(list_copy), ipx_list);
REPROTECT(FList = getListElement(list, "Fleet"), ipx_FList);
REPROTECT(eVar = duplicate(eVar_copy), ipx_eVar);
REPROTECT(fVar = duplicate(fVar_copy), ipx_fVar);

if ((delay<=it) & (gestInd==1) & (it>=1) & isNull(TACbyF)) {
//Rprintf("AA");fichier << "AA" << endl;
////PrintValue(getListElement(getListElement(list, "Fleet"), "nbds_f"));
////PrintValue(getListElement(FList, "nbds_f"));
////PrintValue(mu_nbds);
////PrintValue(mu_nbv);

Gestion(list, it);

//on met � jour les variables sur lesquelles op�re le multiplicateur
    double *mu_nbds_t2 = REAL(mu_nbds);
    double *mu_nbv_t2 = REAL(mu_nbv);
    double *nbdsFM2 = REAL(getListElement(FList, "effort1_f_m"));
    double *nbdsF2 = REAL(getListElement(FList, "effort1_f"));
    double *nbvFM2 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF2 = REAL(getListElement(FList, "nbv_f"));
    double *eff2FM2 = REAL(getListElement(FList, "effort2_f_m"));
    double *eff2F2 = REAL(getListElement(FList, "effort2_f"));


     for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double countEff = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m] + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m]*(1 + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (var==2) {
                if (gestyp==1) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m] + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m]*(1 + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
            }

            if (!ISNA(nbdsFM2[ind_f+nbF*ind_m]) & !ISNA(nbvFM2[ind_f+nbF*ind_m]) & !ISNA(eff2FM2[ind_f+nbF*ind_m])) countEff = countEff + nbdsFM2[ind_f+nbF*ind_m]*nbvFM2[ind_f+nbF*ind_m]*eff2FM2[ind_f+nbF*ind_m];
        }

        if (var==1) nbdsF2[ind_f] = fmax2(countEff/(nbvF2[ind_f]*eff2F2[ind_f]),0.0); //NBDSf
        if (var==2) nbvF2[ind_f] = fmax2(countEff/(nbdsF2[ind_f]*eff2F2[ind_f]),0.0);
     }

      //  for (int e = 0 ; e < nbE ; e++){
                if (Qvec[eTemp]==0) {

                double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

                if (var==1) {
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbds_t2[it]*mpond_oth[eTemp],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbds_t2[it]*mpond_oth[eTemp]),0.0);
                }

                if (var==2){
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbv_t2[it]*mpond_oth[eTemp],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbv_t2[it]*mpond_oth[eTemp]),0.0);
                }


                } else {  //esp�ce SS3


                        double *Fothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 116));
                        double *Fothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 117));
                        double *Fothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 118));
                        double *Fothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 119));
                        double *Fothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 120));
                        double *Fothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 121));
                        double *Fothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 122));
                        double *Fothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 123));
                        double *Fothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 124));
                        double *Fothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 125));
                        double *Fothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 126));
                        double *Fothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 127));
                        double *Fothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 128));
                        double *Fothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 129));
                        double *Fothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 130));
                        double *Fothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 131));

                        double *FRWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 176));
                        double *FRWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 177));
                        double *FRWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 178));
                        double *FRWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 179));
                        double *FRWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 180));
                        double *FRWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 181));
                        double *FRWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 182));
                        double *FRWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 183));
                        double *FRWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 184));
                        double *FRWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 185));
                        double *FRWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 186));
                        double *FRWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 187));
                        double *FRWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 188));
                        double *FRWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 189));
                        double *FRWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 190));
                        double *FRWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 191));

                        double *FDWTothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 208));
                        double *FDWTothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 209));
                        double *FDWTothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 210));
                        double *FDWTothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 211));
                        double *FDWTothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 212));
                        double *FDWTothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 213));
                        double *FDWTothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 214));
                        double *FDWTothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 215));
                        double *FDWTothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 216));
                        double *FDWTothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 217));
                        double *FDWTothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 218));
                        double *FDWTothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 219));
                        double *FDWTothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 220));
                        double *FDWTothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 221));
                        double *FDWTothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 222));
                        double *FDWTothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTemp), 223));

                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"));
                double multip = 0.0;
                if (var==1) multip = mu_nbds_t2[it];
                if (var==2) multip = mu_nbv_t2[it];

                if (gestyp==1)
                    for (int ag = 0; ag < nbi; ag++) {

                            Fothi2_S1M1[ag + it*nbi] = fmax2(Fothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M2[ag + it*nbi] = fmax2(Fothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M3[ag + it*nbi] = fmax2(Fothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S1M4[ag + it*nbi] = fmax2(Fothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M1[ag + it*nbi] = fmax2(Fothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M2[ag + it*nbi] = fmax2(Fothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M3[ag + it*nbi] = fmax2(Fothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S2M4[ag + it*nbi] = fmax2(Fothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M1[ag + it*nbi] = fmax2(Fothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M2[ag + it*nbi] = fmax2(Fothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M3[ag + it*nbi] = fmax2(Fothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S3M4[ag + it*nbi] = fmax2(Fothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M1[ag + it*nbi] = fmax2(Fothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M2[ag + it*nbi] = fmax2(Fothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M3[ag + it*nbi] = fmax2(Fothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            Fothi2_S4M4[ag + it*nbi] = fmax2(Fothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                            FRWTothi2_S1M1[ag + it*nbi] = fmax2(FRWTothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M2[ag + it*nbi] = fmax2(FRWTothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M3[ag + it*nbi] = fmax2(FRWTothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S1M4[ag + it*nbi] = fmax2(FRWTothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M1[ag + it*nbi] = fmax2(FRWTothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M2[ag + it*nbi] = fmax2(FRWTothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M3[ag + it*nbi] = fmax2(FRWTothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S2M4[ag + it*nbi] = fmax2(FRWTothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M1[ag + it*nbi] = fmax2(FRWTothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M2[ag + it*nbi] = fmax2(FRWTothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M3[ag + it*nbi] = fmax2(FRWTothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S3M4[ag + it*nbi] = fmax2(FRWTothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M1[ag + it*nbi] = fmax2(FRWTothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M2[ag + it*nbi] = fmax2(FRWTothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M3[ag + it*nbi] = fmax2(FRWTothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FRWTothi2_S4M4[ag + it*nbi] = fmax2(FRWTothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                            FDWTothi2_S1M1[ag + it*nbi] = fmax2(FDWTothi2_S1M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M2[ag + it*nbi] = fmax2(FDWTothi2_S1M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M3[ag + it*nbi] = fmax2(FDWTothi2_S1M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S1M4[ag + it*nbi] = fmax2(FDWTothi2_S1M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M1[ag + it*nbi] = fmax2(FDWTothi2_S2M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M2[ag + it*nbi] = fmax2(FDWTothi2_S2M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M3[ag + it*nbi] = fmax2(FDWTothi2_S2M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S2M4[ag + it*nbi] = fmax2(FDWTothi2_S2M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M1[ag + it*nbi] = fmax2(FDWTothi2_S3M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M2[ag + it*nbi] = fmax2(FDWTothi2_S3M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M3[ag + it*nbi] = fmax2(FDWTothi2_S3M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S3M4[ag + it*nbi] = fmax2(FDWTothi2_S3M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M1[ag + it*nbi] = fmax2(FDWTothi2_S4M1[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M2[ag + it*nbi] = fmax2(FDWTothi2_S4M2[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M3[ag + it*nbi] = fmax2(FDWTothi2_S4M3[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);
                            FDWTothi2_S4M4[ag + it*nbi] = fmax2(FDWTothi2_S4M4[ag + it*nbi] + multip*mpond_oth[eTemp],0.0);

                    }

                if (gestyp==2)
                    for (int ag = 0; ag < nbi; ag++) {

                            Fothi2_S1M1[ag + it*nbi] = fmax2(Fothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M2[ag + it*nbi] = fmax2(Fothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M3[ag + it*nbi] = fmax2(Fothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S1M4[ag + it*nbi] = fmax2(Fothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M1[ag + it*nbi] = fmax2(Fothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M2[ag + it*nbi] = fmax2(Fothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M3[ag + it*nbi] = fmax2(Fothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S2M4[ag + it*nbi] = fmax2(Fothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M1[ag + it*nbi] = fmax2(Fothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M2[ag + it*nbi] = fmax2(Fothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M3[ag + it*nbi] = fmax2(Fothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S3M4[ag + it*nbi] = fmax2(Fothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M1[ag + it*nbi] = fmax2(Fothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M2[ag + it*nbi] = fmax2(Fothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M3[ag + it*nbi] = fmax2(Fothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            Fothi2_S4M4[ag + it*nbi] = fmax2(Fothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);

                            FRWTothi2_S1M1[ag + it*nbi] = fmax2(FRWTothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M2[ag + it*nbi] = fmax2(FRWTothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M3[ag + it*nbi] = fmax2(FRWTothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S1M4[ag + it*nbi] = fmax2(FRWTothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M1[ag + it*nbi] = fmax2(FRWTothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M2[ag + it*nbi] = fmax2(FRWTothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M3[ag + it*nbi] = fmax2(FRWTothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S2M4[ag + it*nbi] = fmax2(FRWTothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M1[ag + it*nbi] = fmax2(FRWTothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M2[ag + it*nbi] = fmax2(FRWTothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M3[ag + it*nbi] = fmax2(FRWTothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S3M4[ag + it*nbi] = fmax2(FRWTothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M1[ag + it*nbi] = fmax2(FRWTothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M2[ag + it*nbi] = fmax2(FRWTothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M3[ag + it*nbi] = fmax2(FRWTothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FRWTothi2_S4M4[ag + it*nbi] = fmax2(FRWTothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);

                            FDWTothi2_S1M1[ag + it*nbi] = fmax2(FDWTothi2_S1M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M2[ag + it*nbi] = fmax2(FDWTothi2_S1M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M3[ag + it*nbi] = fmax2(FDWTothi2_S1M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S1M4[ag + it*nbi] = fmax2(FDWTothi2_S1M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M1[ag + it*nbi] = fmax2(FDWTothi2_S2M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M2[ag + it*nbi] = fmax2(FDWTothi2_S2M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M3[ag + it*nbi] = fmax2(FDWTothi2_S2M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S2M4[ag + it*nbi] = fmax2(FDWTothi2_S2M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M1[ag + it*nbi] = fmax2(FDWTothi2_S3M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M2[ag + it*nbi] = fmax2(FDWTothi2_S3M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M3[ag + it*nbi] = fmax2(FDWTothi2_S3M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S3M4[ag + it*nbi] = fmax2(FDWTothi2_S3M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M1[ag + it*nbi] = fmax2(FDWTothi2_S4M1[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M2[ag + it*nbi] = fmax2(FDWTothi2_S4M2[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M3[ag + it*nbi] = fmax2(FDWTothi2_S4M3[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                            FDWTothi2_S4M4[ag + it*nbi] = fmax2(FDWTothi2_S4M4[ag + it*nbi]*(1+multip*mpond_oth[eTemp]),0.0);
                    }




                }


        // }

}

////Rprintf("CC2222222");
////PrintValue(getListElement(FList, "nbds_f"));

////PrintValue(NBDSF);

double *nbdsFM4 = REAL(getListElement(FList, "effort1_f_m"));
double *nbdsF4 = REAL(getListElement(FList, "effort1_f"));
double *nbvFM4 = REAL(getListElement(FList, "nbv_f_m"));
double *nbvF4 = REAL(getListElement(FList, "nbv_f"));
double *eff2FM4 = REAL(getListElement(FList, "effort2_f_m"));
double *eff2F4 = REAL(getListElement(FList, "effort2_f"));

for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

    NBVf[ind_f + nbF*it] = nbvF4[ind_f];
    NBDSf[ind_f + nbF*it] = nbdsF4[ind_f];
    EFF2f[ind_f + nbF*it] = eff2F4[ind_f];

    for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

        NBVfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbvFM4[ind_f + nbF*ind_m];
        NBDSfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbdsFM4[ind_f + nbF*ind_m];
        EFF2fm[ind_f + nbF*ind_m + nbF*nbMe*it] = eff2FM4[ind_f + nbF*ind_m];
    }
}

////PrintValue(getListElement(FList, "nbds_f"));


//3 modules avec pas de temps diff�renci� au niveau trimestre  <<-- ANNULE
Mortalite(list, it, eVar);////Rprintf("\nG");////PrintValue(out_Fr_fmi);//PrintValue(VECTOR_ELT(eVar,60));
DynamicPop(list, it, eVar,true);////Rprintf("\nH");//PrintValue(out_Z_eit);//PrintValue(out_N_eitQ);//PrintValue(out_N_eit);
CatchDL(list, it, eVar);////Rprintf("\nI");////PrintValue(out_Y_eit);


//une fois que c'est termin�, il faut remettre certaines choses en place dans le cas n�4
if (conform==4) trgt=2;
if(VERBOSE){Rprintf(" done");}
}



if(VERBOSE){Rprintf("\nMarche");}
Marche(list, it);
if(VERBOSE){Rprintf(" | ");}

if(VERBOSE){Rprintf("EcoDCF");}
EcoDCF(list, it, INTEGER(persCalc)[0], REAL(dr)[0], VERBOSE);
if(VERBOSE){Rprintf(" | ");}

//module gestion : si delay<=it & gestInd==1 & trgt==3 et si Fbar atteint, on rebascule trgt � 2
if (eTemp<nbE) {
  if ((delay<=it) & (gestInd==1) & (trgt==3) & (REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it] <= Fbar_trgt[it]))  trgt=2;    //????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
}
//remise � l'�tat initial des param�tres du module Gestion si n�cessaire (ATTENTION : implique qu'on n'atteint jamais l'effort nul quelle que soit la flottille
//                                                                                 dans le cas contraire, fixer upd=2)

}
////Rprintf("K");
//fichier << "Fin boucle T"  << endl;

SET_VECTOR_ELT(out_effort, 0, NBVF); SET_VECTOR_ELT(out_effort, 1, NBDSF); SET_VECTOR_ELT(out_effort, 2, EFF2F);
SET_VECTOR_ELT(out_effort, 3, NBVFM); SET_VECTOR_ELT(out_effort, 4, NBDSFM); SET_VECTOR_ELT(out_effort, 5, EFF2FM);

const char *nmEf[6] = {"nbv_f","effort1_f","effort2_f","nbv_f_m","effort1_f_m","effort2_f_m"};
PROTECT(nmsEF = allocVector(STRSXP, 6));
for(int k = 0; k < 6; k++) SET_STRING_ELT(nmsEF, k, mkChar(nmEf[k]));
setAttrib(out_effort, R_NamesSymbol, nmsEF);
////Rprintf("L");
//if (eTemp<nbE) {
//  if (Qvec[eTemp]==0) {
//    free_vector(Ztemp,1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI")));
//  } else {
//    free_vector(Ztemp,1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI"))*16);
//  }
//} else {
//  free_vector(Ztemp,1,2);
//}
free_vector(Etemp,1,nbF);


free_vector(Einterm_fm,1,nbF*nbMe);
free_vector(Einterm_fm_copy,1,nbF*nbMe);
free_vector(EffsupTMP_fm,1,nbF*nbMe);
free_vector(multFOTHinterm_e,1,nbE);

//Rprintf("K2\n");fichier << "K2" << endl;
//UNPROTECT(123+nbE+nbE+32+11+1+3+3+2+1+5); //+6 ajout�s apr�s int�gration de 'parOQD'
//if (nbEstat>0) UNPROTECT(nbEstat);
if (pflex){ UNPROTECT(2); } else {UNPROTECT(1);}
UNPROTECT(30-1);
UNPROTECT(17+20+4+16*6+8-1+9);//out_
UNPROTECT(14);
UNPROTECT(10); // PROTECT_WITH_INDEX
Rf_unprotect(2); // DimCst

//fichier.close();
}

//------------------------------------------------------------------------------------
//destructeur de la classe Param
//------------------------------------------------------------------------------------

//BioEcoPar::~BioEcoPar()
//{
//}

