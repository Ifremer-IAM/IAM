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
// #include <Rmath.h>
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
                     SEXP GestInd, SEXP mOth, SEXP bounds, SEXP TACL, SEXP FBAR, SEXP othSpSup, SEXP effSup, SEXP GestParam, SEXP EcoDcf,
                     SEXP persCalc, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyFL, SEXP Ftarg, SEXP W_Ftarg, SEXP MeanRec_Ftarg,
                     SEXP parBHV, SEXP parQEX,
                     SEXP tacCTRL, SEXP stochPrice, SEXP updateE, SEXP parOQD, int VERBOSE)
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

ecodcf = INTEGER(EcoDcf)[0];
drCopy = REAL(dr)[0];

//int conform = 0;
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
pUpdate = !isNull(getListElement(list, "Market"));  //Rprintf("pUpdate = %d \n",pUpdate) ; //
eUpdate = true;    //

if (pUpdate) {
    if(VERBOSE){Rprintf("pUpdate | ");}
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

PROTECT(othSpSupList = othSpSup);
PROTECT(effSupMat = effSup);

gestInd = INTEGER(GestInd)[0];
PROTECT(m_fm = duplicate(mFM)); if (length(m_fm)!=nbF*nbM) error("Check dimension of array 'mFleetMetier'!!\n");
PROTECT(m_oth = duplicate(mOth)); if (length(m_oth)!=nbE) error("Check dimension of array 'mOth'!!\n");
X1 = REAL(bounds)[0];
X2 = REAL(bounds)[1];

//TAC_glob = REAL(TAC);  //� corriger
Fbar_trgt = REAL(FBAR);  //� corriger
PROTECT(TACbyF = TACbyFL);  //� corriger
PROTECT(TAC = TACL); // TODO :pourquoi juste un changement de nom la en fait ?

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
maxIter = INTEGER(getListElement(tacCTRL, "maxIter"))[0];
diffZmax = REAL(getListElement(tacCTRL, "diffZmax"))[0];
lambda = REAL(getListElement(tacCTRL, "lambda"))[0];
t_stop = INTEGER(getListElement(tacCTRL, "t_stop"))[0];
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

//int DCFok = INTEGER(EcoDcf)[0];

SRInd = INTEGER(SRind);

//bool door = true;

//il faut initialiser 'eVar' dans lequel on int�grera toutes les variables interm�diaires � d�cliner par esp�ce
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


if(VERBOSE){Rprintf("Step 0.2 | ");}

PROTECT_WITH_INDEX(fVar = allocVector(VECSXP, 34),&ipx_fVar); //32= rtbsIni_f & 33=rtbsIni_f_m & 33=ETini_f_m
PROTECT_WITH_INDEX(fVar_copy = duplicate(fVar),&ipx_fVar_copy);

//////Rprintf("A");

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
PROTECT(out_Eco = allocVector(VECSXP, 69));
PROTECT(out_EcoDCF = allocVector(VECSXP, 60));
PROTECT(out_effort = allocVector(VECSXP, 6)); //nbv_f, effort1_f, effort2_f, nbv_f_m, effort1_f_m, effort2_f_m

PROTECT(mu_nbds = allocVector(REALSXP, nbT)); //il reste la mise en forme � op�rer
PROTECT(mu_nbv = allocVector(REALSXP, nbT));
PROTECT(out_typeGest = allocVector(INTSXP, nbT));


PROTECT(out_Ytot_fm = NEW_NUMERIC(nbF*nbMe*nbT));
PROTECT(out_DD_efmi = allocVector(VECSXP, nbE));
PROTECT(out_DD_efmc = allocVector(VECSXP, nbE));
PROTECT(out_LD_efmi = allocVector(VECSXP, nbE));
PROTECT(out_LD_efmc = allocVector(VECSXP, nbE));
PROTECT(out_statDD_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLD_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLDst_efm = allocVector(VECSXP, nbEstat));
PROTECT(out_statLDor_efm = allocVector(VECSXP, nbEstat));


if(VERBOSE){Rprintf("Step 0.3 | ");}

//int *typegest = INTEGER(out_typeGest);
double *mu_nbds_t = REAL(mu_nbds); for (int i=0; i<nbT; i++) mu_nbds_t[i] = 0.0; //initialisation
double *mu_nbv_t = REAL(mu_nbv); for (int i=0; i<nbT; i++) mu_nbv_t[i] = 0.0;    //

//double *mpond_fm = REAL(m_fm);
//double *mpond_oth = REAL(m_oth);

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
setAttrib(NBDSF, R_DimNamesSymbol, dnmsF); setAttrib(NBDSFM, R_DimNamesSymbol, dnmsFM);
setAttrib(EFF2F, R_DimNamesSymbol, dnmsF); setAttrib(EFF2FM, R_DimNamesSymbol, dnmsFM);

PROTECT(out_allocEff_fm = alloc3DArray(REALSXP,nbF,nbMe,nbT));
setAttrib(out_allocEff_fm, R_DimNamesSymbol, dnmsFM);
// Intitialisation avec mfm mais modifie dans le module QuotaExchange si appele
for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
for (int ind_m = 0 ; ind_m < nbMe ; ind_m++)
for (int ind_t = 0 ; ind_t < nbT ; ind_t++) {
        REAL(out_allocEff_fm)[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = REAL(m_fm)[ind_f + nbF*ind_m];}
//PrintValue(out_allocEff_fm);

setAttrib(reconcilSPP, R_DimNamesSymbol, dnmsFM);
setAttrib(reconcilSPP_copy, R_DimNamesSymbol, dnmsFM);

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

if(VERBOSE){Rprintf("\nLoop : \n");}

for (int it = 0; it < nbT ; it++) {
    Rprintf("T = %d => ", it);

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
        RecAlea(list, listStochastic, it, 1, recType1); //IMPORTANT : � effectuer AVANT la proc�dure d'optimisation
        RecAlea(list, listStochastic, it, 2, recType2);
        RecAlea(list, listStochastic, it, 3, recType3);
    }

        SRmod(list, listSR, it, TypeSR, SRInd); //important : � envoyer avant 'DynamicPop'

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


//Rprintf("intro2\n");fichier << "intro2" << endl;
// this function is only present in IAM20 and florence work...where does it comes from ?
if ( (delay<=it) & !isNull(Ftarg) & !isNull(W_Ftarg) & (it>=1) & ((t_stop==0) | (t_stop>it)) ) {
    if(VERBOSE){Rprintf(" | EstimationTACfromF");}
   int oooo = EstimationTACfromF(it) ;
   oooo = oooo * 2;
   if(VERBOSE){Rprintf(" | ");}
}

if ((INTEGER(VECTOR_ELT(parQEX,0))[0]==1) & (delay<=it) & !isNull(TACbyF) & !isNull(TAC) & (it>=1) & ((t_stop==0) | (t_stop>it))) {
   if(VERBOSE){Rprintf("QuotaMarket");}
   QuotaMarket(list, VECTOR_ELT(parQEX,1), VECTOR_ELT(parQEX,2), VECTOR_ELT(parQEX,3), REAL(VECTOR_ELT(parQEX,4))[0],REAL(VECTOR_ELT(parQEX,5))[0], REAL(VECTOR_ELT(parQEX,6))[0],INTEGER(VECTOR_ELT(parQEX,7))[0], parBHV, it, INTEGER(persCalc)[0]);
   if(VERBOSE){Rprintf(" | ");}
}


//if ((INTEGER(VECTOR_ELT(parQEX,0))[0]==0) & (delay<=it) & !all_is_na(TACbyF) & !all_is_na(TAC) & (it>=1) & (gestInd==1) & (t_stop==0 | t_stop>it)) {  //optimisation TAC par flottille activ�e si au moins un �l�ment de TACbyF est renseign�
if ( (delay<=it) & !isNull(TACbyF) & !isNull(TAC) & (it>=1) & ((t_stop==0) | (t_stop>it))) {
    if(VERBOSE){Rprintf("GestionF2 | ");}

//Rprintf("adjust\n");

//Rprintf("introOPT\n");fichier << "introOPT" << endl;

        int DELAY = INTEGER(updateE)[0];

        SPPstatOPT = INTEGER(getListElement(tacCTRL, "SPPstatOPT"));
        SPPspictOPT = INTEGER(getListElement(tacCTRL, "SPPspictOPT"));
        SPPdynOPT = INTEGER(getListElement(tacCTRL, "SPPdynOPT"));
        N_SPPstatOPT = length(getListElement(tacCTRL, "SPPstatOPT"));
        N_SPPspictOPT = length(getListElement(tacCTRL, "SPPspictOPT"));
        N_SPPdynOPT = length(getListElement(tacCTRL, "SPPdynOPT"));

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
//            //Rprintf("AA");
//            //PrintValue(getListElement(FList, "effort1_f_m"));

        if (N_SPPspictOPT>0) { //si esp�ce dynamique SPICT

         for (int i = 0; i < N_SPPspictOPT; i++){

            int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPspictOPT[i]))), "modI")); //doit normalement �tre �gal � 1
            double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, SPPspictOPT[i]), 44));
            for (int ag = 0; ag < nbi; ag++) Fothi2[ag + it*nbi] = Fothi2[ag + (DELAY-1)*nbi];

         }

        }

        if (N_SPPdynOPT>0) { //si esp�ce dynamique XSA ou SS3

            for (int i = 0; i < N_SPPdynOPT; i++){

            int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,SPPdynOPT[i]))), "modI"));
            int eTmp = SPPdynOPT[i];

            if (Qvec[eTmp]==0) {

                    double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 44));

                    for (int ag = 0; ag < nbi; ag++) Fothi2[ag + it*nbi] = Fothi2[ag + (DELAY-1)*nbi];


            } else {

                    double *Fothi2_S1M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 116)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + it*nbi] = Fothi2_S1M1[ag + (DELAY-1)*nbi];
                    double *Fothi2_S1M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 117)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + it*nbi] = Fothi2_S1M2[ag + (DELAY-1)*nbi];
                    double *Fothi2_S1M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 118)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + it*nbi] = Fothi2_S1M3[ag + (DELAY-1)*nbi];
                    double *Fothi2_S1M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 119)); for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + it*nbi] = Fothi2_S1M4[ag + (DELAY-1)*nbi];
                    double *Fothi2_S2M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 120)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + it*nbi] = Fothi2_S2M1[ag + (DELAY-1)*nbi];
                    double *Fothi2_S2M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 121)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + it*nbi] = Fothi2_S2M2[ag + (DELAY-1)*nbi];
                    double *Fothi2_S2M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 122)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + it*nbi] = Fothi2_S2M3[ag + (DELAY-1)*nbi];
                    double *Fothi2_S2M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 123)); for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + it*nbi] = Fothi2_S2M4[ag + (DELAY-1)*nbi];
                    double *Fothi2_S3M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 124)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + it*nbi] = Fothi2_S3M1[ag + (DELAY-1)*nbi];
                    double *Fothi2_S3M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 125)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + it*nbi] = Fothi2_S3M2[ag + (DELAY-1)*nbi];
                    double *Fothi2_S3M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 126)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + it*nbi] = Fothi2_S3M3[ag + (DELAY-1)*nbi];
                    double *Fothi2_S3M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 127)); for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + it*nbi] = Fothi2_S3M4[ag + (DELAY-1)*nbi];
                    double *Fothi2_S4M1 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 128)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + it*nbi] = Fothi2_S4M1[ag + (DELAY-1)*nbi];
                    double *Fothi2_S4M2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 129)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + it*nbi] = Fothi2_S4M2[ag + (DELAY-1)*nbi];
                    double *Fothi2_S4M3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 130)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + it*nbi] = Fothi2_S4M3[ag + (DELAY-1)*nbi];
                    double *Fothi2_S4M4 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, eTmp), 131)); for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + it*nbi] = Fothi2_S4M4[ag + (DELAY-1)*nbi];


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
        }}

        }


     //Rprintf("call.GestionF2\n");fichier << "call.GestionF2" << endl;
     int ooo = GestionF2(it);
     ooo = ooo * 2;
     //Rprintf("end.GestionF2\n");fichier << "end.GestionF2" << endl;

  //GestionF(NRmatrix(1,nbF+2,1,nbF+1), NRvector(1,nbF+2), nbF+1, 0.0000001, it);}

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
if(VERBOSE){Rprintf("Mortalite");}
 Mortalite(list, it, eVar);//Rprintf("\nG");fichier << "G" << endl;//if (it>4) error("BBBhh");////PrintValue(out_Fr_fmi);//PrintValue(VECTOR_ELT(eVar,60));
 if(VERBOSE){Rprintf(" | ");}

if(VERBOSE){Rprintf("DynamicPop");}
 DynamicPop(list, it, eVar, true);//Rprintf("\nH");fichier << "H" << endl;////PrintValue(out_Z_eit);//PrintValue(out_N_eitQ);//PrintValue(out_N_eit);
 if(VERBOSE){Rprintf(" | ");}
}
if(VERBOSE){Rprintf("CatchDL");}
CatchDL(list, it, eVar, VERBOSE);//Rprintf("\nI");fichier << "I" << endl;////PrintValue(out_Y_eit);
if(VERBOSE){Rprintf(" | ");}

//if (it>0) error("Unexpectedz scondition occurred");
//if (it<5) {
////Rprintf("Z1\n");
////PrintValue(VECTOR_ELT(out_Z_eit,2));
////Rprintf("N1\n");
////PrintValue(VECTOR_ELT(out_N_eit,2));
////Rprintf("Fbar1\n");
////PrintValue(VECTOR_ELT(out_Fbar_et,2));
//}


if(VERBOSE){Rprintf("Marche | ");}
Marche(list, it);

if(VERBOSE){Rprintf("EcoDCF \n");}
EcoDCF(list, it, INTEGER(persCalc)[0], REAL(dr)[0]);

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
if (pUpdate){ UNPROTECT(2); } else {UNPROTECT(1);}
UNPROTECT(30);
UNPROTECT(17+20+4+16*6+8+9);//out_
UNPROTECT(14);
UNPROTECT(10); // PROTECT_WITH_INDEX

//fichier.close();
}

//------------------------------------------------------------------------------------
//destructeur de la classe Param
//------------------------------------------------------------------------------------

//BioEcoPar::~BioEcoPar()
//{
//}

