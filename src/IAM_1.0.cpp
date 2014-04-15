#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>


using namespace std;

#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED



//classe regroupant les paramètres à intégrer dans le modèle, les variables intermédiaires (discrétisation du processus),
//ainsi que les sorties du modèle
class BioEcoPar
{
public: //normalement, selon les conventions,  les attributs doivent être "private"

typedef  double (BioEcoPar::*BEfn1)(double mult);
typedef  double (BioEcoPar::*BEfn1_F)(double *x);

//   INPUTS   ------------------

SEXP    list;       //liste d'objets R constituant la donnée entrée du modèle. Certaines variables seront remises à jour dans le cadre
                    //de mesures de gestion ou de modélisation de comportement.



//   OUTPUTS  ------------------

//outputs des différents modules (à initialiser pour t=0)
SEXP    out_F_fmi,  //mortalité "captures" par pêche (par espèce)
        out_Fr_fmi,  //mortalité totale (corrigée de la survie) par pêche (par espèce)
        out_Z_eit,  //coefficient de mortalité totale
        out_Fbar_et,  //Fbar par espèce (t)
        out_N_eit,  //effectifs en nombre
        out_B_et,   //biomasse (t)
        out_SSB_et, //biomasse de reproducteurs (t)
        out_C_efmit,//captures en nombres
        out_C_eit,  //captures totales en nombres
        out_Y_efmit,//captures en poids (t)
        out_Y_eit,  //captures totales en poids (t)
        out_D_efmit,//rejets en poids (t)
        out_L_efmit,//débarquements en poids aux âges(t)
        out_L_efmct,//débarquements en poids par catégories(t)
        out_L_efmct2,//débarquements en poids par catégories(t) pour le codage métier eco
        out_P_t,    //prix moyen (en euros) (niveau métier éco si dispo)
        out_CA_eft, //chiffre d'affaires par espèce (en euros)
        out_CAT_ft, //chiffre d'affaires total (en euros)
        out_CA_ft,  //chiffre d'affaires moyen par navire (en euros)
        out_RAP_ft, //reste à partager (en euros)
        out_EBE_ft, //excédent brut d'exploitation (en euros)
        out_ENE_ft, //excédent net d'exploitation (en euros)
        out_SA_ft,  //surplus du capital (en euros)
        out_PS_t,   //surplus producteur (en euros)
        out_ES_t,   //surplus de l'Etat (en euros)
        out_Eco,
        out_EcoDCF,
        out_effort, //variables d'effort utilisées lors de la simulation
        out_SRmod,
        out_N_eitQ,
        out_F_itQ,
        out_SSB_etQ,
        out_PQuot_et;


//    VARIABLES  ---------------

//parties des inputs
SEXP    FList, sppList, fleetList, metierList, metierListEco, namDC, t_init, times, trimInt; //fleetList~~ecoList
SEXP NBVF, NBVFM, NBDSF, NBDSFM, dnmsF, dnmsFM, nmsEF;
//dimensions
int     nbT, nbF, nbM, nbMe, nbE, curQ;


//Module de gestion :
int var, trgt, delay, upd, level, gestInd, gestyp;
SEXP    mu_nbds, mu_nbv;  //mulitplicateurs d'effort
int IND_T, IND_F, eTemp, fTemp;  //indicateurs de temps, d'espèces et de flottilles considérés
double PxQ;
int spQ, tacIT;
double expEff;

SEXP m_f, m_fm, m_oth;
double X1, X2, tacLambda;
double *TAC_glob, *Fbar_trgt, *TAC_byFleet;
int *SRInd, *trim, *EcoIndCopy;
double drCopy;
double *effortIni;
//booléens
double *Zoptim, *FOTHoptim;
bool Zoptim_use, FOTHoptim_use;

int *recType1, *recType2, *recType3; //indicateur conditionnant l'utilisation d'un recrutement aléatoire défini par la méthode implémentée RecAlea
bool boolQ;  //indicateur de présence de données d'effort disponible (calcul capturabilité,...)
bool constMM; //indicateur qui détermine si les niveaux métiers de la partie bio et de la partie éco sont les mêmes (utilisation de l'effort par flottille-métier pour calculer la capturabilité dans le module 'Mortalité')
bool fUpdate, dUpdate, cUpdate, pUpdate, eUpdate; //indicateur de mise à jour des variables de calcul
int scen; //application du scénario??
int bhv_active; //application du module report d'effort
int type, boot, nbBoot;
int ecodcf;
//variables intermédiaires par espèces
SEXP eVar, eVarCopy, eVarQ;

//variables intermédiaires flottilles
SEXP fVar;

//Z fixé pour résoudre l'ajustement par flottille (dimension âge)
double *Ztemp;


//méthodes

	//constructeur
    BioEcoPar(SEXP list, SEXP listSpec, SEXP listStochastic, SEXP listScen,
                SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, SEXP Bootstrp, SEXP nbBootstrp,
                SEXP GestInd, SEXP mF, SEXP mOth, SEXP bounds, SEXP TAC, SEXP FBAR, SEXP GestParam, SEXP EcoDcf,
                SEXP EcoInd, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyF, SEXP parBHV, SEXP parQEX);

	//destructeur
    ~BioEcoPar();

    //accesseur d'un élément de l'input
    SEXP getListElement(SEXP list, const char *str);

    //analyse des NAs dans un objet SEXP
    int all_is_na(SEXP object);

    double finite(double value);

    //indices multipliateurs pour les concordances de dimensions
    SEXP iDim(int *dimInput);

    //fontion d'agrégation d'un objet R accompagné de son attribut 'DimCst'
    SEXP aggregObj(SEXP object, SEXP newDim);

    // fonction de ventilation de la mortalité en fonction d'une matrice de données "capture"
    SEXP allocMortality(SEXP mortality, SEXP capture, SEXP captureTot);

    // fonction de calcul de l'indice de capturabilité en fonction de la mortalité par pêche et d'une variable d'effort quelconque
    SEXP calcCapturabilite(SEXP adjustedMortal, SEXP effortIni);

    void RecAlea(SEXP list, SEXP listSto, int ind_t, int type, int *recTyp);

    void SRmod(SEXP list, SEXP listSR, int ind_t, SEXP TypeSR, int *srind);

    void Scenario(SEXP list, SEXP listScen, int ind_t);

    // MODULES :
    //----------

    // Module 'Mortalité par pêche et survie des rejets'
    void Mortalite(SEXP list, int ind_t, SEXP EVAR, int Qt);

    // Module 'Dynamique de population'
    void DynamicPop(SEXP list, int ind_t, SEXP EVAR, int Qt);

    // Module 'Captures, rejets et débarquements'
	void CatchDL(SEXP list, int ind_t, SEXP EVAR, int Qt);

    // Module 'Marché' : 'modCatch'
    void Marche(SEXP list, int ind_t);

    // Module 'Economie'
    void Economic(SEXP list, int ind_t, int adj, int lev, int ue_choice, int oths, int othsFM, int perscCalc, int report, double dr);

    // Module 'Economie' DCF
    void EcoDCF(SEXP list, int ind_t, int adj, int lev, int ue_choice, int oths, int othsFM, int perscCalc, int report, double dr);

    // Module Gestion
    double fxTAC_glob(double mult);

    void Gestion(SEXP list, int ind_t);

    double zbrent(BEfn1 fx, double x1, double x2, double tol);

    void zbrak(BEfn1 fx, double x1, double x2, int n, double xb1[], double xb2[], int *nb);

    double func(double *x);
    double fxTAC_F(double *x);

    int MinimizeF(double **p, double y[], int ndim, double ftol);
    void amoeba(BEfn1_F funk, double **p, double y[], int ndim, double ftol, int *nfunk);
    double amotry(BEfn1_F funk, double **p, double y[], double psum[], int ndim, int ihi, double fac);

    double *NRvector(long nl, long nh);
    double **NRmatrix(long nrl, long nrh, long ncl, long nch);
    void free_vector(double *v, long nl, long nh);
    void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

    int GestionF(double **p, double y[], int ndim, double ftol, int ind_t);

    void FleetBehav(SEXP list, int ind_t, SEXP paramBehav);

    int QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    double fxMaxProf_FT(double *x);

    double fxTAC_F_customCst(double *x);

    double fxTAC_F_customReport(double *x);

    double fxMaxProf_FT_customCst(double *x);

    double fxMaxProf_FT_customReport(double *x);

    int GestionF2(int spp, int ind_t);

    double fxTAC_F_customCst2(double *x);

    int GestionF2report(int spp, int ind_t);

    double fxTAC_F_customReport2(double *x);

    int QuotaExchV2(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    double fxMaxProf_FT_customCstV2(double *x);

    int QuotaExchV2Report(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    double fxMaxProf_FT_customReportV2(double *x);

};

#endif // PARAM_H_INCLUDED

//------------------------------------------------------------------------------------
//constructeur de la classe Param (voir 'param.h' pour les descriptions de variables)
//------------------------------------------------------------------------------------

BioEcoPar::BioEcoPar(SEXP listInput /* object@input */, SEXP listSpec /* object@specific */, SEXP listStochastic /* object@stochastic */,
                     SEXP listScen /* object@scenario */, SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, SEXP Bootstrp, SEXP nbBootstrp,
                     SEXP GestInd, SEXP mF, SEXP mOth, SEXP bounds, SEXP TAC, SEXP FBAR, SEXP GestParam, SEXP EcoDcf,
                     SEXP EcoInd, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyF, SEXP parBHV, SEXP parQEX)
{

PROTECT(listSpec);
effortIni = REAL(getListElement(getListElement(listInput, "Fleet"), "nbds_f"));
PROTECT(list = duplicate(listInput));
PROTECT(FList = getListElement(list, "Fleet"));
PROTECT(sppList = getListElement(listSpec, "Species"));
PROTECT(fleetList = getListElement(listSpec, "Fleet"));
PROTECT(metierList = getListElement(listSpec, "MetierEco"));
PROTECT(metierListEco = getListElement(listSpec, "MetierEco"));
PROTECT(namDC = getListElement(listSpec, "Ages"));
PROTECT(t_init = getListElement(listSpec, "t_init"));
PROTECT(times = getListElement(listSpec, "times"));
PROTECT(trimInt = getListElement(listSpec, "trimInt"));

Zoptim = effortIni; //REAL(getListElement(getListElement(listInput, "Fleet"), "Zoptim")); //
FOTHoptim = effortIni; //REAL(getListElement(getListElement(listInput, "Fleet"), "FOTHoptim")); //
Zoptim_use = false;
FOTHoptim_use = false;

//Rprintf("A1");

nbT = INTEGER(getListElement(listSpec, "NbSteps"))[0];
nbF = length(fleetList);
nbM = length(metierList);
nbMe = length(metierListEco);
nbE = length(sppList);
trim = INTEGER(getListElement(listSpec, "trim"));

ecodcf = INTEGER(EcoDcf)[0];
EcoIndCopy = INTEGER(EcoInd);
drCopy = REAL(dr)[0];


recType1 = INTEGER(RecType1);//vecteur d'entiers de longueur nbE
recType2 = INTEGER(RecType2);//vecteur d'entiers de longueur nbE
recType3 = INTEGER(RecType3);//vecteur d'entiers de longueur nbE
boot = INTEGER(Bootstrp)[0];
nbBoot = INTEGER(nbBootstrp)[0];
boolQ = true;// paramètre voué à rester fixe -> on calculera toujours la capturabilité afin de moduler la mortalité en fonction de l'effort de pêche
constMM = true; //on calcule la capturabilité via l'effort par flottille (incompatibilité des niveaux métiers entre bio et éco)
fUpdate = true;    // à t=0, on remet à jour
dUpdate = true;    //
cUpdate = true;    //
pUpdate = true;    //
eUpdate = true;    //

scen = INTEGER(Scenarii)[0];//true;
bhv_active = INTEGER(getListElement(parBHV, "active"))[0];

gestInd = INTEGER(GestInd)[0];
PROTECT(m_f = duplicate(mF)); if (length(m_f)!=nbF) error("Check dimension of array 'mFleet'!!\n");
PROTECT(m_fm = duplicate(mFM)); if (length(m_fm)!=nbF*nbM) error("Check dimension of array 'mFleetMetier'!!\n");
PROTECT(m_oth = duplicate(mOth)); if (length(m_oth)!=nbE) error("Check dimension of array 'mOth'!!\n");
X1 = REAL(bounds)[0];
X2 = REAL(bounds)[1];
TAC_glob = REAL(TAC);  //à corriger
Fbar_trgt = REAL(FBAR);  //à corriger
TAC_byFleet = REAL(getListElement(TACbyF, "TACbyF"));
tacLambda = REAL(getListElement(TACbyF, "lambda"))[0];
tacIT = INTEGER(getListElement(TACbyF, "IT"))[0];

eTemp = INTEGER(GestParam)[0];//2;
var = INTEGER(GestParam)[1];//1;
trgt = INTEGER(GestParam)[2];//1;
delay = INTEGER(GestParam)[3];//2;
upd = INTEGER(GestParam)[4];//2;
level = INTEGER(GestParam)[5];//2;

Ztemp = NRvector(1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI")));
expEff = REAL(getListElement(TACbyF, "expEff"))[0]; //1;  //facteur d'expansion de l'effort maximal par flottille autorisé dans le cadre de l'optimisation GestionF2 et QuotaExch

int DCFok = INTEGER(EcoDcf)[0];

if (level==0) gestyp = 2; // paramètre voué à rester fixe (type 1 non pertinent)    //
if (level==1) gestyp = 2;//1;                                                           //A INTEGRER EN TANT QUE PARAMETRE DE L'INTERFACE

SRInd = INTEGER(SRind);

//il faut initialiser 'eVar' dans lequel on intègrera toutes les variables intermédiaires à décliner par espèce

SEXP eltE;
PROTECT(eVar = allocVector(VECSXP, nbE));
setAttrib(eVar, R_NamesSymbol, sppList);
for (int e = 0 ; e < nbE ; e++) {
    PROTECT(eltE = allocVector(VECSXP,62));
    if (trim[e]==1) {
     SET_VECTOR_ELT(eltE,0,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,2,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,3,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,44,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,60,allocVector(VECSXP, 4));      //à compléter
     SET_VECTOR_ELT(eltE,10,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,11,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,23,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,24,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,25,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,32,allocVector(VECSXP, 4));
     SET_VECTOR_ELT(eltE,33,allocVector(VECSXP, 4));
    }
    SET_VECTOR_ELT(eVar, e, eltE);
}



PROTECT(fVar = allocVector(VECSXP, 33)); //32= rtbsIni_f & 33=rtbsIni_f_m

//Rprintf("A");

PROTECT(out_F_fmi = allocVector(VECSXP, nbE));
PROTECT(out_Fr_fmi = allocVector(VECSXP, nbE));
PROTECT(out_Z_eit = allocVector(VECSXP, nbE));
PROTECT(out_SRmod = allocVector(VECSXP, nbE));
PROTECT(out_N_eitQ = allocVector(VECSXP, nbE));
PROTECT(out_F_itQ = allocVector(VECSXP, nbE));
PROTECT(out_SSB_etQ = allocVector(VECSXP, nbE));
PROTECT(out_PQuot_et = allocVector(VECSXP, nbE));
//on assigne 1 élément par trimestre pour les espèces concernées
for (int e = 0 ; e < nbE ; e++) {
 if (trim[e]==1) {
    SET_VECTOR_ELT(out_F_fmi,e,allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_Fr_fmi,e,allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_Z_eit,e,allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_N_eitQ,e,allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_F_itQ,e,allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_SSB_etQ,e,allocVector(VECSXP, 4));
}}
PROTECT(out_Fbar_et = allocVector(VECSXP, nbE));
PROTECT(out_N_eit = allocVector(VECSXP, nbE));
PROTECT(out_B_et = allocVector(VECSXP, nbE));
PROTECT(out_SSB_et = allocVector(VECSXP, nbE));
PROTECT(out_C_efmit = allocVector(VECSXP, nbE));
PROTECT(out_C_eit = allocVector(VECSXP, nbE));
PROTECT(out_Y_efmit = allocVector(VECSXP, nbE));
PROTECT(out_Y_eit = allocVector(VECSXP, nbE));
PROTECT(out_D_efmit = allocVector(VECSXP, nbE));
PROTECT(out_L_efmit = allocVector(VECSXP, nbE));
PROTECT(out_L_efmct = allocVector(VECSXP, nbE));
PROTECT(out_L_efmct2 = allocVector(VECSXP, nbE));
PROTECT(out_P_t = allocVector(VECSXP, nbE));
PROTECT(out_Eco = allocVector(VECSXP, 60));
PROTECT(out_EcoDCF = allocVector(VECSXP, 45));
PROTECT(out_effort = allocVector(VECSXP, 4)); //nbv_f, nbv_f_m, nbds_f, nbds_f_m

PROTECT(mu_nbds = allocVector(REALSXP, nbT)); //il reste la mise en forme à opérer
PROTECT(mu_nbv = allocVector(REALSXP, nbT));

double *mu_nbds_t = REAL(mu_nbds); for (int i=0; i<nbT; i++) mu_nbds_t[i] = 0.0; //initialisation
double *mu_nbv_t = REAL(mu_nbv); for (int i=0; i<nbT; i++) mu_nbv_t[i] = 0.0;    //
double *mpond_f = REAL(m_f);
double *mpond_fm = REAL(m_fm);
double *mpond_oth = REAL(m_oth);

//on n'oublie pas de composer l'objet de sortie décrivant les variables 'nbv' et 'nbds'
//SEXP NBVF, NBVFM, NBDSF, NBDSFM, dnmsF, dnmsFM, nmsEF;
PROTECT(NBVF = allocMatrix(REALSXP,nbF,nbT));
PROTECT(NBVFM = alloc3DArray(REALSXP,nbF,nbMe,nbT));
PROTECT(NBDSF = allocMatrix(REALSXP,nbF,nbT));
PROTECT(NBDSFM = alloc3DArray(REALSXP,nbF,nbMe,nbT));
PROTECT(dnmsF = allocVector(VECSXP,2));
PROTECT(dnmsFM = allocVector(VECSXP,3));

SET_VECTOR_ELT(dnmsF, 0, fleetList); SET_VECTOR_ELT(dnmsF, 1, times);
SET_VECTOR_ELT(dnmsFM, 0, fleetList); SET_VECTOR_ELT(dnmsFM, 1, metierListEco); SET_VECTOR_ELT(dnmsFM, 2, times);
setAttrib(NBVF, R_DimNamesSymbol, dnmsF); setAttrib(NBVFM, R_DimNamesSymbol, dnmsFM);
setAttrib(NBDSF, R_DimNamesSymbol, dnmsF); setAttrib(NBDSFM, R_DimNamesSymbol, dnmsFM);

double *NBVf = REAL(NBVF);
double *NBVfm = REAL(NBVFM);
double *NBDSf = REAL(NBDSF);
double *NBDSfm = REAL(NBDSFM);


SEXP ans_PQuot_et;
setAttrib(out_PQuot_et, R_NamesSymbol, sppList);
for (int e = 0 ; e < nbE ; e++) {
    PROTECT(ans_PQuot_et = NEW_NUMERIC(nbT));
    setAttrib(ans_PQuot_et, R_NamesSymbol, times);
    SET_VECTOR_ELT(out_PQuot_et, e, ans_PQuot_et);
}

//Rprintf("B");
for (int it = 0; it < nbT ; it++) {

//    if (it==0) 	{
////        //float Y=;
////        //float **P=&NRmatrix(1,4,1,3);
//        MinimizeF(NRmatrix(1,4,1,3), NRvector(1,4), 3, 0.0001);
//    }

//Rprintf("A");

if (it>=1) {
    RecAlea(list, listStochastic, it, 1, recType1); //IMPORTANT : à effectuer AVANT la procédure d'optimisation
    RecAlea(list, listStochastic, it, 2, recType2);
    RecAlea(list, listStochastic, it, 3, recType3);
}
//Rprintf("B");
    SRmod(list, listSR, it, TypeSR, SRInd); //important : à envoyer avant 'DynamicPop'
//Rprintf("C");
if (scen & it>=1) Scenario(list, listScen, it); //modif MM 16/01/2012
//Rprintf("C");
if (bhv_active & it>=1) FleetBehav(list, it, parBHV);

if (delay<=it & gestInd==1 & it>=1 & ISNA(TAC_byFleet[0])) {
//Rprintf("D");
Gestion(list, it);
//Rprintf("D");
//on met à jour les variables sur lesquelles opère le multiplicateur
//Rprintf("E");
    double *mu_nbds_t2 = REAL(mu_nbds);
    double *mu_nbv_t2 = REAL(mu_nbv);
    double *nbdsFM2 = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF2 = REAL(getListElement(FList, "nbds_f"));
    double *nbvFM2 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF2 = REAL(getListElement(FList, "nbv_f"));

if (level==0) {//niveau flottille

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            if (gestyp==1) nbdsF2[ind_f] = fmax2(nbdsF2[ind_f] + mu_nbds_t2[it]*mpond_f[ind_f],0.0);
            if (gestyp==2) nbdsF2[ind_f] = fmax2(nbdsF2[ind_f]*(1 + mu_nbds_t2[it]*mpond_f[ind_f]),0.0);
        }

        if (var==2) {
            if (gestyp==1) nbvF2[ind_f] = fmax2(nbvF2[ind_f] + mu_nbv_t2[it]*mpond_f[ind_f],0.0);
            if (gestyp==2) nbvF2[ind_f] = fmax2(nbvF2[ind_f]*(1 + mu_nbv_t2[it]*mpond_f[ind_f]),0.0);
        }

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m] + mu_nbds_t2[it]*mpond_f[ind_f],0.0);
                if (gestyp==2) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m]*(1 + mu_nbds_t2[it]*mpond_f[ind_f]),0.0);
            }

            if (var==2) {
                if (gestyp==1) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m] + mu_nbv_t2[it]*mpond_f[ind_f],0.0);
                if (gestyp==2) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m]*(1 + mu_nbv_t2[it]*mpond_f[ind_f]),0.0);
            }

            for (int e = 0 ; e < nbE ; e++){


                double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));


                if (ind_m==0 & ind_f==0) {

                    if (var==1) {
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbds_t2[it]*mpond_oth[e],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbds_t2[it]*mpond_oth[e]),0.0);
                    }

                    if (var==2){
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbv_t2[it]*mpond_oth[e],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbv_t2[it]*mpond_oth[e]),0.0);
                    }

                }

            }
        }
    }

} else { //niveau flottille-métier


     for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double countEff = 0.0;

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m] + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbdsFM2[ind_f+nbF*ind_m] = fmax2(nbdsFM2[ind_f+nbF*ind_m]*(1 + mu_nbds_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
                if (!ISNA(nbdsFM2[ind_f+nbF*ind_m])) countEff = countEff + nbdsFM2[ind_f+nbF*ind_m] - REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)]; // - NBDSfm à l'instant précédent la mise en action du module
            }

            if (var==2) {
                if (gestyp==1) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m] + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m],0.0);
                if (gestyp==2) nbvFM2[ind_f+nbF*ind_m] = fmax2(nbvFM2[ind_f+nbF*ind_m]*(1 + mu_nbv_t2[it]*mpond_fm[ind_f+nbF*ind_m]),0.0);
                if (!ISNA(nbvFM2[ind_f+nbF*ind_m])) countEff = countEff + nbvFM2[ind_f+nbF*ind_m] - REAL(NBVFM)[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)]; // - NBVfm
            }

//ajout MM 22/11/2013
            for (int e = 0 ; e < nbE ; e++) {


                double *Fothi2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));


                if (ind_m==0 & ind_f==0) {

                    if (var==1) {
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbds_t2[it]*mpond_oth[e],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbds_t2[it]*mpond_oth[e]),0.0);
                    }

                    if (var==2){
                        if (gestyp==1)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi] + mu_nbv_t2[it]*mpond_oth[e],0.0);
                        if (gestyp==2)
                        for (int ag = 0; ag < nbi; ag++)
                            Fothi2[ag + it*nbi] = fmax2(Fothi2[ag + it*nbi]*(1+mu_nbv_t2[it]*mpond_oth[e]),0.0);
                    }

                }

            }

//fin ajout MM 22/11/2013

        }


        if (var==1) {
            if (gestyp==1) nbdsF2[ind_f] = fmax2(countEff + REAL(NBDSF)[ind_f + nbF*(delay-1)],0.0); //NBDSf
            if (gestyp==2) nbdsF2[ind_f] = fmax2(countEff + REAL(NBDSF)[ind_f + nbF*(delay-1)],0.0);
        }

        if (var==2) {
            if (gestyp==1) nbvF2[ind_f] = fmax2(countEff + REAL(NBVF)[ind_f + nbF*(delay-1)],0.0); // NBVf
            if (gestyp==2) nbvF2[ind_f] = fmax2(countEff + REAL(NBVF)[ind_f + nbF*(delay-1)],0.0);
        }
     }

}}

if (INTEGER(VECTOR_ELT(parQEX,0))[0]==0 & delay<=it & !ISNA(TAC_byFleet[0]) & !ISNA(TAC_byFleet[it*(nbF+1)]) & it>=1) {  //optimisation TAC par flottille activée si pas de NA en premier élément de colonne par step

    //if (it==1) {

    GestionF2(eTemp, it);//}

  //GestionF(NRmatrix(1,nbF+2,1,nbF+1), NRvector(1,nbF+2), nbF+1, 0.0000001, it);}

}



if (INTEGER(Bootstrp)[0]==0 & INTEGER(VECTOR_ELT(parQEX,0))[0]==1) {
//QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)

    //QuotaExch(REAL(VECTOR_ELT(parQEX,1))[0],REAL(VECTOR_ELT(parQEX,2))[0],REAL(VECTOR_ELT(parQEX,3))[0],
    //                    REAL(VECTOR_ELT(parQEX,4))[0], eTemp, REAL(VECTOR_ELT(parQEX,5))[0], it);

  //if (it==1) {
    QuotaExchV2(REAL(VECTOR_ELT(parQEX,1))[0],REAL(VECTOR_ELT(parQEX,2))[0],REAL(VECTOR_ELT(parQEX,3))[0],
                        REAL(VECTOR_ELT(parQEX,4))[0], eTemp, REAL(VECTOR_ELT(parQEX,5))[0], it);
  //}

}



//if (scen & it>=1) Scenario(list, listScen, it);
//Rprintf("E");
//on remplit l'objet de sortie décrivant les variables 'nbv' et 'nbds'
double *nbdsFM4 = REAL(getListElement(FList, "nbds_f_m"));
double *nbdsF4 = REAL(getListElement(FList, "nbds_f"));
double *nbvFM4 = REAL(getListElement(FList, "nbv_f_m"));
double *nbvF4 = REAL(getListElement(FList, "nbv_f"));

for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

    NBVf[ind_f + nbF*it] = nbvF4[ind_f];
    NBDSf[ind_f + nbF*it] = nbdsF4[ind_f];

    for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

        NBVfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbvFM4[ind_f + nbF*ind_m];
        NBDSfm[ind_f + nbF*ind_m + nbF*nbMe*it] = nbdsFM4[ind_f + nbF*ind_m];
    }
}

//Rprintf("F\n");
//3 modules avec pas de temps différencié au niveau trimestre
Mortalite(list, it, eVar,0);//Rprintf("\nG");//PrintValue(out_Fr_fmi);PrintValue(VECTOR_ELT(eVar,60));
DynamicPop(list, it, eVar,0);//Rprintf("\nH");PrintValue(out_Z_eit);PrintValue(out_N_eitQ);PrintValue(out_N_eit);
CatchDL(list, it, eVar,0);//Rprintf("\nI");//PrintValue(out_Y_eit);
Mortalite(list, it, eVar,1);//Rprintf("\nG");//PrintValue(out_Fr_fmi);PrintValue(VECTOR_ELT(eVar,60));
DynamicPop(list, it, eVar,1);//Rprintf("\nH");PrintValue(out_Z_eit);PrintValue(out_N_eitQ);PrintValue(out_N_eit);
CatchDL(list, it, eVar,1);//Rprintf("\nI");//PrintValue(out_Y_eit);
Mortalite(list, it, eVar,2);//Rprintf("\nG");//PrintValue(out_Fr_fmi);PrintValue(VECTOR_ELT(eVar,60));
DynamicPop(list, it, eVar,2);//Rprintf("\nH");PrintValue(out_Z_eit);PrintValue(out_N_eitQ);PrintValue(out_N_eit);
CatchDL(list, it, eVar,2);//Rprintf("\nI");//PrintValue(out_Y_eit);
Mortalite(list, it, eVar,3);//Rprintf("\nG");//PrintValue(out_Fr_fmi);PrintValue(VECTOR_ELT(eVar,60));
DynamicPop(list, it, eVar,3);//Rprintf("\nH");PrintValue(out_Z_eit);PrintValue(out_N_eitQ);PrintValue(out_N_eit);
CatchDL(list, it, eVar,3);//Rprintf("\nI");//PrintValue(out_Y_eit);PrintValue(out_C_eit);


Marche(list, it);
//Rprintf("J\n");
if (DCFok==0) {

    Economic(list, it, INTEGER(EcoInd)[0], INTEGER(EcoInd)[1], INTEGER(EcoInd)[2], INTEGER(EcoInd)[3], INTEGER(EcoInd)[4],
                       INTEGER(EcoInd)[5], INTEGER(EcoInd)[6], REAL(dr)[0]);
} else {
    EcoDCF(list, it, INTEGER(EcoInd)[0], INTEGER(EcoInd)[1], INTEGER(EcoInd)[2], INTEGER(EcoInd)[3], INTEGER(EcoInd)[4],
                     INTEGER(EcoInd)[5], INTEGER(EcoInd)[6], REAL(dr)[0]);
}

//module gestion : si delay<=it & gestInd==1 & trgt==3 et si Fbar atteint, on rebascule trgt à 2
if ( (delay<=it) & (gestInd==1) & (trgt==3) & (REAL(VECTOR_ELT(out_Fbar_et, eTemp))[it] <= Fbar_trgt[it]) )  trgt=2;


//remise à l'état initial des paramètres du module Gestion si nécessaire (ATTENTION : implique qu'on n'atteint jamais l'effort nul quelle que soit la flottille
//                                                                                 dans le cas contraire, fixer upd=2)


//cas level=0
if (upd==1 & delay<=it & level==0 & gestInd==1) {

    double *mu_nbds_t3 = REAL(mu_nbds);
    double *mu_nbv_t3 = REAL(mu_nbv);
    double *nbdsFM3 = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF3 = REAL(getListElement(FList, "nbds_f"));
    double *nbvFM3 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF3 = REAL(getListElement(FList, "nbv_f"));


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            if (gestyp==1) nbdsF3[ind_f] = nbdsF3[ind_f] - mu_nbds_t3[it] * mpond_f[ind_f];
            if (gestyp==2) nbdsF3[ind_f] = nbdsF3[ind_f] / (1 + mu_nbds_t3[it] * mpond_f[ind_f]);
        }

        if (var==2) {
            if (gestyp==1) nbvF3[ind_f] = nbvF3[ind_f] - mu_nbv_t3[it] * mpond_f[ind_f];
            if (gestyp==2) nbvF3[ind_f] = nbvF3[ind_f] / (1 + mu_nbv_t3[it] * mpond_f[ind_f]);
        }

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

        if (var==1) {
            if (gestyp==1) nbdsFM3[ind_f+nbF*ind_m] = nbdsFM3[ind_f+nbF*ind_m] - mu_nbds_t3[it] * mpond_f[ind_f];
            if (gestyp==2) nbdsFM3[ind_f+nbF*ind_m] = nbdsFM3[ind_f+nbF*ind_m] / (1 + mu_nbds_t3[it] * mpond_f[ind_f]);
        }

        if (var==2) {
            if (gestyp==1) nbvFM3[ind_f+nbF*ind_m] = nbvFM3[ind_f+nbF*ind_m] - mu_nbv_t3[it] * mpond_f[ind_f];
            if (gestyp==2) nbvFM3[ind_f+nbF*ind_m] = nbvFM3[ind_f+nbF*ind_m] / (1 + mu_nbv_t3[it] * mpond_f[ind_f]);
        }


            for (int e = 0 ; e < nbE ; e++){


                double *Fothi3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));

                //attention : valable que si boolQ = true et si reff indépendant de t


                if (ind_m==0 & ind_f==0) {


                            for (int ag = 0; ag < nbi; ag++)

                                if (delay>=1 & it<(nbT-1)) Fothi3[ag + nbi*(it+1)] = Fothi3[ag + nbi*(delay-1)];

                }

            }
        }
    }


}

//cas level=1
if (delay<=it & level==1 & gestInd==1) {

//on remet au niveau de l'instant précédent la mise en action du module Gestion
    double *mu_nbds_t3 = REAL(mu_nbds);             //ajout MM 22/11/2013
    double *mu_nbv_t3 = REAL(mu_nbv);               //ajout MM 22/11/2013
    double *nbdsFM3 = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF3 = REAL(getListElement(FList, "nbds_f"));
    double *nbvFM3 = REAL(getListElement(FList, "nbv_f_m"));
    double *nbvF3 = REAL(getListElement(FList, "nbv_f"));


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){


        if (var==1) nbdsF3[ind_f] = NBDSf[ind_f + nbF*(delay-1)];
        if (var==2) nbvF3[ind_f] = NBVf[ind_f + nbF*(delay-1)];

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) nbdsFM3[ind_f+nbF*ind_m] = NBDSfm[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)];
            if (var==2) nbvFM3[ind_f+nbF*ind_m] = NBVfm[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)];

//ajout MM 22/11/2013

        for (int e = 0 ; e < nbE ; e++){


                double *Fothi3 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int nbi = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));

                //attention : valable que si boolQ = true et si reff indépendant de t


                if (ind_m==0 & ind_f==0) {


                            for (int ag = 0; ag < nbi; ag++)
                                if (delay>=1 & it<(nbT-1)) Fothi3[ag + nbi*(it+1)] = Fothi3[ag + nbi*(delay-1)];


                }

            }

//fin ajout MM 22/11/2013

        }
    }

}



}
//Rprintf("K");
SET_VECTOR_ELT(out_effort, 0, NBVF); SET_VECTOR_ELT(out_effort, 1, NBDSF);
SET_VECTOR_ELT(out_effort, 2, NBVFM); SET_VECTOR_ELT(out_effort, 3, NBDSFM);

const char *nmEf[4] = {"nbv_f","nbds_f","nbv_f_m","nbds_f_m"};
PROTECT(nmsEF = allocVector(STRSXP, 4));
for(int k = 0; k < 4; k++) SET_STRING_ELT(nmsEF, k, mkChar(nmEf[k]));
setAttrib(out_effort, R_NamesSymbol, nmsEF);
//Rprintf("L");
free_vector(Ztemp,1,length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,eTemp))), "modI")));

UNPROTECT(48+nbE+nbE+1);
}

//------------------------------------------------------------------------------------
//destructeur de la classe Param
//------------------------------------------------------------------------------------

//BioEcoPar::~BioEcoPar()
//{
//}



//extern "C" : pour éviter le "name mangling process" qui renomme les fonctions exportées dans les dll.


//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// UTILITIES
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------




//------------------------------------------
// accesseur à un élément d'une liste donnée (list = liste en question , str {caractère} = intitulé de l'élément de la liste)
//------------------------------------------
extern "C" {

SEXP BioEcoPar::getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }

    return elmt;
}

}


//------------------------------------------
// fonction all.is.na (teste si tous les éléments d'un objet sont à NA ou non)
//------------------------------------------

int BioEcoPar::all_is_na(SEXP object)
{
    int res = 1;
    double *robj = REAL(object);

    for (int i = 0 ; i < length(object) ; i++)
        if (!ISNA(robj[i])) {

            res = 0;
            break;
        }

    return res;
}


//------------------------------------------
// fonction qui rend 0 si la valeur est NA, NaN ou Inf
//------------------------------------------

double BioEcoPar::finite(double value)
{
    if (!R_FINITE(value)) return(0.0); else return(value);

}



//------------------------------------------
// fonction de calcul de multiplicateurs d'indices selon les dimensions d'un objet 'array'
// (permet la généricité des équations en assurant la compatibilité des variables en présence,
//  quelles que soient leurs dimensions respectives)
// INPUT : attribut 'DimCst' de l'objet en question
//------------------------------------------

SEXP BioEcoPar::iDim(int *dimInput) {

    SEXP Tab;
    PROTECT(Tab = allocVector(INTSXP,4));
    int *tab = INTEGER(Tab);

    tab[0] = (dimInput[0]>0);
    tab[1] = (dimInput[1]>0)*(1 + (dimInput[0]-1)*(dimInput[0]>0));
    tab[2] = (dimInput[2]>0)*(1 + (dimInput[1]-1)*(dimInput[1]>0))*(1 + (dimInput[0]-1)*(dimInput[0]>0));
    tab[3] = (dimInput[3]>0)*(1 + (dimInput[2]-1)*(dimInput[2]>0))*(1 + (dimInput[1]-1)*(dimInput[1]>0))*(1 + (dimInput[0]-1)*(dimInput[0]>0));

    UNPROTECT(1);
    return(Tab);

}


//------------------------------------------
// fonction d'agrégation d'un objet attribué type ('object'), en fonction d'un nouveau vecteur dimension DimCst ('newDim')
// NB : toutes les valeurs de 'newDim' doivent être au plus égales aux dimensions correspondantes de l'objet pour que la fonction s'applique
//------------------------------------------

extern "C" {

SEXP BioEcoPar::aggregObj(SEXP object, SEXP newDim)
{
    PROTECT(object=object);
    PROTECT(newDim=newDim);

    SEXP ans, dimObj, dimnames, Dim;

    int *dim, *ndim, *rdim;
    double *rans, *robj;

    PROTECT(dimObj = getAttrib(object, install("DimCst")));

    dim = INTEGER(dimObj); ndim = INTEGER(newDim);

    //tests sur les dimensions
    if (dim[0]==0 & dim[1]==0 & dim[2]==0 & dim[3]==0) {  //c'est terminé, rien à agréger

        return(object);

    } else {

        if ( dim[0]<ndim[0] | dim[1]<ndim[1] | dim[2]<ndim[2] | dim[3]<ndim[3] )
        {
            error("Check input dimensions in 'aggregObj'!!\n");
        }

        //on calcule le nombre de cellules à remplir et le nombre de dimensions nulles
        int nbCell = 1, nbDim = 0, incr = 0, incr2 = 0;
        for (int i = 0 ; i < 4 ; i++) {

            if (ndim[i]>0) {
            nbDim++;
            nbCell = nbCell * ndim[i];
            }

        }

        PROTECT(ans = NEW_NUMERIC(nbCell));

        rans = REAL(ans);
        robj = REAL(object);

        //on initialise
        for (int i = 0 ; i < nbCell ; i++) rans[i] = 0.0;

        if (nbDim>0) {

            //en-têtes
            PROTECT(Dim = allocVector(INTSXP,nbDim));
            rdim = INTEGER(Dim);
            PROTECT(dimnames = allocVector(VECSXP,nbDim));
            for (int i = 0 ; i < 4 ; i++) {

                if (ndim[i]>0) {
                    if (GET_DIMNAMES(object)!=R_NilValue) SET_VECTOR_ELT(dimnames, incr, VECTOR_ELT(GET_DIMNAMES(object), incr2)) ;
                    rdim[incr] = ndim[i] ;
                    incr++;}
                if (dim[i]>0) incr2++;

            }

            setAttrib(ans, R_DimSymbol, Dim);
            if (GET_DIMNAMES(object)!=R_NilValue) setAttrib(ans, R_DimNamesSymbol, dimnames);
        }

        setAttrib(ans, install("DimCst"), newDim);

        //multiplicateurs
        int *index_dim = INTEGER(iDim(dim));
        int *index_ndim = INTEGER(iDim(ndim));

        //il ne reste plus qu'à effectuer l'agrégation
        for (int ind_f = 0 ; ind_f < (1 + (dim[0] - 1)*(dim[0]>0)) ; ind_f++)
        for (int ind_m = 0 ; ind_m< (1 + (dim[1] - 1)*(dim[1]>0)) ; ind_m++)
        for (int ind_i = 0 ; ind_i < (1 + (dim[2] - 1)*(dim[2]>0)) ; ind_i++)
        for (int ind_t = 0 ; ind_t < (1 + (dim[3] - 1)*(dim[3]>0)) ; ind_t++)

            if (!ISNA(robj[ind_f*index_dim[0] + ind_m*index_dim[1] + ind_i*index_dim[2] + ind_t*index_dim[3]])) {
                rans[ind_f*index_ndim[0] + ind_m*index_ndim[1] + ind_i*index_ndim[2] + ind_t*index_ndim[3]] =
                rans[ind_f*index_ndim[0] + ind_m*index_ndim[1] + ind_i*index_ndim[2] + ind_t*index_ndim[3]] +
                robj[ind_f*index_dim[0] + ind_m*index_dim[1] + ind_i*index_dim[2] + ind_t*index_dim[3]];
            }

        if (nbDim>0) {

            UNPROTECT(2);
        }

        UNPROTECT(4);
        return (ans);
}

}
}





//------------------------------------------
// fonction de calcul de l'indice de capturabilité en fonction de la mortalité par pêche et d'une variable d'effort donnée : à opérer à t=0
//------------------------------------------


extern "C" {

SEXP BioEcoPar::calcCapturabilite(SEXP adjustedMortal, SEXP effortIni)
{                                  // adjustedMortal est l'output de la fonction 'allocMortality'
                                   // effortIni est en fait l'objet Effort entier --> la restriction au temps initial se fait en interne

    PROTECT(adjustedMortal=adjustedMortal);
    PROTECT(effortIni=effortIni);

    SEXP ans, formatEff, dimCstEff, dimMort, dimEff;

    int *dimE, *dimM, *dimEffort;
    double *rans, *rEff, *rMort;

    PROTECT(dimMort = getAttrib(adjustedMortal, install("DimCst")));
    PROTECT(dimEff = getAttrib(effortIni, install("DimCst")));

    dimM = INTEGER(dimMort); dimE = INTEGER(dimEff);

    //tests sur les dimensions
    if ((dimE[0]!=0 & dimM[0]!=0 & dimE[0]!=dimM[0]) | (dimE[1]!=0 & dimM[1]!=0 & dimE[1]!=dimM[1]) |
        (dimE[2]!=0 & dimM[2]!=0 & dimE[2]!=dimM[2]) | (dimE[3]!=0 & dimM[3]!=0 & dimE[3]!=dimM[3]))
    {
        error("Non_homogeneous dimensions of 'allocMortality' output object and 'Effort' input object!!\n");
    }
    if (dimM[3]!=0)
    {
        warning("Adjusted 'F_fmi' parameter must be constant within time!! Calculation will be done with initial value! \n");
    }

    PROTECT(dimCstEff = allocVector(INTSXP,4));
    dimEffort = INTEGER(dimCstEff);
    for (int i = 0 ; i < 3 ; i++) dimEffort[i] = imin2(dimM[i], dimE[i]);
    dimEffort[3] = dimE[3]; //on n'agrège pas sur le temps puisque on ne considère ensuite que l'instant initial

        PROTECT(formatEff = aggregObj(effortIni, dimCstEff));
        rEff = REAL(formatEff);
        rMort = REAL(adjustedMortal);

        PROTECT(ans = NEW_NUMERIC(length(adjustedMortal)));
        rans = REAL(ans);


            setAttrib(ans, R_DimSymbol, getAttrib(adjustedMortal,R_DimSymbol));
            if (GET_DIMNAMES(adjustedMortal)!=R_NilValue) setAttrib(ans, R_DimNamesSymbol, getAttrib(adjustedMortal,R_DimNamesSymbol));

            setAttrib(ans, install("DimCst"), dimMort);

        //multiplicateurs
        int *fact1 = INTEGER(iDim(dimM));
        int *fact2 = INTEGER(iDim(dimEffort));

        for (int ind_f = 0 ; ind_f < (1 + (dimM[0] - 1)*(dimM[0]>0)) ; ind_f++)
        for (int ind_m = 0 ; ind_m < (1 + (dimM[1] - 1)*(dimM[1]>0)) ; ind_m++)
        for (int ind_i = 0 ; ind_i < (1 + (dimM[2] - 1)*(dimM[2]>0)) ; ind_i++) {

                rans[ind_f*fact1[0] + ind_m*fact1[1] + ind_i*fact1[2]] =
                rMort[ind_f*fact1[0] + ind_m*fact1[1] + ind_i*fact1[2]] /
                rEff[ind_f*fact2[0] + ind_m*fact2[1] + ind_i*fact2[2]] ;

                if (ISNAN(rans[ind_f*fact1[0] + ind_m*fact1[1] + ind_i*fact1[2]]))
                        rans[ind_f*fact1[0] + ind_m*fact1[1] + ind_i*fact1[2]] = 0.0;
        }

        UNPROTECT(7);
        return (ans);

}}






//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// MODULES
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//------------------------------------------
// Module 'Mortalité par pêche et survie des rejets'
//------------------------------------------

extern "C" {

void BioEcoPar::Mortalite(SEXP list, int ind_t, SEXP EVAR, int Qt)
{

SEXP Flist;
PROTECT(Flist = getListElement(list, "Fleet"));

if (fUpdate) {

SEXP    elmt, dimEff,
        dimCst, Dim, dimCst_Sr_e, dimCst_d_efi, dimCst_doth_ei, dimCst_F_efmi, intAge, //dimCst_Capt_emi, dimCst_Capt_ei,
        v_Sr_e, v_d_efi, v_doth_ei, v_F_efmi, v_F_efmi2, formatEff, dimCstEff, //rDim, v_Capt_emi, v_Capt_ei,
        v_nbNav_f, v_nbds_f, dim_nbNavCst, dim_nbdsCst, dim_Finput,// v_Finput, v_fm, v_ventilMoy_f, dim_fmCst, dim_ventilMoyCst,
        fFACT1, fFACT2, fFACT3, fFACT4, fFACT5, fFACT6, Foth_i, Froth_i, dimI, dimIT, DimIT, fFACTsup1, fFACTsup2; //v_ventil2,

SEXP ans_11 = R_NilValue, ans_11l = R_NilValue, dimnames= R_NilValue, dimnamesIT= R_NilValue, rnames= R_NilValue; //, v_Ffm=R_NilValue;

SEXP effort;

//on intègre la donnée d'effort (qu'on l'utilise ensuite pour le calcul de la capturabilité, ou pas)

    if (level==1) {
        PROTECT(effort = getListElement(Flist, "nbds_f_m_tot"));
    } else {
        PROTECT(effort = getListElement(Flist, "nbds_f_tot"));
    }


PROTECT(dimEff = getAttrib(effort, install("DimCst")));

int *dim_Sr_e, *dim_d_efi, *dim_doth_ei, *dim_F_efmi, *dimC, *dimE, *dimM, *dimEffort, *dimF, //*dim_Capt_emi, *dim_Capt_ei
    *dimNav, *dimNbds; //, *dim_fm, *dim_vMoy,*rdim;
int nbI;

double *rans_11, *rans_11l, *r_Sr_e, *r_d_efi, *r_doth_ei, *r_F_efmi, *rEff, *r_nbNav_f,  //*r_fm, *r_ventilMoy_f,
        *r_nbds_f, *r_Foth_i, *r_Froth_i;

//préparation de l'output
if (ind_t==0) {

    PROTECT(rnames = allocVector(STRSXP, nbE));
    if (Qt==0) setAttrib(out_F_fmi, R_NamesSymbol, rnames);
    if (Qt==0) setAttrib(out_Fr_fmi, R_NamesSymbol, rnames);

}

for (int e = 0 ; e < nbE ; e++) {

    //---------
    // calcul de Fr_efmit
    //---------

if (trim[e]==0){
    if (Qt==0) { //on n'opère dans ce cas que si Qt est le premier trimestre


    //--------------------------------------------------------------------------------------   pas de dimension trimestre



                        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                        PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                        nbI = length(getListElement(elmt, "modI"));

                        PROTECT(v_Sr_e = getListElement(elmt, "sr"));
                        PROTECT(v_d_efi = getListElement(elmt, "d_i"));
                        PROTECT(v_doth_ei = getListElement(elmt, "doth_i"));

                        PROTECT(dimCst_Sr_e = getAttrib(v_Sr_e, install("DimCst")));
                        PROTECT(dimCst_d_efi = getAttrib(v_d_efi, install("DimCst")));
                        PROTECT(dimCst_doth_ei = getAttrib(v_doth_ei, install("DimCst")));

                    //---------------------------------------------------------------------
                    // 1ère étape : on ventile la mortalité par les captures si possible
                    //---------------------------------------------------------------------

                    PROTECT(v_F_efmi = getListElement(elmt, "F_fmi"));
                    PROTECT(dim_Finput = getAttrib(v_F_efmi, install("DimCst")));
                    dimF = INTEGER(dim_Finput);


                    if (level==0) PROTECT(v_nbNav_f = getListElement(Flist, "nbv_f")); else PROTECT(v_nbNav_f = getListElement(Flist, "nbv_f_m"));
                    PROTECT(dim_nbNavCst = getAttrib(v_nbNav_f, install("DimCst")));
                    dimNav = INTEGER(dim_nbNavCst);
                    r_nbNav_f = REAL(v_nbNav_f);

                    if (level==0) PROTECT(v_nbds_f = getListElement(Flist, "nbds_f")); else PROTECT(v_nbds_f = getListElement(Flist, "nbds_f_m"));
                    PROTECT(dim_nbdsCst = getAttrib(v_nbds_f, install("DimCst")));
                    dimNbds = INTEGER(dim_nbdsCst);
                    r_nbds_f = REAL(v_nbds_f);


                    //on calcule la mortalité via la capturabilité

                            PROTECT(v_F_efmi2 = calcCapturabilite(v_F_efmi , effort));
                    //PrintValue(v_F_efmi2);
                            //et dans ce cas, l'effort à appliquer à la capturabilité est...

                        dimE = INTEGER(dimEff);
                        dimM = INTEGER(getAttrib(v_F_efmi2, install("DimCst")));
                        PROTECT(dimCstEff = allocVector(INTSXP,4));
                        dimEffort = INTEGER(dimCstEff);
                        for (int i = 0 ; i < 3 ; i++) dimEffort[i] = imin2( dimM[i] , dimE[i] );

                        //on conserve tout de même la dimension temporelle
                        dimEffort[3] = dimE[3];


                            PROTECT(formatEff = aggregObj(effort, dimCstEff));//PrintValue(formatEff);
                            rEff = REAL(formatEff);
                            r_F_efmi = REAL(v_F_efmi2);
                            PROTECT(dimCst_F_efmi = getAttrib(v_F_efmi2, install("DimCst")));

                        //tests sur les dimensions
                        dim_Sr_e = INTEGER(dimCst_Sr_e);
                        if ((dim_Sr_e[0]!=0 & dim_Sr_e[0]!=nbF) | (dim_Sr_e[1]!=0 & dim_Sr_e[1]!=nbM) |
                            (dim_Sr_e[2]!=0 & dim_Sr_e[2]!=nbI) | (dim_Sr_e[3]!=0 & dim_Sr_e[3]!=nbT))
                        {
                            error("Non_homogeneous dimensions in Sr_e element. Check .ini biological parameters files !!\n");
                        }

                        dim_d_efi = INTEGER(dimCst_d_efi);
                        if ((dim_d_efi[0]!=0 & dim_d_efi[0]!=nbF) | (dim_d_efi[1]!=0 & dim_d_efi[1]!=nbM) |
                            (dim_d_efi[2]!=0 & dim_d_efi[2]!=nbI) | (dim_d_efi[3]!=0 & dim_d_efi[3]!=nbT))
                        {
                            error("Non_homogeneous dimensions in d_efi element. Check .ini biological parameters files !!\n");
                        }

                        dim_doth_ei = INTEGER(dimCst_doth_ei);
                        if ((dim_doth_ei[0]!=0 & dim_doth_ei[0]!=nbF) | (dim_doth_ei[1]!=0 & dim_doth_ei[1]!=nbM) |
                            (dim_doth_ei[2]!=0 & dim_doth_ei[2]!=nbI) | (dim_doth_ei[3]!=0 & dim_doth_ei[3]!=nbT))
                        {
                            error("Non_homogeneous dimensions in doth_ei element. Check .ini biological parameters files !!\n");
                        }

                        dim_F_efmi = INTEGER(dimCst_F_efmi);
                        if ((dim_F_efmi[0]!=0 & dim_F_efmi[0]!=nbF) | (dim_F_efmi[1]!=0 & dim_F_efmi[1]!=nbM) |
                            (dim_F_efmi[2]!=0 & dim_F_efmi[2]!=nbI) | (dim_F_efmi[3]!=0 & dim_F_efmi[3]!=nbT))
                        {
                            error("Non_homogeneous dimensions in F_efmi element. Check .ini biological parameters files !!\n");
                        }

                        //on détermine l'attribut Dimension du tableau résultant -> dimCst (on en profite pour compter les dimensions réelles + nombre de cellules)
                        PROTECT(dimCst = allocVector(INTSXP, 4));
                        dimC = INTEGER(dimCst);
                        int count = 0, prod = 1, count2 = 0, count3 = 0;

                        for (int k = 0 ; k < 4 ; k++) {

                            dimC[k] = imax2( imax2(dim_d_efi[k] , dim_F_efmi[k]) , dimEffort[k]);
                            if (k==3) dimC[3] = nbT; //on considère la donnée temporellement
                            if (dimC[k]>0) {
                                count++;
                                prod = prod * dimC[k];
                            }

                        }

                        PROTECT(Dim = allocVector(INTSXP, count));
                        int *dim = INTEGER(Dim);

                        for (int k = 0 ; k < 4 ; k++) {
                            if (dimC[k]>0) {
                                dim[count2] = dimC[k];
                                count2++;
                                }
                        }


                    if (ind_t==0) {

                        //on crée les tableaux résultat pour l'espèce en question
                        PROTECT(ans_11 = NEW_NUMERIC(prod));
                        PROTECT(ans_11l = NEW_NUMERIC(prod));

                        setAttrib(ans_11, R_DimSymbol, Dim);
                        setAttrib(ans_11l, R_DimSymbol, Dim);

                        PROTECT(dimnames = allocVector(VECSXP,count));
                        if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                        if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                        if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                        if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}

                        rans_11 = REAL(ans_11);
                        rans_11l = REAL(ans_11l);

                    } else {

                        rans_11 = REAL(VECTOR_ELT(out_F_fmi, e));
                        rans_11l = REAL(VECTOR_ELT(out_Fr_fmi, e));

                    }

                        r_Sr_e = REAL(v_Sr_e);
                        r_d_efi = REAL(v_d_efi);
                        r_doth_ei = REAL(v_doth_ei);

                            //facteurs des indices pour genériciser le processus

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


                            //équation

                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                        rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                        rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                            r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                            rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] *
                            (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                              r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                        }

                    if (ind_t==0) {

                        setAttrib(ans_11, R_DimNamesSymbol, dimnames); setAttrib(ans_11l, R_DimNamesSymbol, dimnames);
                        setAttrib(ans_11, install("DimCst"), dimCst); setAttrib(ans_11l, install("DimCst"), dimCst);
                        SET_VECTOR_ELT(out_F_fmi, e, ans_11); SET_VECTOR_ELT(out_Fr_fmi, e, ans_11l);
                        SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
                        UNPROTECT(3);

                    }

                    //il ne reste plus qu'à calculer Foth_i en soutrayant de Ftot_i la somme aux âges de la mortalité ventilée non corrigée, et Froth_i en lui appliquant doth_i

                        PROTECT(Foth_i = NEW_NUMERIC(nbI*nbT)); //attention, on considère la mortalité initiale comme étant définie sans dimension temporelle --> à revoir
                        PROTECT(Froth_i = NEW_NUMERIC(nbI*nbT));
                        PROTECT(dimI = allocVector(INTSXP,4));
                        PROTECT(dimIT = allocVector(INTSXP,4));
                        PROTECT(DimIT = allocVector(INTSXP,2));
                        int *rdimI = INTEGER(dimI); rdimI[0] = 0; rdimI[1] = 0; rdimI[2] = nbI; rdimI[3] = dimF[3];
                        int *rdimIT = INTEGER(dimIT); rdimIT[0] = 0; rdimIT[1] = 0; rdimIT[2] = nbI; rdimIT[3] = nbT;
                        int *rDimIT = INTEGER(DimIT); rDimIT[0] = nbI; rDimIT[1] = nbT;

                        PROTECT(dimnamesIT = allocVector(VECSXP,2));
                        SET_VECTOR_ELT(dimnamesIT, 0, intAge);
                        SET_VECTOR_ELT(dimnamesIT, 1, times);

                        setAttrib(Foth_i, R_DimSymbol, DimIT); setAttrib(Froth_i, R_DimSymbol, DimIT);
                        setAttrib(Foth_i, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i, R_DimNamesSymbol, dimnamesIT);
                        setAttrib(Foth_i, install("DimCst"), dimIT); setAttrib(Froth_i, install("DimCst"), dimIT);

                        r_Foth_i = REAL(Foth_i);
                        r_Froth_i = REAL(Froth_i);

                    //if (ind_t==0) {PrintValue(v_F_efmi); PrintValue(aggregObj(v_F_efmi, dimI)); }

                        double *sumFtot = REAL(getListElement(elmt, "F_i"));

                        rdimI[3] = dimC[3];
                        double *sumFr = REAL(aggregObj(ans_11, dimI));
                    //if (ind_t==0) {PrintValue(aggregObj(ans_11, dimI));}

                        if (ind_t==0) { //on initialise
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                                r_Foth_i[ind_i+ind_t*nbI] = fmax2(0.0 , sumFtot[ind_i+ind_t*nbI*(dimF[3]>0)] - sumFr[ind_i+ind_t*nbI*(dimC[3]>0)]); //ON N'INTEGRE PAS DE MORTALITES NEGATIVES

                                if (FOTHoptim_use & e==eTemp) {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                                } else {
                                    r_Foth_i[ind_i+(ind_t+1)*nbI] = r_Foth_i[ind_i+ind_t*nbI];
                                }
                            }
                            //PrintValue(v_F_efmi);PrintValue(dimI);PrintValue(Foth_i);
                        } else {
                           if (ind_t<(nbT-1)) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
                                if (FOTHoptim_use & e==eTemp) {
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


                        //on n'oublie pas d'archiver dans eVar ce dont on aura besoin dans les itérations suivantes
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 0, v_F_efmi2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 1, formatEff);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 2, v_Sr_e);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 3, v_d_efi);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 4, fFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 5, fFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 6, fFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 7, fFACT4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 8, fFACT5);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 61, fFACT6);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 44, Foth_i);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 60, Froth_i);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 50, fFACTsup1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 51, fFACTsup2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 52, v_nbNav_f);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 53, v_nbds_f);


                    //if (indP==1) UNPROTECT(1);
                    UNPROTECT(34);//43);
    }

} else {

     //----------------------------------------------------------------------------------------    dimension trimestre Qt



                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                    PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                    nbI = length(getListElement(elmt, "modI"));

                    PROTECT(v_Sr_e = getListElement(getListElement(elmt, "sr"),CHAR(STRING_ELT(trimInt,Qt))));
                    PROTECT(v_d_efi = getListElement(getListElement(elmt, "d_i"),CHAR(STRING_ELT(trimInt,Qt))));
                    PROTECT(v_doth_ei = getListElement(getListElement(elmt, "doth_i"),CHAR(STRING_ELT(trimInt,Qt))));

                    PROTECT(dimCst_Sr_e = getAttrib(v_Sr_e, install("DimCst")));
                    PROTECT(dimCst_d_efi = getAttrib(v_d_efi, install("DimCst")));
                    PROTECT(dimCst_doth_ei = getAttrib(v_doth_ei, install("DimCst")));

                //---------------------------------------------------------------------
                // 1ère étape : on ventile la mortalité par les captures si possible
                //---------------------------------------------------------------------

                PROTECT(v_F_efmi = getListElement(getListElement(elmt, "F_fmi"),CHAR(STRING_ELT(trimInt,Qt))));
                PROTECT(dim_Finput = getAttrib(v_F_efmi, install("DimCst")));
                dimF = INTEGER(dim_Finput);


                if (level==0) PROTECT(v_nbNav_f = getListElement(Flist, "nbv_f")); else PROTECT(v_nbNav_f = getListElement(Flist, "nbv_f_m"));
                PROTECT(dim_nbNavCst = getAttrib(v_nbNav_f, install("DimCst")));
                dimNav = INTEGER(dim_nbNavCst);
                r_nbNav_f = REAL(v_nbNav_f);

                if (level==0) PROTECT(v_nbds_f = getListElement(Flist, "nbds_f")); else PROTECT(v_nbds_f = getListElement(Flist, "nbds_f_m"));
                PROTECT(dim_nbdsCst = getAttrib(v_nbds_f, install("DimCst")));
                dimNbds = INTEGER(dim_nbdsCst);
                r_nbds_f = REAL(v_nbds_f);


                //on calcule la mortalité via la capturabilité

                        PROTECT(v_F_efmi2 = calcCapturabilite(v_F_efmi , effort));
                //PrintValue(v_F_efmi2);
                        //et dans ce cas, l'effort à appliquer à la capturabilité est...

                    dimE = INTEGER(dimEff);
                    dimM = INTEGER(getAttrib(v_F_efmi2, install("DimCst")));
                    PROTECT(dimCstEff = allocVector(INTSXP,4));
                    dimEffort = INTEGER(dimCstEff);
                    for (int i = 0 ; i < 3 ; i++) dimEffort[i] = imin2( dimM[i] , dimE[i] );

                    //on conserve tout de même la dimension temporelle
                    dimEffort[3] = dimE[3];


                        PROTECT(formatEff = aggregObj(effort, dimCstEff));//PrintValue(formatEff);
                        rEff = REAL(formatEff);
                        r_F_efmi = REAL(v_F_efmi2);
                        PROTECT(dimCst_F_efmi = getAttrib(v_F_efmi2, install("DimCst")));

                    //tests sur les dimensions
                    dim_Sr_e = INTEGER(dimCst_Sr_e);
                    if ((dim_Sr_e[0]!=0 & dim_Sr_e[0]!=nbF) | (dim_Sr_e[1]!=0 & dim_Sr_e[1]!=nbM) |
                        (dim_Sr_e[2]!=0 & dim_Sr_e[2]!=nbI) | (dim_Sr_e[3]!=0 & dim_Sr_e[3]!=nbT))
                    {
                        error("Non_homogeneous dimensions in Sr_e element. Check .ini biological parameters files !!\n");
                    }

                    dim_d_efi = INTEGER(dimCst_d_efi);
                    if ((dim_d_efi[0]!=0 & dim_d_efi[0]!=nbF) | (dim_d_efi[1]!=0 & dim_d_efi[1]!=nbM) |
                        (dim_d_efi[2]!=0 & dim_d_efi[2]!=nbI) | (dim_d_efi[3]!=0 & dim_d_efi[3]!=nbT))
                    {
                        error("Non_homogeneous dimensions in d_efi element. Check .ini biological parameters files !!\n");
                    }

                    dim_doth_ei = INTEGER(dimCst_doth_ei);
                    if ((dim_doth_ei[0]!=0 & dim_doth_ei[0]!=nbF) | (dim_doth_ei[1]!=0 & dim_doth_ei[1]!=nbM) |
                        (dim_doth_ei[2]!=0 & dim_doth_ei[2]!=nbI) | (dim_doth_ei[3]!=0 & dim_doth_ei[3]!=nbT))
                    {
                        error("Non_homogeneous dimensions in doth_ei element. Check .ini biological parameters files !!\n");
                    }

                    dim_F_efmi = INTEGER(dimCst_F_efmi);
                    if ((dim_F_efmi[0]!=0 & dim_F_efmi[0]!=nbF) | (dim_F_efmi[1]!=0 & dim_F_efmi[1]!=nbM) |
                        (dim_F_efmi[2]!=0 & dim_F_efmi[2]!=nbI) | (dim_F_efmi[3]!=0 & dim_F_efmi[3]!=nbT))
                    {
                        error("Non_homogeneous dimensions in F_efmi element. Check .ini biological parameters files !!\n");
                    }

                    //on détermine l'attribut Dimension du tableau résultant -> dimCst (on en profite pour compter les dimensions réelles + nombre de cellules)
                    PROTECT(dimCst = allocVector(INTSXP, 4));
                    dimC = INTEGER(dimCst);
                    int count = 0, prod = 1, count2 = 0, count3 = 0;

                    for (int k = 0 ; k < 4 ; k++) {

                        dimC[k] = imax2( imax2(dim_d_efi[k] , dim_F_efmi[k]) , dimEffort[k]);
                        if (k==3) dimC[3] = nbT; //on considère la donnée temporellement
                        if (dimC[k]>0) {
                            count++;
                            prod = prod * dimC[k];
                        }

                    }

                    PROTECT(Dim = allocVector(INTSXP, count));
                    int *dim = INTEGER(Dim);

                    for (int k = 0 ; k < 4 ; k++) {
                        if (dimC[k]>0) {
                            dim[count2] = dimC[k];
                            count2++;
                            }
                    }


                if (ind_t==0) {

                    //on crée les tableaux résultat pour l'espèce en question
                    PROTECT(ans_11 = NEW_NUMERIC(prod));
                    PROTECT(ans_11l = NEW_NUMERIC(prod));

                    setAttrib(ans_11, R_DimSymbol, Dim);
                    setAttrib(ans_11l, R_DimSymbol, Dim);

                    PROTECT(dimnames = allocVector(VECSXP,count));
                    if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                    if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                    if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                    if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}

                    rans_11 = REAL(ans_11);
                    rans_11l = REAL(ans_11l);

                } else {

                    rans_11 = REAL(VECTOR_ELT(VECTOR_ELT(out_F_fmi, e),Qt));
                    rans_11l = REAL(VECTOR_ELT(VECTOR_ELT(out_Fr_fmi, e),Qt));

                }

                    r_Sr_e = REAL(v_Sr_e);
                    r_d_efi = REAL(v_d_efi);
                    r_doth_ei = REAL(v_doth_ei);

                        //facteurs des indices pour genériciser le processus

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


                        //équation

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]];

                    rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        rEff[ind_f*fFact5[0] + ind_m*fFact5[1] + ind_i*fFact5[2] + ind_t*fFact5[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                          r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                    }

                if (ind_t==0) {

                    setAttrib(ans_11, R_DimNamesSymbol, dimnames); setAttrib(ans_11l, R_DimNamesSymbol, dimnames);
                    setAttrib(ans_11, install("DimCst"), dimCst); setAttrib(ans_11l, install("DimCst"), dimCst);
                    SET_VECTOR_ELT(VECTOR_ELT(out_F_fmi, e),Qt, ans_11); SET_VECTOR_ELT(VECTOR_ELT(out_Fr_fmi, e),Qt, ans_11l);
                    if (Qt==0) SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
                    UNPROTECT(3);

                }


                //il ne reste plus qu'à calculer Foth_i en soutrayant de Ftot_i la somme aux âges de la mortalité ventilée non corrigée, et Froth_i en lui appliquant doth_i

                    PROTECT(Foth_i = NEW_NUMERIC(nbI*nbT)); //attention, on considère la mortalité initiale comme étant définie sans dimension temporelle --> à revoir
                    PROTECT(Froth_i = NEW_NUMERIC(nbI*nbT));
                    PROTECT(dimI = allocVector(INTSXP,4));
                    PROTECT(dimIT = allocVector(INTSXP,4));
                    PROTECT(DimIT = allocVector(INTSXP,2));
                    int *rdimI = INTEGER(dimI); rdimI[0] = 0; rdimI[1] = 0; rdimI[2] = nbI; rdimI[3] = dimF[3];
                    int *rdimIT = INTEGER(dimIT); rdimIT[0] = 0; rdimIT[1] = 0; rdimIT[2] = nbI; rdimIT[3] = nbT;
                    int *rDimIT = INTEGER(DimIT); rDimIT[0] = nbI; rDimIT[1] = nbT;

                    PROTECT(dimnamesIT = allocVector(VECSXP,2));
                    SET_VECTOR_ELT(dimnamesIT, 0, intAge);
                    SET_VECTOR_ELT(dimnamesIT, 1, times);

                    setAttrib(Foth_i, R_DimSymbol, DimIT); setAttrib(Froth_i, R_DimSymbol, DimIT);
                    setAttrib(Foth_i, R_DimNamesSymbol, dimnamesIT); setAttrib(Froth_i, R_DimNamesSymbol, dimnamesIT);
                    setAttrib(Foth_i, install("DimCst"), dimIT); setAttrib(Froth_i, install("DimCst"), dimIT);

                    r_Foth_i = REAL(Foth_i);
                    r_Froth_i = REAL(Froth_i);

                //if (ind_t==0) {PrintValue(v_F_efmi); PrintValue(aggregObj(v_F_efmi, dimI)); }

                    double *sumFtot = REAL(getListElement(getListElement(elmt, "F_i"),CHAR(STRING_ELT(trimInt,Qt))));

                    rdimI[3] = dimC[3];
                    double *sumFr = REAL(aggregObj(ans_11, dimI));
                //if (ind_t==0) {PrintValue(aggregObj(ans_11, dimI));}

                    if (ind_t==0) { //on initialise
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {
                            r_Foth_i[ind_i+ind_t*nbI] = fmax2(0.0 , sumFtot[ind_i+ind_t*nbI*(dimF[3]>0)] - sumFr[ind_i+ind_t*nbI*(dimC[3]>0)]); //ON N'INTEGRE PAS DE MORTALITES NEGATIVES
                            r_Foth_i[ind_i+(ind_t+1)*nbI] = r_Foth_i[ind_i+ind_t*nbI];
                        }
                        //PrintValue(v_F_efmi);PrintValue(dimI);PrintValue(Foth_i);
                    } else {
                       if (ind_t<(nbT-1)) {
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_i[ind_i+(ind_t+1)*nbI] = r_Foth_i[ind_i+ind_t*nbI];
                       }
                    }


                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                            r_Froth_i[ind_i+ind_t*nbI] = r_Foth_i[ind_i+ind_t*nbI] *
                                (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                                r_doth_ei[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);


                    //on n'oublie pas d'archiver dans eVar ce dont on aura besoin dans les itérations suivantes
                    SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 0),Qt, v_F_efmi2);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 1, formatEff);
                    SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 2),Qt, v_Sr_e);
                    SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 3),Qt, v_d_efi);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 4, fFACT1);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 5, fFACT2);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 6, fFACT3);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 7, fFACT4);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 8, fFACT5);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 61, fFACT6);
                    SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44),Qt, Foth_i);
                    SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60),Qt, Froth_i);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 50, fFACTsup1);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 51, fFACTsup2);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 52, v_nbNav_f);
                    SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 53, v_nbds_f);


                UNPROTECT(34);//43);



    }
}

if (Qt==3) fUpdate = false;

if (ind_t==0) UNPROTECT(1);
UNPROTECT(2);

} else {

for (int e = 0 ; e < nbE ; e++) {

    if (trim[e]==0) { //-----------------------------------------------------------------------    pas de dimension trimestre

            if (Qt==0) {
                    int nbI = length(VECTOR_ELT(namDC,e));

                    SEXP elmt;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                    double *rans_11 = REAL(VECTOR_ELT(out_F_fmi, e));
                    double *rans_11l = REAL(VECTOR_ELT(out_Fr_fmi, e));
                    double *r_F_efmi = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 0));
                    double *r_Sr_e = REAL(getListElement(elmt, "sr"));
                    double *r_d_efi = REAL(getListElement(elmt, "d_i"));
                    double *r_doth_ei = REAL(getListElement(elmt, "doth_i"));
                    double *r_nbv_f, *r_nbds_f;
                    if (level==0) r_nbv_f = REAL(getListElement(Flist, "nbv_f")); else r_nbv_f = REAL(getListElement(Flist, "nbv_f_m"));
                    if (level==0) r_nbds_f = REAL(getListElement(Flist, "nbds_f")); else r_nbds_f = REAL(getListElement(Flist, "nbds_f_m"));
                    double *r_Foth_it = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));
                    double *r_Froth_it = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));



                    int *fFact1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 4)),
                        *fFact2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 5)),
                        *fFact3 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 6)),
                        *fFact4 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 7)),
                        *fFact6 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 61)),
                        *fFactSup1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 50)),
                        *fFactSup2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 51));



                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];


                    rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                        r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                    }

                    if (ind_t<(nbT-1)) {
                        if (FOTHoptim_use & e==eTemp) {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_it[ind_i+(ind_t+1)*nbI] = FOTHoptim[ind_i+(ind_t+1)*nbI];
                        } else {
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_it[ind_i+(ind_t+1)*nbI] = r_Foth_it[ind_i+ind_t*nbI];   //à modifier quand on considèrera une mortalité "autres" variable
                        }
                    }

                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                        r_Froth_it[ind_i+ind_t*nbI] = r_Foth_it[ind_i+ind_t*nbI] *
                            (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            r_doth_ei[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);


                UNPROTECT(1);
            }

} else {


     //--------------------------------------------------------------------------------------   dimension trimestre



                    int nbI = length(VECTOR_ELT(namDC,e));

                    SEXP elmt;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                    double *rans_11 = REAL(VECTOR_ELT(VECTOR_ELT(out_F_fmi, e),Qt));
                    double *rans_11l = REAL(VECTOR_ELT(VECTOR_ELT(out_Fr_fmi, e),Qt));
                    double *r_F_efmi = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 0),Qt));
                    double *r_Sr_e = REAL(getListElement(getListElement(elmt, "sr"),CHAR(STRING_ELT(trimInt,Qt))));
                    double *r_d_efi = REAL(getListElement(getListElement(elmt, "d_i"),CHAR(STRING_ELT(trimInt,Qt))));
                    double *r_doth_ei = REAL(getListElement(getListElement(elmt, "doth_i"),CHAR(STRING_ELT(trimInt,Qt))));
                    double *r_nbv_f, *r_nbds_f;
                    if (level==0) r_nbv_f = REAL(getListElement(Flist, "nbv_f")); else r_nbv_f = REAL(getListElement(Flist, "nbv_f_m"));
                    if (level==0) r_nbds_f = REAL(getListElement(Flist, "nbds_f")); else r_nbds_f = REAL(getListElement(Flist, "nbds_f_m"));
                    double *r_Foth_it = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44),Qt));
                    double *r_Froth_it = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60),Qt));

                    int *fFact1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 4)),
                        *fFact2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 5)),
                        *fFact3 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 6)),
                        *fFact4 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 7)),
                        *fFact6 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 61)),
                        *fFactSup1 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 50)),
                        *fFactSup2 = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 51));

                    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                    for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                    rans_11[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]];

                    rans_11l[ind_f*fFact1[0] + ind_m*fFact1[1] + ind_i*fFact1[2] + ind_t*fFact1[3]] =
                        r_F_efmi[ind_f*fFact4[0] + ind_m*fFact4[1] + ind_i*fFact4[2] + ind_t*fFact4[3]] *
                        r_nbv_f[ind_f*fFactSup1[0] + ind_m*fFactSup1[1] + ind_i*fFactSup1[2] + ind_t*fFactSup1[3]] *
                        r_nbds_f[ind_f*fFactSup2[0] + ind_m*fFactSup2[1] + ind_i*fFactSup2[2] + ind_t*fFactSup2[3]] *
                        (1 - r_Sr_e[ind_f*fFact3[0] + ind_m*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                        r_d_efi[ind_f*fFact2[0] + ind_m*fFact2[1] + ind_i*fFact2[2] + ind_t*fFact2[3]]);

                    }

                    if (ind_t<(nbT-1)) {

                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) r_Foth_it[ind_i+(ind_t+1)*nbI] = r_Foth_it[ind_i+ind_t*nbI];   //à modifier quand on considèrera une mortalité "autres" variable

                    }

                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)
                        r_Froth_it[ind_i+ind_t*nbI] = r_Foth_it[ind_i+ind_t*nbI] *
                            (1 - r_Sr_e[0*fFact3[0] + 0*fFact3[1] + ind_i*fFact3[2] + ind_t*fFact3[3]] *
                            r_doth_ei[0*fFact6[0] + 0*fFact6[1] + ind_i*fFact6[2] + ind_t*fFact6[3]]);


                UNPROTECT(1);

    }


}

}
UNPROTECT(1);

}}










//------------------------------------------
// Module 'Dynamique de population' : 'out_F_fmi' = output de la fonction 'Mortalite' MAJ 27/09/2010 ajout de l'output SSB_et
//------------------------------------------

extern "C" {

void BioEcoPar::DynamicPop(SEXP list, int ind_t, SEXP EVAR, int Qt)
{

if (dUpdate) {

    SEXP    elmt, dFACT1, dFACT2, dFACT3, dFACT4, dFACT5, dFACT6, dFACT7, dFACT8, dFACT9, dFACT10,
            dimCst1, dimCst2, dimCst3, dimCst4, Dim1, Dim2, Dim3, Dim4,
            dimCst_Fr_efmit, dimCst_M_ei, dimCst_w_ei, dimCst_N_ei0, dimCst_N_e0t, dimCst_mat_ei,
            intAge, v_Fr_efmit, v_F_efmit, v_M_ei, v_w_ei, v_N_ei0, v_N_e0t, v_mat_ei, v_Fbar;

    SEXP dimnames1=R_NilValue, dimnames2=R_NilValue, dimnames3=R_NilValue, dimnames4=R_NilValue, rnames_Esp=R_NilValue;

    int *dim_Fr_efmit, *dim_M_ei, *dim_w_ei, *dim_N_ei0, *dim_N_e0t, *dim_mat_ei, *dimC1, *dimC2, *dimC3, *dimC4, *dim1, *dim2, *dim3, *dim4;
    int nbI;

    double *rans_Z_eit, *rans_N_eit, *rans_B_et, *rans_Fbar_et, *rans_SSB_et, *r_Fr_efmit,  *r_F_efmit, *r_Fbar, *r_M_ei, *r_w_ei,
                *r_N_ei0, *r_N_e0t, *r_mat_ei, *r_N_eitQ, *r_F_itQ, *r_SSB_etQ;

    SEXP ans_Z_eit=R_NilValue, ans_N_eit=R_NilValue, ans_Fbar_et=R_NilValue, ans_B_et=R_NilValue, ans_SSB_et=R_NilValue,
                                ans_N_eitQ=R_NilValue, ans_F_itQ=R_NilValue, ans_SSB_etQ=R_NilValue;

if (ind_t==0) {


    //à t=0, préparation des outputs

    PROTECT(rnames_Esp = allocVector(STRSXP, nbE));

    if (Qt==0) setAttrib(out_Z_eit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_N_eit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_Fbar_et, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_B_et, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_SSB_et, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_N_eitQ, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_F_itQ, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_SSB_etQ, R_NamesSymbol, rnames_Esp);
}


    for (int e = 0 ; e < nbE ; e++) {//Rprintf("G1");

     if (trim[e]==0) { //-------------------------------------------------------- pas de dimension trimestre
//Rprintf("G2");
        if (Qt==0) {
                                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                                    nbI = length(getListElement(elmt, "modI"));

                                    PROTECT(v_M_ei = getListElement(elmt, "M_i"));
                                    PROTECT(v_w_ei = getListElement(elmt, "wStock_i"));
                                    PROTECT(v_mat_ei = getListElement(elmt, "mat_i"));
                                    PROTECT(v_Fr_efmit = getListElement(out_Fr_fmi, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_F_efmit = getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                                    PROTECT(v_N_ei0 = getListElement(elmt, "N_it0"));
                                    PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));

                                    PROTECT(v_Fbar = getListElement(elmt, "Fbar"));
                                    PROTECT(dimCst_M_ei = getAttrib(v_M_ei, install("DimCst")));
                                    PROTECT(dimCst_w_ei = getAttrib(v_w_ei, install("DimCst")));
                                    PROTECT(dimCst_mat_ei = getAttrib(v_mat_ei, install("DimCst")));
                                    PROTECT(dimCst_Fr_efmit = getAttrib(v_Fr_efmit, install("DimCst")));
                                    PROTECT(dimCst_N_ei0 = getAttrib(v_N_ei0, install("DimCst")));
                                    PROTECT(dimCst_N_e0t = getAttrib(v_N_e0t, install("DimCst")));

                                    //tests sur les dimensions :
                                    dim_M_ei = INTEGER(dimCst_M_ei);
                                    if ((dim_M_ei[0]!=0) | (dim_M_ei[1]!=0) |
                                        (dim_M_ei[2]!=0 & dim_M_ei[2]!=nbI) | (dim_M_ei[3]!=0 & dim_M_ei[3]!=nbT)) //on laisse une ouverture pour un indice temporel
                                    {
                                        error("Non_homogeneous dimensions in M_ei element. Check .ini biological parameters files !!\n");
                                    }

                                    dim_w_ei = INTEGER(dimCst_w_ei);
                                    if ((dim_w_ei[0]!=0) | (dim_w_ei[1]!=0) |
                                        (dim_w_ei[2]!=0 & dim_w_ei[2]!=nbI) | (dim_w_ei[3]!=0))
                                    {
                                        error("Non_homogeneous dimensions in w_ei element. Check .ini biological parameters files !!\n");
                                    }

                                    dim_mat_ei = INTEGER(dimCst_mat_ei);
                                    if ((dim_mat_ei[0]!=0) | (dim_mat_ei[1]!=0) |
                                        (dim_mat_ei[2]!=0 & dim_mat_ei[2]!=nbI) | (dim_mat_ei[3]!=0))
                                    {
                                        error("Non_homogeneous dimensions in mat_ei element. Check .ini biological parameters files !!\n");
                                    }

                                    dim_Fr_efmit = INTEGER(dimCst_Fr_efmit);
                                    if ((dim_Fr_efmit[0]!=0 & dim_Fr_efmit[0]!=nbF) | (dim_Fr_efmit[1]!=0 & dim_Fr_efmit[1]!=nbM) |
                                        (dim_Fr_efmit[2]!=0 & dim_Fr_efmit[2]!=nbI) | (dim_Fr_efmit[3]!=0 & dim_Fr_efmit[3]!=nbT))
                                    {
                                        error("Non_homogeneous dimensions in Fr_efmit element. Check .ini biological parameters files !!\n");
                                    }

                                    dim_N_ei0 = INTEGER(dimCst_N_ei0);
                                    if ((dim_N_ei0[0]!=0) | (dim_N_ei0[1]!=0) |
                                        (dim_N_ei0[2]!=0 & dim_N_ei0[2]!=nbI)) // | (dim_N_ei0[3]!=0)) --> peu importe, on ne prendra de toute façon que la donnée à t0
                                    {
                                        error("Non_homogeneous dimensions in N_ei0 element. Check .ini biological parameters files !!\n");
                                    }

                                    dim_N_e0t = INTEGER(dimCst_N_e0t);
                                    if ((dim_N_e0t[0]!=0) | (dim_N_e0t[1]!=0) |
                                        (dim_N_e0t[2]!=0) | (dim_N_e0t[3]!=0 & dim_N_e0t[3]!=nbT))
                                    {
                                        error("Non_homogeneous dimensions in N_e0t element. Check .ini biological parameters files !!\n");
                                    }
//Rprintf("G3");
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // initialisation de out_Z_eit
                                    //---------

                                    //on détermine l'attribut Dimension de Z_eit
                                    PROTECT(dimCst1 = allocVector(INTSXP, 4));
                                    dimC1 = INTEGER(dimCst1);
                                    dimC1[0] = 0 ; dimC1[1] = 0 ; dimC1[2] = imax2(dim_M_ei[2] , dim_Fr_efmit[2]);
                                    dimC1[3] = imax2(dim_M_ei[3] , dim_Fr_efmit[3]);

                                    int count = 0, prod = 1, count2 = 0, count3 = 0, count4 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC1[k]>0) {
                                            count++;
                                            prod = prod * dimC1[k];
                                        }
                                    }

                                    PROTECT(Dim1 = allocVector(INTSXP, count));
                                    dim1 = INTEGER(Dim1);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC1[k]>0) {
                                            dim1[count2] = dimC1[k];
                                            count2++;
                                            }
                                    }

                            if (ind_t==0){
                                    //on crée le tableau résultat pour l'espèce en question
                                    PROTECT(ans_Z_eit = NEW_NUMERIC(prod));
                                    setAttrib(ans_Z_eit, R_DimSymbol, Dim1);

                                    PROTECT(dimnames1 = allocVector(VECSXP,count));
                                    if (dimC1[0]>0) {SET_VECTOR_ELT(dimnames1, count3, fleetList) ; count3++;}
                                    if (dimC1[1]>0) {SET_VECTOR_ELT(dimnames1, count3, metierList) ; count3++;}
                                    if (dimC1[2]>0) {SET_VECTOR_ELT(dimnames1, count3, intAge) ; count3++;}
                                    if (dimC1[3]>0) {SET_VECTOR_ELT(dimnames1, count3, times) ; count3++;}

                                    rans_Z_eit = REAL(ans_Z_eit);

                            } else {

                                    rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit, e));

                            }
//Rprintf("G4");
                                    r_Fr_efmit = REAL(v_Fr_efmit);
                                    r_F_efmit = REAL(v_F_efmit);
                                    r_M_ei = REAL(v_M_ei);
                                    r_Fbar = REAL(v_Fbar);

                                    //facteurs des indices
                                    PROTECT(dFACT1 = iDim(dimC1));
                                    PROTECT(dFACT2 = iDim(dim_Fr_efmit));
                                    PROTECT(dFACT3 = iDim(dim_M_ei));

                                    int *fact1_D = INTEGER(dFACT1);
                                    int *fact2_D = INTEGER(dFACT2);
                                    int *fact3_D = INTEGER(dFACT3);

                                    double *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));

                                    double fmax = 0.0, sumWt = 0.0;

                                    //équation
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                        double temp = 0.0, tempCap = 0.0;

                                        for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                        for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                        if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                            tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                        }

                                    if (Zoptim_use & e==eTemp) {

                                     rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                        Zoptim[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];


                                    } else {

                                        if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                        } else {
                                          rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                            temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                            r_Froth_i[ind_i + nbI*ind_t];

                                           //on initialise aussi Ztemp (attention : indexé à partir de 1)

                                           Ztemp[ind_i+1] = rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];
                                        }
                                    }

                                    //on en profite pour calculer Fbar

                                    fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                                    sumWt = sumWt + r_Fbar[ind_i];

                                    }

//Rprintf("G5");
                            if (ind_t==0) {

                                    setAttrib(ans_Z_eit, R_DimNamesSymbol, dimnames1);
                                    setAttrib(ans_Z_eit, install("DimCst"), dimCst1);

                                    SET_VECTOR_ELT(out_Z_eit, e, ans_Z_eit);
                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 9, dimCst_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 10, v_Fr_efmit);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 11, v_M_ei);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 12, dFACT1);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 13, dFACT2);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 14, dFACT3);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 59, v_Fbar);



                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de N_eit
                                    //---------

                                    //on détermine l'attribut Dimension de N_eit
                                    PROTECT(dimCst2 = allocVector(INTSXP, 4));
                                    dimC2 = INTEGER(dimCst2);

                                    dimC2[0] = 0 ; dimC2[1] = 0 ; dimC2[2] = imax2(dimC1[2] , dim_N_ei0[2]) ; dimC2[3] = nbT;

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            count++;
                                            prod = prod * dimC2[k];
                                        }
                                    }

                                    PROTECT(Dim2 = allocVector(INTSXP, count));
                                    dim2 = INTEGER(Dim2);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC2[k]>0) {
                                            dim2[count2] = dimC2[k];
                                            count2++;
                                            }
                                    }
//Rprintf("G6");

                            if (ind_t==0) {

                                    //on crée le tableau résultat pour l'espèce en question
                                    PROTECT(ans_N_eit = NEW_NUMERIC(prod));
                                    setAttrib(ans_N_eit, R_DimSymbol, Dim2);

                                    PROTECT(dimnames2 = allocVector(VECSXP,count));
                                    if (dimC2[0]>0) {SET_VECTOR_ELT(dimnames2, count3, fleetList) ; count3++;}
                                    if (dimC2[1]>0) {SET_VECTOR_ELT(dimnames2, count3, metierList) ; count3++;}
                                    if (dimC2[2]>0) {SET_VECTOR_ELT(dimnames2, count3, intAge) ; count3++;}
                                    if (dimC2[3]>0) {SET_VECTOR_ELT(dimnames2, count3, times) ; count3++;}

                                    rans_N_eit = REAL(ans_N_eit);

                            } else {

                                    rans_N_eit = REAL(VECTOR_ELT(out_N_eit, e));

                            }

                                    r_N_ei0 = REAL(v_N_ei0);
                                    r_N_e0t = REAL(v_N_e0t);

                                    //facteurs des indices
                                    PROTECT(dFACT4 = iDim(dimC2));
                                    PROTECT(dFACT5 = iDim(dim_N_ei0));
                                    PROTECT(dFACT6 = iDim(dim_N_e0t));

                                    int *fact4_D = INTEGER(dFACT4);
                                    int *fact5_D = INTEGER(dFACT5);
                                    int *fact6_D = INTEGER(dFACT6);

                                    //équation

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                            if (ind_i == 0) { //recrutement

                                                if (SRInd[e]==1 & ind_t>0) {

                                                    if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                            REAL(VECTOR_ELT(out_SRmod,e))[ind_t];

                                                    } else {

                                                     if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial défini

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                                    }}


                                                } else {

                                                    if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial défini

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                                    }
                                                }

                                            } else {

                                                if (ind_t == 0) {

                                                    rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                                } else {

                                                    if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                    } else {

                                                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                          rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                          exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                    }
                                                }
                                            }
                                    }


                            if (ind_t==0) {

                                    setAttrib(ans_N_eit, R_DimNamesSymbol, dimnames2);
                                    setAttrib(ans_N_eit, install("DimCst"), dimCst2);

                                    SET_VECTOR_ELT(out_N_eit, e, ans_N_eit);

                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 15, v_N_ei0);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 16, v_N_e0t);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 17, dFACT4);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 18, dFACT5);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 19, dFACT6);
//Rprintf("G7");

                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de SSB_et
                                    //---------

                                    //on détermine l'attribut Dimension de SSB_et
                                    PROTECT(dimCst4 = allocVector(INTSXP, 4));
                                    dimC4 = INTEGER(dimCst4);

                                    dimC4[0] = 0 ; dimC4[1] = 0 ; dimC4[2] = 0 ; dimC4[3] = dimC2[3];

                                    count = 0 ; prod = 1 ; count2 = 0 ; count4 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC4[k]>0) {
                                            count++;
                                            prod = prod * dimC4[k];
                                        }
                                    }


                                    PROTECT(Dim4 = allocVector(INTSXP, count));
                                    dim4 = INTEGER(Dim4);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC4[k]>0) {
                                            dim4[count2] = dimC4[k];
                                            count2++;
                                            }
                                    }

                            if (ind_t==0) {
                                    //on crée le tableau résultat pour l'espèce en question (on en profite pour faire de même avec Fbar --> même dimension)
                                    PROTECT(ans_SSB_et = NEW_NUMERIC(prod));
                                    PROTECT(ans_Fbar_et = NEW_NUMERIC(prod));

                                    if (count>0) { //valable seulement si SSB_et n'est pas seulement un scalaire

                                    setAttrib(ans_SSB_et, R_DimSymbol, Dim4);
                                    setAttrib(ans_Fbar_et, R_DimSymbol, Dim4);

                                    }

                                    PROTECT(dimnames4 = allocVector(VECSXP,count));
                                    if (dimC4[0]>0) {SET_VECTOR_ELT(dimnames4, count4, fleetList) ; count4++;}
                                    if (dimC4[1]>0) {SET_VECTOR_ELT(dimnames4, count4, metierList) ; count4++;}
                                    if (dimC4[2]>0) {SET_VECTOR_ELT(dimnames4, count4, intAge) ; count4++;}
                                    if (dimC4[3]>0) {SET_VECTOR_ELT(dimnames4, count4, times) ; count4++;}


                                    rans_SSB_et = REAL(ans_SSB_et);
                                    rans_Fbar_et = REAL(ans_Fbar_et);

                            } else {


                                    rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et, e));
                                    rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et, e));

                            }

                                    r_w_ei = REAL(v_w_ei);
                                    r_mat_ei = REAL(v_mat_ei);

                                    //facteurs des indices
                                    PROTECT(dFACT9 = iDim(dimC4));
                                    PROTECT(dFACT8 = iDim(dim_w_ei));
                                    PROTECT(dFACT10 = iDim(dim_mat_ei));

                                    int *fact9_D = INTEGER(dFACT9);
                                    int *fact8_D = INTEGER(dFACT8);
                                    int *fact10_D = INTEGER(dFACT10);

                                    //équation

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                                temp = temp +
                                                 rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                                 r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                            rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                            rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;

                                        }

                            if(ind_t==0) {

                                    if (count>0) setAttrib(ans_SSB_et, R_DimNamesSymbol, dimnames4);
                                    setAttrib(ans_SSB_et, install("DimCst"), dimCst4);

                                    SET_VECTOR_ELT(out_SSB_et, e, ans_SSB_et);

                                    if (count>0) setAttrib(ans_Fbar_et, R_DimNamesSymbol, dimnames4);
                                    setAttrib(ans_Fbar_et, install("DimCst"), dimCst4);

                                    SET_VECTOR_ELT(out_Fbar_et, e, ans_Fbar_et);

                                    SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 48, dFACT9);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 49, dFACT10);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 47, v_mat_ei);

//Rprintf("G8");


                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////
                                /////////////////////////////////////////////////////////////////////////////////////////////////

                                    //---------
                                    // calcul de B_et
                                    //---------

                                    //on détermine l'attribut Dimension de B_et
                                    PROTECT(dimCst3 = allocVector(INTSXP, 4));
                                    dimC3 = INTEGER(dimCst3);

                                    dimC3[0] = 0 ; dimC3[1] = 0 ; dimC3[2] = 0 ; dimC3[3] = dimC2[3];

                                    count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC3[k]>0) {
                                            count++;
                                            prod = prod * dimC3[k];
                                        }
                                    }

                                    PROTECT(Dim3 = allocVector(INTSXP, count));
                                    dim3 = INTEGER(Dim3);

                                    for (int k = 0 ; k < 4 ; k++) {

                                        if (dimC3[k]>0) {
                                            dim3[count2] = dimC3[k];
                                            count2++;
                                            }
                                    }

                            if (ind_t==0) {
                                    //on crée le tableau résultat pour l'espèce en question
                                    PROTECT(ans_B_et = NEW_NUMERIC(prod));

                                    if (count>0) { //valable seulement si B_et n'est pas seulement un scalaire

                                    setAttrib(ans_B_et, R_DimSymbol, Dim3);

                                    }

                                    PROTECT(dimnames3 = allocVector(VECSXP,count));
                                    if (dimC3[0]>0) {SET_VECTOR_ELT(dimnames3, count3, fleetList) ; count3++;}
                                    if (dimC3[1]>0) {SET_VECTOR_ELT(dimnames3, count3, metierList) ; count3++;}
                                    if (dimC3[2]>0) {SET_VECTOR_ELT(dimnames3, count3, intAge) ; count3++;}
                                    if (dimC3[3]>0) {SET_VECTOR_ELT(dimnames3, count3, times) ; count3++;}

                                    rans_B_et = REAL(ans_B_et);

                            } else {

                                    rans_B_et = REAL(VECTOR_ELT(out_B_et, e));

                            }

                                    //facteurs des indices
                                    PROTECT(dFACT7 = iDim(dimC3));

                                    int *fact7_D = INTEGER(dFACT7);

                                    //équation

                                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                            double temp = 0.0;
                                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                                temp = temp +
                                                 rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                 r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                            rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                        }

                            if(ind_t==0) {

                                    if (count>0) setAttrib(ans_B_et, R_DimNamesSymbol, dimnames3);
                                    setAttrib(ans_B_et, install("DimCst"), dimCst3);

                                    SET_VECTOR_ELT(out_B_et, e, ans_B_et);
                                    SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                                    UNPROTECT(9);
                            }

                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 20, v_w_ei);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 21, dFACT7);
                                SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 22, dFACT8);

                                UNPROTECT(34);


            }
            } else {

     //-------------------------------------------------------------------------------------------------  dimension trimestre


//Rprintf("G9");
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                            PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                            nbI = length(getListElement(elmt, "modI"));

                            PROTECT(v_M_ei = getListElement(getListElement(elmt, "M_i"),CHAR(STRING_ELT(trimInt,Qt))));
                            PROTECT(v_w_ei = getListElement(getListElement(elmt, "wStock_i"),CHAR(STRING_ELT(trimInt,0))));
                            PROTECT(v_mat_ei = getListElement(getListElement(elmt, "mat_i"),CHAR(STRING_ELT(trimInt,0))));
                            PROTECT(v_Fr_efmit = VECTOR_ELT(getListElement(out_Fr_fmi, CHAR(STRING_ELT(sppList,e))),Qt));
                            PROTECT(v_F_efmit = VECTOR_ELT(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))),Qt));
                            PROTECT(v_N_ei0 = getListElement(elmt, "N_it0"));
                            PROTECT(v_N_e0t = getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,Qt))));

                            PROTECT(v_Fbar = getListElement(elmt, "Fbar"));
                            PROTECT(dimCst_M_ei = getAttrib(v_M_ei, install("DimCst")));
                            PROTECT(dimCst_w_ei = getAttrib(v_w_ei, install("DimCst")));
                            PROTECT(dimCst_mat_ei = getAttrib(v_mat_ei, install("DimCst")));
                            PROTECT(dimCst_Fr_efmit = getAttrib(v_Fr_efmit, install("DimCst")));
                            PROTECT(dimCst_N_ei0 = getAttrib(v_N_ei0, install("DimCst")));
                            PROTECT(dimCst_N_e0t = getAttrib(v_N_e0t, install("DimCst")));

                            double NtotIniRec = REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,0))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,1))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,2))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,3))))[0];


                            //tests sur les dimensions :
                            dim_M_ei = INTEGER(dimCst_M_ei);
                            if ((dim_M_ei[0]!=0) | (dim_M_ei[1]!=0) |
                                (dim_M_ei[2]!=0 & dim_M_ei[2]!=nbI) | (dim_M_ei[3]!=0 & dim_M_ei[3]!=nbT)) //on laisse une ouverture pour un indice temporel
                            {
                                error("Non_homogeneous dimensions in M_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_w_ei = INTEGER(dimCst_w_ei);
                            if ((dim_w_ei[0]!=0) | (dim_w_ei[1]!=0) |
                                (dim_w_ei[2]!=0 & dim_w_ei[2]!=nbI) | (dim_w_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in w_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_mat_ei = INTEGER(dimCst_mat_ei);
                            if ((dim_mat_ei[0]!=0) | (dim_mat_ei[1]!=0) |
                                (dim_mat_ei[2]!=0 & dim_mat_ei[2]!=nbI) | (dim_mat_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in mat_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_Fr_efmit = INTEGER(dimCst_Fr_efmit);
                            if ((dim_Fr_efmit[0]!=0 & dim_Fr_efmit[0]!=nbF) | (dim_Fr_efmit[1]!=0 & dim_Fr_efmit[1]!=nbM) |
                                (dim_Fr_efmit[2]!=0 & dim_Fr_efmit[2]!=nbI) | (dim_Fr_efmit[3]!=0 & dim_Fr_efmit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in Fr_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_N_ei0 = INTEGER(dimCst_N_ei0);
                            if ((dim_N_ei0[0]!=0) | (dim_N_ei0[1]!=0) |
                                (dim_N_ei0[2]!=0 & dim_N_ei0[2]!=nbI)) // | (dim_N_ei0[3]!=0)) --> peu importe, on ne prendra de toute façon que la donnée à t0
                            {
                                error("Non_homogeneous dimensions in N_ei0 element. Check .ini biological parameters files !!\n");
                            }

                            dim_N_e0t = INTEGER(dimCst_N_e0t);
                            if ((dim_N_e0t[0]!=0) | (dim_N_e0t[1]!=0) |
                                (dim_N_e0t[2]!=0) | (dim_N_e0t[3]!=0 & dim_N_e0t[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in N_e0t element. Check .ini biological parameters files !!\n");
                            }

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // initialisation de out_Z_eit
                            //---------

                            //on détermine l'attribut Dimension de Z_eit
                            PROTECT(dimCst1 = allocVector(INTSXP, 4));
                            dimC1 = INTEGER(dimCst1);
                            dimC1[0] = 0 ; dimC1[1] = 0 ; dimC1[2] = imax2(dim_M_ei[2] , dim_Fr_efmit[2]);
                            dimC1[3] = imax2(dim_M_ei[3] , dim_Fr_efmit[3]);

                            int count = 0, prod = 1, count2 = 0, count3 = 0, count4 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC1[k]>0) {
                                    count++;
                                    prod = prod * dimC1[k];
                                }
                            }

                            PROTECT(Dim1 = allocVector(INTSXP, count));
                            dim1 = INTEGER(Dim1);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC1[k]>0) {
                                    dim1[count2] = dimC1[k];
                                    count2++;
                                    }
                            }

                    if (ind_t==0){
                            //on crée le tableau résultat pour l'espèce en question
                            PROTECT(ans_Z_eit = NEW_NUMERIC(prod));
                            setAttrib(ans_Z_eit, R_DimSymbol, Dim1);

                            PROTECT(dimnames1 = allocVector(VECSXP,count));
                            if (dimC1[0]>0) {SET_VECTOR_ELT(dimnames1, count3, fleetList) ; count3++;}
                            if (dimC1[1]>0) {SET_VECTOR_ELT(dimnames1, count3, metierList) ; count3++;}
                            if (dimC1[2]>0) {SET_VECTOR_ELT(dimnames1, count3, intAge) ; count3++;}
                            if (dimC1[3]>0) {SET_VECTOR_ELT(dimnames1, count3, times) ; count3++;}

                            rans_Z_eit = REAL(ans_Z_eit);

                            PROTECT(ans_F_itQ = NEW_NUMERIC(prod)); //mêmes dimensions que Z
                            setAttrib(ans_F_itQ, R_DimSymbol, Dim1);
                            r_F_itQ = REAL(ans_F_itQ);

                    } else {

                            rans_Z_eit = REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt));
                            r_F_itQ = REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),Qt));

                    }
//Rprintf("G10");
                            r_Fr_efmit = REAL(v_Fr_efmit);
                            r_F_efmit = REAL(v_F_efmit);
                            r_M_ei = REAL(v_M_ei);
                            r_Fbar = REAL(v_Fbar);

                            //facteurs des indices
                            PROTECT(dFACT1 = iDim(dimC1));
                            PROTECT(dFACT2 = iDim(dim_Fr_efmit));
                            PROTECT(dFACT3 = iDim(dim_M_ei));

                            int *fact1_D = INTEGER(dFACT1);
                            int *fact2_D = INTEGER(dFACT2);
                            int *fact3_D = INTEGER(dFACT3);

                            double *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60),Qt));

                            //double fmax = 0.0, sumWt = 0.0;

                            //équation
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0, tempCap = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                }

                                if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i[ind_i + nbI*ind_t];
                                }

                            //on en profite pour mettre à jour F_itQ

                            r_F_itQ[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = temp + r_Froth_i[ind_i + nbI*ind_t];

                            //fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                            //sumWt = sumWt + r_Fbar[ind_i];

                            }

//Rprintf("G11");
                    if (ind_t==0) {

                            setAttrib(ans_Z_eit, R_DimNamesSymbol, dimnames1);//Rprintf("G11_1");
                            setAttrib(ans_Z_eit, install("DimCst"), dimCst1);//Rprintf("G11_2");
                            setAttrib(ans_F_itQ, R_DimNamesSymbol, dimnames1);//Rprintf("G11_3");
                            setAttrib(ans_F_itQ, install("DimCst"), dimCst1);//Rprintf("G11_4");

                            SET_VECTOR_ELT(VECTOR_ELT(out_Z_eit, e), Qt, ans_Z_eit);
                            SET_VECTOR_ELT(VECTOR_ELT(out_F_itQ, e), Qt, ans_F_itQ);
                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 9, dimCst_Fr_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 10), Qt, v_Fr_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 11), Qt, v_M_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 12, dFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 13, dFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 14, dFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 59, v_Fbar);



                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////


                            //---------
                            // calcul de N_eit
                            //---------

                            //on détermine l'attribut Dimension de N_eit et des composantes trimestrielles N_eitQ
                            PROTECT(dimCst2 = allocVector(INTSXP, 4));
                            dimC2 = INTEGER(dimCst2);

                            dimC2[0] = 0 ; dimC2[1] = 0 ; dimC2[2] = imax2(dimC1[2] , dim_N_ei0[2]) ; dimC2[3] = nbT;

                            count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC2[k]>0) {
                                    count++;
                                    prod = prod * dimC2[k];
                                }
                            }

                            PROTECT(Dim2 = allocVector(INTSXP, count));
                            dim2 = INTEGER(Dim2);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC2[k]>0) {
                                    dim2[count2] = dimC2[k];
                                    count2++;
                                    }
                            }


                    if (ind_t==0) {

                            //on crée le tableau résultat pour l'espèce en question
                            if (Qt==0) {
                                PROTECT(ans_N_eit = NEW_NUMERIC(prod));
                                setAttrib(ans_N_eit, R_DimSymbol, Dim2);
                                rans_N_eit = REAL(ans_N_eit);
                            } else {
                                rans_N_eit = REAL(VECTOR_ELT(out_N_eit, e));
                            }


                            PROTECT(dimnames2 = allocVector(VECSXP,count));
                            if (dimC2[0]>0) {SET_VECTOR_ELT(dimnames2, count3, fleetList) ; count3++;}
                            if (dimC2[1]>0) {SET_VECTOR_ELT(dimnames2, count3, metierList) ; count3++;}
                            if (dimC2[2]>0) {SET_VECTOR_ELT(dimnames2, count3, intAge) ; count3++;}
                            if (dimC2[3]>0) {SET_VECTOR_ELT(dimnames2, count3, times) ; count3++;}


                            PROTECT(ans_N_eitQ = NEW_NUMERIC(prod));                    //+1
                            setAttrib(ans_N_eitQ, R_DimSymbol, Dim2);
                            r_N_eitQ = REAL(ans_N_eitQ);

                    } else {

                            rans_N_eit = REAL(VECTOR_ELT(out_N_eit, e));
                            r_N_eitQ = REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt));

                    }

                            r_N_ei0 = REAL(v_N_ei0);
                            r_N_e0t = REAL(v_N_e0t);

                            //facteurs des indices
                            PROTECT(dFACT4 = iDim(dimC2));
                            PROTECT(dFACT5 = iDim(dim_N_ei0));
                            PROTECT(dFACT6 = iDim(dim_N_e0t));

                            int *fact4_D = INTEGER(dFACT4);
                            int *fact5_D = INTEGER(dFACT5);
                            int *fact6_D = INTEGER(dFACT6);
//Rprintf("G12");
                            //équation

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ind_i == 0) { //recrutement

                                        if (SRInd[e]==1 & ind_t>0) {//Rprintf("G12a");

                                            if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                    REAL(VECTOR_ELT(out_SRmod,e))[ind_t] * r_N_e0t[0] / NtotIniRec;

                                            } else {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }


                                        } else {//Rprintf("G12b");

                                            if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial défini

                                            } else {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }
                                        }

                                    } else {//Rprintf("G12c");

                                        if (ind_t==0 & Qt==0) {

                                              r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                        } else {//Rprintf("G12d");

                                            if (Qt>0) { //trimestre 2, 3 ou 4

                                                if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + ind_t*fact1_D[3]]) +
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]);

                                                } else {

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + ind_t*fact1_D[3]]);

                                                }

                                            } else {  //trimestre 1

                                               if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                } else {

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                }

                                            }

                                        }
                                    }
                            }

                    //on transfère tout dans out_N_eit si Qt=0

                    if (Qt==0) {

                        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                        rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]];

                    }


                    if (ind_t==0) {

                            if (Qt==0) setAttrib(ans_N_eit, R_DimNamesSymbol, dimnames2);
                            if (Qt==0) setAttrib(ans_N_eit, install("DimCst"), dimCst2);
                            setAttrib(ans_N_eitQ, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_N_eitQ, install("DimCst"), dimCst2);

                            if (Qt==0) SET_VECTOR_ELT(out_N_eit, e, ans_N_eit);
                            SET_VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e), Qt, ans_N_eitQ);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 15, v_N_ei0);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 16, v_N_e0t);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 17, dFACT4);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 18, dFACT5);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 19, dFACT6);

//Rprintf("G13");
                //on calcule Fbar si Qt=3 (condition nbI%%4=0)   r_Fbar[ind_i]

                    double CNTnum = 0.0, CNTdenom = 0.0;

                    if (Qt==3) {

                        double Cnt0 = 0.0, Cnt1 = 0.0, Cnt2 = 0.0, Cnt3 = 0.0;
                        double CNT0 = 0.0, CNT1 = 0.0, CNT2 = 0.0, CNT3 = 0.0;

                        for (int IND_i = 0 ; IND_i < nbI ; IND_i++) {

                            if ((IND_i%4)==0) {Cnt1 = 0.0 ; Cnt2 = 0.0 ; Cnt3 = 0.0 ; Cnt0 = 0.0;
                                                CNT0 = 0.0, CNT1 = 0.0, CNT2 = 0.0, CNT3 = 0.0;}

                            Cnt0 = Cnt0 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt1 = Cnt1 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt2 = Cnt2 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt3 = Cnt3 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT0 = CNT0 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT1 = CNT1 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT2 = CNT2 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT3 = CNT3 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];

                            if ((IND_i%4)==3) {
                               CNTnum = CNTnum + r_Fbar[IND_i] * ((CNT0/Cnt0) + (CNT1/Cnt1) + (CNT2/Cnt2) + (CNT3/Cnt3));
                               CNTdenom = CNTdenom + r_Fbar[IND_i];
                            }

                        }

                    }



                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de SSB_et
                            //---------

                            //on détermine l'attribut Dimension de SSB_et
                            PROTECT(dimCst4 = allocVector(INTSXP, 4));
                            dimC4 = INTEGER(dimCst4);

                            dimC4[0] = 0 ; dimC4[1] = 0 ; dimC4[2] = 0 ; dimC4[3] = dimC2[3];

                            count = 0 ; prod = 1 ; count2 = 0 ; count4 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC4[k]>0) {
                                    count++;
                                    prod = prod * dimC4[k];
                                }
                            }


                            PROTECT(Dim4 = allocVector(INTSXP, count));
                            dim4 = INTEGER(Dim4);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC4[k]>0) {
                                    dim4[count2] = dimC4[k];
                                    count2++;
                                    }
                            }
//Rprintf("G14");
                    if (ind_t==0) {
                            //on crée le tableau résultat pour l'espèce en question (on en profite pour faire de même avec Fbar --> même dimension)
                            if (Qt==0) {
                                PROTECT(ans_SSB_et = NEW_NUMERIC(prod));
                                PROTECT(ans_Fbar_et = NEW_NUMERIC(prod));
                                rans_SSB_et = REAL(ans_SSB_et);
                                rans_Fbar_et = REAL(ans_Fbar_et);
                            } else {
                                rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et, e));
                                rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et, e));
                            }

                            PROTECT(ans_SSB_etQ = NEW_NUMERIC(prod));                    //+1

                            if (count>0) { //valable seulement si SSB_et n'est pas seulement un scalaire

                            if (Qt==0) setAttrib(ans_SSB_et, R_DimSymbol, Dim4);
                            if (Qt==0) setAttrib(ans_Fbar_et, R_DimSymbol, Dim4);
                            setAttrib(ans_SSB_etQ, R_DimSymbol, Dim4);

                            }

                            PROTECT(dimnames4 = allocVector(VECSXP,count));
                            if (dimC4[0]>0) {SET_VECTOR_ELT(dimnames4, count4, fleetList) ; count4++;}
                            if (dimC4[1]>0) {SET_VECTOR_ELT(dimnames4, count4, metierList) ; count4++;}
                            if (dimC4[2]>0) {SET_VECTOR_ELT(dimnames4, count4, intAge) ; count4++;}
                            if (dimC4[3]>0) {SET_VECTOR_ELT(dimnames4, count4, times) ; count4++;}


                            r_SSB_etQ = REAL(ans_SSB_etQ);

                    } else {


                            rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et, e));
                            rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et, e));
                            r_SSB_etQ = REAL(VECTOR_ELT(VECTOR_ELT(out_SSB_etQ, e),Qt));

                    }

                            r_w_ei = REAL(v_w_ei);
                            r_mat_ei = REAL(v_mat_ei);

                            //facteurs des indices
                            PROTECT(dFACT9 = iDim(dimC4));
                            PROTECT(dFACT8 = iDim(dim_w_ei));
                            PROTECT(dFACT10 = iDim(dim_mat_ei));

                            int *fact9_D = INTEGER(dFACT9);
                            int *fact8_D = INTEGER(dFACT8);
                            int *fact10_D = INTEGER(dFACT10);

                            //équation

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                    double temp = 0.0;
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                        temp = temp +
                                         r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    r_SSB_etQ[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;

                                    if (Qt==0) rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;

                                    if (Qt==3) rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = CNTnum/(4*CNTdenom);

                                }

                    if(ind_t==0) {

                            if (count>0) {
                                if (Qt==0) setAttrib(ans_SSB_et, R_DimNamesSymbol, dimnames4);
                                setAttrib(ans_SSB_etQ, R_DimNamesSymbol, dimnames4);
                            }
                            if (Qt==0) setAttrib(ans_SSB_et, install("DimCst"), dimCst4);
                            setAttrib(ans_SSB_etQ, install("DimCst"), dimCst4);

                            if (Qt==0) SET_VECTOR_ELT(out_SSB_et, e, ans_SSB_et);
                            SET_VECTOR_ELT(VECTOR_ELT(out_SSB_etQ, e), Qt, ans_SSB_etQ);

                            if (count>0 & Qt==0) setAttrib(ans_Fbar_et, R_DimNamesSymbol, dimnames4);
                            if (Qt==0) setAttrib(ans_Fbar_et, install("DimCst"), dimCst4);

                            if (Qt==0) SET_VECTOR_ELT(out_Fbar_et, e, ans_Fbar_et);

                            if (Qt==0) SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 48, dFACT9);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 49, dFACT10);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 47, v_mat_ei);


//Rprintf("G15");

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de B_et
                            //---------

                            //on détermine l'attribut Dimension de B_et
                            PROTECT(dimCst3 = allocVector(INTSXP, 4));
                            dimC3 = INTEGER(dimCst3);

                            dimC3[0] = 0 ; dimC3[1] = 0 ; dimC3[2] = 0 ; dimC3[3] = dimC2[3];

                            count = 0 ; prod = 1 ; count2 = 0 ; count3 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC3[k]>0) {
                                    count++;
                                    prod = prod * dimC3[k];
                                }
                            }

                            PROTECT(Dim3 = allocVector(INTSXP, count));
                            dim3 = INTEGER(Dim3);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC3[k]>0) {
                                    dim3[count2] = dimC3[k];
                                    count2++;
                                    }
                            }

                    if (ind_t==0 & Qt==0) {
                            //on crée le tableau résultat pour l'espèce en question
                            PROTECT(ans_B_et = NEW_NUMERIC(prod));

                            if (count>0) { //valable seulement si B_et n'est pas seulement un scalaire

                            setAttrib(ans_B_et, R_DimSymbol, Dim3);

                            }

                            PROTECT(dimnames3 = allocVector(VECSXP,count));
                            if (dimC3[0]>0) {SET_VECTOR_ELT(dimnames3, count3, fleetList) ; count3++;}
                            if (dimC3[1]>0) {SET_VECTOR_ELT(dimnames3, count3, metierList) ; count3++;}
                            if (dimC3[2]>0) {SET_VECTOR_ELT(dimnames3, count3, intAge) ; count3++;}
                            if (dimC3[3]>0) {SET_VECTOR_ELT(dimnames3, count3, times) ; count3++;}

                            rans_B_et = REAL(ans_B_et);

                    } else {

                            rans_B_et = REAL(VECTOR_ELT(out_B_et, e));

                    }

                            //facteurs des indices
                            PROTECT(dFACT7 = iDim(dimC3));

                            int *fact7_D = INTEGER(dFACT7);

                            //équation

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                    double temp = 0.0;
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    if (Qt==0) rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                }

                    if(ind_t==0 & Qt==0) {

                            if (count>0) setAttrib(ans_B_et, R_DimNamesSymbol, dimnames3);
                            setAttrib(ans_B_et, install("DimCst"), dimCst3);

                            SET_VECTOR_ELT(out_B_et, e, ans_B_et);
                            if (Qt==0) SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                            UNPROTECT(5);
                    }

                    if (ind_t==0) UNPROTECT(7);

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 20, v_w_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 21, dFACT7);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 22, dFACT8);

                        UNPROTECT(34);


            }

}

if (Qt==3) dUpdate = false;

if (ind_t==0) UNPROTECT(1);


} else {


for (int e = 0 ; e < nbE ; e++) {

    if (trim[e]==0) {
//Rprintf("G16");
        if (Qt==0) {

                    SEXP elmt;
                    PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                    int nbI = length(VECTOR_ELT(namDC, e));
                    SEXP v_N_e0t;

                    PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));

                    double  *rans_Z_eit = REAL(VECTOR_ELT(out_Z_eit,e));
                    double  *rans_N_eit = REAL(VECTOR_ELT(out_N_eit,e));
                    double  *rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et,e));
                    double  *rans_B_et = REAL(VECTOR_ELT(out_B_et,e));
                    double  *rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et,e));
                    double  *r_Fr_efmit = REAL(VECTOR_ELT(out_Fr_fmi, e));
                    double  *r_F_efmit = REAL(VECTOR_ELT(out_F_fmi, e));
                    double  *r_M_ei = REAL(getListElement(elmt, "M_i"));
                    double  *r_N_ei0 = REAL(getListElement(elmt, "N_it0"));
                    double  *r_N_e0t = REAL(v_N_e0t);
                    double  *r_w_ei = REAL(getListElement(elmt, "wStock_i"));
                    double  *r_mat_ei = REAL(getListElement(elmt, "mat_i"));
                    double  *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60));
                    double  *r_Fbar = REAL(getListElement(elmt, "Fbar"));

                    int *dim_Fr_efmit = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 9)),
                        *fact1_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 12)),
                        *fact2_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 13)),
                        *fact3_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 14)),
                        *fact4_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 17)),
                        *fact5_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 18)),
                        *fact6_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 19)),
                        *fact7_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 21)),
                        *fact8_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 22)),
                        *fact9_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 48)),
                        *fact10_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 49));

                    double fmax = 0.0, sumWt = 0.0;

                    //équation n°1 : out_Z_eit

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0, tempCap = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                }

                            if (Zoptim_use & e==eTemp) {

                               rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    Zoptim[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]];


                            } else {

                                if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i[ind_i + ind_t*nbI];
                                }
                            }
                            //on en profite pour calculer Fbar

                            fmax = fmax + (temp + r_Froth_i[ind_i + ind_t*nbI])*r_Fbar[ind_i];
                            sumWt = sumWt + r_Fbar[ind_i];

                            }
//Rprintf("G17");

                    //équation n°2 : out_N_eit

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ind_i == 0) {

                                        if (SRInd[e]==1 & ind_t>0) {

                                            if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                    REAL(VECTOR_ELT(out_SRmod,e))[ind_t];

                                            } else {

                                             if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }}




                                        } else {

                                            if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }
                                        }

                                    } else {

                                        if (ind_t == 0) {

                                            rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                        } else {

                                            if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            } else {

                                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                  exp(-rans_Z_eit[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                            }
                                        }
                                    }
                            }


                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                    rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = fmax/sumWt;
                                }

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                    double temp = 0.0;

                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                                }

                    UNPROTECT(2);
    }
    } else {

//Rprintf("G18");
                SEXP elmt;
                PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                int nbI = length(VECTOR_ELT(namDC, e));
                SEXP v_N_e0t;

                PROTECT(v_N_e0t = getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,Qt))));

                double  *rans_Z_eit = REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit,e),Qt));
                double  *rans_N_eit = REAL(VECTOR_ELT(out_N_eit,e));
                double  *r_N_eitQ = REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ,e),Qt));
                double  *r_F_itQ = REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),Qt));
                double  *r_SSB_etQ = REAL(VECTOR_ELT(VECTOR_ELT(out_SSB_etQ, e),Qt));
                double  *rans_Fbar_et = REAL(VECTOR_ELT(out_Fbar_et,e));
                double  *rans_B_et = REAL(VECTOR_ELT(out_B_et,e));
                double  *rans_SSB_et = REAL(VECTOR_ELT(out_SSB_et,e));
                double  *r_Fr_efmit = REAL(VECTOR_ELT(getListElement(out_Fr_fmi, CHAR(STRING_ELT(sppList,e))),Qt));
                double  *r_F_efmit = REAL(VECTOR_ELT(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))),Qt));
                double  *r_M_ei = REAL(getListElement(getListElement(elmt, "M_i"),CHAR(STRING_ELT(trimInt,Qt))));
                double  *r_N_ei0 = REAL(getListElement(elmt, "N_it0"));
                double  *r_N_e0t = REAL(v_N_e0t);
                double  *r_w_ei = REAL(getListElement(getListElement(elmt, "wStock_i"),CHAR(STRING_ELT(trimInt,0))));
                double  *r_mat_ei = REAL(getListElement(getListElement(elmt, "mat_i"),CHAR(STRING_ELT(trimInt,0))));
                double  *r_Froth_i = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 60),Qt));
                double  *r_Fbar = REAL(getListElement(elmt, "Fbar"));

                double NtotIniRec = REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,0))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,1))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,2))))[0] +
                                REAL(getListElement(getListElement(elmt, "N_i0t"),CHAR(STRING_ELT(trimInt,3))))[0];

                int *dim_Fr_efmit = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 9)),
                    *fact1_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 12)),
                    *fact2_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 13)),
                    *fact3_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 14)),
                    *fact4_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 17)),
                    *fact5_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 18)),
                    *fact6_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 19)),
                    *fact7_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 21)),
                    *fact8_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 22)),
                    *fact9_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 48)),
                    *fact10_D = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 49));

                //double fmax = 0.0, sumWt = 0.0;

                //équation n°1 : out_Z_eit

                        for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0, tempCap = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF - 1)*(dim_Fr_efmit[0]>0)) ; ind_f++) //attention aux sommations multiples!!!
                                for (int ind_m = 0 ; ind_m < (1 + (nbM - 1)*(dim_Fr_efmit[1]>0)) ; ind_m++) {

                                if (!ISNA(r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    temp = temp +  r_Fr_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];
                                if (!ISNA(r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]]))
                                    tempCap = tempCap +  r_F_efmit[ind_f*fact2_D[0] + ind_m*fact2_D[1] + ind_i*fact2_D[2] + ind_t*fact2_D[3]];

                                }

                                if (ISNA(r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]])) {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = NA_REAL;
                                } else {
                                  rans_Z_eit[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] =
                                    temp + r_M_ei[0*fact3_D[0] + 0*fact3_D[1] + ind_i*fact3_D[2] + ind_t*fact3_D[3]] +
                                    r_Froth_i[ind_i + nbI*ind_t];
                                }

                            //on en profite pour mettre à jour F_itQ

                            r_F_itQ[0*fact1_D[0] + 0*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]] = temp + r_Froth_i[ind_i + nbI*ind_t];

                            //fmax = fmax + (temp + r_Froth_i[ind_i + nbI*ind_t])*r_Fbar[ind_i];
                            //sumWt = sumWt + r_Fbar[ind_i];

                            }

//Rprintf("G19");
                //équation n°2 : out_N_eit

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ind_i == 0) { //recrutement

                                        if (SRInd[e]==1 & ind_t>0) {

                                            if (!ISNA(REAL(VECTOR_ELT(out_SRmod,e))[ind_t])) {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                    REAL(VECTOR_ELT(out_SRmod,e))[ind_t] * r_N_e0t[0] / NtotIniRec;;

                                            } else {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }


                                        } else {

                                            if (ISNA(r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]])) {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + 0*fact5_D[3]]; //seul instant initial défini

                                            } else {

                                                r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                  r_N_e0t[ind_f*fact6_D[0] + ind_m*fact6_D[1] + ind_i*fact6_D[2] + ind_t*fact6_D[3]];

                                            }
                                        }

                                    } else {

                                        if (ind_t == 0) {

                                            r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                              r_N_ei0[ind_f*fact5_D[0] + ind_m*fact5_D[1] + ind_i*fact5_D[2] + ind_t*fact5_D[3]];

                                        } else {

                                            if (Qt>0) { //trimestre 2, 3 ou 4

                                                if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + ind_t*fact1_D[3]]) +
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + ind_t*fact1_D[3]]);

                                                } else {

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt-1))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + ind_t*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt-1))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + ind_t*fact1_D[3]]);

                                                }

                                            } else {  //trimestre 1

                                               if (ind_i == (nbI-1)) {  //groupe d'âge +

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]) +
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + ind_i*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                } else {

                                                    r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                                      REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[ind_f*fact4_D[0] + ind_m*fact4_D[1] + (ind_i-1)*fact4_D[2] + (ind_t-1)*fact4_D[3]] *
                                                      exp(-REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),3))[ind_f*fact1_D[0] + ind_m*fact1_D[1] + (ind_i-1)*fact1_D[2] + (ind_t-1)*fact1_D[3]]);

                                                }

                                            }

                                        }
                                    }
                            }


                            if (Qt==0) {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] =
                                        r_N_eitQ[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]];

                            }




               //on calcule Fbar si Qt=3 (condition nbI%%4=0)   r_Fbar[ind_i]
//Rprintf("G20");
                    double CNTnum = 0.0, CNTdenom = 0.0;

                    if (Qt==3) {

                        double Cnt0 = 0.0, Cnt1 = 0.0, Cnt2 = 0.0, Cnt3 = 0.0;
                        double CNT0 = 0.0, CNT1 = 0.0, CNT2 = 0.0, CNT3 = 0.0;

                        for (int IND_i = 0 ; IND_i < nbI ; IND_i++) {

                            if ((IND_i%4)==0) {
                                Cnt1 = 0.0 ; Cnt2 = 0.0 ; Cnt3 = 0.0 ; Cnt0 = 0.0;
                                CNT0 = 0.0, CNT1 = 0.0, CNT2 = 0.0, CNT3 = 0.0;}

                            Cnt0 = Cnt0 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt1 = Cnt1 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt2 = Cnt2 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            Cnt3 = Cnt3 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT0 = CNT0 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),0))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT1 = CNT1 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),1))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT2 = CNT2 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),2))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];
                            CNT3 = CNT3 + REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]]*
                                REAL(VECTOR_ELT(VECTOR_ELT(out_F_itQ, e),3))[IND_i*fact4_D[2] + ind_t*fact4_D[3]];

                            if ((IND_i%4)==3) {
                               CNTnum = CNTnum + r_Fbar[IND_i] * ((CNT0/Cnt0) + (CNT1/Cnt1) + (CNT2/Cnt2) + (CNT3/Cnt3));
                               CNTdenom = CNTdenom + r_Fbar[IND_i];
                            }

                        }

                    }



                // SSB et finalisation Fbar

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                    double temp = 0.0;
                                    for (int ind_i = 0 ; ind_i < nbI ; ind_i++) //on fait ici l'hypothèse que la dimension âge est toujours présente
                                        temp = temp +
                                         rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                         r_mat_ei[ind_f*fact10_D[0] + ind_m*fact10_D[1] + ind_i*fact10_D[2] + ind_t*fact10_D[3]] *
                                         r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                    r_SSB_etQ[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;

                                    rans_SSB_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = temp;
                                    if (Qt==3) rans_Fbar_et[ind_f*fact9_D[0] + ind_m*fact9_D[1] + 0*fact9_D[2] + ind_t*fact9_D[3]] = CNTnum/(4*CNTdenom);

                            }

                // B
                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++){

                                double temp = 0.0;

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                    temp = temp +
                                     rans_N_eit[ind_f*fact4_D[0] + ind_m*fact4_D[1] + ind_i*fact4_D[2] + ind_t*fact4_D[3]] *
                                     r_w_ei[ind_f*fact8_D[0] + ind_m*fact8_D[1] + ind_i*fact8_D[2] + ind_t*fact8_D[3]] / 1000;

                                rans_B_et[ind_f*fact7_D[0] + ind_m*fact7_D[1] + 0*fact7_D[2] + ind_t*fact7_D[3]] = temp;

                            }

                UNPROTECT(2);


    }



}}

}}












//------------------------------------------
// Module 'Captures, rejets et débarquements'
//------------------------------------------

extern "C" {

void BioEcoPar::CatchDL(SEXP list, int ind_t, SEXP EVAR, int Qt)
{

if (cUpdate) {


    SEXP    elmt, dimCst, Dim, dimCst_F_efmit, dimCst_N_eit, dimCst_Z_eit, dimCst_wL_ei, dimCst_wD_ei, dimCst_d_efmit,
            intAge, v_F_efmit, v_N_eit, v_Z_eit, v_wL_ei, v_wD_ei, v_d_efmit, dimCst2, Dim2,
            cFACT1, cFACT2, cFACT3, cFACT4, cFACT5, cFACT6, cFACT7;

    SEXP ans_C_efmit=R_NilValue, ans_Y_efmit=R_NilValue, ans_D_efmit=R_NilValue, ans_L_efmit=R_NilValue,
         dimnames=R_NilValue, rnames_Esp=R_NilValue, ans_C_eit=R_NilValue, ans_Y_eit=R_NilValue, dimnames2=R_NilValue;

    int *dim_F_efmit, *dim_N_eit, *dim_Z_eit, *dim_wL_ei, *dim_wD_ei, *dim_d_efmit, *dimC, *dim, *dim2, *dimcst2;
    int nbI;

    double *rans_C_efmit, *rans_Y_efmit, *rans_D_efmit, *rans_L_efmit, *r_F_efmit, *r_N_eit, *r_Z_eit, *r_wL_ei, *r_wD_ei, *r_d_efmit,
            *rans_C_eit, *rans_Y_eit;

if (ind_t==0) {

    PROTECT(rnames_Esp = allocVector(STRSXP, nbE));
    if (Qt==0) setAttrib(out_C_efmit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_C_eit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_Y_eit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_Y_efmit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_D_efmit, R_NamesSymbol, rnames_Esp);

    if (Qt==0) setAttrib(out_L_efmit, R_NamesSymbol, rnames_Esp);
}


    for (int e = 0 ; e < nbE ; e++) {

         if (trim[e]==0) { //-------------------------------------------------------------- pas de dimension trimestre

            if (Qt==0) {
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                            PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                            nbI = length(getListElement(elmt, "modI"));

                            PROTECT(v_wL_ei = getListElement(elmt, "wL_i"));
                            PROTECT(v_wD_ei = getListElement(elmt, "wD_i"));
                            PROTECT(v_d_efmit = getListElement(elmt, "d_i"));
                            PROTECT(v_F_efmit = getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))));
                            PROTECT(v_Z_eit = getListElement(out_Z_eit , CHAR(STRING_ELT(sppList,e))));
                            PROTECT(v_N_eit = getListElement(out_N_eit , CHAR(STRING_ELT(sppList,e))));

                            PROTECT(dimCst_wL_ei = getAttrib(v_wL_ei, install("DimCst")));
                            PROTECT(dimCst_wD_ei = getAttrib(v_wD_ei, install("DimCst")));
                            PROTECT(dimCst_d_efmit = getAttrib(v_d_efmit, install("DimCst")));
                            PROTECT(dimCst_F_efmit = getAttrib(v_F_efmit, install("DimCst")));
                            PROTECT(dimCst_N_eit = getAttrib(v_N_eit, install("DimCst")));
                            PROTECT(dimCst_Z_eit = getAttrib(v_Z_eit, install("DimCst")));

                            //tests sur les dimensions :
                            dim_d_efmit = INTEGER(dimCst_d_efmit);
                            if ((dim_d_efmit[0]!=0 & dim_d_efmit[0]!=nbF) | (dim_d_efmit[1]!=0 & dim_d_efmit[1]!=nbM) |
                                (dim_d_efmit[2]!=0 & dim_d_efmit[2]!=nbI) | (dim_d_efmit[3]!=0 & dim_d_efmit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in d_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_wL_ei = INTEGER(dimCst_wL_ei);
                            if ((dim_wL_ei[0]!=0) | (dim_wL_ei[1]!=0) |
                                (dim_wL_ei[2]!=0 & dim_wL_ei[2]!=nbI) | (dim_wL_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wL_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_wD_ei = INTEGER(dimCst_wD_ei);
                            if ((dim_wD_ei[0]!=0) | (dim_wD_ei[1]!=0) |
                                (dim_wD_ei[2]!=0 & dim_wD_ei[2]!=nbI) | (dim_wD_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wD_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_F_efmit = INTEGER(dimCst_F_efmit);
                            if ((dim_F_efmit[0]!=0 & dim_F_efmit[0]!=nbF) | (dim_F_efmit[1]!=0 & dim_F_efmit[1]!=nbM) |
                                (dim_F_efmit[2]!=0 & dim_F_efmit[2]!=nbI) | (dim_F_efmit[3]!=0 & dim_F_efmit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in F_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_N_eit = INTEGER(dimCst_N_eit);
                            if ((dim_N_eit[0]!=0) | (dim_N_eit[1]!=0) |
                                (dim_N_eit[2]!=0 & dim_N_eit[2]!=nbI) | (dim_N_eit[3]!=0 & dim_N_eit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in N_eit element. Check .ini biological parameters files !!\n");
                            }

                            dim_Z_eit = INTEGER(dimCst_Z_eit);
                            if ((dim_Z_eit[0]!=0) | (dim_Z_eit[1]!=0) |
                                (dim_Z_eit[2]!=0 & dim_Z_eit[2]!=nbI) | (dim_Z_eit[3]!=0 & dim_Z_eit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in Z_eit element. Check .ini biological parameters files !!\n");
                            }

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de C_efmit
                            //---------

                            //on détermine l'attribut Dimension de C_efmit
                            PROTECT(dimCst = allocVector(INTSXP, 4));
                            dimC = INTEGER(dimCst);
                            dimC[0] = dim_F_efmit[0] ; dimC[1] = dim_F_efmit[1] ; dimC[2] = dim_F_efmit[2];
                            dimC[3] = imax2(dim_N_eit[3] , dim_F_efmit[3]);

                            int count = 0, prod = 1, count2 = 0, count3 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    count++;
                                    prod = prod * dimC[k];
                                }

                            }

                            PROTECT(Dim = allocVector(INTSXP, count));
                            dim = INTEGER(Dim);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    dim[count2] = dimC[k];
                                    count2++;
                                    }
                            }

                            PROTECT(dimCst2 = allocVector(INTSXP, 4));
                            dimcst2 = INTEGER(dimCst2); dimcst2[0] = 0; dimcst2[1] = 0; dimcst2[2] = nbI; dimcst2[3] = nbT;
                            PROTECT(Dim2 = allocVector(INTSXP, 2));
                            dim2 = INTEGER(Dim2); dim2[0] = nbI; dim2[1] = nbT;


                    if (ind_t==0){

                            //on crée le tableau résultat pour l'espèce en question
                            PROTECT(ans_C_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_C_efmit, R_DimSymbol, Dim);

                            PROTECT(ans_C_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_C_eit, R_DimSymbol, Dim2);

                            PROTECT(dimnames = allocVector(VECSXP,count));
                            if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                            if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                            if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                            if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}

                            PROTECT(dimnames2 = allocVector(VECSXP,2));
                            SET_VECTOR_ELT(dimnames2, 0, intAge);
                            SET_VECTOR_ELT(dimnames2, 1, times);

                            rans_C_efmit = REAL(ans_C_efmit);
                            rans_C_eit = REAL(ans_C_eit);

                    } else {

                            rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit, e));
                            rans_C_eit = REAL(VECTOR_ELT(out_C_eit, e));

                    }


                            r_F_efmit = REAL(v_F_efmit);
                            r_N_eit = REAL(v_N_eit);
                            r_Z_eit = REAL(v_Z_eit);

                            //facteurs des indices
                            PROTECT(cFACT1 = iDim(dimC));
                            PROTECT(cFACT2 = iDim(dim_F_efmit));
                            PROTECT(cFACT3 = iDim(dim_N_eit));
                            PROTECT(cFACT4 = iDim(dim_Z_eit));

                            int *fact1_C = INTEGER(cFACT1);
                            int *fact2_C = INTEGER(cFACT2);
                            int *fact3_C = INTEGER(cFACT3);
                            int *fact4_C = INTEGER(cFACT4);

                            double *r_Fot_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                  rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                    if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                    temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                                rans_C_eit[ind_i + ind_t*nbI] =
                                    (temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                        }

                    if (ind_t==0){

                            setAttrib(ans_C_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_C_efmit, install("DimCst"), dimCst);

                            setAttrib(ans_C_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_C_eit, install("DimCst"), dimCst2);


                            SET_VECTOR_ELT(out_C_efmit, e, ans_C_efmit);
                            SET_VECTOR_ELT(out_C_eit, e, ans_C_eit);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 23, v_F_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 24, v_N_eit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 25, v_Z_eit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 26, cFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 27, cFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 28, cFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 29, cFACT4);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de Y_efmit
                            //---------

                        //on considère les dimensions de C, Y, D et L homogènes sur tout le module --> pas besoin de les redéfinir

                    if (ind_t==0){

                            PROTECT(ans_Y_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_Y_efmit, R_DimSymbol, Dim);

                            rans_Y_efmit = REAL(ans_Y_efmit);

                            PROTECT(ans_Y_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_Y_eit, R_DimSymbol, Dim2);

                            rans_Y_eit = REAL(ans_Y_eit);


                    } else {

                            rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e));
                            rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e));
                    }


                            r_wL_ei = REAL(v_wL_ei);

                            //facteurs des indices
                            PROTECT(cFACT5 = iDim(dim_wL_ei));

                            int *fact5_C = INTEGER(cFACT5);

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                  rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;


                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                  rans_Y_eit[ind_i + ind_t*nbI] =
                                    r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    rans_C_eit[ind_i + ind_t*nbI] / 1000;

                    if (ind_t==0) {

                            setAttrib(ans_Y_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_Y_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_Y_efmit, e, ans_Y_efmit);

                            setAttrib(ans_Y_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_Y_eit, install("DimCst"), dimCst2);

                            SET_VECTOR_ELT(out_Y_eit, e, ans_Y_eit);


                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 30, v_wL_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 31, cFACT5);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de D_efmit
                            //---------

                    if (ind_t==0) {

                            PROTECT(ans_D_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_D_efmit, R_DimSymbol, Dim);

                            rans_D_efmit = REAL(ans_D_efmit);

                    } else {

                            rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e));

                    }

                            r_wD_ei = REAL(v_wD_ei);
                            r_d_efmit = REAL(v_d_efmit);

                            //facteurs des indices
                            PROTECT(cFACT6 = iDim(dim_d_efmit));
                            PROTECT(cFACT7 = iDim(dim_wD_ei));

                            int *fact6_C = INTEGER(cFACT6);
                            int *fact7_C = INTEGER(cFACT7);


                            //équation : 2 manières de calculer selon la disponibilité de wD_i

                    if (all_is_na(v_wD_ei)) { //1ère méthode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                            }

                    } else {                 //2ème méthode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;

                            }

                    }


                    if (ind_t==0) {

                            setAttrib(ans_D_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_D_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_D_efmit, e, ans_D_efmit);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 32, v_wD_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 33, v_d_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 34, cFACT6);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 35, cFACT7);

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de L_efmit
                            //---------

                    if (ind_t==0) {

                            PROTECT(ans_L_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_L_efmit, R_DimSymbol, Dim);

                            rans_L_efmit = REAL(ans_L_efmit);

                    } else {

                            rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e));

                    }

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                  rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                    if (ind_t==0) {

                            setAttrib(ans_L_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_L_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_L_efmit, e, ans_L_efmit);
                            SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                            UNPROTECT(8);
                    }


                        UNPROTECT(25);
         }
         } else {    // ----------------------------------------------------------------------------- dimension trimestre


                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                            PROTECT(intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e))));

                            nbI = length(getListElement(elmt, "modI"));

                            PROTECT(v_wL_ei = getListElement(getListElement(elmt, "wL_i"),CHAR(STRING_ELT(trimInt,Qt))));
                            PROTECT(v_wD_ei = getListElement(getListElement(elmt, "wD_i"),CHAR(STRING_ELT(trimInt,Qt))));
                            PROTECT(v_d_efmit = getListElement(getListElement(elmt, "d_i"),CHAR(STRING_ELT(trimInt,Qt))));
                            PROTECT(v_F_efmit = VECTOR_ELT(getListElement(out_F_fmi, CHAR(STRING_ELT(sppList,e))),Qt));
                            PROTECT(v_Z_eit = VECTOR_ELT(getListElement(out_Z_eit , CHAR(STRING_ELT(sppList,e))),Qt));
                            PROTECT(v_N_eit = VECTOR_ELT(getListElement(out_N_eitQ , CHAR(STRING_ELT(sppList,e))),Qt));

                            PROTECT(dimCst_wL_ei = getAttrib(v_wL_ei, install("DimCst")));
                            PROTECT(dimCst_wD_ei = getAttrib(v_wD_ei, install("DimCst")));
                            PROTECT(dimCst_d_efmit = getAttrib(v_d_efmit, install("DimCst")));
                            PROTECT(dimCst_F_efmit = getAttrib(v_F_efmit, install("DimCst")));
                            PROTECT(dimCst_N_eit = getAttrib(v_N_eit, install("DimCst")));
                            PROTECT(dimCst_Z_eit = getAttrib(v_Z_eit, install("DimCst")));

                            //tests sur les dimensions :
                            dim_d_efmit = INTEGER(dimCst_d_efmit);
                            if ((dim_d_efmit[0]!=0 & dim_d_efmit[0]!=nbF) | (dim_d_efmit[1]!=0 & dim_d_efmit[1]!=nbM) |
                                (dim_d_efmit[2]!=0 & dim_d_efmit[2]!=nbI) | (dim_d_efmit[3]!=0 & dim_d_efmit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in d_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_wL_ei = INTEGER(dimCst_wL_ei);
                            if ((dim_wL_ei[0]!=0) | (dim_wL_ei[1]!=0) |
                                (dim_wL_ei[2]!=0 & dim_wL_ei[2]!=nbI) | (dim_wL_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wL_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_wD_ei = INTEGER(dimCst_wD_ei);
                            if ((dim_wD_ei[0]!=0) | (dim_wD_ei[1]!=0) |
                                (dim_wD_ei[2]!=0 & dim_wD_ei[2]!=nbI) | (dim_wD_ei[3]!=0))
                            {
                                error("Non_homogeneous dimensions in wD_ei element. Check .ini biological parameters files !!\n");
                            }

                            dim_F_efmit = INTEGER(dimCst_F_efmit);
                            if ((dim_F_efmit[0]!=0 & dim_F_efmit[0]!=nbF) | (dim_F_efmit[1]!=0 & dim_F_efmit[1]!=nbM) |
                                (dim_F_efmit[2]!=0 & dim_F_efmit[2]!=nbI) | (dim_F_efmit[3]!=0 & dim_F_efmit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in F_efmit element. Check .ini biological parameters files !!\n");
                            }

                            dim_N_eit = INTEGER(dimCst_N_eit);
                            if ((dim_N_eit[0]!=0) | (dim_N_eit[1]!=0) |
                                (dim_N_eit[2]!=0 & dim_N_eit[2]!=nbI) | (dim_N_eit[3]!=0 & dim_N_eit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in N_eit element. Check .ini biological parameters files !!\n");
                            }

                            dim_Z_eit = INTEGER(dimCst_Z_eit);
                            if ((dim_Z_eit[0]!=0) | (dim_Z_eit[1]!=0) |
                                (dim_Z_eit[2]!=0 & dim_Z_eit[2]!=nbI) | (dim_Z_eit[3]!=0 & dim_Z_eit[3]!=nbT))
                            {
                                error("Non_homogeneous dimensions in Z_eit element. Check .ini biological parameters files !!\n");
                            }

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                            //---------
                            // calcul de C_efmit
                            //---------

                            //on détermine l'attribut Dimension de C_efmit
                            PROTECT(dimCst = allocVector(INTSXP, 4));
                            dimC = INTEGER(dimCst);
                            dimC[0] = dim_F_efmit[0] ; dimC[1] = dim_F_efmit[1] ; dimC[2] = dim_F_efmit[2];
                            dimC[3] = imax2(dim_N_eit[3] , dim_F_efmit[3]);

                            int count = 0, prod = 1, count2 = 0, count3 = 0;

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    count++;
                                    prod = prod * dimC[k];
                                }

                            }

                            PROTECT(Dim = allocVector(INTSXP, count));
                            dim = INTEGER(Dim);

                            for (int k = 0 ; k < 4 ; k++) {

                                if (dimC[k]>0) {
                                    dim[count2] = dimC[k];
                                    count2++;
                                    }
                            }

                            PROTECT(dimCst2 = allocVector(INTSXP, 4));
                            dimcst2 = INTEGER(dimCst2); dimcst2[0] = 0; dimcst2[1] = 0; dimcst2[2] = nbI; dimcst2[3] = nbT;
                            PROTECT(Dim2 = allocVector(INTSXP, 2));
                            dim2 = INTEGER(Dim2); dim2[0] = nbI; dim2[1] = nbT;


                    if (ind_t==0 & Qt==0){   //(ind_t==0 & Qt==0)

                            //on crée le tableau résultat pour l'espèce en question
                            PROTECT(ans_C_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_C_efmit, R_DimSymbol, Dim);

                            PROTECT(ans_C_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_C_eit, R_DimSymbol, Dim2);

                            PROTECT(dimnames = allocVector(VECSXP,count));
                            if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, fleetList) ; count3++;}
                            if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, metierList) ; count3++;}
                            if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
                            if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, times) ; count3++;}

                            PROTECT(dimnames2 = allocVector(VECSXP,2));
                            SET_VECTOR_ELT(dimnames2, 0, intAge);
                            SET_VECTOR_ELT(dimnames2, 1, times);

                            rans_C_efmit = REAL(ans_C_efmit);//Rprintf("JJ7");
                            rans_C_eit = REAL(ans_C_eit);

                    } else {

                            rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit, e));//Rprintf("JJ8");
                            rans_C_eit = REAL(VECTOR_ELT(out_C_eit, e));

                    }


                            r_F_efmit = REAL(v_F_efmit);
                            r_N_eit = REAL(v_N_eit);
                            r_Z_eit = REAL(v_Z_eit);

                            //facteurs des indices
                            PROTECT(cFACT1 = iDim(dimC));
                            PROTECT(cFACT2 = iDim(dim_F_efmit));
                            PROTECT(cFACT3 = iDim(dim_N_eit));
                            PROTECT(cFACT4 = iDim(dim_Z_eit));

                            int *fact1_C = INTEGER(cFACT1);
                            int *fact2_C = INTEGER(cFACT2);
                            int *fact3_C = INTEGER(cFACT3);
                            int *fact4_C = INTEGER(cFACT4);

                            double *r_Fot_i = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44),Qt));

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {


                             if (Qt==0) rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                  rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]);

                            }

                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                    if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                    temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                            if (Qt==0) rans_C_eit[ind_i + ind_t*nbI] = 0.0;

                            rans_C_eit[ind_i + ind_t*nbI] = rans_C_eit[ind_i + ind_t*nbI] +
                                    ((temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]);

                        }

                    if (ind_t==0 & Qt==0){

                            setAttrib(ans_C_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_C_efmit, install("DimCst"), dimCst);

                            setAttrib(ans_C_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_C_eit, install("DimCst"), dimCst2);


                            SET_VECTOR_ELT(out_C_efmit, e, ans_C_efmit);
                            SET_VECTOR_ELT(out_C_eit, e, ans_C_eit);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 23), Qt, v_F_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 24), Qt, v_N_eit);
                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 25), Qt, v_Z_eit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 26, cFACT1);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 27, cFACT2);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 28, cFACT3);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 29, cFACT4);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de Y_efmit
                            //---------

                        //on considère les dimensions de C, Y, D et L homogènes sur tout le module --> pas besoin de les redéfinir

                    if (ind_t==0 & Qt==0){

                            PROTECT(ans_Y_efmit = NEW_NUMERIC(prod));//Rprintf("JJ5");
                            setAttrib(ans_Y_efmit, R_DimSymbol, Dim);

                            rans_Y_efmit = REAL(ans_Y_efmit);

                            PROTECT(ans_Y_eit = NEW_NUMERIC(nbI*nbT));
                            setAttrib(ans_Y_eit, R_DimSymbol, Dim2);

                            rans_Y_eit = REAL(ans_Y_eit);


                    } else {

                            rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e));//Rprintf("JJ6");
                            rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e));
                    }


                            r_wL_ei = REAL(v_wL_ei);

                            //facteurs des indices
                            PROTECT(cFACT5 = iDim(dim_wL_ei));

                            int *fact5_C = INTEGER(cFACT5);

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                if (Qt==0) rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    //rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000);
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);

                            }


                           for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                               double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                    if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                    temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                               if (Qt==0) rans_Y_eit[ind_i + ind_t*nbI] = 0.0;

                                  rans_Y_eit[ind_i + ind_t*nbI] = rans_Y_eit[ind_i + ind_t*nbI] +
                                    (r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    //rans_C_eit[ind_i + ind_t*nbI] / 1000;
                                    ((temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);
                           }

                    if (ind_t==0 & Qt==0) {

                            setAttrib(ans_Y_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_Y_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_Y_efmit, e, ans_Y_efmit);

                            setAttrib(ans_Y_eit, R_DimNamesSymbol, dimnames2);
                            setAttrib(ans_Y_eit, install("DimCst"), dimCst2);

                            SET_VECTOR_ELT(out_Y_eit, e, ans_Y_eit);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 30, v_wL_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 31, cFACT5);


                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de D_efmit
                            //---------

                    if (ind_t==0 & Qt==0) {//Rprintf("JJ2");

                            PROTECT(ans_D_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_D_efmit, R_DimSymbol, Dim);

                            rans_D_efmit = REAL(ans_D_efmit);

                    } else {

                            rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e));//Rprintf("JJ3");

                    }

                            r_wD_ei = REAL(v_wD_ei);
                            r_d_efmit = REAL(v_d_efmit);

                            //facteurs des indices
                            PROTECT(cFACT6 = iDim(dim_d_efmit));
                            PROTECT(cFACT7 = iDim(dim_wD_ei));

                            int *fact6_C = INTEGER(cFACT6);
                            int *fact7_C = INTEGER(cFACT7);


                            //équation : 2 manières de calculer selon la disponibilité de wD_i

                    if (all_is_na(v_wD_ei)) { //1ère méthode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (Qt==0) rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    //rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                    (r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000));
                            }

                    } else {                 //2ème méthode

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (Qt==0) rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    //rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);
                            }

                    }


                    if (ind_t==0 & Qt==0) {

                            setAttrib(ans_D_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_D_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_D_efmit, e, ans_D_efmit);

                    }

                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 32), Qt, v_wD_ei);
                        SET_VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 33), Qt, v_d_efmit);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 34, cFACT6);
                        SET_VECTOR_ELT(VECTOR_ELT(EVAR, e), 35, cFACT7);

                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////
                        /////////////////////////////////////////////////////////////////////////////////////////////////

                             //---------
                            // calcul de L_efmit
                            //---------

                    if (ind_t==0 & Qt==0) { //Rprintf("JJ1");

                            PROTECT(ans_L_efmit = NEW_NUMERIC(prod));
                            setAttrib(ans_L_efmit, R_DimSymbol, Dim);

                            rans_L_efmit = REAL(ans_L_efmit);

                    } else {

                            rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e));

                    }

                            //équation

                            for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                            for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                  rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                    if (ind_t==0 & Qt==0) {

                            setAttrib(ans_L_efmit, R_DimNamesSymbol, dimnames);
                            setAttrib(ans_L_efmit, install("DimCst"), dimCst);

                            SET_VECTOR_ELT(out_L_efmit, e, ans_L_efmit);
                            SET_STRING_ELT(rnames_Esp, e, STRING_ELT(sppList,e));

                    }

                        if (ind_t==0 & Qt==0) UNPROTECT(8);
                        UNPROTECT(25);

         }
}

if (Qt==3) cUpdate = false;

if (ind_t==0) UNPROTECT(1);


} else {

    for (int e = 0 ; e < nbE ; e++) {

        if (trim[e]==0) { //--------------------------------------------------------------------------- pas de dimension trimestre

            if (Qt==0) {
                                int nbI = length(VECTOR_ELT(namDC,e));

                            SEXP elmt;
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                                SEXP v_wD_ei;

                                PROTECT(v_wD_ei = getListElement(elmt, "wD_i"));

                                double  *rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit,e)),
                                        *rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e)),
                                        *rans_C_eit = REAL(VECTOR_ELT(out_C_eit,e)),
                                        *rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e)),
                                        *rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e)),
                                        *rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e)),
                                        *r_F_efmit = REAL(VECTOR_ELT(out_F_fmi,e)),
                                        *r_N_eit = REAL(VECTOR_ELT(out_N_eit, e)),
                                        *r_Z_eit = REAL(VECTOR_ELT(out_Z_eit, e)),
                                        *r_wL_ei = REAL(getListElement(elmt, "wL_i")),
                                        *r_wD_ei = REAL(getListElement(elmt, "wD_i")),
                                        *r_d_efmit = REAL(getListElement(elmt, "d_i"));

                                double *r_Fot_i = REAL(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44));

                                int     *fact1_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 26)),
                                        *fact2_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 27)),
                                        *fact3_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 28)),
                                        *fact4_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 29)),
                                        *fact5_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 31)),
                                        *fact6_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 34)),
                                        *fact7_C  = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 35));

                                //équation n°1

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                     rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                        r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];



                                //équation

                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                    double temp = 0.0;

                                    for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                    for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                        if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                        temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                    }

                                    rans_C_eit[ind_i + ind_t*nbI] =
                                        (temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                        (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                        r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]];

                                }


                               //équation n°2

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                      rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                        rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;


                              for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                      rans_Y_eit[ind_i + ind_t*nbI] =
                                        r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                        rans_C_eit[ind_i + ind_t*nbI] / 1000;




                               //équation n°3

                            if (all_is_na(v_wD_ei)) {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                }

                            } else {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                      r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;


                                      rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                        r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                        rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                }

                            }

                               //équation n°4

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                      rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            UNPROTECT(2);
        }
    } else { //--------------------------------------------------------------------------- dimension trimestre


                            int nbI = length(VECTOR_ELT(namDC,e));

                            SEXP elmt;
                            PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

                                SEXP v_wD_ei;

                                PROTECT(v_wD_ei = getListElement(getListElement(elmt, "wD_i"),CHAR(STRING_ELT(trimInt,Qt))));

                                double  *rans_C_efmit = REAL(VECTOR_ELT(out_C_efmit,e)),
                                        *rans_Y_efmit = REAL(VECTOR_ELT(out_Y_efmit,e)),
                                        *rans_C_eit = REAL(VECTOR_ELT(out_C_eit,e)),
                                        *rans_Y_eit = REAL(VECTOR_ELT(out_Y_eit,e)),
                                        *rans_D_efmit = REAL(VECTOR_ELT(out_D_efmit,e)),
                                        *rans_L_efmit = REAL(VECTOR_ELT(out_L_efmit,e)),
                                        *r_F_efmit = REAL(VECTOR_ELT(VECTOR_ELT(out_F_fmi,e),Qt)),
                                        *r_N_eit = REAL(VECTOR_ELT(VECTOR_ELT(out_N_eitQ, e),Qt)),
                                        *r_Z_eit = REAL(VECTOR_ELT(VECTOR_ELT(out_Z_eit, e),Qt)),
                                        *r_wL_ei = REAL(getListElement(getListElement(elmt, "wL_i"),CHAR(STRING_ELT(trimInt,Qt)))),
                                        *r_wD_ei = REAL(getListElement(getListElement(elmt, "wD_i"),CHAR(STRING_ELT(trimInt,Qt)))),
                                        *r_d_efmit = REAL(getListElement(getListElement(elmt, "d_i"),CHAR(STRING_ELT(trimInt,Qt))));

                                double *r_Fot_i = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(EVAR, e), 44),Qt));

                                int     *fact1_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 26)),
                                        *fact2_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 27)),
                                        *fact3_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 28)),
                                        *fact4_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 29)),
                                        *fact5_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 31)),
                                        *fact6_C = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 34)),
                                        *fact7_C  = INTEGER(VECTOR_ELT(VECTOR_ELT(EVAR, e), 35));

                                //équation n°1

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                 if (Qt==0) rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                  rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]);

                                }




                                //équation


                            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                    if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                    temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                            if (Qt==0) rans_C_eit[ind_i + ind_t*nbI] = 0.0;

                            rans_C_eit[ind_i + ind_t*nbI] = rans_C_eit[ind_i + ind_t*nbI] +
                                    ((temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]);

                            }


                               //équation n°2

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (Qt==0) rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1]  + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    //rans_C_efmit[ind_f*fact1_C[0]  + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000);
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);

                            }

                              for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                                double temp = 0.0;

                                for (int ind_f = 0 ; ind_f < (1 + (nbF-1)*(fact2_C[0]>0)) ; ind_f++)
                                for (int ind_m = 0 ; ind_m < (1 + (nbM-1)*(fact2_C[1]>0)) ; ind_m++) {

                                    if (!ISNA(r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]]))
                                    temp = temp + r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]];

                                }

                               if (Qt==0) rans_Y_eit[ind_i + ind_t*nbI] = 0.0;

                                  rans_Y_eit[ind_i + ind_t*nbI] = rans_Y_eit[ind_i + ind_t*nbI] +
                                    (r_wL_ei[0*fact5_C[0]  + 0*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *   //attention : poids individuels en kg
                                    //rans_C_eit[ind_i + ind_t*nbI] / 1000;
                                    ((temp + r_Fot_i[ind_i + ind_t*nbI]) * r_N_eit[0*fact3_C[0] + 0*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[0*fact4_C[0] + 0*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);


                              }


                               //équation n°3

                            if (all_is_na(v_wD_ei)) {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                    if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                    if (Qt==0) rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                  rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    //rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];
                                    (r_wL_ei[ind_f*fact5_C[0]  + ind_m*fact5_C[1] + ind_i*fact5_C[2] + ind_t*fact5_C[3]] *
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000));
                                }

                            } else {

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++){

                                if (ISNA(r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]]))  //si NA, alors pas de rejets
                                  r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] = 0.0;

                                if (Qt==0) rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] = 0.0;

                                rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                    rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] +
                                    (r_wD_ei[ind_f*fact7_C[0]  + ind_m*fact7_C[1] + ind_i*fact7_C[2] + ind_t*fact7_C[3]] *
                                    r_d_efmit[ind_f*fact6_C[0] + ind_m*fact6_C[1] + ind_i*fact6_C[2] + ind_t*fact6_C[3]] *
                                    //rans_C_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] / 1000;
                                    (r_F_efmit[ind_f*fact2_C[0] + ind_m*fact2_C[1] + ind_i*fact2_C[2] + ind_t*fact2_C[3]] *
                                    r_N_eit[ind_f*fact3_C[0] + ind_m*fact3_C[1] + ind_i*fact3_C[2] + ind_t*fact3_C[3]] *
                                    (1 - exp( -r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]])) /
                                    r_Z_eit[ind_f*fact4_C[0] + ind_m*fact4_C[1] + ind_i*fact4_C[2] + ind_t*fact4_C[3]]) / 1000);
                                }

                            }

                               //équation n°4

                                for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
                                for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
                                for (int ind_i = 0 ; ind_i < nbI ; ind_i++)

                                      rans_L_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] =
                                        rans_Y_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]] -
                                        rans_D_efmit[ind_f*fact1_C[0] + ind_m*fact1_C[1] + ind_i*fact1_C[2] + ind_t*fact1_C[3]];

                            UNPROTECT(2);

    }
    }
}

}}







//
////------------------------------------------
//// Module 'Marché'
////------------------------------------------
//
extern "C" {

void BioEcoPar::Marche(SEXP list, int ind_t)
{


    SEXP    elmt, intC, v_P_fmce, v_icat, v_L_efmit, dimCst_P_fmce, dimCst_L_efmit, dimCst_L_efmct, Dim_L_efmct,
            ans_L_efmct = R_NilValue, dimnames_Lc = R_NilValue, rnames_Esp, cFACTc, cFACTi; //dimnames_Lc2 = R_NilValue,
            //v_mme, dimCst_mme, cFACTc2, cFACTmm;

    //SEXP ans_L_efmct2 = R_NilValue,  dimCst_L_efmct2, Dim_L_efmct2;

    int *dim_P_fmce, *dim_L_efmit, *dim_icat, *dim_L_efmct, *dimLc;
    //int *dim_L_efmct2, *dim_mme, *dimLc2, *r_mme;

    int nbI, nbC;

    double *rans_L_efmct, *r_L_efmit, *r_P_fmce, *r_icat;
    //double *rans_L_efmct2;
//Rprintf("CCC1");

if (ind_t==0) {


    PROTECT(rnames_Esp = getAttrib(out_L_efmit,R_NamesSymbol));//PrintValue(rnames_Esp);
    setAttrib(out_L_efmct, R_NamesSymbol, rnames_Esp);

    setAttrib(out_P_t, R_NamesSymbol, rnames_Esp);

    setAttrib(out_L_efmct2, R_NamesSymbol, rnames_Esp);

}


    for (int e = 0 ; e < nbE ; e++) {

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        intC = getListElement(elmt, "modC");
        nbC = length(intC);

        PROTECT(v_P_fmce = getListElement(elmt, "P_fmce"));//Rprintf("CCC1");
        PROTECT(v_icat = getListElement(elmt, "icat"));   //qqsoit i, sum_c icat = 1
        PROTECT(v_L_efmit = getListElement(out_L_efmit, CHAR(STRING_ELT(sppList,e))));//Rprintf("BBB1");
//        PROTECT(v_mme = getListElement(elmt, "mm"));

        PROTECT(dimCst_P_fmce = getAttrib(v_P_fmce, install("DimCst")));
        PROTECT(dimCst_L_efmit = getAttrib(v_L_efmit, install("DimCst")));
//        PROTECT(dimCst_mme = getAttrib(v_mme, install("DimCst")));

        //tests sur les dimensions :
        dim_P_fmce = INTEGER(dimCst_P_fmce);//Rprintf("AAA1");
        if ((dim_P_fmce[0]!=0 & dim_P_fmce[0]!=nbF) | (dim_P_fmce[1]!=0 & dim_P_fmce[1]!=nbMe) |
            (dim_P_fmce[2]!=0 & dim_P_fmce[2]!=nbC) | (dim_P_fmce[3]!=0 & dim_P_fmce[3]!=nbT))
        {
            error("Non_homogeneous dimensions in P_fmce element. Check .ini biological parameters files !!\n");
        }

        dim_L_efmit = INTEGER(dimCst_L_efmit);//Rprintf("AAA2");
        if ((dim_L_efmit[0]!=0 & dim_L_efmit[0]!=nbF) | (dim_L_efmit[1]!=0 & dim_L_efmit[1]!=nbM) |
            (dim_L_efmit[2]!=0 & dim_L_efmit[2]!=nbI) | (dim_L_efmit[3]!=0 & dim_L_efmit[3]!=nbT))
        {
            error("Non_homogeneous dimensions in L_efmit element. Check .ini biological parameters files !!\n");
        }

        dim_icat = INTEGER(getAttrib(v_icat, R_DimSymbol));//Rprintf("AAA3");
        if ((dim_icat[0]!=nbI) & (dim_icat[1]!=nbC))
        {
            error("Non_homogeneous dimensions in icat element. Check .ini biological parameters files !!\n");
        }

//        dim_mme = INTEGER(getAttrib(v_mme, R_DimSymbol));
//        if ((dim_mme[0]!=nbF) & (dim_mme[1]!=nbM))
//        {
//            error("Non_homogeneous dimensions in mm element. Check .ini biological parameters files !!\n");
//        }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

        //---------
        // calcul de L_efmct
        //---------

        PROTECT(dimCst_L_efmct = allocVector(INTSXP, 4));
//        PROTECT(dimCst_L_efmct2 = allocVector(INTSXP, 4));
        dim_L_efmct = INTEGER(dimCst_L_efmct);//Rprintf("AAA4");
//        dim_L_efmct2 = INTEGER(dimCst_L_efmct2);
        dim_L_efmct[0] = dim_L_efmit[0] ; dim_L_efmct[1] = dim_L_efmit[1] ; dim_L_efmct[2] = nbC; dim_L_efmct[3] = dim_L_efmit[3];
//        dim_L_efmct2[0] = dim_L_efmit[0] ; dim_L_efmct2[1] = nbMe*(dim_L_efmit[1]>0) ; dim_L_efmct2[2] = nbC; dim_L_efmct2[3] = dim_L_efmit[3];

        int count = 0, prod = 1, count2 = 0, count3 = 0; //prod2 = 1, count22 = 0,
        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                count++;
                prod = prod * dim_L_efmct[k];
//                prod2 = prod2 * dim_L_efmct2[k];
            }

        }

        PROTECT(Dim_L_efmct = allocVector(INTSXP, count));
        dimLc = INTEGER(Dim_L_efmct);//Rprintf("AAA5");
//        PROTECT(Dim_L_efmct2 = allocVector(INTSXP, count));
//        dimLc2 = INTEGER(Dim_L_efmct2);


        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                dimLc[count2] = dim_L_efmct[k];
                count2++;
            }

//            if (dim_L_efmct2[k]>0) {
//                dimLc2[count22] = dim_L_efmct2[k];
//                count22++;
//            }

        }


if (ind_t==0){

        //on crée le tableau résultat pour l'espèce en question
        PROTECT(ans_L_efmct = NEW_NUMERIC(prod));
        setAttrib(ans_L_efmct, R_DimSymbol, Dim_L_efmct);
//        PROTECT(ans_L_efmct2 = NEW_NUMERIC(prod2));
//        setAttrib(ans_L_efmct2, R_DimSymbol, Dim_L_efmct2);

        PROTECT(dimnames_Lc = allocVector(VECSXP,count));
        if (dim_L_efmct[0]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, fleetList) ; count3++;}
        if (dim_L_efmct[1]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, metierList) ; count3++;}
        if (dim_L_efmct[2]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, intC) ; count3++;}
        if (dim_L_efmct[3]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, times) ; count3++;}

//        count3 = 0;
//        PROTECT(dimnames_Lc2 = allocVector(VECSXP,count));
//        if (dim_L_efmct2[0]>0) {SET_VECTOR_ELT(dimnames_Lc2, count3, fleetList) ; count3++;}
//        if (dim_L_efmct2[1]>0) {SET_VECTOR_ELT(dimnames_Lc2, count3, metierListEco) ; count3++;}
//        if (dim_L_efmct2[2]>0) {SET_VECTOR_ELT(dimnames_Lc2, count3, intC) ; count3++;}
//        if (dim_L_efmct2[3]>0) {SET_VECTOR_ELT(dimnames_Lc2, count3, times) ; count3++;}

        rans_L_efmct = REAL(ans_L_efmct);
//        rans_L_efmct2 = REAL(ans_L_efmct2);

} else {

        rans_L_efmct = REAL(VECTOR_ELT(out_L_efmct, e));
//        rans_L_efmct2 = REAL(VECTOR_ELT(out_L_efmct2, e));

}

        r_L_efmit = REAL(v_L_efmit);
        r_P_fmce = REAL(v_P_fmce);
        r_icat = REAL(v_icat);
//        r_mme = INTEGER(AS_INTEGER(v_mme));
//        double *r_mme2 = REAL(v_mme);

        //facteurs des indices
        PROTECT(cFACTc = iDim(dim_L_efmct));
//        PROTECT(cFACTc2 = iDim(dim_L_efmct2));
        PROTECT(cFACTi = iDim(dim_L_efmit));
//        PROTECT(cFACTmm = iDim(dim_mme));

        int *fact_Cc = INTEGER(cFACTc);
        int *fact_Ci = INTEGER(cFACTi);//Rprintf("AAA6");
//        int *fact_Cc2 = INTEGER(cFACTc2);

        //équation n°1 : conversion âge/catgégorie

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
        for (int ind_c = 0 ; ind_c < nbC ; ind_c++) {

            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                if (ind_i ==0) {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                } else {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                }
            }
        }




if (ind_t==0) {

        setAttrib(ans_L_efmct, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_L_efmct, install("DimCst"), dimCst_L_efmct);

//        setAttrib(ans_L_efmct2, R_DimNamesSymbol, dimnames_Lc2);
//        setAttrib(ans_L_efmct2, install("DimCst"), dimCst_L_efmct2);

        SET_VECTOR_ELT(out_L_efmct, e, ans_L_efmct);
        SET_VECTOR_ELT(out_L_efmct2, e, ans_L_efmct);
        SET_VECTOR_ELT(out_P_t, e, v_P_fmce);

}

if (ind_t==0) UNPROTECT(2);//4);
UNPROTECT(10);//6);

    }

if (ind_t==0) UNPROTECT(1);

}}




//extern "C" {
//
//void BioEcoPar::Marche(SEXP list, int t)
//{
//
//if (t==0){
//
//    SEXP    ans_1, elmt,
//            dimCst, Dim, dimnames, dimCst_L_efmit, dimCst_cat_i, dimCst_alpha_i, dimCst_beta_i, dimCst_gamma_i, dimCst_P_it, intAge,
//            v_L_efmit, v_cat_i, v_alpha_i, v_beta_i, v_gamma_i, v_P_it, tab_sum_i, tab_sum_not_i, dimCoeff;
//
//    SEXP rnames;
//
//    int *dim_L_efmit, *dim_cat_i, *dim_alpha_i, *dim_beta_i, *dim_gamma_i, *dim_P_it, *dimC, *dimCo;
//    int nbI;
//
//    double *rans_1, *r_L_efmit, *r_alpha_i, *r_beta_i, *r_gamma_i, *r_P_it, *sum_i, *sum_not_i, *rtab_sum_i, *rtab_sum_not_i;
//
//    PROTECT(out_P_t = allocVector(VECSXP, nbE));
//    PROTECT(rnames = allocVector(STRSXP, nbE));
//    setAttrib(out_P_t , R_NamesSymbol, rnames);
//
//    for (int e = 0 ; e < nbE ; e++) {
//
//        //---------
//        // calcul de P_eit
//        //---------
//
//        elmt = getListElement(bioList, CHAR(STRING_ELT(sppList,e)));
//        intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e)));
//
//        nbI = length(getListElement(elmt, "age"));
//
//        v_cat_i = getListElement(elmt, "cat_i");
//        v_alpha_i = getListElement(elmt, "alpha_c");    //attention : nom de variable à remettre à jour
//        v_beta_i = getListElement(elmt, "beta_c");
//        v_gamma_i = getListElement(elmt, "gamma_c");
//        v_P_it = getListElement(elmt, "P_ct");
//        v_L_efmit = getListElement( out_L_efmit , CHAR(STRING_ELT(sppList,e))) ;
//
//        dimCst_cat_i = getAttrib(v_cat_i, install("DimCst"));
//        dimCst_alpha_i = getAttrib(v_alpha_i, install("DimCst"));
//        dimCst_beta_i = getAttrib(v_beta_i, install("DimCst"));
//        dimCst_gamma_i = getAttrib(v_gamma_i, install("DimCst"));
//        dimCst_P_it = getAttrib(v_P_it, install("DimCst"));
//        dimCst_L_efmit = getAttrib(v_L_efmit, install("DimCst"));
//
//        //tests sur les dimensions
//        dim_cat_i = INTEGER(dimCst_cat_i);
//        if ((dim_cat_i[0]!=0) | (dim_cat_i[1]!=0) |
//            (dim_cat_i[2]!=0 & dim_cat_i[2]!=nbI) | (dim_cat_i[3]!=0))
//        {
//            error("Non_homogeneous dimensions in cat_i element. Check .ini biological parameters files !!\n");
//        }
//
//        dim_alpha_i = INTEGER(dimCst_alpha_i);
//        if ((dim_alpha_i[0]!=0 & dim_alpha_i[0]!=nbF) | (dim_alpha_i[1]!=0 & dim_alpha_i[1]!=nbM) |
//            (dim_alpha_i[2]!=0 & dim_alpha_i[2]!=nbI) | (dim_alpha_i[3]!=0 & dim_alpha_i[3]!=nbT))
//        {
//            error("Non_homogeneous dimensions in alpha_mi element. Check .ini biological parameters files !!\n");
//        }
//
//        dim_beta_i = INTEGER(dimCst_beta_i); //les facteurs alpha, beta et gamma doivent avoir même dimension
//        if ((dim_beta_i[0]!=dim_alpha_i[0]) | (dim_beta_i[1]!=dim_alpha_i[1]) |
//            (dim_beta_i[2]!=dim_alpha_i[2]) | (dim_beta_i[3]!=dim_alpha_i[3]))
//        {
//            error("Non_homogeneous dimensions in beta_mi element. Check .ini biological parameters files !!\n");
//        }
//
//        dim_gamma_i = INTEGER(dimCst_gamma_i);
//        if ((dim_gamma_i[0]!=dim_alpha_i[0]) | (dim_gamma_i[1]!=dim_alpha_i[1]) |
//            (dim_gamma_i[2]!=dim_alpha_i[2]) | (dim_gamma_i[3]!=dim_alpha_i[3]))
//        {
//            error("Non_homogeneous dimensions in gamma_mi element. Check .ini biological parameters files !!\n");
//        }
//
//        dim_P_it = INTEGER(dimCst_P_it);
//        if ((dim_P_it[0]!=0 & dim_P_it[0]!=nbF) | (dim_P_it[1]!=0 & dim_P_it[1]!=nbM) |
//            (dim_P_it[2]!=0 & dim_P_it[2]!=nbI) | (dim_P_it[3]!=0 & dim_P_it[3]!=nbT))
//        {
//            error("Non_homogeneous dimensions in P_mit element. Check .ini biological parameters files !!\n");
//        }
//
//        dim_L_efmit = INTEGER(dimCst_L_efmit);
//        if ((dim_L_efmit[0]!=0 & dim_L_efmit[0]!=nbF) | (dim_L_efmit[1]!=0 & dim_L_efmit[1]!=nbM) |
//            (dim_L_efmit[2]!=nbI) | (dim_L_efmit[3]!=0 & dim_L_efmit[3]!=nbT))
//        {
//            error("Non_homogeneous dimensions in L_efmit element. Check .ini biological parameters files !!\n");
//        }
//
//        //on détermine l'attribut Dimension du tableau résultant -> dimCst (on en profite pour compter les dimensions réelles + nombre de cellules)
//        PROTECT(dimCst = allocVector(INTSXP, 4));
//        dimC = INTEGER(dimCst);
//        dimC[0] = nbF; dimC[1] = nbM; dimC[2] = nbI; dimC[3] =nbT;
//        int count = 0, prod = 1, count2 = 0, count3 = 0;
//
//        for (int k = 0 ; k < 4 ; k++) {
//
//            if (dimC[k]>0) {
//                count++;
//                prod = prod * dimC[k];
//            }
//
//        }
//
//        PROTECT(Dim = allocVector(INTSXP, count));
//        int *dim = INTEGER(Dim);
//
//        for (int k = 0 ; k < 4 ; k++) {
//            if (dimC[k]>0) {
//                dim[count2] = dimC[k];
//                count2++;
//                }
//        }
//
//
//        //on crée le tableau résultat pour l'espèce en question -> ans_1
//        ans_1 = PROTECT(NEW_NUMERIC(prod));
//        setAttrib(ans_1, R_DimSymbol, Dim);
//
//        PROTECT(dimnames = allocVector(VECSXP,count));
//        if (dimC[0]>0) {SET_VECTOR_ELT(dimnames, count3, getListElement(paramList, "Fleet")) ; count3++;}
//        if (dimC[1]>0) {SET_VECTOR_ELT(dimnames, count3, getListElement(paramList, "Metier")) ; count3++;}
//        if (dimC[2]>0) {SET_VECTOR_ELT(dimnames, count3, intAge) ; count3++;}
//        if (dimC[3]>0) {SET_VECTOR_ELT(dimnames, count3, getListElement(paramList, "times")) ; count3++;}
//
//        rans_1 = REAL(ans_1);
//        r_L_efmit = REAL(v_L_efmit);
//        r_alpha_i = REAL(v_alpha_i);
//        r_beta_i = REAL(v_beta_i);
//        r_gamma_i = REAL(v_gamma_i);
//        r_P_it = REAL(v_P_it);
//
//        //facteurs des indices
//        fact1_P = iDim(dimC);
//        fact2_P = iDim(dim_alpha_i);
//    //   fact3_P = iDim(dim_beta_i);
//    //   fact4_P = iDim(dim_gamma_i);
//        fact5_P = iDim(dim_P_it);
//        fact6_P = iDim(dim_L_efmit);
//
//
//        //il faut avant tout créer les tableaux sum_L_fmeit et sumNot_L_fmeit
//            //1ère étape : somme sur les âges de chaque classe
//        tab_sum_i = PROTECT(NEW_NUMERIC(prod));
//        tab_sum_not_i = PROTECT(NEW_NUMERIC(prod));
//        rtab_sum_i = REAL(tab_sum_i);
//        rtab_sum_not_i = REAL(tab_sum_not_i);
//
//                //initialisation
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//            rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] = 0.0;
//            rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] = 0.0;
//
//        }
//
//                //somme sur i
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//        //classe associée à l'âge
//
//        for (int I = 0 ; I < nbI ; I++) {
//
//            if (CHAR(STRING_ELT(v_cat_i, I))==CHAR(STRING_ELT(v_cat_i, ind_i))){   //ATTENTION : ceci implique que cat_i reste un vecteur par âge --> pas d'autres déclinaisons
//
//                rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] =
//                 rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] +
//                 r_L_efmit[ind_f*fact6_P[0] + ind_m*fact6_P[1] + I*fact6_P[2] + 0*fact6_P[3]];
//
//            } else {
//
//                rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] =
//                 rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] +
//                 r_L_efmit[ind_f*fact6_P[0] + ind_m*fact6_P[1] + I*fact6_P[2] + 0*fact6_P[3]];
//
//            }
//
//        }
//        }
//            //2ème étape : on agrège en fonction des dimensions des coefficients
//
//        PROTECT(dimCoeff = allocVector(INTSXP, 4));
//        dimCo = INTEGER(dimCoeff);
//        dimCo[0] = dim_alpha_i[0]; dimCo[1] = dim_alpha_i[1]; dimCo[2] = nbI; dimCo[3] = nbT;
//        setAttrib(tab_sum_i, install("DimCst"), dimCst);
//        setAttrib(tab_sum_not_i, install("DimCst"), dimCst);
//
//        sum_i = REAL(aggregObj(tab_sum_i,dimCoeff));
//        sum_not_i = REAL(aggregObj(tab_sum_not_i,dimCoeff));
//
//        //équation
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//            if (!ISNA(r_P_it[ind_f*fact5_P[0] + ind_m*fact5_P[1] + ind_i*fact5_P[2] + 0*fact5_P[3]])) {
//
//                rans_1[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] =
//                    r_P_it[ind_f*fact5_P[0] + ind_m*fact5_P[1] + ind_i*fact5_P[2] + 0*fact5_P[3]];
//
//            } else {
//
//                rans_1[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + 0*fact1_P[3]] =
//                 exp(r_alpha_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + 0*fact2_P[3]] +
//                  r_beta_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + 0*fact2_P[3]] *
//                  log(sum_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + 0*fact2_P[3]]) +
//                  r_gamma_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + 0*fact2_P[3]] *
//                  log(sum_not_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + 0*fact2_P[3]]));
//            }
//        }
//
//        setAttrib(ans_1, R_DimNamesSymbol, dimnames);
//        setAttrib(ans_1, install("DimCst"), dimCst);
//
//        SET_VECTOR_ELT(out_P_t, e, ans_1);
//        SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
//
//        UNPROTECT(7);
//    }
//
//    UNPROTECT(2);
//
//
//
//} else {
//
//
//
//
//    SEXP    ans_1, elmt,
//            dimCst, dimCst_L_efmit, dimCst_cat_i, dimCst_alpha_i, dimCst_beta_i, dimCst_gamma_i, dimCst_P_it, intAge,
//            v_L_efmit, v_cat_i, v_alpha_i, v_beta_i, v_gamma_i, v_P_it, tab_sum_i, tab_sum_not_i, dimCoeff;
//
//    int *dim_alpha_i, *dimC, *dimCo;
//    int nbI;
//
//    double *rans_1, *r_L_efmit, *r_alpha_i, *r_beta_i, *r_gamma_i, *r_P_it, *sum_i, *sum_not_i, *rtab_sum_i, *rtab_sum_not_i;
//
//    for (int e = 0 ; e < nbE ; e++) {
//
//        //---------
//        // calcul de P_eit
//        //---------
//
//        elmt = getListElement(bioList, CHAR(STRING_ELT(sppList,e)));
//        intAge = getListElement(namDC, CHAR(STRING_ELT(sppList,e)));
//
//        nbI = length(getListElement(elmt, "age"));
//
//        v_cat_i = getListElement(elmt, "cat_i");
//        v_alpha_i = getListElement(elmt, "alpha_c");    //attention : nom de variable à remettre à jour
//        v_beta_i = getListElement(elmt, "beta_c");
//        v_gamma_i = getListElement(elmt, "gamma_c");
//        v_P_it = getListElement(elmt, "P_ct");
//        v_L_efmit = getListElement( out_L_efmit, CHAR(STRING_ELT(sppList,e))) ;
//
//        dimCst_cat_i = getAttrib(v_cat_i, install("DimCst"));
//        dimCst_alpha_i = getAttrib(v_alpha_i, install("DimCst"));
//        dimCst_beta_i = getAttrib(v_beta_i, install("DimCst"));
//        dimCst_gamma_i = getAttrib(v_gamma_i, install("DimCst"));
//        dimCst_P_it = getAttrib(v_P_it, install("DimCst"));
//        dimCst_L_efmit = getAttrib(v_L_efmit, install("DimCst"));
//
//        dim_alpha_i = INTEGER(dimCst_alpha_i);
//
//        //on crée le tableau résultat pour l'espèce en question -> ans_1
//        ans_1 = getListElement(out_P_t, CHAR(STRING_ELT(sppList,e)));
//
//        rans_1 = REAL(ans_1);
//        r_L_efmit = REAL(v_L_efmit);
//        r_alpha_i = REAL(v_alpha_i);
//        r_beta_i = REAL(v_beta_i);
//        r_gamma_i = REAL(v_gamma_i);
//        r_P_it = REAL(v_P_it);
//
//        PROTECT(dimCst = allocVector(INTSXP, 4));
//        dimC = INTEGER(dimCst);
//        dimC[0] = nbF; dimC[1] = nbM; dimC[2] = nbI; dimC[3] =nbT;
//        int count = 0, prod = 1;
//
//        for (int k = 0 ; k < 4 ; k++) {
//
//            if (dimC[k]>0) {
//                count++;
//                prod = prod * dimC[k];
//            }
//
//        }
//
//      //il faut avant tout créer les tableaux sum_L_fmeit et sumNot_L_fmeit
//            //1ère étape : somme sur les ages de chaque classe
//        tab_sum_i = PROTECT(NEW_NUMERIC(prod));
//        tab_sum_not_i = PROTECT(NEW_NUMERIC(prod));
//        rtab_sum_i = REAL(tab_sum_i);
//        rtab_sum_not_i = REAL(tab_sum_not_i);
//
//                //initialisation
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//            rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] = 0.0;
//            rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] = 0.0;
//
//        }
//
//                //somme sur i
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//        //classe associée à l'âge
//
//        for (int I = 0 ; I < nbI ; I++) {
//
//            if (CHAR(STRING_ELT(v_cat_i, I))==CHAR(STRING_ELT(v_cat_i, ind_i))){   //ATTENTION : ceci implique que cat_i reste un vecteur par âge --> pas d'autres déclinaisons
//
//                rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] =
//                 rtab_sum_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] +
//                 r_L_efmit[ind_f*fact6_P[0] + ind_m*fact6_P[1] + I*fact6_P[2] + t*fact6_P[3]];
//
//            } else {
//
//                rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] =
//                 rtab_sum_not_i[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] +
//                 r_L_efmit[ind_f*fact6_P[0] + ind_m*fact6_P[1] + I*fact6_P[2] + t*fact6_P[3]];
//
//            }
//
//        }
//        }
//            //2ème étape : on agrège en fonction des dimensions des coefficients
//
//        PROTECT(dimCoeff = allocVector(INTSXP, 4));
//        dimCo = INTEGER(dimCoeff);
//        dimCo[0] = dim_alpha_i[0]; dimCo[1] = dim_alpha_i[1]; dimCo[2] = nbI; dimCo[3] = nbT;
//        setAttrib(tab_sum_i, install("DimCst"), dimCst);
//        setAttrib(tab_sum_not_i, install("DimCst"), dimCst);
//
//        sum_i = REAL(aggregObj(tab_sum_i,dimCoeff));
//        sum_not_i = REAL(aggregObj(tab_sum_not_i,dimCoeff));
//
//        //équation
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
//        for (int ind_i = 0 ; ind_i < nbI ; ind_i++){
//
//            if (!ISNA(r_P_it[ind_f*fact5_P[0] + ind_m*fact5_P[1] + ind_i*fact5_P[2] + t*fact5_P[3]])) {
//
//                rans_1[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] =
//                    r_P_it[ind_f*fact5_P[0] + ind_m*fact5_P[1] + ind_i*fact5_P[2] + t*fact5_P[3]];
//
//            } else {
//
//                rans_1[ind_f*fact1_P[0] + ind_m*fact1_P[1] + ind_i*fact1_P[2] + t*fact1_P[3]] =
//                 exp(r_alpha_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + t*fact2_P[3]] +
//                  r_beta_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + t*fact2_P[3]] *
//                  log(sum_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + t*fact2_P[3]]) +
//                  r_gamma_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + t*fact2_P[3]] *
//                  log(sum_not_i[ind_f*fact2_P[0] + ind_m*fact2_P[1] + ind_i*fact2_P[2] + t*fact2_P[3]]));
//            }
//        }
//
//        UNPROTECT(4);
//    }
//}
//
//}
//}
//
//
//
//




//------------------------------------------
// Module 'Economie'
//------------------------------------------

extern "C" {

void BioEcoPar::Economic(SEXP list, int ind_t, int adj, int lev, int ue_choice, int oths, int othsFM, int perscCalc, int report, double dr)
{

    SEXP Flist;
    PROTECT(Flist = getListElement(list, "Fleet"));

    PROTECT(out_Eco);

//2 protect
    SEXP dimCstF, DimF, dimnamesF, dimCstFM, dimCstFini, dimCstFMini, DimFM, DimFMini, dimnamesFM, dimnamesFMini; //formatage des objets résultats

    SEXP eFACTf, eFACTfm, elmt;

    SEXP    Lref_f, GVLref_f, GVLref_f_m, nbv_f, nbv_f_m, lc_f, gc_f, gc_f_m, nbh_f, ue_f, ue_f_m, fc_f, fc_f_m, vf_f, vf_f_m, ovc_f, ovc_f_m,
            oilc_f, oilc_f_m, bc_f, bc_f_m, foc_f, foc_f_m, cnb_f, cnb_f_m, icec_f, icec_f_m, cshr_f, cshr_f_m, eec_f, mwh_f, altwh_f,
            rep_f, onvc_f, insp_f, ownc_f, mngc_f, licc_f, comc_f, finc_f, dep_f, ic_f, K_f, vc_f, persc_f, ecc_f, pl_f;

    SEXP    dc_Lref_f, dc_GVLref_f, dc_GVLref_f_m, dc_nbv_f, dc_nbv_f_m, dc_lc_f, dc_gc_f, dc_gc_f_m, dc_nbh_f, dc_ue_f, dc_ue_f_m, dc_fc_f,
            dc_fc_f_m, dc_vf_f, dc_vf_f_m, dc_ovc_f, dc_ovc_f_m, dc_oilc_f, dc_oilc_f_m, dc_bc_f, dc_bc_f_m, dc_foc_f, dc_foc_f_m,
            dc_cnb_f, dc_cnb_f_m, dc_icec_f, dc_icec_f_m, dc_cshr_f, dc_cshr_f_m, dc_eec_f, dc_mwh_f, dc_altwh_f, dc_rep_f, dc_onvc_f,
            dc_insp_f, dc_ownc_f, dc_mngc_f, dc_licc_f, dc_comc_f, dc_finc_f, dc_dep_f, dc_ic_f, dc_K_f, dc_vc_f, dc_persc_f, dc_ecc_f,
            dc_pl_f;

    int *dCF,*dCFM,*dCFini,*dCFMini,*DF,*DFM, *DFMini;

    int *dim_Lref_f, *dim_GVLref_f, *dim_GVLref_f_m, *dim_nbv_f, *dim_nbv_f_m, *dim_lc_f, *dim_gc_f, *dim_gc_f_m, *dim_nbh_f, *dim_ue_f, *dim_ue_f_m,
        *dim_fc_f, *dim_fc_f_m, *dim_vf_f, *dim_vf_f_m, *dim_ovc_f, *dim_ovc_f_m, *dim_oilc_f, *dim_oilc_f_m, *dim_bc_f, *dim_bc_f_m,
        *dim_foc_f, *dim_foc_f_m, *dim_cnb_f, *dim_cnb_f_m, *dim_icec_f, *dim_icec_f_m, *dim_cshr_f, *dim_cshr_f_m, *dim_eec_f, *dim_mwh_f,
        *dim_altwh_f, *dim_rep_f, *dim_onvc_f, *dim_insp_f, *dim_ownc_f, *dim_mngc_f, *dim_licc_f, *dim_comc_f, *dim_finc_f,
        *dim_dep_f, *dim_ic_f, *dim_K_f, *dim_vc_f, *dim_persc_f, *dim_ecc_f, *dim_pl_f;

    double  *r_Lref_f, *r_GVLref_f, *r_GVLref_f_m, *r_nbv_f, *r_nbv_f_m, *r_lc_f, *r_gc_f, *r_gc_f_m, *r_nbh_f, *r_ue_f, *r_ue_f_m, *r_fc_f,
            *r_fc_f_m, *r_vf_f, *r_vf_f_m, *r_ovc_f, *r_ovc_f_m, *r_oilc_f, *r_oilc_f_m, *r_bc_f, *r_bc_f_m, *r_foc_f, *r_foc_f_m,
            *r_cnb_f, *r_cnb_f_m, *r_icec_f, *r_icec_f_m, *r_cshr_f, *r_cshr_f_m, *r_eec_f, *r_mwh_f, *r_altwh_f, *r_rep_f, *r_onvc_f,
            *r_insp_f, *r_ownc_f, *r_mngc_f, *r_licc_f, *r_comc_f, *r_finc_f, *r_dep_f, *r_ic_f, *r_K_f, *r_vc_f, *r_persc_f, *r_ecc_f,
            *r_pl_f;

    double  *r_Lbio_f_out, *r_GVLtot_f_m_out, *r_GVLav_f_m_out, *r_GVLtot_f_out, *r_GVLav_f_out, *r_NGVLav_f_m_out, *r_NGVLav_f_out,
            *r_vcst_f_m_out, *r_vcst_f_out, *r_rtbs_f_m_out, *r_rtbs_f_out, *r_cshrT_f_m_out, *r_cshrT_f_out, *r_sshr_f_m_out,
            *r_sshr_f_out, *r_ncshr_f_out, *r_ocl_f_out, *r_cs_f_out, *r_csTot_f_out, *r_gva_f_out, *r_ccw_f_out, *r_ccwCr_f_out,
            *r_wageg_f_out, *r_wagen_f_out, *r_gcf_f_out, *r_ngcf_f_out, *r_gp_f_out, *r_ssTot_f_out, *r_ps_f_out, *r_sts_f_out,
            *r_ber_f_out, *r_ratio_gva_GVL_f_out, *r_ratio_gcf_GVL_f_out, *r_ratio_fc_GVL_f_out,
            *r_ratio_oilc_GVL_f_out, *r_ratio_bc_GVL_f_out, *r_ratio_foc_GVL_f_out, *r_ratio_icec_GVL_f_out, *r_ratio_gc_GVL_f_out,
            *r_ratio_vc_GVL_f_out, *r_ratio_rep_GVL_f_out, *r_ratio_mngc_GVL_f_out, *r_ratio_licc_GVL_f_out, *r_ratio_fvol_GVL_f_out,
            *r_ratio_fvol_Lbio_f_out, *r_ratio_fvol_gva_f_out, *r_ratio_gcf_gva_f_out, *r_ratio_K_cnb_f_out, *r_ratio_GVL_K_f_out,
            *r_ratio_gcf_K_f_out, *r_ratio_ngcf_K_f_out, *r_ratio_gp_K_f_out, *r_ratio_GVL_cnb_ue_f_out,
            *r_rtbsAct_f_out, *r_csAct_f_out, *r_gvaAct_f_out, *r_gcfAct_f_out, *r_psAct_f_out, *r_stsAct_f_out;


//définition des dimensions
    PROTECT(dimCstF = allocVector(INTSXP, 4));
    PROTECT(DimF = allocVector(INTSXP, 2));
    PROTECT(dimnamesF = allocVector(VECSXP,2));
    PROTECT(dimCstFM = allocVector(INTSXP, 4));
    PROTECT(DimFM = allocVector(INTSXP, 3));
    PROTECT(dimnamesFM = allocVector(VECSXP,3));

    PROTECT(dimCstFini = allocVector(INTSXP, 4));
    PROTECT(dimCstFMini = allocVector(INTSXP, 4));
    PROTECT(DimFMini = allocVector(INTSXP, 2));
    PROTECT(dimnamesFMini = allocVector(VECSXP,2));

    SET_VECTOR_ELT(dimnamesF, 0, fleetList); SET_VECTOR_ELT(dimnamesF, 1, times);
    SET_VECTOR_ELT(dimnamesFM, 0, fleetList); SET_VECTOR_ELT(dimnamesFM, 1, metierListEco); SET_VECTOR_ELT(dimnamesFM, 2, times);
    SET_VECTOR_ELT(dimnamesFMini, 0, fleetList); SET_VECTOR_ELT(dimnamesFMini, 1, metierListEco);

    dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
    dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;
    dCFini = INTEGER(dimCstFini) ; dCFini[0] = nbF; dCFini[1] = 0; dCFini[2] = 0; dCFini[3] = 0;
    dCFMini = INTEGER(dimCstFMini) ; dCFMini[0] = nbF; dCFMini[1] = nbMe; dCFMini[2] = 0; dCFMini[3] = 0;

    DF = INTEGER(DimF) ; DF[0] = nbF; DF[1] = nbT;
    DFM = INTEGER(DimFM) ; DFM[0] = nbF; DFM[1] = nbMe; DFM[2] = nbT;
    DFMini = INTEGER(DimFMini) ; DFMini[0] = nbF; DFMini[1] = nbMe;

    //facteurs des indices génériques F/FM
    PROTECT(eFACTf = iDim(dCF));
    PROTECT(eFACTfm = iDim(dCFM));

//12 protect -> 14

    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    PROTECT(Lref_f = getListElement(Flist, "Lref_f"));      PROTECT(dc_Lref_f = iDim(INTEGER(getAttrib(Lref_f, install("DimCst")))));
    PROTECT(GVLref_f = getListElement(Flist, "GVLref_f"));    PROTECT(dc_GVLref_f = iDim(INTEGER(getAttrib(GVLref_f, install("DimCst")))));
    PROTECT(GVLref_f_m = getListElement(Flist, "GVLref_f_m")); PROTECT(dc_GVLref_f_m = iDim(INTEGER(getAttrib(GVLref_f_m, install("DimCst")))));
    PROTECT(nbv_f = getListElement(Flist, "nbv_f"));        PROTECT(dc_nbv_f = iDim(INTEGER(getAttrib(nbv_f, install("DimCst")))));
    PROTECT(nbv_f_m = getListElement(Flist, "nbv_f_m"));    PROTECT(dc_nbv_f_m = iDim(INTEGER(getAttrib(nbv_f_m, install("DimCst")))));
    PROTECT(lc_f = getListElement(Flist, "lc_f"));          PROTECT(dc_lc_f = iDim(INTEGER(getAttrib(lc_f, install("DimCst")))));
    PROTECT(gc_f = getListElement(Flist, "gc_f"));          PROTECT(dc_gc_f = iDim(INTEGER(getAttrib(gc_f, install("DimCst")))));
    PROTECT(gc_f_m = getListElement(Flist, "gc_f_m"));      PROTECT(dc_gc_f_m = iDim(INTEGER(getAttrib(gc_f_m, install("DimCst")))));
    PROTECT(nbh_f = getListElement(Flist, "nbh_f"));        PROTECT(dc_nbh_f = iDim(INTEGER(getAttrib(nbh_f, install("DimCst")))));
    PROTECT(fc_f = getListElement(Flist, "fc_f"));          PROTECT(dc_fc_f = iDim(INTEGER(getAttrib(fc_f, install("DimCst")))));
    PROTECT(fc_f_m = getListElement(Flist, "fc_f_m"));      PROTECT(dc_fc_f_m = iDim(INTEGER(getAttrib(fc_f_m, install("DimCst")))));
    PROTECT(vf_f = getListElement(Flist, "vf_f"));          PROTECT(dc_vf_f = iDim(INTEGER(getAttrib(vf_f, install("DimCst")))));
    PROTECT(vf_f_m = getListElement(Flist, "vf_f_m"));      PROTECT(dc_vf_f_m = iDim(INTEGER(getAttrib(vf_f_m, install("DimCst")))));
    PROTECT(ovc_f = getListElement(Flist, "ovc_f"));        PROTECT(dc_ovc_f = iDim(INTEGER(getAttrib(ovc_f, install("DimCst")))));
    PROTECT(ovc_f_m = getListElement(Flist, "ovc_f_m"));    PROTECT(dc_ovc_f_m = iDim(INTEGER(getAttrib(ovc_f_m, install("DimCst")))));
    PROTECT(oilc_f = getListElement(Flist, "oilc_f"));      PROTECT(dc_oilc_f = iDim(INTEGER(getAttrib(oilc_f, install("DimCst")))));
    PROTECT(oilc_f_m = getListElement(Flist, "oilc_f_m"));  PROTECT(dc_oilc_f_m = iDim(INTEGER(getAttrib(oilc_f_m, install("DimCst")))));
    PROTECT(bc_f = getListElement(Flist, "bc_f"));          PROTECT(dc_bc_f = iDim(INTEGER(getAttrib(bc_f, install("DimCst")))));
    PROTECT(bc_f_m = getListElement(Flist, "bc_f_m"));      PROTECT(dc_bc_f_m = iDim(INTEGER(getAttrib(bc_f_m, install("DimCst")))));
    PROTECT(foc_f = getListElement(Flist, "foc_f"));        PROTECT(dc_foc_f = iDim(INTEGER(getAttrib(foc_f, install("DimCst")))));
    PROTECT(foc_f_m = getListElement(Flist, "foc_f_m"));    PROTECT(dc_foc_f_m = iDim(INTEGER(getAttrib(foc_f_m, install("DimCst")))));
    PROTECT(cnb_f = getListElement(Flist, "cnb_f"));        PROTECT(dc_cnb_f = iDim(INTEGER(getAttrib(cnb_f, install("DimCst")))));
    PROTECT(cnb_f_m = getListElement(Flist, "cnb_f_m"));    PROTECT(dc_cnb_f_m = iDim(INTEGER(getAttrib(cnb_f_m, install("DimCst")))));
    PROTECT(icec_f = getListElement(Flist, "icec_f"));      PROTECT(dc_icec_f = iDim(INTEGER(getAttrib(icec_f, install("DimCst")))));
    PROTECT(icec_f_m = getListElement(Flist, "icec_f_m"));  PROTECT(dc_icec_f_m = iDim(INTEGER(getAttrib(icec_f_m, install("DimCst")))));
    PROTECT(cshr_f = getListElement(Flist, "cshr_f"));      PROTECT(dc_cshr_f = iDim(INTEGER(getAttrib(cshr_f, install("DimCst")))));
    PROTECT(cshr_f_m = getListElement(Flist, "cshr_f_m"));  PROTECT(dc_cshr_f_m = iDim(INTEGER(getAttrib(cshr_f_m, install("DimCst")))));
    PROTECT(eec_f = getListElement(Flist, "eec_f"));        PROTECT(dc_eec_f = iDim(INTEGER(getAttrib(eec_f, install("DimCst")))));
    PROTECT(mwh_f = getListElement(Flist, "mwh_f"));        PROTECT(dc_mwh_f = iDim(INTEGER(getAttrib(mwh_f, install("DimCst")))));
    PROTECT(altwh_f = getListElement(Flist, "altwh_f"));    PROTECT(dc_altwh_f = iDim(INTEGER(getAttrib(altwh_f, install("DimCst")))));
    PROTECT(rep_f = getListElement(Flist, "rep_f"));        PROTECT(dc_rep_f = iDim(INTEGER(getAttrib(rep_f, install("DimCst")))));
    PROTECT(onvc_f = getListElement(Flist, "onvc_f"));      PROTECT(dc_onvc_f = iDim(INTEGER(getAttrib(onvc_f, install("DimCst")))));
    PROTECT(insp_f = getListElement(Flist, "insp_f"));      PROTECT(dc_insp_f = iDim(INTEGER(getAttrib(insp_f, install("DimCst")))));
    PROTECT(ownc_f = getListElement(Flist, "ownc_f"));      PROTECT(dc_ownc_f = iDim(INTEGER(getAttrib(ownc_f, install("DimCst")))));
    PROTECT(mngc_f = getListElement(Flist, "mngc_f"));      PROTECT(dc_mngc_f = iDim(INTEGER(getAttrib(mngc_f, install("DimCst")))));
    PROTECT(licc_f = getListElement(Flist, "licc_f"));      PROTECT(dc_licc_f = iDim(INTEGER(getAttrib(licc_f, install("DimCst")))));
    PROTECT(comc_f = getListElement(Flist, "comc_f"));      PROTECT(dc_comc_f = iDim(INTEGER(getAttrib(comc_f, install("DimCst")))));
    PROTECT(finc_f = getListElement(Flist, "finc_f"));      PROTECT(dc_finc_f = iDim(INTEGER(getAttrib(finc_f, install("DimCst")))));
    PROTECT(dep_f = getListElement(Flist, "dep_f"));        PROTECT(dc_dep_f = iDim(INTEGER(getAttrib(dep_f, install("DimCst")))));
    PROTECT(ic_f = getListElement(Flist, "ic_f"));          PROTECT(dc_ic_f = iDim(INTEGER(getAttrib(ic_f, install("DimCst")))));
    PROTECT(K_f = getListElement(Flist, "K_f"));            PROTECT(dc_K_f = iDim(INTEGER(getAttrib(K_f, install("DimCst")))));
    PROTECT(vc_f = getListElement(Flist, "vc_f"));          PROTECT(dc_vc_f = iDim(INTEGER(getAttrib(vc_f, install("DimCst")))));
    PROTECT(persc_f = getListElement(Flist, "persc_f"));    PROTECT(dc_persc_f = iDim(INTEGER(getAttrib(persc_f, install("DimCst")))));
    PROTECT(ecc_f = getListElement(Flist, "ecc_f"));        PROTECT(dc_ecc_f = iDim(INTEGER(getAttrib(ecc_f, install("DimCst")))));
    PROTECT(pl_f = getListElement(Flist, "pl_f"));          PROTECT(dc_pl_f = iDim(INTEGER(getAttrib(pl_f, install("DimCst")))));
//90 protect  ->104

if (ue_choice == 1) {

    PROTECT(ue_f = getListElement(Flist, "nbds_f"));
    PROTECT(ue_f_m = getListElement(Flist, "nbds_f_m"));

} else {

    if (ue_choice == 2) {

        PROTECT(ue_f = getListElement(Flist, "nbh_f"));
        PROTECT(ue_f_m = getListElement(Flist, "nbh_f_m"));

    } else {

        PROTECT(ue_f = getListElement(Flist, "nbtrip_f"));
        PROTECT(ue_f_m = getListElement(Flist, "nbtrip_f_m"));

    }
}

    PROTECT(dc_ue_f = iDim(INTEGER(getAttrib(ue_f, install("DimCst")))));
    PROTECT(dc_ue_f_m = iDim(INTEGER(getAttrib(ue_f_m, install("DimCst")))));
//4 protect -> 108

    dim_Lref_f = INTEGER(dc_Lref_f);                        r_Lref_f = REAL(Lref_f);
    dim_GVLref_f = INTEGER(dc_GVLref_f);                    r_GVLref_f = REAL(GVLref_f);
    dim_GVLref_f_m = INTEGER(dc_GVLref_f_m);                r_GVLref_f_m = REAL(GVLref_f_m);
    dim_nbv_f = INTEGER(dc_nbv_f);                          r_nbv_f = REAL(nbv_f);
    dim_nbv_f_m = INTEGER(dc_nbv_f_m);                      r_nbv_f_m = REAL(nbv_f_m);
    dim_lc_f = INTEGER(dc_lc_f);                            r_lc_f = REAL(lc_f);
    dim_gc_f = INTEGER(dc_gc_f);                            r_gc_f = REAL(gc_f);
    dim_gc_f_m = INTEGER(dc_gc_f_m);                        r_gc_f_m = REAL(gc_f_m);
    dim_nbh_f = INTEGER(dc_nbh_f);                          r_nbh_f = REAL(nbh_f);
    dim_ue_f = INTEGER(dc_ue_f);                            r_ue_f = REAL(ue_f);
    dim_ue_f_m = INTEGER(dc_ue_f_m);                        r_ue_f_m = REAL(ue_f_m);
    dim_fc_f = INTEGER(dc_fc_f);                            r_fc_f = REAL(fc_f);
    dim_fc_f_m = INTEGER(dc_fc_f_m);                        r_fc_f_m = REAL(fc_f_m);
    dim_vf_f = INTEGER(dc_vf_f);                            r_vf_f = REAL(vf_f);
    dim_vf_f_m = INTEGER(dc_vf_f_m);                        r_vf_f_m = REAL(vf_f_m);
    dim_ovc_f = INTEGER(dc_ovc_f);                          r_ovc_f = REAL(ovc_f);
    dim_ovc_f_m = INTEGER(dc_ovc_f_m);                      r_ovc_f_m = REAL(ovc_f_m);
    dim_oilc_f = INTEGER(dc_oilc_f);                        r_oilc_f = REAL(oilc_f);
    dim_oilc_f_m = INTEGER(dc_oilc_f_m);                    r_oilc_f_m = REAL(oilc_f_m);
    dim_bc_f = INTEGER(dc_bc_f);                            r_bc_f = REAL(bc_f);
    dim_bc_f_m = INTEGER(dc_bc_f_m);                        r_bc_f_m = REAL(bc_f_m);
    dim_foc_f = INTEGER(dc_foc_f);                          r_foc_f = REAL(foc_f);
    dim_foc_f_m = INTEGER(dc_foc_f_m);                      r_foc_f_m = REAL(foc_f_m);
    dim_cnb_f = INTEGER(dc_cnb_f);                          r_cnb_f = REAL(cnb_f);
    dim_cnb_f_m = INTEGER(dc_cnb_f_m);                      r_cnb_f_m = REAL(cnb_f_m);
    dim_icec_f = INTEGER(dc_icec_f);                        r_icec_f = REAL(icec_f);
    dim_icec_f_m = INTEGER(dc_icec_f_m);                    r_icec_f_m = REAL(icec_f_m);
    dim_cshr_f = INTEGER(dc_cshr_f);                        r_cshr_f = REAL(cshr_f);
    dim_cshr_f_m = INTEGER(dc_cshr_f_m);                    r_cshr_f_m = REAL(cshr_f_m);
    dim_eec_f = INTEGER(dc_eec_f);                          r_eec_f = REAL(eec_f);
    dim_mwh_f = INTEGER(dc_mwh_f);                          r_mwh_f = REAL(mwh_f);
    dim_altwh_f = INTEGER(dc_altwh_f);                      r_altwh_f = REAL(altwh_f);
    dim_rep_f = INTEGER(dc_rep_f);                          r_rep_f = REAL(rep_f);
    dim_onvc_f = INTEGER(dc_onvc_f);                        r_onvc_f = REAL(onvc_f);
    dim_insp_f = INTEGER(dc_insp_f);                        r_insp_f = REAL(insp_f);
    dim_ownc_f = INTEGER(dc_ownc_f);                        r_ownc_f = REAL(ownc_f);
    dim_mngc_f = INTEGER(dc_mngc_f);                        r_mngc_f = REAL(mngc_f);
    dim_licc_f = INTEGER(dc_licc_f);                        r_licc_f = REAL(licc_f);
    dim_comc_f = INTEGER(dc_comc_f);                        r_comc_f = REAL(comc_f);
    dim_finc_f = INTEGER(dc_finc_f);                        r_finc_f = REAL(finc_f);
    dim_dep_f = INTEGER(dc_dep_f);                          r_dep_f = REAL(dep_f);
    dim_ic_f = INTEGER(dc_ic_f);                            r_ic_f = REAL(ic_f);
    dim_K_f = INTEGER(dc_K_f);                              r_K_f = REAL(K_f);
    dim_vc_f = INTEGER(dc_vc_f);                            r_vc_f = REAL(vc_f);
    dim_persc_f = INTEGER(dc_persc_f);                      r_persc_f = REAL(persc_f);
    dim_ecc_f = INTEGER(dc_ecc_f);                          r_ecc_f = REAL(ecc_f);
    dim_pl_f = INTEGER(dc_pl_f);                            r_pl_f = REAL(pl_f);

    int nbI, nbC;

if (ind_t==0) {

    SEXP  Lbio_f_out, GVL_f_m_e_out, GVLtot_f_m_out, GVLav_f_m_out, GVLtot_f_out, GVLav_f_out, NGVLav_f_m_out, NGVLav_f_out, vcst_f_m_out,
          vcst_f_out, rtbs_f_m_out, rtbs_f_out,
          cshrT_f_m_out, cshrT_f_out, sshr_f_m_out, sshr_f_out, ncshr_f_out, ocl_f_out, cs_f_out, csTot_f_out, gva_f_out,
          ccw_f_out, ccwCr_f_out, wageg_f_out, wagen_f_out, gcf_f_out,
          ngcf_f_out, gp_f_out, ssTot_f_out, ps_f_out, sts_f_out, ber_f_out, ratio_gva_GVL_f_out, ratio_gcf_GVL_f_out, ratio_fc_GVL_f_out,
          ratio_oilc_GVL_f_out, ratio_bc_GVL_f_out, ratio_foc_GVL_f_out, ratio_icec_GVL_f_out, ratio_gc_GVL_f_out, ratio_vc_GVL_f_out,
          ratio_rep_GVL_f_out, ratio_mngc_GVL_f_out, ratio_licc_GVL_f_out, ratio_fvol_GVL_f_out, ratio_fvol_Lbio_f_out, ratio_fvol_gva_f_out,
          ratio_gcf_gva_f_out, ratio_K_cnb_f_out, ratio_GVL_K_f_out, ratio_gcf_K_f_out, ratio_ngcf_K_f_out, ratio_gp_K_f_out, ratio_GVL_cnb_ue_f_out,
          rtbsAct_f_out, csAct_f_out, gvaAct_f_out, gcfAct_f_out, psAct_f_out, stsAct_f_out;

    SEXP  Loths_f, Lothm_f_e, GVLtot_f_m_e, GVLreftot_f_m_e, GVLreftot_f_m, GVLoths_f_m, GVLothsref_f_m, GVLothsue_f_m, GVLothsrefue_f_m,
    GVLtot_f_e, GVLreftot_f, GVLoths_f, GVLothsue_f, GVLothmet_f, GVLothmetue_f,
    gcue_f_m, fvolue_f_m, ovcue_f_m, oilcue_f_m, bcue_f_m, focue_f_m, focuecnb_f_m, icecue_f_m, gcue_f, fvolue_f, ovcue_f, oilcue_f,
    bcue_f, focue_f, focuecnb_f, icecue_f, onvcr_f, fvol_f, NGVLav_f, NGVLav_f_m, vcst_f_m, vcst_f, vcstOthm_f, rtbs_f,
    rtbs_f_m, rtbsOthm_f, cshrT_f_m, cshrT_f, cshrTothm_f, sshr_f_m, sshr_f, sshrOthm_f, perscr_f, ccwr_f, opersc_f, eco_names;

    double *r_Loths_f, *r_GVLreftot_f_m, *r_GVLoths_f_m, *r_GVLothsref_f_m,
    *r_GVLothsue_f_m, *r_GVLothsrefue_f_m, *r_GVLreftot_f, *r_GVLoths_f, *r_GVLothsue_f, *r_GVLothmet_f, *r_GVLothmetue_f,
    *r_gcue_f_m, *r_fvolue_f_m, *r_ovcue_f_m, *r_oilcue_f_m, *r_bcue_f_m, *r_focue_f_m, *r_focuecnb_f_m, *r_icecue_f_m, *r_gcue_f,
    *r_fvolue_f, *r_ovcue_f, *r_oilcue_f, *r_bcue_f, *r_focue_f, *r_focuecnb_f, *r_icecue_f, *r_onvcr_f, *r_fvol_f, *r_NGVLav_f,
    *r_NGVLav_f_m, *r_vcst_f_m, *r_vcst_f, *r_vcstOthm_f, *r_rtbs_f, *r_rtbs_f_m, *r_rtbsOthm_f, *r_cshrT_f_m, *r_cshrT_f,
    *r_cshrTothm_f, *r_sshr_f_m, *r_sshr_f, *r_sshrOthm_f, *r_perscr_f, *r_ccwr_f, *r_opersc_f;


//-------------------------
// Stade préliminaire (temps initial)
//-------------------------

    PROTECT(Loths_f = NEW_NUMERIC(nbF));                     r_Loths_f = REAL(Loths_f);
    PROTECT(GVLreftot_f_m = NEW_NUMERIC(nbF*nbMe));          r_GVLreftot_f_m = REAL(GVLreftot_f_m);
    PROTECT(GVLoths_f_m = NEW_NUMERIC(nbF*nbMe));            r_GVLoths_f_m = REAL(GVLoths_f_m);
    PROTECT(GVLothsref_f_m = NEW_NUMERIC(nbF*nbMe));         r_GVLothsref_f_m = REAL(GVLothsref_f_m);
    PROTECT(GVLothsue_f_m = NEW_NUMERIC(nbF*nbMe));          r_GVLothsue_f_m = REAL(GVLothsue_f_m);
    PROTECT(GVLothsrefue_f_m = NEW_NUMERIC(nbF*nbMe));       r_GVLothsrefue_f_m = REAL(GVLothsrefue_f_m);
    PROTECT(GVLreftot_f = NEW_NUMERIC(nbF));                 r_GVLreftot_f = REAL(GVLreftot_f);
    PROTECT(GVLoths_f = NEW_NUMERIC(nbF));                   r_GVLoths_f = REAL(GVLoths_f);
    PROTECT(GVLothsue_f = NEW_NUMERIC(nbF));                 r_GVLothsue_f = REAL(GVLothsue_f);
    PROTECT(GVLothmet_f = NEW_NUMERIC(nbF));                 r_GVLothmet_f = REAL(GVLothmet_f);
    PROTECT(GVLothmetue_f = NEW_NUMERIC(nbF));               r_GVLothmetue_f = REAL(GVLothmetue_f);
    PROTECT(gcue_f_m = NEW_NUMERIC(nbF*nbMe));               r_gcue_f_m = REAL(gcue_f_m);
    PROTECT(fvolue_f_m = NEW_NUMERIC(nbF*nbMe));             r_fvolue_f_m = REAL(fvolue_f_m);
    PROTECT(ovcue_f_m = NEW_NUMERIC(nbF*nbMe));              r_ovcue_f_m = REAL(ovcue_f_m);
    PROTECT(oilcue_f_m = NEW_NUMERIC(nbF*nbMe));             r_oilcue_f_m = REAL(oilcue_f_m);
    PROTECT(bcue_f_m = NEW_NUMERIC(nbF*nbMe));               r_bcue_f_m = REAL(bcue_f_m);
    PROTECT(focue_f_m = NEW_NUMERIC(nbF*nbMe));              r_focue_f_m = REAL(focue_f_m);
    PROTECT(focuecnb_f_m = NEW_NUMERIC(nbF*nbMe));           r_focuecnb_f_m = REAL(focuecnb_f_m);
    PROTECT(icecue_f_m = NEW_NUMERIC(nbF*nbMe));             r_icecue_f_m = REAL(icecue_f_m);
    PROTECT(gcue_f = NEW_NUMERIC(nbF));                     r_gcue_f = REAL(gcue_f);
    PROTECT(fvolue_f = NEW_NUMERIC(nbF));                   r_fvolue_f = REAL(fvolue_f);
    PROTECT(ovcue_f = NEW_NUMERIC(nbF));                    r_ovcue_f = REAL(ovcue_f);
    PROTECT(oilcue_f = NEW_NUMERIC(nbF));                   r_oilcue_f = REAL(oilcue_f);
    PROTECT(bcue_f = NEW_NUMERIC(nbF));                     r_bcue_f = REAL(bcue_f);
    PROTECT(focue_f = NEW_NUMERIC(nbF));                    r_focue_f = REAL(focue_f);
    PROTECT(focuecnb_f = NEW_NUMERIC(nbF));                 r_focuecnb_f = REAL(focuecnb_f);
    PROTECT(icecue_f = NEW_NUMERIC(nbF));                   r_icecue_f = REAL(icecue_f);
    PROTECT(onvcr_f = NEW_NUMERIC(nbF));                    r_onvcr_f = REAL(onvcr_f);
    PROTECT(fvol_f = NEW_NUMERIC(nbF));                     r_fvol_f = REAL(fvol_f);
    PROTECT(NGVLav_f = NEW_NUMERIC(nbF));                   r_NGVLav_f = REAL(NGVLav_f);
    PROTECT(NGVLav_f_m = NEW_NUMERIC(nbF*nbMe));             r_NGVLav_f_m = REAL(NGVLav_f_m);
    PROTECT(vcst_f_m = NEW_NUMERIC(nbF*nbMe));               r_vcst_f_m = REAL(vcst_f_m);
    PROTECT(vcst_f = NEW_NUMERIC(nbF));                     r_vcst_f = REAL(vcst_f );
    PROTECT(vcstOthm_f = NEW_NUMERIC(nbF));                 r_vcstOthm_f = REAL(vcstOthm_f);
    PROTECT(rtbs_f = NEW_NUMERIC(nbF));                     r_rtbs_f = REAL(rtbs_f);
    PROTECT(rtbs_f_m = NEW_NUMERIC(nbF*nbMe));               r_rtbs_f_m = REAL(rtbs_f_m);
    PROTECT(rtbsOthm_f = NEW_NUMERIC(nbF));                 r_rtbsOthm_f = REAL(rtbsOthm_f);
    PROTECT(cshrT_f_m = NEW_NUMERIC(nbF*nbMe));              r_cshrT_f_m = REAL(cshrT_f_m);
    PROTECT(cshrT_f = NEW_NUMERIC(nbF));                    r_cshrT_f = REAL(cshrT_f);
    PROTECT(cshrTothm_f = NEW_NUMERIC(nbF));                r_cshrTothm_f = REAL(cshrTothm_f);
    PROTECT(sshr_f_m = NEW_NUMERIC(nbF*nbMe));               r_sshr_f_m = REAL(sshr_f_m);
    PROTECT(sshr_f = NEW_NUMERIC(nbF));                     r_sshr_f = REAL(sshr_f);
    PROTECT(sshrOthm_f = NEW_NUMERIC(nbF));                 r_sshrOthm_f = REAL(sshrOthm_f);
    PROTECT(perscr_f = NEW_NUMERIC(nbF));                   r_perscr_f = REAL(perscr_f);
    PROTECT(ccwr_f = NEW_NUMERIC(nbF));                   r_ccwr_f = REAL(ccwr_f);
    PROTECT(opersc_f = NEW_NUMERIC(nbF));                   r_opersc_f = REAL(opersc_f);
//46 protect  (+2)  -> 46 (t0)


    for (int e = 0 ; e < nbE ; e++) {

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        nbC = length(getListElement(elmt, "modC"));

        PROTECT(Lothm_f_e = NEW_NUMERIC(nbF));
        PROTECT(GVLtot_f_m_e = NEW_NUMERIC(nbF*nbMe*nbT));
        PROTECT(GVLreftot_f_m_e = NEW_NUMERIC(nbF*nbMe));
        PROTECT(GVLtot_f_e = NEW_NUMERIC(nbF));

        double *r_Lothm_f_e = REAL(Lothm_f_e);
        double *r_GVLtot_f_m_e = REAL(GVLtot_f_m_e);
        double *r_GVLreftot_f_m_e = REAL(GVLreftot_f_m_e);
        double *r_GVLtot_f_e = REAL(GVLtot_f_e);
        double *r_Lref_f_e = REAL(getListElement(elmt, "Lref_f_e"));
        double *r_Lref_f_sum_e = REAL(aggregObj(getListElement(elmt, "Lref_f_m_e"),dimCstFini));
        double *r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct2, e));
        double *r_Lbio_f_sum_e = REAL(aggregObj(VECTOR_ELT(out_L_efmct2, e),dimCstF));
        double *r_GVLref_f_m_e = REAL(getListElement(elmt, "GVLref_f_m_e"));
        double *r_GVLref_f_e = REAL(getListElement(elmt, "GVLref_f_e"));
        double *r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
        int *dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct2, e), install("DimCst")))));
        int *dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));


        //------------------------------
        //équations de la table "p"
        //------------------------------

        //-- 1. Loths_f

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){   //on rappelle ici que ind_t est en fait égal à 0

            if (e==0) {

                r_Loths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_Lref_f[ind_f*dim_Lref_f[0] + 0*dim_Lref_f[1] + 0*dim_Lref_f[2] + ind_t*dim_Lref_f[3]] -
                    r_Lref_f_e[ind_f*dim_Lref_f[0] + 0*dim_Lref_f[1] + 0*dim_Lref_f[2] + ind_t*dim_Lref_f[3]]; //on suppose que la donnée par espèce est de même dimension que la donnée totale

            } else {

                r_Loths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_Loths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_Lref_f_e[ind_f*dim_Lref_f[0] + 0*dim_Lref_f[1] + 0*dim_Lref_f[2] + ind_t*dim_Lref_f[3]];

            }

        //-- 4. Lothm_f_e

    if (adj==2) { //alors on cale la variable sur les valeurs issues des bases activité

        r_Lothm_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_Lref_f_e[ind_f*dim_Lref_f[0] + 0*dim_Lref_f[1] + 0*dim_Lref_f[2] + ind_t*dim_Lref_f[3]] -
                r_Lref_f_sum_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];

    } else {        //autrement, on se rapporte à la valeur estimée à t initial

        r_Lothm_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_Lref_f_e[ind_f*dim_Lref_f[0] + 0*dim_Lref_f[1] + 0*dim_Lref_f[2] + ind_t*dim_Lref_f[3]] -
                r_Lbio_f_sum_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];

    }


    double countGVLtotf = 0.0; //pour sommer GVLtot_f_m_e sur les métiers

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

         //-- 5. GVLtot_f_m_e

             double count = 0.0;

             for (int ind_c = 0 ; ind_c < nbC ; ind_c++){

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]))

                count = count +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] ;

             }

            countGVLtotf = countGVLtotf + count;
            r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = count;


         //-- 6. GVLreftot_f_m_e

                r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLref_f_m_e[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                  r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]; //rappel : ind_t = 0 ici




        //-- 7. GVLreftot_f_m

            if (e==0) {

                r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLref_f_m[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                  r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]; //rappel : ind_t = 0 ici

            }

         //-- 8. GVLoths_f_m

            if (e==0) {

                if (!ISNA(r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }

            } else {

                if (!ISNA(r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }



        //-- 9. GVLothsref_f_m

            if (e==0) {

                if (!ISNA(r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }

            } else {

                if (!ISNA(r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }

            }

    //-- 12. GVLtot_f_e

    if (lev==2 & adj==1) {

            r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = countGVLtotf;

    } else {

            r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_GVLref_f_e[ind_f*dim_GVLref_f[0] + 0*dim_GVLref_f[1] + 0*dim_GVLref_f[2] + ind_t*dim_GVLref_f[3]] *
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];
    }


    //-- 13. GVLreftot_f

    if (e==0) {

            r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLref_f[ind_f*dim_GVLref_f[0] + 0*dim_GVLref_f[1] + 0*dim_GVLref_f[2] + ind_t*dim_GVLref_f[3]] *
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

    }

    //-- 14. GVLoths_f

            if (e==0) {

                if (!ISNA(r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }

            } else {

                if (!ISNA(r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_e[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }



        }

        //on formatte le(s) résultat(s) et on les intègre à 'eVar'

        setAttrib(Lothm_f_e, R_NamesSymbol, fleetList);
        setAttrib(Lothm_f_e, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 40, Lothm_f_e);
        setAttrib(GVLtot_f_m_e, R_DimSymbol, DimFM);
        setAttrib(GVLtot_f_m_e, R_DimNamesSymbol, dimnamesFM);
        setAttrib(GVLtot_f_m_e, install("DimCst"), dimCstFM);
        SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 41, GVLtot_f_m_e);

        UNPROTECT(5);

}

    // à ce stade, plus de considération d'espèce pour les variables à initialiser

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

 double countEff = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){


        //-- 16. GVLothmet_f

            if (lev==2) {

                r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

            } else {

                if (ind_m==0) {

                    if (!ISNA(r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                    } else {

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                    }

                } else {

                    if (!ISNA(r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }
            }

           //-- terme : SUM_m (ue_f_m * nbv_f_m)

            countEff = countEff +
                r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

            }

        //-- 15. GVLothsue_f

        if (report==0) {

                r_GVLothsue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    (r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) );


        } else {

                r_GVLothsue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLoths_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    ((r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) - countEff) );

        }


        //-- 17. GVLothmetue_f

        if (lev==2) {

                r_GVLothmetue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

        } else {

                r_GVLothmetue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLothmet_f[ind_f*eF_fm[0] + 0*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    ((r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) - countEff) );

        }


    for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 10. GVLothsue_f_m

                r_GVLothsue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    (r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) );


        //-- 11. GVLothsrefue_f_m

                r_GVLothsrefue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    (r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) );





        //-- 18. gcue_f_m

                r_gcue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_gc_f_m[ind_f*dim_gc_f_m[0] + ind_m*dim_gc_f_m[1] + 0*dim_gc_f_m[2] + ind_t*dim_gc_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        //-- 19. fvolue_f_m

                r_fvolue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]] /
                    (r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]] *
                     r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]]) );

        //-- 20. ovcue_f_m

                r_ovcue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_ovc_f_m[ind_f*dim_ovc_f_m[0] + ind_m*dim_ovc_f_m[1] + 0*dim_ovc_f_m[2] + ind_t*dim_ovc_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        //-- 21. oilcue_f_m

                r_oilcue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_oilc_f_m[ind_f*dim_oilc_f_m[0] + ind_m*dim_oilc_f_m[1] + 0*dim_oilc_f_m[2] + ind_t*dim_oilc_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        //-- 22. bcue_f_m

                r_bcue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_bc_f_m[ind_f*dim_bc_f_m[0] + ind_m*dim_bc_f_m[1] + 0*dim_bc_f_m[2] + ind_t*dim_bc_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        //-- 23. focue_f_m

                r_focue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite(r_foc_f_m[ind_f*dim_foc_f_m[0] + ind_m*dim_foc_f_m[1] + 0*dim_foc_f_m[2] + ind_t*dim_foc_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        //-- 24. focuecnb_f_m

                r_focuecnb_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_focue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    r_cnb_f_m[ind_f*dim_cnb_f_m[0] + ind_m*dim_cnb_f_m[1] + 0*dim_cnb_f_m[2] + ind_t*dim_cnb_f_m[3]] );

        //-- 25. icecue_f_m

                r_icecue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_icec_f_m[ind_f*dim_icec_f_m[0] + ind_m*dim_icec_f_m[1] + 0*dim_icec_f_m[2] + ind_t*dim_icec_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

        }



        //-- 26. gcue_f

                r_gcue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 27. fvolue_f

                r_fvolue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_fc_f[ind_f*dim_fc_f_m[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]] /
                    (r_vf_f[ind_f*dim_vf_f_m[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] *
                     r_ue_f[ind_f*dim_ue_f_m[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) );

        //-- 28. ovcue_f

                r_ovcue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_ovc_f[ind_f*dim_ovc_f[0] + 0*dim_ovc_f[1] + 0*dim_ovc_f[2] + ind_t*dim_ovc_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 29. oilcue_f

                r_oilcue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_oilc_f[ind_f*dim_oilc_f[0] + 0*dim_oilc_f[1] + 0*dim_oilc_f[2] + ind_t*dim_oilc_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 30. bcue_f

                r_bcue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_bc_f[ind_f*dim_bc_f[0] + 0*dim_bc_f[1] + 0*dim_bc_f[2] + ind_t*dim_bc_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 31. focue_f

                r_focue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_foc_f[ind_f*dim_foc_f[0] + 0*dim_foc_f[1] + 0*dim_foc_f[2] + ind_t*dim_foc_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 32. focuecnb_f

                r_focuecnb_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_focue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] );

        //-- 33. icecue_f

                r_icecue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_icec_f[ind_f*dim_icec_f[0] + 0*dim_icec_f[1] + 0*dim_icec_f[2] + ind_t*dim_icec_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );

        //-- 34. onvcr_f

                r_onvcr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_onvc_f[ind_f*dim_onvc_f[0] + 0*dim_onvc_f[1] + 0*dim_onvc_f[2] + ind_t*dim_onvc_f[3]] -
                    r_insp_f[ind_f*dim_insp_f[0] + 0*dim_insp_f[1] + 0*dim_insp_f[2] + ind_t*dim_insp_f[3]] -
                    r_ownc_f[ind_f*dim_ownc_f[0] + 0*dim_ownc_f[1] + 0*dim_ownc_f[2] + ind_t*dim_ownc_f[3]] -
                    r_mngc_f[ind_f*dim_mngc_f[0] + 0*dim_mngc_f[1] + 0*dim_mngc_f[2] + ind_t*dim_mngc_f[3]] -
                    r_licc_f[ind_f*dim_licc_f[0] + 0*dim_licc_f[1] + 0*dim_licc_f[2] + ind_t*dim_licc_f[3]] -
                    r_comc_f[ind_f*dim_comc_f[0] + 0*dim_comc_f[1] + 0*dim_comc_f[2] + ind_t*dim_comc_f[3]];

        //-- 35. fvol_f

            r_fvol_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                finite( r_fc_f[ind_f*dim_fc_f[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]] /
                r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] );

        //-- 36. NGVLav_f

            r_NGVLav_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]])/
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]] ;


       //-- 39. vcst_f

                if (!all_is_na(oilc_f) & !(all_is_na(bc_f)) & !(all_is_na(foc_f)) & !(all_is_na(icec_f))) {

                    r_vcst_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_fc_f[ind_f*dim_fc_f[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]] +
                        r_oilc_f[ind_f*dim_oilc_f[0] + 0*dim_oilc_f[1] + 0*dim_oilc_f[2] + ind_t*dim_oilc_f[3]] +
                        r_bc_f[ind_f*dim_bc_f[0] + 0*dim_bc_f[1] + 0*dim_bc_f[2] + ind_t*dim_bc_f[3]] +
                        r_foc_f[ind_f*dim_foc_f[0] + 0*dim_foc_f[1] + 0*dim_foc_f[2] + ind_t*dim_foc_f[3]] +
                        r_icec_f[ind_f*dim_icec_f[0] + 0*dim_icec_f[1] + 0*dim_icec_f[2] + ind_t*dim_icec_f[3]];

                } else {

                    r_vcst_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_ovc_f[ind_f*dim_ovc_f[0] + 0*dim_ovc_f[1] + 0*dim_ovc_f[2] + ind_t*dim_ovc_f[3]] +
                        r_fc_f[ind_f*dim_fc_f[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]];

                }


        //-- 41. rtbs_f

            r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_NGVLav_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                r_vcst_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 45. cshrT_f

            r_cshrT_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 48. sshr_f

            r_sshr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                r_cshrT_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 37. NGVLav_f_m

            r_NGVLav_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]]) /
                r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

        //-- 38. vcst_f_m

                if (!all_is_na(oilc_f_m) & !(all_is_na(bc_f_m)) & !(all_is_na(foc_f_m)) & !(all_is_na(icec_f_m))) {

                    r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]] +
                        r_oilc_f_m[ind_f*dim_oilc_f_m[0] + ind_m*dim_oilc_f_m[1] + 0*dim_oilc_f_m[2] + ind_t*dim_oilc_f_m[3]] +
                        r_bc_f_m[ind_f*dim_bc_f_m[0] + ind_m*dim_bc_f_m[1] + 0*dim_bc_f_m[2] + ind_t*dim_bc_f_m[3]] +
                        r_foc_f_m[ind_f*dim_foc_f_m[0] + ind_m*dim_foc_f_m[1] + 0*dim_foc_f_m[2] + ind_t*dim_foc_f_m[3]] +
                        r_icec_f_m[ind_f*dim_icec_f_m[0] + ind_m*dim_icec_f_m[1] + 0*dim_icec_f_m[2] + ind_t*dim_icec_f_m[3]];

                } else {

                    r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        r_ovc_f_m[ind_f*dim_ovc_f_m[0] + ind_m*dim_ovc_f_m[1] + 0*dim_ovc_f_m[2] + ind_t*dim_ovc_f_m[3]] +
                        r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]];

                }

        //-- 40. vcstOthm_f

                if (ind_m==0) {

                    if (!ISNA(r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])) {

                    r_vcstOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_vcst_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) -
                        (r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                    } else {

                    r_vcstOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_vcst_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                    }

                } else {

                if (!ISNA(r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_vcstOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_vcstOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        (r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                }

        //-- 42. rtbs_f_m

                r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_NGVLav_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_vcst_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

        //-- 43. rtbsOthm_f

                if (ind_m==0) {

                    if (!ISNA(r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])) {


                    r_rtbsOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) -
                        (r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                    } else {

                    r_rtbsOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                    }

                } else {

                    if (!ISNA(r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_rtbsOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_rtbsOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        (r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                }


        //-- 44. cshrT_f_m

                r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    0.01*r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                    r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

        //-- 46. cshrTothm_f

                if (ind_m==0) {

                    if (!ISNA(r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])) {

                    r_cshrTothm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_cshrT_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) -
                        (r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                    } else {

                    r_cshrTothm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_cshrT_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                    }

                } else {

                    if (!ISNA(r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_cshrTothm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_cshrTothm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        (r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                }


        //-- 47. sshr_f_m

                r_sshr_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_rtbs_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_cshrT_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

        //-- 49. sshrOthm_f

                if (ind_m==0) {

                    if (!ISNA(r_sshr_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])) {

                    r_sshrOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_sshr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) -
                        (r_sshr_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                    } else {

                    r_sshrOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_sshr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                         r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                    }

                } else {

                    if (!ISNA(r_sshr_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_sshrOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_sshrOthm_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        (r_sshr_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]);

                }

}

        //-- 50. perscr_f

            r_perscr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] -
                r_eec_f[ind_f*dim_eec_f[0] + 0*dim_eec_f[1] + 0*dim_eec_f[2] + ind_t*dim_eec_f[3]] -
                r_ecc_f[ind_f*dim_ecc_f[0] + 0*dim_ecc_f[1] + 0*dim_ecc_f[2] + ind_t*dim_ecc_f[3]] -
                r_pl_f[ind_f*dim_pl_f[0] + 0*dim_pl_f[1] + 0*dim_pl_f[2] + ind_t*dim_pl_f[3]];



        //-- 51. ccwr_f

            r_ccwr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] /
                r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 52. opersc_f

            r_opersc_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] -
                r_cshrT_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

}



//on formatte le(s) résultat(s) et on intègre à fVar

        setAttrib(Loths_f, R_NamesSymbol, fleetList);
        setAttrib(Loths_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 0, Loths_f);

        setAttrib(GVLoths_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLoths_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLoths_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 1, GVLoths_f_m);

        setAttrib(GVLothsref_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsref_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsref_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 2, GVLothsref_f_m);

        setAttrib(GVLothsue_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 23, GVLothsue_f_m);

        setAttrib(GVLothsrefue_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsrefue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsrefue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 24, GVLothsrefue_f_m);

        setAttrib(GVLothmet_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothmet_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 3, GVLothmet_f);

        setAttrib(GVLothmetue_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothmetue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 25, GVLothmetue_f);

        setAttrib(GVLothsue_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothsue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 26, GVLothsue_f);

        setAttrib(fvolue_f_m, R_DimSymbol, DimFMini);
        setAttrib(fvolue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(fvolue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 4, fvolue_f_m);

        setAttrib(oilcue_f_m, R_DimSymbol, DimFMini);
        setAttrib(oilcue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(oilcue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 5, oilcue_f_m);

        setAttrib(bcue_f_m, R_DimSymbol, DimFMini);
        setAttrib(bcue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(bcue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 6, bcue_f_m);

        setAttrib(focuecnb_f_m, R_DimSymbol, DimFMini);
        setAttrib(focuecnb_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(focuecnb_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 7, focuecnb_f_m);

        setAttrib(icecue_f_m, R_DimSymbol, DimFMini);
        setAttrib(icecue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(icecue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 8, icecue_f_m);

        setAttrib(gcue_f_m, R_DimSymbol, DimFMini);
        setAttrib(gcue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(gcue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 9, gcue_f_m);

        setAttrib(ovcue_f_m, R_DimSymbol, DimFMini);
        setAttrib(ovcue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(ovcue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 10, ovcue_f_m);

        setAttrib(vcstOthm_f, R_NamesSymbol, fleetList);
        setAttrib(vcstOthm_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 11, vcstOthm_f);

        setAttrib(rtbsOthm_f, R_NamesSymbol, fleetList);
        setAttrib(rtbsOthm_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 12, rtbsOthm_f);

        setAttrib(cshrTothm_f, R_NamesSymbol, fleetList);
        setAttrib(cshrTothm_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 13, cshrTothm_f);

        setAttrib(sshrOthm_f, R_NamesSymbol, fleetList);
        setAttrib(sshrOthm_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 14, sshrOthm_f);

        setAttrib(fvol_f, R_NamesSymbol, fleetList);
        setAttrib(fvol_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 15, fvol_f);

        setAttrib(gcue_f, R_NamesSymbol, fleetList);
        setAttrib(gcue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 16, gcue_f);

        setAttrib(fvolue_f, R_NamesSymbol, fleetList);
        setAttrib(fvolue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 17, fvolue_f);

        setAttrib(ovcue_f, R_NamesSymbol, fleetList);
        setAttrib(ovcue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 18, ovcue_f);

        setAttrib(oilcue_f, R_NamesSymbol, fleetList);
        setAttrib(oilcue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 19, oilcue_f);

        setAttrib(bcue_f, R_NamesSymbol, fleetList);
        setAttrib(bcue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 20, bcue_f);

        setAttrib(focuecnb_f, R_NamesSymbol, fleetList);
        setAttrib(focuecnb_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 21, focuecnb_f);

        setAttrib(icecue_f, R_NamesSymbol, fleetList);
        setAttrib(icecue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 22, icecue_f);

        setAttrib(ccwr_f, R_NamesSymbol, fleetList);
        setAttrib(ccwr_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 27, ccwr_f);

        setAttrib(opersc_f, R_NamesSymbol, fleetList);
        setAttrib(opersc_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 28, opersc_f);

        setAttrib(GVLoths_f, R_NamesSymbol, fleetList);
        setAttrib(GVLoths_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 29, GVLoths_f);

        setAttrib(onvcr_f, R_NamesSymbol, fleetList);
        setAttrib(onvcr_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 30, onvcr_f);

        SET_VECTOR_ELT(fVar, 31, rtbs_f);
        SET_VECTOR_ELT(fVar, 32, rtbs_f_m);


//enfin, on initialise l'output


    PROTECT(Lbio_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(Lbio_f_out, R_DimSymbol, DimF);
    setAttrib(Lbio_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(Lbio_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 0, Lbio_f_out);

    PROTECT(GVL_f_m_e_out = allocVector(VECSXP, nbE));
    setAttrib(GVL_f_m_e_out, R_NamesSymbol, sppList);
    SET_VECTOR_ELT(out_Eco, 1, GVL_f_m_e_out);

    PROTECT(GVLtot_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLtot_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLtot_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLtot_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 2, GVLtot_f_m_out);

    PROTECT(GVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 3, GVLav_f_m_out);

    PROTECT(GVLtot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLtot_f_out, R_DimSymbol, DimF);
    setAttrib(GVLtot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLtot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 4, GVLtot_f_out);

    PROTECT(GVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLav_f_out, R_DimSymbol, DimF);
    setAttrib(GVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 5, GVLav_f_out);

    PROTECT(NGVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(NGVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(NGVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(NGVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 6, NGVLav_f_m_out);

    PROTECT(NGVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(NGVLav_f_out, R_DimSymbol, DimF);
    setAttrib(NGVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(NGVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 7, NGVLav_f_out);

    PROTECT(vcst_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(vcst_f_m_out, R_DimSymbol, DimFM);
    setAttrib(vcst_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(vcst_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 8, vcst_f_m_out);

    PROTECT(vcst_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(vcst_f_out, R_DimSymbol, DimF);
    setAttrib(vcst_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(vcst_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 9, vcst_f_out);

    PROTECT(rtbs_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(rtbs_f_m_out, R_DimSymbol, DimFM);
    setAttrib(rtbs_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(rtbs_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 10, rtbs_f_m_out);

    PROTECT(rtbs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbs_f_out, R_DimSymbol, DimF);
    setAttrib(rtbs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 11, rtbs_f_out);

    PROTECT(cshrT_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(cshrT_f_m_out, R_DimSymbol, DimFM);
    setAttrib(cshrT_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(cshrT_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 12, cshrT_f_m_out);

    PROTECT(cshrT_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cshrT_f_out, R_DimSymbol, DimF);
    setAttrib(cshrT_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cshrT_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 13, cshrT_f_out);

    PROTECT(sshr_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(sshr_f_m_out, R_DimSymbol, DimFM);
    setAttrib(sshr_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(sshr_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_Eco, 14, sshr_f_m_out);

    PROTECT(sshr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(sshr_f_out, R_DimSymbol, DimF);
    setAttrib(sshr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(sshr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 15, sshr_f_out);

    PROTECT(ncshr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ncshr_f_out, R_DimSymbol, DimF);
    setAttrib(ncshr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ncshr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 16, ncshr_f_out);

    PROTECT(ocl_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ocl_f_out, R_DimSymbol, DimF);
    setAttrib(ocl_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ocl_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 17, ocl_f_out);

    PROTECT(cs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cs_f_out, R_DimSymbol, DimF);
    setAttrib(cs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 18, cs_f_out);

    PROTECT(csTot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csTot_f_out, R_DimSymbol, DimF);
    setAttrib(csTot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csTot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 19, csTot_f_out);

    PROTECT(gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gva_f_out, R_DimSymbol, DimF);
    setAttrib(gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 20, gva_f_out);

    PROTECT(ccw_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccw_f_out, R_DimSymbol, DimF);
    setAttrib(ccw_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccw_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 21, ccw_f_out);

    PROTECT(ccwCr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccwCr_f_out, R_DimSymbol, DimF);
    setAttrib(ccwCr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccwCr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 22, ccwCr_f_out);

    PROTECT(wageg_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wageg_f_out, R_DimSymbol, DimF);
    setAttrib(wageg_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wageg_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 23, wageg_f_out);

    PROTECT(wagen_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wagen_f_out, R_DimSymbol, DimF);
    setAttrib(wagen_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wagen_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 24, wagen_f_out);

    PROTECT(gcf_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gcf_f_out, R_DimSymbol, DimF);
    setAttrib(gcf_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gcf_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 25, gcf_f_out);

    PROTECT(ngcf_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ngcf_f_out, R_DimSymbol, DimF);
    setAttrib(ngcf_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ngcf_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 26, ngcf_f_out);

    PROTECT(gp_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gp_f_out, R_DimSymbol, DimF);
    setAttrib(gp_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gp_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 27, gp_f_out);

    PROTECT(ssTot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ssTot_f_out, R_DimSymbol, DimF);
    setAttrib(ssTot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ssTot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 28, ssTot_f_out);

    PROTECT(ps_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ps_f_out, R_DimSymbol, DimF);
    setAttrib(ps_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ps_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 29, ps_f_out);

    PROTECT(sts_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(sts_f_out, R_DimSymbol, DimF);
    setAttrib(sts_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(sts_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 30, sts_f_out);

    PROTECT(ber_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ber_f_out, R_DimSymbol, DimF);
    setAttrib(ber_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ber_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 31, ber_f_out);

    PROTECT(ratio_gva_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gva_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gva_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gva_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 32, ratio_gva_GVL_f_out);

    PROTECT(ratio_gcf_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 33, ratio_gcf_GVL_f_out);

    PROTECT(ratio_fc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 34, ratio_fc_GVL_f_out);

    PROTECT(ratio_oilc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_oilc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_oilc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_oilc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 35, ratio_oilc_GVL_f_out);

    PROTECT(ratio_bc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_bc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_bc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_bc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 36, ratio_bc_GVL_f_out);

    PROTECT(ratio_foc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_foc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_foc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_foc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 37, ratio_foc_GVL_f_out);

    PROTECT(ratio_icec_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_icec_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_icec_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_icec_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 38, ratio_icec_GVL_f_out);

    PROTECT(ratio_gc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 39, ratio_gc_GVL_f_out);

    PROTECT(ratio_vc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_vc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_vc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_vc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 40, ratio_vc_GVL_f_out);

    PROTECT(ratio_rep_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_rep_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_rep_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_rep_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 41, ratio_rep_GVL_f_out);

    PROTECT(ratio_mngc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_mngc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_mngc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_mngc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 42, ratio_mngc_GVL_f_out);

    PROTECT(ratio_licc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_licc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_licc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_licc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 43, ratio_licc_GVL_f_out);

    PROTECT(ratio_fvol_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 44, ratio_fvol_GVL_f_out);

    PROTECT(ratio_fvol_Lbio_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_Lbio_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_Lbio_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_Lbio_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 45, ratio_fvol_Lbio_f_out);

    PROTECT(ratio_fvol_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 46, ratio_fvol_gva_f_out);

    PROTECT(ratio_gcf_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 47, ratio_gcf_gva_f_out);

    PROTECT(ratio_K_cnb_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_K_cnb_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_K_cnb_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_K_cnb_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 48, ratio_K_cnb_f_out);

    PROTECT(ratio_GVL_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 49, ratio_GVL_K_f_out);

    PROTECT(ratio_gcf_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 50, ratio_gcf_K_f_out);

    PROTECT(ratio_ngcf_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_ngcf_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_ngcf_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_ngcf_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 51, ratio_ngcf_K_f_out);

    PROTECT(ratio_gp_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gp_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gp_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gp_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 52, ratio_gp_K_f_out);

    PROTECT(ratio_GVL_cnb_ue_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_cnb_ue_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 53, ratio_GVL_cnb_ue_f_out);

    PROTECT(rtbsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbsAct_f_out, R_DimSymbol, DimF);
    setAttrib(rtbsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 54, rtbsAct_f_out);

    PROTECT(csAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csAct_f_out, R_DimSymbol, DimF);
    setAttrib(csAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 55, csAct_f_out);

    PROTECT(gvaAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gvaAct_f_out, R_DimSymbol, DimF);
    setAttrib(gvaAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gvaAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 56, gvaAct_f_out);

    PROTECT(gcfAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gcfAct_f_out, R_DimSymbol, DimF);
    setAttrib(gcfAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gcfAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 57, gcfAct_f_out);

    PROTECT(psAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(psAct_f_out, R_DimSymbol, DimF);
    setAttrib(psAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(psAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 58, psAct_f_out);

    PROTECT(stsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(stsAct_f_out, R_DimSymbol, DimF);
    setAttrib(stsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(stsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_Eco, 59, stsAct_f_out);


    //on nomme les éléments de out_Eco
    const char *namesEco[60] = {"Lbio_f","GVL_f_m_e","GVLtot_f_m","GVLav_f_m","GVLtot_f","GVLav_f","NGVLav_f_m","NGVLav_f","vcst_f_m",
                          "vcst_f","rtbs_f_m","rtbs_f","cshrT_f_m","cshrT_f","sshr_f_m","sshr_f","ncshr_f","ocl_f","cs_f","csTot_f",
                          "gva_f","ccw_f","ccwCr_f","wageg_f","wagen_f",
                          "gcf_f","ngcf_f","gp_f","ssTot_f","ps_f","sts_f","ber_f","ratio_gva_GVL_f","ratio_gcf_GVL_f","ratio_fc_GVL_f",
                          "ratio_oilc_GVL_f","ratio_bc_GVL_f","ratio_foc_GVL_f","ratio_icec_GVL_f","ratio_gc_GVL_f","ratio_vc_GVL_f",
                          "ratio_rep_GVL_f","ratio_mngc_GVL_f","ratio_licc_GVL_f","ratio_fvol_GVL_f","ratio_fvol_Lbio_f",
                          "ratio_fvol_gva_f","ratio_gcf_gva_f","ratio_K_cnb_f","ratio_GVL_K_f","ratio_gcf_K_f","ratio_ngcf_K_f",
                          "ratio_gp_K_f","ratio_GVL_cnb_ue_f","rtbsAct_f","csAct_f","gvaAct_f","gcfAct_f","psAct_f",
                          "stsAct_f"};

    PROTECT(eco_names = allocVector(STRSXP, 60));

    for(int ct = 0; ct < 60; ct++) SET_STRING_ELT(eco_names, ct, mkChar(namesEco[ct]));

    setAttrib(out_Eco, R_NamesSymbol, eco_names);

//55 protect    --> 101
}

//on importe les outputs afin de les mettre à jour à l'instant ind_t

    r_Lbio_f_out = REAL(VECTOR_ELT(out_Eco,0));
    r_GVLtot_f_m_out = REAL(VECTOR_ELT(out_Eco,2));
    r_GVLav_f_m_out = REAL(VECTOR_ELT(out_Eco,3));
    r_GVLtot_f_out = REAL(VECTOR_ELT(out_Eco,4));
    r_GVLav_f_out = REAL(VECTOR_ELT(out_Eco,5));
    r_NGVLav_f_m_out = REAL(VECTOR_ELT(out_Eco,6));
    r_NGVLav_f_out = REAL(VECTOR_ELT(out_Eco,7));
    r_vcst_f_m_out = REAL(VECTOR_ELT(out_Eco,8));
    r_vcst_f_out = REAL(VECTOR_ELT(out_Eco,9));
    r_rtbs_f_m_out = REAL(VECTOR_ELT(out_Eco,10));
    r_rtbs_f_out = REAL(VECTOR_ELT(out_Eco,11));
    r_cshrT_f_m_out = REAL(VECTOR_ELT(out_Eco,12));
    r_cshrT_f_out = REAL(VECTOR_ELT(out_Eco,13));
    r_sshr_f_m_out = REAL(VECTOR_ELT(out_Eco,14));
    r_sshr_f_out = REAL(VECTOR_ELT(out_Eco,15));
    r_ncshr_f_out = REAL(VECTOR_ELT(out_Eco,16));
    r_ocl_f_out = REAL(VECTOR_ELT(out_Eco,17));
    r_cs_f_out = REAL(VECTOR_ELT(out_Eco,18));
    r_csTot_f_out = REAL(VECTOR_ELT(out_Eco,19));
    r_gva_f_out = REAL(VECTOR_ELT(out_Eco,20));
    r_ccw_f_out = REAL(VECTOR_ELT(out_Eco,21));
    r_ccwCr_f_out = REAL(VECTOR_ELT(out_Eco,22));
    r_wageg_f_out = REAL(VECTOR_ELT(out_Eco,23));
    r_wagen_f_out = REAL(VECTOR_ELT(out_Eco,24));
    r_gcf_f_out = REAL(VECTOR_ELT(out_Eco,25));
    r_ngcf_f_out = REAL(VECTOR_ELT(out_Eco,26));
    r_gp_f_out = REAL(VECTOR_ELT(out_Eco,27));
    r_ssTot_f_out = REAL(VECTOR_ELT(out_Eco,28));
    r_ps_f_out = REAL(VECTOR_ELT(out_Eco,29));
    r_sts_f_out = REAL(VECTOR_ELT(out_Eco,30));
    r_ber_f_out = REAL(VECTOR_ELT(out_Eco,31));
    r_ratio_gva_GVL_f_out = REAL(VECTOR_ELT(out_Eco,32));
    r_ratio_gcf_GVL_f_out = REAL(VECTOR_ELT(out_Eco,33));
    r_ratio_fc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,34));
    r_ratio_oilc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,35));
    r_ratio_bc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,36));
    r_ratio_foc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,37));
    r_ratio_icec_GVL_f_out = REAL(VECTOR_ELT(out_Eco,38));
    r_ratio_gc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,39));
    r_ratio_vc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,40));
    r_ratio_rep_GVL_f_out = REAL(VECTOR_ELT(out_Eco,41));
    r_ratio_mngc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,42));
    r_ratio_licc_GVL_f_out = REAL(VECTOR_ELT(out_Eco,43));
    r_ratio_fvol_GVL_f_out = REAL(VECTOR_ELT(out_Eco,44));
    r_ratio_fvol_Lbio_f_out = REAL(VECTOR_ELT(out_Eco,45));
    r_ratio_fvol_gva_f_out = REAL(VECTOR_ELT(out_Eco,46));
    r_ratio_gcf_gva_f_out = REAL(VECTOR_ELT(out_Eco,47));
    r_ratio_K_cnb_f_out = REAL(VECTOR_ELT(out_Eco,48));
    r_ratio_GVL_K_f_out = REAL(VECTOR_ELT(out_Eco,49));
    r_ratio_gcf_K_f_out = REAL(VECTOR_ELT(out_Eco,50));
    r_ratio_ngcf_K_f_out = REAL(VECTOR_ELT(out_Eco,51));
    r_ratio_gp_K_f_out = REAL(VECTOR_ELT(out_Eco,52));
    r_ratio_GVL_cnb_ue_f_out = REAL(VECTOR_ELT(out_Eco,53));
    r_rtbsAct_f_out = REAL(VECTOR_ELT(out_Eco,54));
    r_csAct_f_out = REAL(VECTOR_ELT(out_Eco,55));
    r_gvaAct_f_out = REAL(VECTOR_ELT(out_Eco,56));
    r_gcfAct_f_out = REAL(VECTOR_ELT(out_Eco,57));
    r_psAct_f_out = REAL(VECTOR_ELT(out_Eco,58));
    r_stsAct_f_out = REAL(VECTOR_ELT(out_Eco,59));


    double *r_Loths_f2 = REAL(VECTOR_ELT(fVar,0));
    double *r_GVLoths_f_m2 = REAL(VECTOR_ELT(fVar,1));
    double *r_GVLothsref_f_m2 = REAL(VECTOR_ELT(fVar,2));
    double *r_GVLothmet_f2 = REAL(VECTOR_ELT(fVar,3));
    double *r_fvolue_f_m2 = REAL(VECTOR_ELT(fVar,4));
    double *r_oilcue_f_m2 = REAL(VECTOR_ELT(fVar,5));
    double *r_bcue_f_m2 = REAL(VECTOR_ELT(fVar,6));
    double *r_focuecnb_f_m2 = REAL(VECTOR_ELT(fVar,7));
    double *r_icecue_f_m2 = REAL(VECTOR_ELT(fVar,8));
    double *r_ovcue_f_m2 = REAL(VECTOR_ELT(fVar,10));
    double *r_vcstOthm_f2 = REAL(VECTOR_ELT(fVar,11));
    double *r_rtbsOthm_f2 = REAL(VECTOR_ELT(fVar,12));
    double *r_cshrTothm_f2 = REAL(VECTOR_ELT(fVar,13));
    double *r_sshrOthm_f2 = REAL(VECTOR_ELT(fVar,14));
    double *r_gcue_f2 = REAL(VECTOR_ELT(fVar,16));
    double *r_fvolue_f2 = REAL(VECTOR_ELT(fVar,17));
    double *r_ovcue_f2 = REAL(VECTOR_ELT(fVar,18));
    double *r_oilcue_f2 = REAL(VECTOR_ELT(fVar,19));
    double *r_bcue_f2 = REAL(VECTOR_ELT(fVar,20));
    double *r_focuecnb_f2 = REAL(VECTOR_ELT(fVar,21));
    double *r_icecue_f2 = REAL(VECTOR_ELT(fVar,22));
    double *r_GVLothsue_f_m2 = REAL(VECTOR_ELT(fVar,23));
    double *r_GVLothsrefue_f_m2 = REAL(VECTOR_ELT(fVar,24));
    double *r_GVLothmetue_f2 = REAL(VECTOR_ELT(fVar,25));
    double *r_GVLothsue_f2 = REAL(VECTOR_ELT(fVar,26));
    double *r_ccwr_f2 = REAL(VECTOR_ELT(fVar,27));
    //double *r_opersc_f2 = REAL(VECTOR_ELT(fVar,28));
    double *r_GVLoths_f2 = REAL(VECTOR_ELT(fVar,29));
    double *r_onvcr_f2 = REAL(VECTOR_ELT(fVar,30));
    double *r_rtbs_f2 = REAL(VECTOR_ELT(fVar,31));
    double *r_rtbs_f_m2 = REAL(VECTOR_ELT(fVar,32));


   SEXP countGVLf;
   PROTECT(countGVLf = NEW_NUMERIC(nbF)); // --> 109
   double *r_countGVLf = REAL(countGVLf);


   for (int e = 0 ; e < nbE ; e++) {

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        nbC = length(getListElement(elmt, "modC"));

        double *r_Lbio_f_sum_t_e = REAL(aggregObj(VECTOR_ELT(out_L_efmct2, e),dimCstF)); //on conserve la dimension temporelle
        double *r_Lothm_f_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),40));
        double *r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
        double *r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct2, e));
        double *r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
        int *dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct2, e), install("DimCst")))));
        int *dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));

        //---------------------
        //équations de la table "t"
        //---------------------

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        //-- 1. Lbio_f

            if (e==0) {

                r_Lbio_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_Loths_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] + //Loths_f
                    r_Lothm_f_e2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] + //Lothm_f_e
                    r_Lbio_f_sum_t_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            } else {

                r_Lbio_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_Lbio_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                    r_Lothm_f_e2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] +
                    r_Lbio_f_sum_t_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            }

        r_countGVLf[ind_f] = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 2. GVL_f_m_e

             double count = 0.0;

             for (int ind_c = 0 ; ind_c < nbC ; ind_c++){

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]))

                count = count +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

             }

            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = count / pow(1+0.0,ind_t) ;
            r_countGVLf[ind_f] = r_countGVLf[ind_f] + (count / pow(1+0.0,ind_t));



        //-- 3. GVLtot_f_m


            if (e==0) {

                if (othsFM==1) {

                    if (adj==2) {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLothsrefue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                            r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                            r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]
                             / pow(1+0.0,ind_t)) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] ;

                    } else {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLothsue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                            r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                            r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]
                            / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                    }

                } else {

                    if (adj==2) {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            ( r_GVLothsref_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] ;

                    } else {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLoths_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                    }

                }

            } else {

                r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
                  r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }
        }
   }

  SET_VECTOR_ELT(VECTOR_ELT(out_Eco,1), e, VECTOR_ELT(VECTOR_ELT(eVar, e),41));

  UNPROTECT(1);

}


    // à ce stade, plus de considération d'espèce pour les indicateurs

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            double countEff = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++)

                countEff = countEff +
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];


            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            //-- 4. GVLav_f_m

                r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + 0*dim_nbv_f_m[3]];

             //-- 5. GVLtot_f

                if (ind_m==0) {

                    if (lev==1) {

                        if (report==0) {

                            if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmet_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                            } else {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmet_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) );
                            }

                        } else {

                            if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmetue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                    ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff) +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                            } else {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmetue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                    ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff);

                            }
                        }

                    } else { //lev = 2

                        if (oths==0) {

                            r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                r_countGVLf[ind_f] +
                                (r_GVLoths_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));

                        } else {

                           if (report==0) {

                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                        r_countGVLf[ind_f] +
                                        (r_GVLothsue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                        ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff);

                            } else {

                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                        r_countGVLf[ind_f] +
                                        ((r_GVLothsue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]);

                            }
                        }
                    }

                } else { //ind_m>0

                    if (lev==1) {

                        if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                        }
                    }
                }



            //-- 7. NGVLav_f_m

                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                    (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]]);


            } //on sort de la boucle sur les niveaux métiers


            //-- 6. GVLav_f

                r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]];


            //-- 8. NGVLav_f

                r_NGVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]]);

            //---------------------------------------------------------------------------------------------


            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            //-- 9. vcst_f_m

                if (!all_is_na(VECTOR_ELT(fVar,5)) & !(all_is_na(VECTOR_ELT(fVar,6))) &
                    !(all_is_na(VECTOR_ELT(fVar,7))) & !(all_is_na(VECTOR_ELT(fVar,8)))) {

                    r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        (r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                        r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]] +
                        r_oilcue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] +
                        r_bcue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] +
                        r_focuecnb_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                        r_cnb_f_m[ind_f*dim_cnb_f_m[0] + ind_m*dim_cnb_f_m[1] + ind_t*dim_cnb_f_m[2] + ind_t*dim_cnb_f_m[3]] +
                        r_icecue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) *
                        r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] / pow(1+0.0,ind_t);

                } else {

                    r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        (r_ovcue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] +
                        r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                        r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]) *
                        r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] / pow(1+0.0,ind_t);

                }


             //-- 10. vcst_f

            if (!all_is_na(VECTOR_ELT(out_Eco,7))) { //alias 'vcst_f_m_out'

                if (ind_m==0) {

                    if (!ISNA(r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])){

                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        ((r_vcstOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                        r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    } else {

                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_vcstOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    }

                } else {

                    if (!ISNA(r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        ((r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                }

            } else {

                if (!all_is_na(VECTOR_ELT(fVar,19)) & !(all_is_na(VECTOR_ELT(fVar,20))) &
                    !(all_is_na(VECTOR_ELT(fVar,21))) & !(all_is_na(VECTOR_ELT(fVar,22)))) {

                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                        r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] +
                        r_oilcue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] +
                        r_bcue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] +
                        r_focuecnb_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                        r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] +
                        r_icecue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]) *
                        r_ue_f[ind_f*dim_ue_f[0] + ind_m*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t);

                } else {

                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_ovcue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] +
                        r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                        r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]]) *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t);

                }
            }


             //-- 11. rtbs_f_m

                r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_vcst_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];


             //-- 12. rtbs_f

            if (lev==1) {

                if (ind_m==0) {

                    if (!ISNA(r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])){


                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        ((r_rtbsOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                        r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    } else {

                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_rtbsOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                        }

                } else {

                    if (!ISNA(r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        ((r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                }

            } else {


                 r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_NGVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            }

            //version actualisée
            r_rtbsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 13. cshrT_f_m


             if (perscCalc==0) {

               r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    0.01*r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                    r_rtbs_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]];

             } else {

                r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    0.01*r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
             }


            //-- 14. cshrT_f


            if (perscCalc==0) {

               r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];

            } else {

            if (lev==1) {

                if (ind_m==0) {

                    if (!ISNA(r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])){


                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        ((r_cshrTothm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                        r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    } else {

                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_cshrTothm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    }

                } else {

                    if (!ISNA(r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        ((r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                }

            } else {


                 r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            }
            }

             //-- 15. sshr_f_m

                r_sshr_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];


            //-- 16. sshr_f

            if (lev==1) {

                if (ind_m==0) {

                    if (!ISNA(r_sshr_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]])){

                    r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        ((r_sshrOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                        r_sshr_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    } else {

                    r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        (r_sshrOthm_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                    }

                } else {

                    if (!ISNA(r_sshr_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]))

                    r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        ((r_sshr_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) /
                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]);

                }

            } else {

                 r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            }

            } //on sort de la boucle sur les niveaux métiers


             //-- 17. ncshr_f

                r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_eec_f[ind_f*dim_eec_f[0] + 0*dim_eec_f[1] + 0*dim_eec_f[2] + ind_t*dim_eec_f[3]] / pow(1+0.0,ind_t) );

             //-- 18. ocl_f

                r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_mwh_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                    r_nbh_f[ind_f*dim_nbh_f[0] + 0*dim_nbh_f[1] + 0*dim_nbh_f[2] + ind_t*dim_nbh_f[3]] / pow(1+0.0,ind_t);

             //-- 19. cs_f

                r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                        //version actualisée
            r_csAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 20. csTot_f

                r_csTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

             //-- 21. gva_f

                r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] +
                    r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]] +
                    r_insp_f[ind_f*dim_insp_f[0] + 0*dim_insp_f[1] + 0*dim_insp_f[2] + ind_t*dim_insp_f[3]] +
                    r_ownc_f[ind_f*dim_ownc_f[0] + 0*dim_ownc_f[1] + 0*dim_ownc_f[2] + ind_t*dim_ownc_f[3]] +
                    r_mngc_f[ind_f*dim_mngc_f[0] + 0*dim_mngc_f[1] + 0*dim_mngc_f[2] + ind_t*dim_mngc_f[3]] +
                    r_licc_f[ind_f*dim_licc_f[0] + 0*dim_licc_f[1] + 0*dim_licc_f[2] + ind_t*dim_licc_f[3]] +
                    r_comc_f[ind_f*dim_comc_f[0] + 0*dim_comc_f[1] + 0*dim_comc_f[2] + ind_t*dim_comc_f[3]] +
                    r_onvcr_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]) / pow(1+0.0,ind_t) ;

                //version actualisée
            r_gvaAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 22. ccw_f


            if (perscCalc==1) {

                r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ; //+
                    //(r_opersc_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));

            } else {

                if (perscCalc==2) {

                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ccwr_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                } else {

                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] ;//+
                    //(r_opersc_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));
                }
            }
             //-- 23. ccwCr_f

             r_ccwCr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];


            //-- 24. wageg_f

             r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];

            //-- 25. wagen_f

             r_wagen_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];


             //-- 26. gcf_f

                r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

               //version actualisée
            r_gcfAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



             //-- 27. ngcf_f

                r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]] / pow(1+0.0,ind_t) );


             //-- 28. gp_f

                r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]] / pow(1+0.0,ind_t) );

             //-- 29. ssTot_f

                r_ssTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

             //-- 30. ps_f

                r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                     r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]) *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

                //version actualisée
            r_psAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                         r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 31. sts_f

                r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]] *
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

            //version actualisée
            r_stsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                         r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 32. ber_f

                r_ber_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_vcst_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                    ((r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] +
                    r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]] +
                    r_insp_f[ind_f*dim_insp_f[0] + 0*dim_insp_f[1] + 0*dim_insp_f[2] + ind_t*dim_insp_f[3]] +
                    r_ownc_f[ind_f*dim_ownc_f[0] + 0*dim_ownc_f[1] + 0*dim_ownc_f[2] + ind_t*dim_ownc_f[3]] +
                    r_mngc_f[ind_f*dim_mngc_f[0] + 0*dim_mngc_f[1] + 0*dim_mngc_f[2] + ind_t*dim_mngc_f[3]] +
                    r_licc_f[ind_f*dim_licc_f[0] + 0*dim_licc_f[1] + 0*dim_licc_f[2] + ind_t*dim_licc_f[3]] +
                    r_comc_f[ind_f*dim_comc_f[0] + 0*dim_comc_f[1] + 0*dim_comc_f[2] + ind_t*dim_comc_f[3]] +
                    r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]] +
                    r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]] +
                    r_onvcr_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]) / pow(1+0.0,ind_t)) ;


             //-- 33. ratio_gva_GVL_f

                r_ratio_gva_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 34. ratio_gcf_GVL_f

                r_ratio_gcf_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 35. ratio_fc_GVL_f

                r_ratio_fc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                        r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t)) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 36. ratio_oilc_GVL_f

                r_ratio_oilc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_oilcue_f2[ind_f*dim_oilc_f[0] + 0*dim_oilc_f[1] + 0*dim_oilc_f[2] + ind_t*dim_oilc_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t) ) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 37. ratio_bc_GVL_f

                r_ratio_bc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_bcue_f2[ind_f*dim_bc_f[0] + 0*dim_bc_f[1] + 0*dim_bc_f[2] + ind_t*dim_bc_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t) ) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 38. ratio_foc_GVL_f

                r_ratio_foc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_focuecnb_f2[ind_f*dim_foc_f[0] + 0*dim_foc_f[1] + 0*dim_foc_f[2] + ind_t*dim_foc_f[3]] *
                        r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t) ) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 39. ratio_icec_GVL_f

                r_ratio_icec_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_icecue_f2[ind_f*dim_icec_f[0] + 0*dim_icec_f[1] + 0*dim_icec_f[2] + ind_t*dim_icec_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 40. ratio_gc_GVL_f

                r_ratio_gc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_gcue_f2[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t) ) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 41. ratio_vc_GVL_f

                r_ratio_vc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_vc_f[ind_f*dim_vc_f[0] + 0*dim_vc_f[1] + 0*dim_vc_f[2] + ind_t*dim_vc_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 42. ratio_rep_GVL_f

                r_ratio_rep_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 43. ratio_mngc_GVL_f

                r_ratio_mngc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_mngc_f[ind_f*dim_mngc_f[0] + 0*dim_mngc_f[1] + 0*dim_mngc_f[2] + ind_t*dim_mngc_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 44. ratio_licc_GVL_f

                r_ratio_licc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_licc_f[ind_f*dim_licc_f[0] + 0*dim_licc_f[1] + 0*dim_licc_f[2] + ind_t*dim_licc_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 45. ratio_fvol_GVL_f

                r_ratio_fvol_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 46. ratio_fvol_Lbio_f

                r_ratio_fvol_Lbio_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] /
                    r_Lbio_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 47. ratio_fvol_gva_f

                r_ratio_fvol_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] /
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 48. ratio_gcf_gva_f

                r_ratio_gcf_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 49. ratio_K_cnb_f

                r_ratio_K_cnb_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) ) /
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]];

            //-- 50. ratio_GVL_K_f

                r_ratio_GVL_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

           //-- 51. ratio_gcf_K_f

                r_ratio_gcf_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

           //-- 52. ratio_ngcf_K_f

                r_ratio_ngcf_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

          //-- 53. ratio_gp_K_f

                r_ratio_gp_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

          //-- 54. ratio_GVL_cnb_ue_f

                r_ratio_GVL_cnb_ue_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                     r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]);

        }

if (ind_t==0) UNPROTECT(107);
UNPROTECT(109);

}}














//---------------------------------
//
// Module de recrutement aléatoire
//
//---------------------------------



extern "C" {

void BioEcoPar::RecAlea(SEXP list, SEXP listSto, int ind_t, int type, int *recTyp) //list : liste des paramètres d'entrée ; listSto : liste des variables d'opérations stochastiques ; type : 1 -> samples sur l'historique (temps variable), 2 -> samples sur l'historique (temps constant), 3 -> loi de distribution
{

if (type<3) {

       SEXP elmtIn, elmtMeanSto, elmtResSto, MeanSto, ResSto, Rec, dimRec;
    //on tire au sort pour chacune des espèces modélisées un résidu et on l'ajoute à la moyenne géométrique pré-calculée
       int index = 0;

    for (int e = 0 ; e < nbE ; e++) {

        if (trim[e]==0) { //------------------------------------------------------------------ pas de dimension trimestre

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

                        if (type==1) { //multiple tirage d'indice (1 par espèce)

                            index = ll;
                            while (index >= ll) index = (int)(rand() / (((double)RAND_MAX + 1)/ ll));

                        } else {        //unique tirage d'indice pour les espèces considérées (historiques de même taille)

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

        } else { // ------------------------------------------------------------------------ dimension trimestre

             for (int Qt=0 ; Qt<4; Qt++){

                        if (recTyp[e]==1) {

                        PROTECT(elmtIn = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                        if (type==1) {

                            PROTECT(elmtMeanSto = getListElement(listSto, "GeoMeanRec"));
                            PROTECT(elmtResSto = getListElement(listSto, "RecResiduals"));

                        } else {

                            PROTECT(elmtMeanSto = getListElement(listSto, "GeoMeanRecLink"));
                            PROTECT(elmtResSto = getListElement(listSto, "RecResidualsLink"));

                        }

                        PROTECT(MeanSto = getListElement(getListElement(elmtMeanSto, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));
                        PROTECT(ResSto = getListElement(getListElement(elmtResSto, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));

                        int ll = length(ResSto);

                    if (ll > 0) {

                        if (type==1) { //multiple tirage d'indice (1 par espèce)

                            index = ll;
                            while (index >= ll) index = (int)(rand() / (((double)RAND_MAX + 1)/ ll));

                        } else {        //unique tirage d'indice pour les espèces considérées (historiques de même taille)

                            if (e==0 & Qt==0) {

                                index = ll;
                                while (index >= ll) index = (int)(rand() / (((double)RAND_MAX + 1)/ ll));

                            }

                        }

                        PROTECT(Rec = getListElement(getListElement(elmtIn, "N_i0t"),CHAR(STRING_ELT(trimInt,Qt))));
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

        }

    }



} else {



    SEXP elmtIn, elmtDist, elmtDistParOne, elmtDistParTwo, elmtDistParThree,
         elmtDistSp, elmtDistParOneSp, elmtDistParTwoSp, elmtDistParThreeSp, Rec, dimRec;
    //on génère une variable aléatoire suivant une loi log-normale de paramètres spécifiés

    for (int e = 0 ; e < nbE ; e++) {

        if (trim[e]==0) { //------------------------------------------------------------------ pas de dimension trimestre

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

                //if ... pour les autres lois --> à compléter

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

    } else { // ------------------------------------------------------------------- dimension trimestrielle

        for (int Qt=0 ; Qt<4; Qt++){

                    if (recTyp[e]==1) {

                    PROTECT(elmtIn = getListElement(list, CHAR(STRING_ELT(sppList,e))));
                    PROTECT(elmtDist = getListElement(listSto, "RecDist"));
                    PROTECT(elmtDistParOne = getListElement(listSto, "RecDistPar1"));
                    PROTECT(elmtDistParTwo = getListElement(listSto, "RecDistPar2"));
                    PROTECT(elmtDistParThree = getListElement(listSto, "RecDistPar3"));

                    PROTECT(elmtDistSp = getListElement(getListElement(elmtDist, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));
                    PROTECT(elmtDistParOneSp = getListElement(getListElement(elmtDistParOne, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));
                    PROTECT(elmtDistParTwoSp = getListElement(getListElement(elmtDistParTwo, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));
                    PROTECT(elmtDistParThreeSp = getListElement(getListElement(elmtDistParThree, CHAR(STRING_ELT(sppList,e))),CHAR(STRING_ELT(trimInt,Qt))));

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

                    //if ... pour les autres lois --> à compléter

                    PROTECT(Rec = getListElement(getListElement(elmtIn, "N_i0t"),CHAR(STRING_ELT(trimInt,Qt))));
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

    }
}

}}





//---------------------------------
//
// Module de modélisation de relations S/R
//
//---------------------------------



extern "C" {

void BioEcoPar::SRmod(SEXP list, SEXP listSR, int ind_t, SEXP TypeSR, int *srind)
        //list : liste des paramètres d'entrée;
        //listSR : liste des paramètres a,b&c du modèle SR + e.t bruit normal ou lognormal + type de bruit (1=normal, 2=lognormal) (un vecteur de longueur 5 par espèce modélisée contenant des "doubles")
        //type : type de relation Stock-Recrutement : (liste de longueur "nb d'espèces modélisées" contenant des entiers)
        //                                              1 -> recrutement constant moyen (rec~a)
        //                                              2 -> Hockey stick (rec ~ (si (ssb<=b) a*ssb sinon a*b))
        //                                              3 -> Beverton & Holt (rec ~ a*ssb/(b+ssb))
        //                                              4 -> Ricker (rec ~ a*ssb*exp(-b*ssb))
        //                                              5 -> Shepherd (rec ~ a*ssb/(1+ (ssb/b)^c))
        //                                              6 -> Hockey Stick Quadratic (rec ~ (si (ssb<=b*(1-c)) a*ssb ; si (b*(1-c)<ssb<b*(1+c)) a*(ssb-((ssb-b*(1-c))^2)/(4*b*c)) ; sinon a*b))
{

SEXP ans, rnames=R_NilValue;
double *rans, *paramet, *ssb;
int typeSR, fstAge;

if (ind_t==0) {  //on formatte l'objet en sortie

    PROTECT(rnames = allocVector(STRSXP, nbE));
    setAttrib(out_SRmod, R_NamesSymbol, rnames);

}


for (int e = 0 ; e < nbE ; e++) {

    if (srind[e]==1) {  //activation du module

    if (ind_t==0) { //deuxième étape d'initialisation (niveau espèce)

        PROTECT(ans = NEW_NUMERIC(nbT));
        setAttrib(ans, R_NamesSymbol, times);
        SET_VECTOR_ELT(out_SRmod, e, ans);
        SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
    }

    rans = REAL(VECTOR_ELT(out_SRmod, e));
    paramet = REAL(VECTOR_ELT(listSR, e));
    typeSR = INTEGER(VECTOR_ELT(TypeSR, e))[0];
    if (ind_t>0) ssb = REAL(VECTOR_ELT(out_SSB_et, e)); else ssb = &NA_REAL;

    //il nous faut aussi le décalage temporel dû au premier âge modélisé -> un SSB(t) générera un R(t+age0+1)
    fstAge = CHAR(STRING_ELT(VECTOR_ELT(namDC, e),0))[0] - '0'; ++fstAge;

    //on en profite pour initialiser l'objet pour les premières années pour lesquelles on devra aller chercher l'info dans Ni0
    if (ind_t<fstAge) {

        rans[ind_t] = NA_REAL;

    } else {

    switch (typeSR) {

        case 1 :

        rans[ind_t] = paramet[0]; break;

        case 2 :

        if (ssb[ind_t-fstAge] <= paramet[1])
            rans[ind_t] = paramet[0] * ssb[ind_t-fstAge];
        else
            rans[ind_t] = paramet[0] * paramet[1]; break;

        case 3 :

        rans[ind_t] = paramet[0] * ssb[ind_t-fstAge] / (paramet[1] + ssb[ind_t-fstAge]); break;

        case 4 :

        rans[ind_t] = paramet[0] * ssb[ind_t-fstAge] * exp(-1.0 * paramet[1] * ssb[ind_t-fstAge]); break;

        case 5 :

        rans[ind_t] = paramet[0] * ssb[ind_t-fstAge] / (1 + pow(ssb[ind_t-fstAge] / paramet[1] , paramet[2])); break;

        case 6 :

        if (ssb[ind_t-fstAge] <= (paramet[1]*(1-paramet[2])))

            rans[ind_t] = paramet[0] * ssb[ind_t-fstAge];

        else {

            if (ssb[ind_t-fstAge] >= (paramet[1]*(1+paramet[2])))

                rans[ind_t] = paramet[0] * paramet[1];

            else

                rans[ind_t] = paramet[0]*(ssb[ind_t-fstAge] - (pow(ssb[ind_t-fstAge] - paramet[1]*(1-paramet[2]),2.0)/(4*paramet[1]*paramet[2])));

            } break;

        default :

        rans[ind_t] = NA_REAL;

    }
    //il ne reste plus qu'à ajouter le bruit blanc issue de N(0,sigma) avec sigma = paramet[3]
    double v_alea = 0.0;
GetRNGstate();
        if (!ISNA(paramet[3])) v_alea = rnorm(0.0,paramet[3]); ////Rprintf("%f ",v_alea);Rprintf("%f ",rnorm(0.0,0.157));
PutRNGstate();
    if (paramet[4]==1)  //bruit de loi normale

        rans[ind_t] = rans[ind_t] + v_alea;

    else                //bruit de loi lognormale

        rans[ind_t] = rans[ind_t] * exp(v_alea);

    }

      if (ind_t==0) UNPROTECT(1);

    }}

    if (ind_t==0) UNPROTECT(1);

}}







//---------------------------------
//
// Module de gestion des scénarios
//
//---------------------------------



extern "C" {

void BioEcoPar::Scenario(SEXP list, SEXP listScen, int ind_t) //list : liste des paramètres d'entrée ; listScen : liste des multiplicateurs pour un scénario donné
{

//1er niveau de la liste de multiplicateurs : Fleet ou Species --> on cible la partie de "list" correspondante

SEXP mult_lvl_1, target_lvl_1, mult_lvl_2, target_lvl_2, namVar, namElt, dimMult, dimTarget, fTarg, fMult;

int nbElt = length(listScen);

for (int elt = 0 ; elt < nbElt ; elt++) {

    PROTECT(namElt = STRING_ELT(getAttrib(listScen, R_NamesSymbol), elt));

    PROTECT(mult_lvl_1 = getListElement(listScen, CHAR(namElt)));

    if (mult_lvl_1 != NULL) {

        PROTECT(target_lvl_1 = getListElement(list, CHAR(namElt)));
        int nbVar = length(mult_lvl_1);

        for (int i = 0 ; i < nbVar ; i++) {

            PROTECT(namVar = STRING_ELT(getAttrib(mult_lvl_1, R_NamesSymbol), i));
            PROTECT(mult_lvl_2 = getListElement(mult_lvl_1, CHAR(namVar)));

            //ici, selon que la variable considérée est un input ou une variable interne (ex : Foth_i), on agit différemment

            if (strcmp(CHAR(namVar), "Foth_i") == 0) {

                //Rprintf("%i",IS_NUMERIC(VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 44)));
                PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 44));

            } else {

                if (strcmp(CHAR(namVar), "F_fmi") == 0) {

                        PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 0));

                } else {

                    if (strcmp(CHAR(namVar), "GVLoths_fm") == 0) {

                        PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 1));

                    } else {

                        if (strcmp(CHAR(namVar), "GVLothsref_fm") == 0) {

                            PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 2));

                        } else {

                            if (strcmp(CHAR(namVar), "GVLothsue_fm") == 0) {

                                PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 23));

                            } else {

                                if (strcmp(CHAR(namVar), "GVLothsrefue_fm") == 0) {

                                    PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 24));

                                } else {

                                    if (strcmp(CHAR(namVar), "GVLothsue_f") == 0) {

                                        PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 26));

                                    } else {

                                        if (strcmp(CHAR(namVar), "GVLoths_f") == 0) {

                                            PROTECT(target_lvl_2 = VECTOR_ELT(fVar, 29));

                                        } else {

                                            PROTECT(target_lvl_2 = getListElement(target_lvl_1, CHAR(namVar)));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (IS_NUMERIC(target_lvl_2)) { //alors, pas de niveau trimestre

                int *dimM, *dimT;

                PROTECT(dimMult = getAttrib(mult_lvl_2, install("DimCst")));
                //si 'target_lvl_2' est un élément de eVar, s'assurer au préalable de l'existence de l'attribut DimCst
                PROTECT(dimTarget = getAttrib(target_lvl_2, install("DimCst")));

                dimM = INTEGER(dimMult); dimT = INTEGER(dimTarget);

            //tests sur les dimensions
                if (dimM[0]>dimT[0] | dimM[1]>dimT[1] | dimM[2]>dimT[2]) error("Wrong dimensions specification in 'Scenario' !!\n");

                PROTECT(fTarg = iDim(dimT));
                PROTECT(fMult = iDim(dimM));

                int *ftarg = INTEGER(fTarg);
                int *fmult = INTEGER(fMult);

                double *target = REAL(target_lvl_2), *mult = REAL(mult_lvl_2);

            //et on applique la mise à jour

                for (int ind_f = 0 ; ind_f < imax2(1,dimT[0]) ; ind_f++)
                for (int ind_m = 0 ; ind_m < imax2(1,dimT[1]) ; ind_m++)
                for (int ind_i = 0 ; ind_i < imax2(1,dimT[2]) ; ind_i++) {

                    target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] =
                    target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] *
                    mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]];

                }

                UNPROTECT(4);

            } else { //niveau trimestre

              for (int Qt=0 ; Qt<4; Qt++){

                int *dimM, *dimT;

                if (getListElement(mult_lvl_2,CHAR(STRING_ELT(trimInt,Qt))) != NULL) {

                PROTECT(dimMult = getAttrib(getListElement(mult_lvl_2,CHAR(STRING_ELT(trimInt,Qt))), install("DimCst")));
                //si 'target_lvl_2' est un élément de eVar, s'assurer au préalable de l'existence de l'attribut DimCst
                PROTECT(dimTarget = getAttrib(VECTOR_ELT(target_lvl_2, Qt), install("DimCst")));

                dimM = INTEGER(dimMult); dimT = INTEGER(dimTarget);

            //tests sur les dimensions
                if (dimM[0]>dimT[0] | dimM[1]>dimT[1] | dimM[2]>dimT[2]) error("Wrong dimensions specification in 'Scenario' !!\n");

                PROTECT(fTarg = iDim(dimT));
                PROTECT(fMult = iDim(dimM));

                int *ftarg = INTEGER(fTarg);
                int *fmult = INTEGER(fMult);

                double *target = REAL(VECTOR_ELT(target_lvl_2,Qt)), *mult = REAL(getListElement(mult_lvl_2,CHAR(STRING_ELT(trimInt,Qt))));

            //et on applique la mise à jour

                for (int ind_f = 0 ; ind_f < imax2(1,dimT[0]) ; ind_f++)
                for (int ind_m = 0 ; ind_m < imax2(1,dimT[1]) ; ind_m++)
                for (int ind_i = 0 ; ind_i < imax2(1,dimT[2]) ; ind_i++) {

                    target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] =
                    target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] *
                    mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]];

                }

                UNPROTECT(4);

              }
              }


            }

            UNPROTECT(3);

        }
        UNPROTECT(1);
    }

    UNPROTECT(2);
}

}}



//---------------------------------
//
//   Module de gestion
//
//---------------------------------



extern "C" {

double BioEcoPar::fxTAC_glob(double mult) //par temps IND_T pour une espèce donnée // m_f : vecteur de pondération par flottille
{
    SEXP listTemp;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));
//Rprintf("3");
    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));//Rprintf("4");
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));//Rprintf("5");
    double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));//Rprintf("6");
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));//Rprintf("7");
    double *mpond_f = REAL(m_f);//Rprintf("8");
    double *mpond_fm = REAL(m_fm);//Rprintf("9");
    double *mpond_oth = REAL(m_oth);//Rprintf("10");

    double result;


if (level==0) {

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            if (gestyp==1) g_nbdsF[ind_f] = fmax2(g_nbdsF[ind_f] + mpond_f[ind_f]*mult,0.0);
            if (gestyp==2) g_nbdsF[ind_f] = fmax2(g_nbdsF[ind_f] * (1 + mpond_f[ind_f]*mult),0.0);
        }
        if (var==2) {
            if (gestyp==1) g_nbvF[ind_f] = fmax2(g_nbvF[ind_f] + mpond_f[ind_f]*mult,0.0);
            if (gestyp==2) g_nbvF[ind_f] = fmax2(g_nbvF[ind_f] * (1 + mpond_f[ind_f]*mult),0.0);
        }

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m] + mpond_f[ind_f]*mult,0.0);
                if (gestyp==2) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m] * (1 + mpond_f[ind_f]*mult),0.0);
            }

            if (var==2) {
                if (gestyp==1) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m] + mpond_f[ind_f]*mult,0.0);
                if (gestyp==2) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m] * (1 + mpond_f[ind_f]*mult),0.0);
            }

         if (ind_m==0 & ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){//Rprintf("11");

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));//Rprintf("12");
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));

                 if (gestyp==1)
                    for (int ag = 0; ag < ni; ag++)
                        g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T] + mpond_oth[e]*mult,0.0);
                 if (gestyp==2)
                    for (int ag = 0; ag < ni; ag++)
                        g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T] * (1 + mpond_oth[e]*mult),0.0);

                }

        }


    }}


} else {


   for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        double countEff = 0.0;

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) {
                if (gestyp==1) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m] + mpond_fm[ind_f+nbF*ind_m]*mult,0.0);
                if (gestyp==2) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m] * (1 + mpond_fm[ind_f+nbF*ind_m]*mult),0.0);//Rprintf("13");
                if (!ISNA(g_nbdsFM[ind_f+nbF*ind_m])) countEff = countEff + g_nbdsFM[ind_f+nbF*ind_m] - REAL(NBDSFM)[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)];//Rprintf("14");
            }

            if (var==2) {
                if (gestyp==1) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m] + mpond_fm[ind_f+nbF*ind_m]*mult,0.0);
                if (gestyp==2) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m] * (1 + mpond_fm[ind_f+nbF*ind_m]*mult),0.0);
                if (!ISNA(g_nbvFM[ind_f+nbF*ind_m])) countEff = countEff + g_nbvFM[ind_f+nbF*ind_m] - REAL(NBVFM)[ind_f + nbF*ind_m + nbF*nbMe*(delay-1)];
            }

//ajout MM 22/11/2013

//         if (ind_m==0 & ind_f==0) {
//
//            for (int e = 0 ; e < nbE ; e++){//Rprintf("11");
//
//                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));//Rprintf("12");
//                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
//
//                 if (gestyp==1)
//                    for (int ag = 0; ag < ni; ag++)
//                        g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T] + mpond_oth[e]*mult,0.0);
//                 if (gestyp==2)
//                    for (int ag = 0; ag < ni; ag++)
//                        g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T] * (1 + mpond_oth[e]*mult),0.0);
//
//                }
//
//        }
//
//fin ajout MM 22/11/2013

    }
//Rprintf("1");
        if (var==1) {
            if (gestyp==1) g_nbdsF[ind_f] = fmax2(countEff + REAL(NBDSF)[ind_f + nbF*(delay-1)],0.0);
            if (gestyp==2) g_nbdsF[ind_f] = fmax2(countEff + REAL(NBDSF)[ind_f + nbF*(delay-1)],0.0);
        }
        if (var==2) {
            if (gestyp==1) g_nbvF[ind_f] = fmax2(countEff + REAL(NBVF)[ind_f + nbF*(delay-1)],0.0);
            if (gestyp==2) g_nbvF[ind_f] = fmax2(countEff + REAL(NBVF)[ind_f + nbF*(delay-1)],0.0);
        }

    }
//Rprintf("15");
}

if (trgt==1 | trgt==3) {//on vise un TAC

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

//Rprintf("16");
    SEXP nDim = allocVector(INTSXP,4);
    int *nd = INTEGER(nDim); for (int i = 0; i<3; i++) nd[i] = 0; nd[3] = nbT;
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    result = TAC_glob[IND_T]-tot[IND_T];

} else { //on vise un Fbar

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    double *tot = REAL(VECTOR_ELT(out_Fbar_et, eTemp));
    result = Fbar_trgt[IND_T]-tot[IND_T];

}


    UNPROTECT(2);

    return result;

}
}



//------------------------------------------
// Module de gestion : ajustement des variables d'effort (nbds (paramètre "var" = 1) ou nbv (paramètre "var" = 2))
// avec objectif d'atteinte du TAC (paramètre "trgt = 1") OU du Fbar (paramètre "trgt = 2")
// Un 3ème paramètre "delay" spécifie le délai de première applicaton de l'ajustement (valeur par défaut et minimale = 1).
// Enfin, un 4ème paramètre "upd" (update) permet de spécifier si le multiplicateur s'applique à la donnée initiale à chaque pas de temps ("upd" = 1),
// ou si elle s'applique à la donnée à l'instant précédent ("upd" = 2).
//------------------------------------------



extern "C" {

void BioEcoPar::Gestion(SEXP list, int ind_t) //paramètres en entrée pas forcément utiles dans la mesure où ils doivent rester constant tout au long de la simulation
{

//on teste la validité des paramètres d'entrée

    if (var!=1 & var!=2) error("Wrong 'var' parameter in 'Gestion' module!!\n");
    if (trgt!=1 & trgt!=2 & trgt!=3) error("Wrong 'trgt' parameter in 'Gestion' module!!\n");
    if (delay<1 | delay>nbT) error("Wrong 'delay' parameter in 'Gestion' module!!\n");
    if (upd!=1 & upd!=2) error("Wrong 'upd' parameter in 'Gestion' module!!\n");
    if (level!=0 & level!=1) error("Wrong 'level' parameter in 'Gestion' module!!\n");

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;
    double *mu_;//Rprintf("1");
    if (var==1) mu_ = REAL(mu_nbds); else mu_ = REAL(mu_nbv);

    int NBMAX=1;
	int nb=NBMAX;
	float tol;
	int NbInter = 10;

    BEfn1 p = &BioEcoPar::fxTAC_glob;

    double *xb1,*xb2;
    xb1 = new double[NBMAX+1];
    xb2 = new double[NBMAX+1];
//Rprintf("2");
    zbrak(p,X1,X2,NbInter,xb1,xb2,&nb);
    for (int i=1;i<=nb;i++) {
        tol=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
        double result=zbrent(p,xb1[i],xb2[i],tol);
        mu_[IND_T] = result;
    }
    delete xb1;
    delete xb2;

   }

}}







// Numerical Recipes //----------------------------------------------------------------------------------------

// --------  détermination racine (unidimensionnel)

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


#define ITMAX 10000
#define EPS 3.0e-6
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double  BioEcoPar::zbrent(BEfn1 fx, double x1, double x2, double tol)
{
	int iter;
	double a=x1,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
	double fa=(this->*fx)(a),fb=(this->*fx)(b),fc,p,q,r,s,tol1,xm;

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
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
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
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
Rprintf("Maximum number of iterations exceeded in zbrent \n");
	return b;//0.0;  //modif 17/05/2013
}

#undef ITMAX
#undef EPS
#undef SIGN

// --------  simplex multi-dimensionnel

//int MinimizeF(void);
//float func(float x[]);
//float *NRvector(long nl, long nh);
//float **NRmatrix(long nrl, long nrh, long ncl, long nch);
//void free_vector(float *v, long nl, long nh);
//void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
//void amoeba(float **p, float y[], int ndim, float ftol, float (*funk)(float []), int *nfunk);

#define NR_END 1
#define FREE_ARG char*


double *BioEcoPar::NRvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) Rprintf("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **BioEcoPar::NRmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) Rprintf("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) Rprintf("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void BioEcoPar::free_vector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void BioEcoPar::free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


#define TINY 0.001
#define NMAX 10000
#define GET_PSUM \
                    for (j=1;j<=ndim;j++) {\
                    for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
                    psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}



void BioEcoPar::amoeba(BEfn1_F funk, double **p, double y[], int ndim, double ftol, int *nfunk) {

    //float amotry(float **p, float y[], float psum[], int ndim, float (*funk)(float []), int ihi, float fac);
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    double rtol,sum,swap,ysave,ytry,*psum;

    psum=NRvector(1,ndim);
    *nfunk=0;
    GET_PSUM
    for (;;) {
        ilo=1;
        ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
        for (i=1;i<=mpts;i++) {
            if (y[i] <= y[ilo]) ilo=i;
            if (y[i] > y[ihi]) {
                inhi=ihi;
                ihi=i;
            } else if (y[i] > y[inhi] && i != ihi) inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        //Rprintf("rtol %f \n",rtol);
        if (rtol < ftol) {
            SWAP(y[1],y[ilo])
            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
            break;
        }
        if (*nfunk >= NMAX) Rprintf("NMAX exceeded : rtol = %f\n",rtol); //break;}//'Rprintf' remplace 'nrerror'
        *nfunk += 2;

        //BEfn1_F FUNK = &BioEcoPar::fxTAC_glob;

        ytry=amotry(funk,p,y,psum,ndim,ihi,-1.0);
        if (ytry <= y[ilo]) {
            ytry=amotry(funk,p,y,psum,ndim,ihi,2.0);
        } else if (ytry >= y[inhi]) {
                ysave=y[ihi];
                ytry=amotry(funk,p,y,psum,ndim,ihi,0.5);
                if (ytry >= ysave) {
                    for (i=1;i<=mpts;i++) {
                        if (i != ilo) {
                            for (j=1;j<=ndim;j++)
                                p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
                            y[i]=(this->*funk)(psum);
                        }
                    }
                    *nfunk += ndim;
                    GET_PSUM
                }
               } else --(*nfunk);
    }
    free_vector(psum,1,ndim);
}


double BioEcoPar::amotry(BEfn1_F funk, double **p, double y[], double psum[], int ndim, int ihi, double fac) {
    int j;
    double fac1,fac2,ytry,*ptry;
    ptry=NRvector(1,ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=(this->*funk)(ptry);
    if (ytry < y[ihi]) {
        y[ihi]=ytry;
        for (j=1;j<=ndim;j++) {
            psum[j] += ptry[j]-p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    free_vector(ptry,1,ndim);
    return ytry;
}


#define MP 4
#define NP 3
#define FTOL 0.0000001

extern "C" {

double BioEcoPar::func(double *x)
{
	return ((x[1]-23.14)*(x[1]-23.14) + (x[2]-0.256)*(x[2]-0.256) + (x[3]+17.45)*(x[3]+17.45));
}

}




//------------------------------

extern "C" {

int BioEcoPar::QuotaExchV2(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
{

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;
    PxQ = pxQuIni;
    spQ = spp;

    int nfunc;
	int ITquotaExch = 20;
	int ITTOT = 7;

    Rprintf("TIME %i \n",IND_T);

    double **q = NRmatrix(1,2,1,1);
    double *multF = NRvector(1,nbF+1);
    double *z = NRvector(1,2);
    double *x = NRvector(1,1);

    bool GoOn = true;

	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT_customCstV2; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)
	BEfn1_F foo3 = &BioEcoPar::fxTAC_F_customCst2;

    double DIFF = 0.0;
    int IT = 0;

    while (GoOn) {

        for (int IT2 = 0 ; IT2 < ITTOT ; IT2++){

        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
    //Rprintf("IT %i \n",IT);


        for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

        Rprintf("T %i F %i \n",ind_t,ind_f);

        IND_F = ind_f;

        if (ind_f==nbF) { //optimisation de Foth --> x multiplicateur

            q[1][1]=x[1]=0.85;
            z[1]=(this->*foo3)(x);
            q[2][1]=x[1]=1;
            z[2]=(this->*foo3)(x);

        amoeba(foo3, q,z,1,ftol,&nfunc);

        multF[ind_f+1] = q[2][1];

        } else {

        x[1]=1.0;
        z[1]=(this->*foo2)(x);//Rprintf("z1 %f \n",z[1]);
        x[1]=350.0;
        z[2]=(this->*foo2)(x);//Rprintf("z2 %f \n",z[2]);

        if (z[1]>z[2]) multF[ind_f+1] = 350.0; else multF[ind_f+1] = 0.0;

        //if (IT==0)  Rprintf("Mult %f \n",multF[ind_f+1]);

        }
        }

                //2ème étape : redéfinition de Ztemp

        SEXP listTemp;

        PROTECT(listTemp = duplicate(list));
        PROTECT(eVarCopy = duplicate(eVar));

        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),fmin2(effortIni[ind_f]*expEff,350.0)-g_nbdsFM[ind_f+nbF*1]);
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

        }

        int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,spQ))), "modI"));

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        for (int i = 0 ; i < NBI ; i++) {

            Rprintf("diffZZ %f ", REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);
            Ztemp[i+1] = Ztemp[i+1] + lambda*(REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);

        }

        //... et on recommence jusqu'à convergence de Z vers la valeur coïncidant avec les mortalités marginales

        UNPROTECT(2);

        }

        //3ème étape : redéfinition de PxQ

    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),fmin2(effortIni[ind_f]*expEff,350.0)-g_nbdsFM[ind_f+nbF*1]);
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

           if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }


    }

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        PROTECT(nDimF = allocVector(INTSXP,4));
        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

            //calcul de diffLQ
        double diffLQ = 0.0;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);

        }

        if ((DIFF*diffLQ)<0) lambda = lambda/3; //ie si DIFF et diffLQ de signe différent

        if ((diffLQ<=0) & (((diffLQ*lambda)*(diffLQ*lambda)<0.25) | (IT>ITquotaExch))) GoOn = false; //on ne s'arrête que si diffLQ<=0 (Quota respecté)

        IT++;

        DIFF = diffLQ;

//        Rprintf("\n");
        Rprintf("diffLQ %f \n",diffLQ);
//
        Rprintf("PxQ_1 %f ",PxQ);

        if (diffLQ<0) {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmin2(lambda*diffLQ,-0.1),pxQuMax));
        } else {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmax2(lambda*diffLQ,0.1),pxQuMax));
        }

        Rprintf("PxQ_2 %f \n",PxQ);
        Rprintf("lambda %f \n",lambda);
        //... et on recommence

        UNPROTECT(3);

    }


    double *pxQuot = REAL(VECTOR_ELT(out_PQuot_et,spQ));
    pxQuot[ind_t] = PxQ;

	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),fmin2(effortIni[ind_f]*expEff,350.0)-nbdsFM_G[ind_f+nbF*1]);
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

    }

    free_matrix(q,1,2,1,1);
    free_vector(multF,1,nbF+1);
    free_vector(z,1,2);
    free_vector(x,1,1);

    }

	return 0;
}

}

extern "C" {

double BioEcoPar::fxMaxProf_FT_customCstV2(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));

    double *gcfF;

    //for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[IND_F+nbF*0] = fmin2(fmax2(x[1],0.0),fmin2(effortIni[IND_F]*expEff,350.0)-g_nbdsFM[IND_F+nbF*1]);
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

//    Rprintf("O %f \n",x[ind_f+1]);
//    Rprintf("A %f \n",g_nbdsFM[ind_f+nbF*0]);
//    Rprintf("B %f \n",g_nbdsFM[ind_f+nbF*1]);
//    Rprintf("C %f \n",g_nbdsF[ind_f]);
//    Rprintf("D %f \n",effortIni[ind_f]);
//    Rprintf("F %i \n",ind_f);

    //}

    int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,3);

    Marche(listTemp, IND_T);



    if (ecodcf==0) {

        Economic(listTemp,IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                           EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_Eco, 25));

    } else {

        EcoDCF(listTemp, IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                         EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_EcoDCF, 20));

    }

    //Rprintf("GCF %f \n",gcfF[IND_F + nbF*IND_T]);

    //calcul de l'indicateur de profit à optimiser
    //OUT@output$gcf_f[indF,indT]-pxQ[indT]*(sum(OUT@outputSp$Li$Sole_commune[indF,,,indT],na.rm=TRUE)-TAC_f_t[indF,1])

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

    double result = 0.0;
    result = gcfF[IND_F + nbF*IND_T]*g_nbvF[IND_F] - PxQ * (totF[IND_F + nbF*IND_T] - TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("A %f \n",TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("B %f \n",gcfF[IND_F + nbF*IND_T]);
    //Rprintf("C %f \n",g_nbvF[IND_F]);
    //Rprintf("D %f \n",totF[IND_F + nbF*IND_T]);

    UNPROTECT(3);

    return (-1*result);

}
}

//Report

extern "C" {

int BioEcoPar::QuotaExchV2Report(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
{

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;
    PxQ = pxQuIni;
    spQ = spp;

    int nfunc;
	int ITquotaExch = 20;
	int ITTOT = 6;

    Rprintf("TIME %i \n",IND_T);

    double **q = NRmatrix(1,2,1,1);
    double *multF = NRvector(1,nbF+1);
    double *z = NRvector(1,2);
    double *x = NRvector(1,1);

    bool GoOn = true;

	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT_customReportV2; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)
	BEfn1_F foo3 = &BioEcoPar::fxTAC_F_customReport2;

    double DIFF = 0.0;
    int IT = 0;

    while (GoOn) {

        for (int IT2 = 0 ; IT2 < ITTOT ; IT2++){

        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
    //Rprintf("IT %i \n",IT);


        for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

        Rprintf("T %i F %i \n",ind_t,ind_f);

        IND_F = ind_f;

        if (ind_f==nbF) { //optimisation de Foth --> x multiplicateur

            q[1][1]=x[1]=0.85;
            z[1]=(this->*foo3)(x);
            q[2][1]=x[1]=1;
            z[2]=(this->*foo3)(x);

        amoeba(foo3, q,z,1,ftol,&nfunc);

        multF[ind_f+1] = q[2][1];

        } else {

        x[1]=1.0;
        z[1]=(this->*foo2)(x);//Rprintf("z1 %f \n",z[1]);
        x[1]=350.0;
        z[2]=(this->*foo2)(x);//Rprintf("z2 %f \n",z[2]);

        if (z[1]>z[2]) multF[ind_f+1] = 350.0; else multF[ind_f+1] = 0.0;

        //if (IT==0)  Rprintf("Mult %f \n",multF[ind_f+1]);

        }
        }

                //2ème étape : redéfinition de Ztemp

        SEXP listTemp;

        PROTECT(listTemp = duplicate(list));
        PROTECT(eVarCopy = duplicate(eVar));

        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
            g_nbdsFM[ind_f+nbF*1] = effortIni[ind_f] - g_nbdsFM[ind_f+nbF*0];
        }
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

        }

        int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,spQ))), "modI"));

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        for (int i = 0 ; i < NBI ; i++) {

            Rprintf("diffZZ %f ", REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);
            Ztemp[i+1] = Ztemp[i+1] + lambda*(REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);

        }

        //... et on recommence jusqu'à convergence de Z vers la valeur coïncidant avec les mortalités marginales

        UNPROTECT(2);

        }

        //3ème étape : redéfinition de PxQ

    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
            g_nbdsFM[ind_f+nbF*1] = effortIni[ind_f] - g_nbdsFM[ind_f+nbF*0];
        }
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

           if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }


    }

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        PROTECT(nDimF = allocVector(INTSXP,4));
        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

            //calcul de diffLQ
        double diffLQ = 0.0;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);

        }

        if ((DIFF*diffLQ)<0) lambda = lambda/3; //ie si DIFF et diffLQ de signe différent

        if ((diffLQ<=0) & (((diffLQ*lambda)*(diffLQ*lambda)<0.25) | (IT>ITquotaExch))) GoOn = false; //on ne s'arrête que si diffLQ<=0 (Quota respecté)

        IT++;

        DIFF = diffLQ;

//        Rprintf("\n");
        Rprintf("diffLQ %f \n",diffLQ);
//
        Rprintf("PxQ_1 %f ",PxQ);

        if (diffLQ<0) {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmin2(lambda*diffLQ,-0.1),pxQuMax));
        } else {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmax2(lambda*diffLQ,0.1),pxQuMax));
        }

        Rprintf("PxQ_2 %f \n",PxQ);
        Rprintf("lambda %f \n",lambda);
        //... et on recommence

        UNPROTECT(3);

    }


    double *pxQuot = REAL(VECTOR_ELT(out_PQuot_et,spQ));
    pxQuot[ind_t] = PxQ;

	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
            nbdsFM_G[ind_f+nbF*1] = effortIni[ind_f] - nbdsFM_G[ind_f+nbF*0];
        }
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

    }

    free_matrix(q,1,2,1,1);
    free_vector(multF,1,nbF+1);
    free_vector(z,1,2);
    free_vector(x,1,1);

    }

	return 0;
}

}

extern "C" {

double BioEcoPar::fxMaxProf_FT_customReportV2(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));

    double *gcfF;

    //for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            g_nbdsFM[IND_F+nbF*0] = fmin2(fmax2(x[1],0.0),effortIni[IND_F]);
            g_nbdsFM[IND_F+nbF*1] = effortIni[IND_F] - g_nbdsFM[IND_F+nbF*0];
        }
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

//    Rprintf("O %f \n",x[ind_f+1]);
//    Rprintf("A %f \n",g_nbdsFM[ind_f+nbF*0]);
//    Rprintf("B %f \n",g_nbdsFM[ind_f+nbF*1]);
//    Rprintf("C %f \n",g_nbdsF[ind_f]);
//    Rprintf("D %f \n",effortIni[ind_f]);
//    Rprintf("F %i \n",ind_f);

    //}

    int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,3);

    Marche(listTemp, IND_T);



    if (ecodcf==0) {

        Economic(listTemp,IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                           EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_Eco, 25));

    } else {

        EcoDCF(listTemp, IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                         EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_EcoDCF, 20));

    }

    //Rprintf("GCF %f \n",gcfF[IND_F + nbF*IND_T]);

    //calcul de l'indicateur de profit à optimiser
    //OUT@output$gcf_f[indF,indT]-pxQ[indT]*(sum(OUT@outputSp$Li$Sole_commune[indF,,,indT],na.rm=TRUE)-TAC_f_t[indF,1])

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

    double result = 0.0;
    result = gcfF[IND_F + nbF*IND_T]*g_nbvF[IND_F] - PxQ * (totF[IND_F + nbF*IND_T] - TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("A %f \n",TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("B %f \n",gcfF[IND_F + nbF*IND_T]);
    //Rprintf("C %f \n",g_nbvF[IND_F]);
    //Rprintf("D %f \n",totF[IND_F + nbF*IND_T]);

    UNPROTECT(3);

    return (-1*result);

}
}




//hypothese effort métier autre constant

extern "C" {

double BioEcoPar::fxTAC_F_customCst2(double *x) //cas métier Sole des flottilles modélisées seulement impacté
{
    SEXP listTemp, nDimF, nDim, nDimFM;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    if (IND_F < nbF) {

        if (var==1) g_nbdsFM[IND_F+nbF*0] = x[1];//fmin2(fmax2(x[1],0.0),fmin2(effortIni[IND_F]*expEff,350.0)-g_nbdsFM[IND_F+nbF*1]);
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

    } else {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*x[1],0.0);
            }
    }


    int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,0);

    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,3);

    PROTECT(nDimF = allocVector(INTSXP,4)); PROTECT(nDimFM = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    int *ndFM = INTEGER(nDimFM); ndFM[0] = nbF; ndFM[1] = nbMe; ndFM[2] = 0; ndFM[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimF));
    double *totFM = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimFM));
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    double *totMod = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDim));
    //PrintValue(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));

    //if (IND_F < nbF) Rprintf("totF %f TAC %f\n",totF[IND_F + nbF*IND_T],TAC_byFleet[IND_F + (nbF+1)*IND_T]);

    double result = 0.0;
    if (IND_F < nbF) {
            //result = result + fabs(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            result = (TAC_byFleet[IND_F + (nbF+1)*IND_T] - totFM[IND_F + nbF*1 + nbF*nbMe*IND_T])/totFM[IND_F + nbF*0 + nbF*nbMe*IND_T];
            result = fmin2(fmax2(result,0.0),fmin2(effortIni[IND_F]*expEff,350.0)-g_nbdsFM[IND_F+nbF*1]);
            //Rprintf("%12.6f ",result);
    } else {

    //result = result + fabs(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
            double TACoth = TAC_byFleet[nbF + (nbF+1)*IND_T];
            for (int ind_f = 0 ; ind_f < nbF ; ind_f++) TACoth = TACoth - TAC_byFleet[ind_f + (nbF+1)*IND_T];
            result = (tot[IND_T] - totMod[IND_T] - TACoth);//*(tot[IND_T] - totMod[IND_T] - TACoth);

    }
    //Rprintf("ccc");
    //Rprintf("%12.6f \n",result);
    //Rprintf("result %f x %f\n",result,x[1]);


    UNPROTECT(5);

    if (IND_F < nbF) return(result); else return (result*result);

}
}



extern "C" {

int BioEcoPar::GestionF2(int spp, int ind_t)
{

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;
    spQ = spp;

	int nfunc;
	int ITtot = tacIT;
	double lambda = tacLambda;
	double ftol = 0.00000001;


	//on déclare q et z nécessaire à la procédure d'optimisation

	double **q = NRmatrix(1,2,1,1);
	double *z = NRvector(1,2);
	double *x = NRvector(1,1);
    double *multF = NRvector(1,nbF+1);

	BEfn1_F foo2 = &BioEcoPar::fxTAC_F_customCst2;

    for (int IT = 0 ; IT < ITtot ; IT++){

        for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

        Rprintf("T %i F %i IT %i\n",ind_t,ind_f,IT);

        IND_F = ind_f;

        if (ind_f==nbF) { //optimisation de Foth --> x multiplicateur
            q[1][1]=x[1]=0.85;
            z[1]=(this->*foo2)(x);
            q[2][1]=x[1]=1;
            z[2]=(this->*foo2)(x);

            amoeba(foo2, q,z,1,ftol,&nfunc);

            multF[ind_f+1] = q[2][1];

        } else {            //optimisation de nbds --> x effort
            x[1]=1.0;
            multF[ind_f+1] = fxTAC_F_customCst2(x);
        }



//        Rprintf("Mult %f \n",multF[ind_f+1]);

        }

        //2ème étape : redéfinition de Ztemp

        SEXP listTemp;

        PROTECT(listTemp = duplicate(list));
        PROTECT(eVarCopy = duplicate(eVar));

        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),fmin2(effortIni[ind_f]*expEff,350.0)-g_nbdsFM[ind_f+nbF*1]);
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

        }

        int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,spQ))), "modI"));

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        for (int i = 0 ; i < NBI ; i++) {

            Rprintf("diffZZ %f ", REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);
            Ztemp[i+1] = Ztemp[i+1] + lambda*(REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);

        }

        //... et on recommence

        UNPROTECT(2);
    }


	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),fmin2(effortIni[ind_f]*expEff,350.0)-nbdsFM_G[ind_f+nbF*1]);
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }
    }


    free_matrix(q,1,2,1,1);
	free_vector(z,1,2);
	free_vector(x,1,1);
    free_vector(multF,1,nbF+1);


    }

	return 0;
}

}



//hypothese report d'effort

extern "C" {

double BioEcoPar::fxTAC_F_customReport2(double *x) //cas métier Sole des flottilles modélisées seulement impacté
{
    SEXP listTemp, nDimF, nDim;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    if (IND_F < nbF) {

        if (var==1) {
            g_nbdsFM[IND_F+nbF*0] = fmin2(fmax2(x[1],0.0),effortIni[IND_F]);
            g_nbdsFM[IND_F+nbF*1] = effortIni[IND_F] - g_nbdsFM[IND_F+nbF*0];
        }
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

    } else {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*x[1],0.0);
            }

    }



    int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,3);

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimF));
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    double *totMod = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDim));
    //PrintValue(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));

    double result = 0.0;
    if (IND_F < nbF) {
            //result = result + fabs(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            result = (totF[IND_F + nbF*IND_T]-TAC_byFleet[IND_F + (nbF+1)*IND_T])*(totF[IND_F + nbF*IND_T]-TAC_byFleet[IND_F + (nbF+1)*IND_T]);
            //Rprintf("%12.6f ",result);
    } else {

    //result = result + fabs(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
            double TACoth = TAC_byFleet[nbF + (nbF+1)*IND_T];
            for (int ind_f = 0 ; ind_f < nbF ; ind_f++) TACoth = TACoth - TAC_byFleet[ind_f + (nbF+1)*IND_T];
            result = (tot[IND_T] - totMod[IND_T] - TACoth)*(tot[IND_T] - totMod[IND_T] - TACoth);

    }
    //Rprintf("ccc");
    //Rprintf("%12.6f \n",result);
    //Rprintf("result %f \n",result);
    Rprintf("result %f x %f\n",result,x[1]);


    UNPROTECT(4);

    return result;

}
}



extern "C" {

int BioEcoPar::GestionF2report(int spp, int ind_t)
{

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;
    spQ = spp;

	int nfunc;
	int ITtot = 6;
	double lambda = 0.9;
	double ftol = 0.00000001;

    Rprintf("TIME %i \n",IND_T);
	//on déclare q et z nécessaire à la procédure d'optimisation

	double **q = NRmatrix(1,2,1,1);
	double *z = NRvector(1,2);
	double *x = NRvector(1,1);
    double *multF = NRvector(1,nbF+1);

	BEfn1_F foo2 = &BioEcoPar::fxTAC_F_customReport2;

    for (int IT = 0 ; IT < ITtot ; IT++){

        for (int ind_f = 0 ; ind_f <= nbF ; ind_f++){

        Rprintf("F %i \n",ind_f);

        IND_F = ind_f;

        if (ind_f==nbF) { //optimisation de Foth --> x multiplicateur
            q[1][1]=x[1]=0.85;
            z[1]=(this->*foo2)(x);
            q[2][1]=x[1]=1;
            z[2]=(this->*foo2)(x);
        } else {            //optimisation de nbds --> x effort
            q[1][1]=x[1]=REAL(getListElement(getListElement(list, "Fleet"), "nbds_f_m"))[IND_F+nbF*0] - 0.5;//0.0;
            z[1]=(this->*foo2)(x);
            q[2][1]=x[1]=REAL(getListElement(getListElement(list, "Fleet"), "nbds_f_m"))[IND_F+nbF*0] + 0.5;//1.0;
            z[2]=(this->*foo2)(x);
        }

        amoeba(foo2, q,z,1,ftol,&nfunc);

        multF[ind_f+1] = q[2][1];

//        Rprintf("Mult %f \n",multF[ind_f+1]);

        }

        //2ème étape : redéfinition de Ztemp

        SEXP listTemp;

        PROTECT(listTemp = duplicate(list));
        PROTECT(eVarCopy = duplicate(eVar));

        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
        if (var==1) g_nbdsFM[ind_f+nbF*1] = effortIni[ind_f] - g_nbdsFM[ind_f+nbF*0];
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];


        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

        }

        int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,spQ))), "modI"));

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        for (int i = 0 ; i < NBI ; i++) {

            Rprintf("diffZZ %f ", REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);
            Ztemp[i+1] = Ztemp[i+1] + lambda*(REAL(VECTOR_ELT(out_Z_eit,spQ))[i+NBI*IND_T] - Ztemp[i+1]);

        }

        //... et on recommence

        UNPROTECT(2);
    }


	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {

        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
        if (var==1) nbdsFM_G[ind_f+nbF*1] = effortIni[ind_f] - nbdsFM_G[ind_f+nbF*0];
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];


        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++) {

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*multF[nbF+1],0.0);
            }
         }

    }

    free_matrix(q,1,2,1,1);
	free_vector(z,1,2);
	free_vector(x,1,1);
    free_vector(multF,1,nbF+1);


    }

	return 0;
}

}




//------------------------------




extern "C" {

double BioEcoPar::fxTAC_F(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF, nDim;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsF[ind_f] = fmax2(g_nbdsF[ind_f]*x[ind_f+1],0.0);
        if (var==2) g_nbvF[ind_f] = fmax2(g_nbvF[ind_f]*x[ind_f+1],0.0);

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1) g_nbdsFM[ind_f+nbF*ind_m] = fmax2(g_nbdsFM[ind_f+nbF*ind_m]*x[ind_f+1],0.0);
            if (var==2) g_nbvFM[ind_f+nbF*ind_m] = fmax2(g_nbvFM[ind_f+nbF*ind_m]*x[ind_f+1],0.0);

         if (ind_m==0 & ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*x[nbF+1],0.0);
            }
         }
        }
    }


    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimF));
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    //PrintValue(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));

    double result = 0.0;
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
            //result = result + fabs(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            result = result + (totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T])*(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            //Rprintf("%12.6f ",result);
    }
    //result = result + fabs(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    result = result + (tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T])*(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    //Rprintf("result %f \n",result);
    //Rprintf("%12.6f \n",result);


    UNPROTECT(4);

    return result;

}
}

//Attention : fonctions customisées seulement valides pour le cas modèle individuel Quota Sole (2 métiers)

extern "C" {

double BioEcoPar::fxTAC_F_customCst(double *x) //cas métier Sole des flottilles modélisées seulement impacté
{
    SEXP listTemp, nDimF, nDim;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(g_nbdsFM[ind_f+nbF*0]*x[ind_f+1],0.0),effortIni[ind_f]-g_nbdsFM[ind_f+nbF*1]);
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*x[nbF+1],0.0);
            }
         }

    }

    int NBI = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,eTemp))), "modI"));

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);

    for (int i = 0 ; i < NBI ; i++)
        REAL(VECTOR_ELT(out_Z_eit,eTemp))[i+NBI*IND_T] = Ztemp[i+1];

    CatchDL(listTemp, IND_T, eVarCopy,3);

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimF));
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    //PrintValue(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));

    double result = 0.0;
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
            //result = result + fabs(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            result = result + (totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T])*(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            //Rprintf("%12.6f ",result);
    }
    //result = result + fabs(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    result = result + (tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T])*(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    //Rprintf("ccc");
    //Rprintf("%12.6f \n",result);
    Rprintf("result %f \n",result);


    UNPROTECT(4);

    return result;

}
}


extern "C" {

double BioEcoPar::fxTAC_F_customReport(double *x) //cas report d'effort
{
    SEXP listTemp, nDimF, nDim;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(g_nbdsFM[ind_f+nbF*0]*x[ind_f+1],0.0),effortIni[ind_f]);
            g_nbdsFM[ind_f+nbF*1] = effortIni[ind_f] - g_nbdsFM[ind_f+nbF*0];
        }
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];


        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *g_Fothi = REAL(VECTOR_ELT(VECTOR_ELT(eVarCopy, e), 44));
                int ni = length(getListElement(getListElement(listTemp, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) g_Fothi[ag + ni*IND_T] = fmax2(g_Fothi[ag + ni*IND_T]*x[nbF+1],0.0);
            }
         }

    }


    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, eTemp),nDimF));
    double *tot = REAL(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));
    //PrintValue(aggregObj(VECTOR_ELT(out_Y_eit, eTemp),nDim));

    double result = 0.0;
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
            //result = result + fabs(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            result = result + (totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T])*(totF[ind_f + nbF*IND_T]-TAC_byFleet[ind_f + (nbF+1)*IND_T]);
            //Rprintf("%12.6f ",result);
    }
    //result = result + fabs(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    result = result + (tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T])*(tot[IND_T]-TAC_byFleet[nbF + (nbF+1)*IND_T]);
    //Rprintf("ccc");
    //Rprintf("%12.6f \n",result);


    UNPROTECT(4);

    return result;

}
}


extern "C" {

int BioEcoPar::MinimizeF(double **p, double y[], int ndim, double ftol)
{
	int i,nfunc,j;//,ndim=3;
	//float *x,*y,**p;

	double *x=NRvector(1,NP);
	//y=NRvector(1,4);
	//p=NRmatrix(1,4,1,3);

	BEfn1_F foo = &BioEcoPar::func;

	for (i=1;i<=MP;i++) {
		for (j=1;j<=NP;j++)
			x[j]=p[i][j]=(i == (j+1) ? 1.0 : 0.0);
		y[i]=(this->*foo)(x);
	}

	amoeba(foo, p,y,ndim,ftol,&nfunc);
	Rprintf("\nNumber of function evaluations: %3d\n",nfunc);
	Rprintf("Vertices of final 3-d simplex and\n");
	Rprintf("function values at the vertices:\n\n");
	Rprintf("%3s %10s %12s %12s %14s\n\n",
		"i","x[i]","y[i]","z[i]","function");
	for (i=1;i<=MP;i++) {
		Rprintf("%3d ",i);
		for (j=1;j<=NP;j++) Rprintf("%12.6f ",p[i][j]);
		Rprintf("%12.6f\n",y[i]);
	}
	Rprintf("\nTrue minimum is at (0.5,0.6,0.7)\n");
	//free_matrix(p,1,MP,1,NP);
	//free_vector(y,1,MP);
	free_vector(x,1,NP);
	return 0;
}

}

//original 'GestionF' à conserver

//extern "C" {
//
//int BioEcoPar::GestionF(double **p, double y[], int ndim, double ftol, int ind_t)
//{
//
//    if (ind_t<delay) {
//
//    } else {
//
//    IND_T = ind_t;
//
//	int i,nfunc,j;
//
//	double *x=NRvector(1,ndim);
//	//y=NRvector(1,4);
//	//p=NRmatrix(1,4,1,3);
//
//	BEfn1_F foo = &BioEcoPar::fxTAC_F;
//
//	for (i=1;i<=(ndim+1);i++) {
//		for (j=1;j<=ndim;j++)
//			x[j]=p[i][j]=1+(i == (j+1) ? 0.1 : -0.1);
//		y[i]=(this->*foo)(x);
//	}
//
//	amoeba(foo, p,y,ndim,ftol,&nfunc);
//	//Rprintf("\nNumber of function evaluations: %3d\n",nfunc);
//	//Rprintf("Vertices of final 3-d simplex and\n");
//	//Rprintf("function values at the vertices:\n\n");
//	//Rprintf("%3s %10s %12s %12s %14s\n\n",
//	//	"i","x[i]","y[i]","z[i]","function");
//	//for (i=1;i<=(ndim+1);i++) {
//	//	Rprintf("%3d ",i);
//	//	for (j=1;j<=ndim;j++) Rprintf("%12.6f ",p[i][j]);
//	//	Rprintf("%12.6f\n",y[i]);
//	//}
//	//Rprintf("\nTrue minimum is at (0.5,0.6,0.7)\n");
//
//	free_vector(x,1,ndim);
//
//	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
//    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));
//    double *nbvFM_G = REAL(getListElement(FList, "nbv_f_m"));
//    double *nbvF_G = REAL(getListElement(FList, "nbv_f"));
//
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        if (var==1) nbdsF_G[ind_f] = fmax2(nbdsF_G[ind_f]*p[nbF+2][ind_f+1],0.0);
//        if (var==2) nbvF_G[ind_f] = fmax2(nbvF_G[ind_f]*p[nbF+2][ind_f+1],0.0);
//
//        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//
//            if (var==1) nbdsFM_G[ind_f+nbF*ind_m] = fmax2(nbdsFM_G[ind_f+nbF*ind_m]*p[nbF+2][ind_f+1],0.0);
//            if (var==2) nbvFM_G[ind_f+nbF*ind_m] = fmax2(nbvFM_G[ind_f+nbF*ind_m]*p[nbF+2][ind_f+1],0.0);
//
//         if (ind_m==0 & ind_f==0) {
//
//            for (int e = 0 ; e < nbE ; e++){
//
//                double *Fothi_G = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
//                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
//                for (int ag = 0; ag < ni; ag++) Fothi_G[ag + ni*IND_T] = fmax2(Fothi_G[ag + ni*IND_T]*p[nbF+2][nbF+1],0.0);
//            }
//         }
//        }
//    }
//    }
//	return 0;
//}
//
//}


//custom Cst

extern "C" {

int BioEcoPar::GestionF(double **p, double y[], int ndim, double ftol, int ind_t)
{

    if (ind_t<delay) {

    } else {

    IND_T = ind_t;

	int i,nfunc,j;

	double *x=NRvector(1,ndim);
	//y=NRvector(1,4);
	//p=NRmatrix(1,4,1,3);

	BEfn1_F foo = &BioEcoPar::fxTAC_F_customCst;

	for (i=1;i<=(ndim+1);i++) {
		for (j=1;j<=ndim;j++)
			x[j]=p[i][j]=1+(i == (j+1) ? 0.1 : -0.1);
		y[i]=(this->*foo)(x);
	}

	amoeba(foo, p,y,ndim,ftol,&nfunc);

    Rprintf("Time %i\n",ind_t);
	Rprintf("END--------------------------------\n");
	//Rprintf("\nNumber of function evaluations: %3d\n",nfunc);
	//Rprintf("Vertices of final 3-d simplex and\n");
	//Rprintf("function values at the vertices:\n\n");
	//Rprintf("%3s %10s %12s %12s %14s\n\n",
	//	"i","x[i]","y[i]","z[i]","function");
	//for (i=1;i<=(ndim+1);i++) {
	//	Rprintf("%3d ",i);
	//	for (j=1;j<=ndim;j++) Rprintf("%12.6f ",p[i][j]);
	//	Rprintf("%12.6f\n",y[i]);
	//}
	//Rprintf("\nTrue minimum is at (0.5,0.6,0.7)\n");

	free_vector(x,1,ndim);

    double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(nbdsFM_G[ind_f+nbF*0]*p[nbF+2][ind_f+1],0.0),effortIni[ind_f]-nbdsFM_G[ind_f+nbF*1]);
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];

        if (ind_f==0) {

            for (int e = 0 ; e < nbE ; e++){

                double *Fothi_G = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
                for (int ag = 0; ag < ni; ag++) Fothi_G[ag + ni*IND_T] = fmax2(Fothi_G[ag + ni*IND_T]*p[nbF+2][nbF+1],0.0);
            }
        }

    }

    }

	return 0;
}

}



//custom Report

//extern "C" {
//
//int BioEcoPar::GestionF(double **p, double y[], int ndim, double ftol, int ind_t)
//{
//
//    if (ind_t<delay) {
//
//    } else {
//
//    IND_T = ind_t;
//
//	int i,nfunc,j;
//
//	double *x=NRvector(1,ndim);
//	//y=NRvector(1,4);
//	//p=NRmatrix(1,4,1,3);
//
//	BEfn1_F foo = &BioEcoPar::fxTAC_F_customReport;
//
//	for (i=1;i<=(ndim+1);i++) {
//		for (j=1;j<=ndim;j++)
//			x[j]=p[i][j]=1+(i == (j+1) ? 0.1 : -0.1);
//		y[i]=(this->*foo)(x);
//	}
//
//	amoeba(foo, p,y,ndim,ftol,&nfunc);
//	//Rprintf("\nNumber of function evaluations: %3d\n",nfunc);
//	//Rprintf("Vertices of final 3-d simplex and\n");
//	//Rprintf("function values at the vertices:\n\n");
//	//Rprintf("%3s %10s %12s %12s %14s\n\n",
//	//	"i","x[i]","y[i]","z[i]","function");
//	//for (i=1;i<=(ndim+1);i++) {
//	//	Rprintf("%3d ",i);
//	//	for (j=1;j<=ndim;j++) Rprintf("%12.6f ",p[i][j]);
//	//	Rprintf("%12.6f\n",y[i]);
//	//}
//	//Rprintf("\nTrue minimum is at (0.5,0.6,0.7)\n");
//
//	free_vector(x,1,ndim);
//
//    double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
//    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        if (var==1) {
//            nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(nbdsFM_G[ind_f+nbF*0]*p[nbF+2][ind_f+1],0.0),effortIni[ind_f]);
//            nbdsFM_G[ind_f+nbF*1] = effortIni[ind_f] - nbdsFM_G[ind_f+nbF*0];
//        }
//        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];
//
//
//        if (ind_f==0) {
//
//            for (int e = 0 ; e < nbE ; e++){
//
//                double *Fothi_G = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e), 44));
//                int ni = length(getListElement(getListElement(list, CHAR(STRING_ELT(sppList,e))), "modI"));
//                for (int ag = 0; ag < ni; ag++) Fothi_G[ag + ni*IND_T] = fmax2(Fothi_G[ag + ni*IND_T]*p[nbF+2][nbF+1],0.0);
//            }
//        }
//    }
//
//    }
//
//	return 0;
//}
//
//}


/* (C) Copr. 1986-92 Numerical Recipes Software *pA24. */
//------------------------------------------------------------------------------------------------------




//------------------------------------------
// Module 'Report d'effort' selon une pondération des ratio profit par métier et effort par métier anticipés
//------------------------------------------

extern "C" {

void BioEcoPar::FleetBehav(SEXP list, int ind_t, SEXP paramBehav) //ind_t>0
{

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

 //type n°1 : pas de report d'effort. Intervention sur l'effort au niveau flottille-métier via la matrice FMT
 // qui opère additivement, avec redressement en cas d'effort résultant négatif ou supérieur à 365 sommé sur les métiers
 // L'effort au niveau flottille est ensuite réévalué par agrégation du niveau flottille-métier

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


 //type n°2 : reports d'effort pilotés. Intervention sur les métiers par flottille avec report conditionné par une matrice FMT
 // de type :   | xx  xx   1 -0.5 -0.5   xx |
 //             | xx 0.7 0.3   xx -0.2 -0.8 |
 //             | ...                       |
 //
 // La quantité brute de report par flottille-métier est ensuite évaluée par multiplication de FMT par un vecteur MU de dimension nbF
 // MU est contraint pour que les reports soient cohérents


    if (typeBhv==2 & (fmtBhv != NULL) & (muBhv != NULL)) {

       double *r_fmtBhv = REAL(fmtBhv), *r_muBhv = REAL(muBhv);

       for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

           //détermination de la validité de MU_f et correction le cas échéant

        double mu_limSup=-1.0, mu_limInf=0.0, finalMu=0.0;

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (!ISNA(r_nbds_f_m[ind_f + nbF*ind_m])) {

                if (mu_limSup<0) { //première évaluation
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

        //normalement, pas besoin de réévaluer nbds_f car la conservation de l'effort est assurée par la méthodo
       }
    }




 //type n°3 : report d'effort orienté par pondération des ratio de profit et d'effort de l'année précédente (cf P. Marchal).

    if (typeBhv==3 & ind_t>0 & (alphaBhv != NULL)) {

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





//------------------------------------------
// Module 'Echange de quotas'
//------------------------------------------

extern "C" {

double BioEcoPar::fxMaxProf_FT(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));

    double initVal = 0.0;
    double *gcfF;

        if (var==1 & g_nbdsF[IND_F]>0) {initVal = g_nbdsF[IND_F]; g_nbdsF[IND_F] = fmin2(fmax2(x[1],1),350);}
        if (var==2 & g_nbvF[IND_F]>0) {initVal = g_nbvF[IND_F]; g_nbvF[IND_F] = fmin2(fmax2(x[1],1),350);}

        for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

            if (var==1 & initVal>0) g_nbdsFM[IND_F+nbF*ind_m] = g_nbdsFM[IND_F+nbF*ind_m]*g_nbdsF[IND_F]/initVal;
            if (var==2 & initVal>0) g_nbvFM[IND_F+nbF*ind_m] = g_nbvFM[IND_F+nbF*ind_m]*g_nbvF[IND_F]/initVal;

        }

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

    Marche(listTemp, IND_T);



    if (ecodcf==0) {

        Economic(listTemp,IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                           EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_Eco, 25));

    } else {

        EcoDCF(listTemp, IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                         EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_EcoDCF, 20));

    }

    //Rprintf("GCF %f \n",gcfF[IND_F + nbF*IND_T]);

    //calcul de l'indicateur de profit à optimiser
    //OUT@output$gcf_f[indF,indT]-pxQ[indT]*(sum(OUT@outputSp$Li$Sole_commune[indF,,,indT],na.rm=TRUE)-TAC_f_t[indF,1])

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

    double result = 0.0;
    result = gcfF[IND_F + nbF*IND_T]*g_nbvF[IND_F] - PxQ * (totF[IND_F + nbF*IND_T] - TAC_byFleet[IND_F + (nbF+1)*IND_T]);


    UNPROTECT(3);

    return (-1*result);

}
}


extern "C" {

double BioEcoPar::fxMaxProf_FT_customCst(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));

    double *gcfF;

    //for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[IND_F+nbF*0] = fmin2(fmax2(x[1],0.0),effortIni[IND_F]-g_nbdsFM[IND_F+nbF*1]);
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

//    Rprintf("O %f \n",x[ind_f+1]);
//    Rprintf("A %f \n",g_nbdsFM[ind_f+nbF*0]);
//    Rprintf("B %f \n",g_nbdsFM[ind_f+nbF*1]);
//    Rprintf("C %f \n",g_nbdsF[ind_f]);
//    Rprintf("D %f \n",effortIni[ind_f]);
//    Rprintf("F %i \n",ind_f);

    //}

    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

    Marche(listTemp, IND_T);



    if (ecodcf==0) {

        Economic(listTemp,IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                           EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_Eco, 25));

    } else {

        EcoDCF(listTemp, IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                         EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_EcoDCF, 20));

    }

    //Rprintf("GCF %f \n",gcfF[IND_F + nbF*IND_T]);

    //calcul de l'indicateur de profit à optimiser
    //OUT@output$gcf_f[indF,indT]-pxQ[indT]*(sum(OUT@outputSp$Li$Sole_commune[indF,,,indT],na.rm=TRUE)-TAC_f_t[indF,1])

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

    double result = 0.0;
    result = gcfF[IND_F + nbF*IND_T]*g_nbvF[IND_F] - PxQ * (totF[IND_F + nbF*IND_T] - TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("A %f \n",TAC_byFleet[IND_F + (nbF+1)*IND_T]);
    //Rprintf("B %f \n",gcfF[IND_F + nbF*IND_T]);
    //Rprintf("C %f \n",g_nbvF[IND_F]);
    //Rprintf("D %f \n",totF[IND_F + nbF*IND_T]);

    UNPROTECT(3);

    return (-1*result);

}
}




extern "C" {

double BioEcoPar::fxMaxProf_FT_customReport(double *x) //attention : l'indexation de x commence à 1 et non 0
{
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
    double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));

    double *gcfF;

    //for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) {
            g_nbdsFM[IND_F+nbF*0] = fmin2(fmax2(x[1],0.0),effortIni[IND_F]);
            g_nbdsFM[IND_F+nbF*1] = effortIni[IND_F] - g_nbdsFM[IND_F+nbF*0];
        }
        if (var==1) g_nbdsF[IND_F] = g_nbdsFM[IND_F+nbF*0] + g_nbdsFM[IND_F+nbF*1];

    //}


    Mortalite(listTemp, IND_T, eVarCopy,0);
    DynamicPop(listTemp, IND_T, eVarCopy,0);
    CatchDL(listTemp, IND_T, eVarCopy,0);
    Mortalite(listTemp, IND_T, eVarCopy,1);
    DynamicPop(listTemp, IND_T, eVarCopy,1);
    CatchDL(listTemp, IND_T, eVarCopy,1);
    Mortalite(listTemp, IND_T, eVarCopy,2);
    DynamicPop(listTemp, IND_T, eVarCopy,2);
    CatchDL(listTemp, IND_T, eVarCopy,2);
    Mortalite(listTemp, IND_T, eVarCopy,3);
    DynamicPop(listTemp, IND_T, eVarCopy,3);
    CatchDL(listTemp, IND_T, eVarCopy,3);

    Marche(listTemp, IND_T);



    if (ecodcf==0) {

        Economic(listTemp,IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                           EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_Eco, 25));

    } else {

        EcoDCF(listTemp, IND_T, EcoIndCopy[0], EcoIndCopy[1], EcoIndCopy[2], EcoIndCopy[3], EcoIndCopy[4],
                         EcoIndCopy[5], EcoIndCopy[6], drCopy);

        gcfF = REAL(VECTOR_ELT(out_EcoDCF, 20));

    }

    //Rprintf("GCF %f \n",gcfF[IND_F + nbF*IND_T]);

    //calcul de l'indicateur de profit à optimiser
    //OUT@output$gcf_f[indF,indT]-pxQ[indT]*(sum(OUT@outputSp$Li$Sole_commune[indF,,,indT],na.rm=TRUE)-TAC_f_t[indF,1])

    PROTECT(nDimF = allocVector(INTSXP,4));
    int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
    double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

    double result = 0.0;
    result = gcfF[IND_F + nbF*IND_T]*g_nbvF[IND_F] - PxQ * (totF[IND_F + nbF*IND_T] - TAC_byFleet[IND_F + (nbF+1)*IND_T]);


    UNPROTECT(3);

    return (-1*result);

}
}




//extern "C" {
//
//int BioEcoPar::QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
//{
//
//    if (ind_t<delay) {
//
//    } else {
//
//    IND_T = ind_t;
//    PxQ = pxQuIni;
//    spQ = spp;
//
//	int nfunc;
//	int ITquotaExch = 20;
//
//    Rprintf("TIME %i \n",IND_T);
//	//on déclare q et z nécessaire à la procédure d'optimisation
//
//	double **q = NRmatrix(1,2,1,1);
//	double *z = NRvector(1,2);
//	double *x = NRvector(1,1);
//    double *multF = NRvector(1,nbF);
//
//	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)
//
//    for (int IT = 0 ; IT < ITquotaExch ; IT++){
//
//        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
//    //Rprintf("IT %i \n",IT);
//
//    if (IT==5) lambda = lambda/6;
//    if (IT==10) lambda = lambda/6;
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        IND_F = ind_f;
//
////        x[1]=1;
////        Rprintf("Val : %f ",(this->*foo2)(x));
////        x[1]=50;
////        Rprintf("%f ",(this->*foo2)(x));
////        x[1]=125;
////        Rprintf("%f ",(this->*foo2)(x));
////        x[1]=200;
////        Rprintf("%f ",(this->*foo2)(x));
////        x[1]=350;
////        Rprintf("%f \n",(this->*foo2)(x));
//
//        q[1][1]=x[1]=175;
//        z[1]=(this->*foo2)(x);
//        q[2][1]=x[1]=120;
//        z[2]=(this->*foo2)(x);
//
//        amoeba(foo2, q,z,1,ftol,&nfunc);
//
//        multF[ind_f+1] = q[2][1];
//
////        Rprintf("Mult %f \n",multF[ind_f+1]);
//
//        }
//
//        //2ème étape : redéfinition de PxQ
//
//        SEXP listTemp, nDimF;
//
//        PROTECT(listTemp = duplicate(list));
//        PROTECT(eVarCopy = duplicate(eVar));
//
//        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
//        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
//        double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));
//        double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));
//
//        double initVal = 0.0;
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//            if (var==1 & g_nbdsF[ind_f]>0) {initVal = g_nbdsF[ind_f] ; g_nbdsF[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//            if (var==2 & g_nbvF[ind_f]>0) {initVal = g_nbvF[ind_f] ; g_nbvF[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//
//            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//
//                if (var==1 & initVal>0) g_nbdsFM[ind_f+nbF*ind_m] = g_nbdsFM[ind_f+nbF*ind_m]*g_nbdsF[ind_f]/initVal;
//                if (var==2 & initVal>0) g_nbvFM[ind_f+nbF*ind_m] = g_nbvFM[ind_f+nbF*ind_m]*g_nbvF[ind_f]/initVal;
//
//            }
//        }
//
//        Mortalite(listTemp, IND_T, eVarCopy,0);
//        DynamicPop(listTemp, IND_T, eVarCopy,0);
//        CatchDL(listTemp, IND_T, eVarCopy,0);
//        Mortalite(listTemp, IND_T, eVarCopy,1);
//        DynamicPop(listTemp, IND_T, eVarCopy,1);
//        CatchDL(listTemp, IND_T, eVarCopy,1);
//        Mortalite(listTemp, IND_T, eVarCopy,2);
//        DynamicPop(listTemp, IND_T, eVarCopy,2);
//        CatchDL(listTemp, IND_T, eVarCopy,2);
//        Mortalite(listTemp, IND_T, eVarCopy,3);
//        DynamicPop(listTemp, IND_T, eVarCopy,3);
//        CatchDL(listTemp, IND_T, eVarCopy,3);
//
//        PROTECT(nDimF = allocVector(INTSXP,4));
//        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
//        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));
//
//            //calcul de diffLQ
//        double diffLQ = 0.0;
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);
//
//        }
//
////        Rprintf("\n");
//        Rprintf("diffLQ %f \n",diffLQ);
////
//        Rprintf("PxQ_1 %f ",PxQ);
//
//        PxQ = fmax2(pxQuMin,fmin2(PxQ + lambda*diffLQ,pxQuMax));
//
//        Rprintf("PxQ_2 %f \n",PxQ);
//        //... et on recommence
//
//        UNPROTECT(3);
//    }
//
//
//	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
//    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));
//    double *nbvFM_G = REAL(getListElement(FList, "nbv_f_m"));
//    double *nbvF_G = REAL(getListElement(FList, "nbv_f"));
//
//    double initVal2 = 0.0;
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//            if (var==1 & nbdsF_G[ind_f]>0) {initVal2 = nbdsF_G[ind_f] ; nbdsF_G[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//            if (var==2 & nbvF_G[ind_f]>0) {initVal2 = nbvF_G[ind_f] ; nbvF_G[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//
//            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//
//                if (var==1 & initVal2>0) nbdsFM_G[ind_f+nbF*ind_m] = nbdsFM_G[ind_f+nbF*ind_m]*nbdsF_G[ind_f]/initVal2;
//                if (var==2 & initVal2>0) nbvFM_G[ind_f+nbF*ind_m] = nbvFM_G[ind_f+nbF*ind_m]*nbvF_G[ind_f]/initVal2;
//
//            }
//        }
//
//    free_matrix(q,1,2,1,1);
//	free_vector(z,1,2);
//	free_vector(x,1,1);
//    free_vector(multF,1,nbF);
//
//
//    }
//
//	return 0;
//}
//
//}



//
////version sans algorithme d'optimisation (hypothèse de fonction de cout linéaire --> solutions en coin)
//
//extern "C" {
//
//int BioEcoPar::QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
//{
//
//    if (ind_t<delay) {
//
//    } else {
//
//    IND_T = ind_t;
//    PxQ = pxQuIni;
//    spQ = spp;
//
//	int ITquotaExch = 30;
//
//    Rprintf("TIME %i \n",IND_T);
//
//    double *multF = NRvector(1,nbF);
//    double *z = NRvector(1,2);
//    double *x = NRvector(1,1);
//
//    bool GoOn = true;
//
//	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)
//
//    double DIFF = 0.0;
//    int IT = 0;
//
//    while (GoOn) {
//
//    //for (int IT = 0 ; IT < ITquotaExch ; IT++){
//
//        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
//    //Rprintf("IT %i \n",IT);
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        IND_F = ind_f;
//
//        x[1]=1.0;
//        z[1]=(this->*foo2)(x);Rprintf("z1 %f \n",z[1]);
//        x[1]=350.0;
//        z[2]=(this->*foo2)(x);Rprintf("z2 %f \n",z[2]);
//
//        if (z[1]>z[2]) multF[ind_f+1] = 350.0; else multF[ind_f+1] = 1.0;
//
//        //if (IT==0)  Rprintf("Mult %f \n",multF[ind_f+1]);
//
//        }
//
//        //2ème étape : redéfinition de PxQ
//
//        SEXP listTemp, nDimF;
//
//        PROTECT(listTemp = duplicate(list));
//        PROTECT(eVarCopy = duplicate(eVar));
//
//        double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
//        double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
//        double *g_nbvFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f_m"));
//        double *g_nbvF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbv_f"));
//
//        double initVal = 0.0;
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
//
//            if (var==1 & g_nbdsF[ind_f]>0) {initVal = g_nbdsF[ind_f] ; g_nbdsF[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//            if (var==2 & g_nbvF[ind_f]>0) {initVal = g_nbvF[ind_f] ; g_nbvF[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//
//            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//
//                if (var==1 & initVal>0) g_nbdsFM[ind_f+nbF*ind_m] = g_nbdsFM[ind_f+nbF*ind_m]*g_nbdsF[ind_f]/initVal;
//                if (var==2 & initVal>0) g_nbvFM[ind_f+nbF*ind_m] = g_nbvFM[ind_f+nbF*ind_m]*g_nbvF[ind_f]/initVal;
//
//            }
//        }
//
//        Mortalite(listTemp, IND_T, eVarCopy,0);
//        DynamicPop(listTemp, IND_T, eVarCopy,0);
//        CatchDL(listTemp, IND_T, eVarCopy,0);
//        Mortalite(listTemp, IND_T, eVarCopy,1);
//        DynamicPop(listTemp, IND_T, eVarCopy,1);
//        CatchDL(listTemp, IND_T, eVarCopy,1);
//        Mortalite(listTemp, IND_T, eVarCopy,2);
//        DynamicPop(listTemp, IND_T, eVarCopy,2);
//        CatchDL(listTemp, IND_T, eVarCopy,2);
//        Mortalite(listTemp, IND_T, eVarCopy,3);
//        DynamicPop(listTemp, IND_T, eVarCopy,3);
//        CatchDL(listTemp, IND_T, eVarCopy,3);
//
//        PROTECT(nDimF = allocVector(INTSXP,4));
//        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
//        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));
//
//            //calcul de diffLQ
//        double diffLQ = 0.0;
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);
//
//        }
//
//        if ((DIFF*diffLQ)<0) lambda = lambda/3; //ie si DIFF et diffLQ de signe différent
//
//        if ((diffLQ<=0) & ((lambda<0.002) | (IT>ITquotaExch))) GoOn = false; //on ne s'arrête que si diffLQ<=0 (Quota respecté)
//
//        IT++;
//
//        DIFF = diffLQ;
//
////        Rprintf("\n");
//        Rprintf("diffLQ %f \n",diffLQ);
////
//        Rprintf("PxQ_1 %f ",PxQ);
//
//        PxQ = fmax2(pxQuMin,fmin2(PxQ + lambda*diffLQ,pxQuMax));
//
//        Rprintf("PxQ_2 %f \n",PxQ);
//        Rprintf("lambda %f \n",lambda);
//        //... et on recommence
//
//        UNPROTECT(3);
//
//    }
//
//
//    double *pxQuot = REAL(VECTOR_ELT(out_PQuot_et,spQ));
//    pxQuot[ind_t] = PxQ;
//
//	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
//    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));
//    double *nbvFM_G = REAL(getListElement(FList, "nbv_f_m"));
//    double *nbvF_G = REAL(getListElement(FList, "nbv_f"));
//
//    double initVal2 = 0.0;
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//            if (var==1 & nbdsF_G[ind_f]>0) {initVal2 = nbdsF_G[ind_f] ; nbdsF_G[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//            if (var==2 & nbvF_G[ind_f]>0) {initVal2 = nbvF_G[ind_f] ; nbvF_G[ind_f] = fmin2(fmax2(multF[ind_f+1],1),350);}
//
//            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//
//                if (var==1 & initVal2>0) nbdsFM_G[ind_f+nbF*ind_m] = nbdsFM_G[ind_f+nbF*ind_m]*nbdsF_G[ind_f]/initVal2;
//                if (var==2 & initVal2>0) nbvFM_G[ind_f+nbF*ind_m] = nbvFM_G[ind_f+nbF*ind_m]*nbvF_G[ind_f]/initVal2;
//
//            }
//        }
//
//    free_vector(multF,1,nbF);
//    free_vector(z,1,2);
//    free_vector(x,1,1);
//
//    }
//
//	return 0;
//}
//
//}

//---------------------------------
//
//
//
//---------------------------------

//Version customCst

extern "C" {

int BioEcoPar::QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
{

    if (ind_t<delay) {

    } else {

    FOTHoptim_use = true;

    IND_T = ind_t;
    PxQ = pxQuIni;
    spQ = spp;

	int ITquotaExch = 20;

    Rprintf("TIME %i \n",IND_T);

    double *multF = NRvector(1,nbF);
    double *z = NRvector(1,2);
    double *x = NRvector(1,1);

    bool GoOn = true;

	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT_customCst; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)

    double DIFF = 0.0;
    int IT = 0;

    while (GoOn) {

    //for (int IT = 0 ; IT < ITquotaExch ; IT++){

        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
    //Rprintf("IT %i \n",IT);

        Zoptim_use = true;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        IND_F = ind_f;

        x[1]=1.0;
        z[1]=(this->*foo2)(x);//Rprintf("z1 %f \n",z[1]);
        x[1]=350.0;
        z[2]=(this->*foo2)(x);//Rprintf("z2 %f \n",z[2]);

        if (z[1]>z[2]) multF[ind_f+1] = 350.0; else multF[ind_f+1] = 0.0;

        //if (IT==0)  Rprintf("Mult %f \n",multF[ind_f+1]);

        }

        //2ème étape : redéfinition de PxQ
    Zoptim_use = false;
    SEXP listTemp, nDimF;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]-g_nbdsFM[ind_f+nbF*1]);
        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];

    }

        Mortalite(listTemp, IND_T, eVarCopy,0);
        DynamicPop(listTemp, IND_T, eVarCopy,0);
        CatchDL(listTemp, IND_T, eVarCopy,0);
        Mortalite(listTemp, IND_T, eVarCopy,1);
        DynamicPop(listTemp, IND_T, eVarCopy,1);
        CatchDL(listTemp, IND_T, eVarCopy,1);
        Mortalite(listTemp, IND_T, eVarCopy,2);
        DynamicPop(listTemp, IND_T, eVarCopy,2);
        CatchDL(listTemp, IND_T, eVarCopy,2);
        Mortalite(listTemp, IND_T, eVarCopy,3);
        DynamicPop(listTemp, IND_T, eVarCopy,3);
        CatchDL(listTemp, IND_T, eVarCopy,3);

        PROTECT(nDimF = allocVector(INTSXP,4));
        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));

            //calcul de diffLQ
        double diffLQ = 0.0;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);

        }

        if ((DIFF*diffLQ)<0) lambda = lambda/3; //ie si DIFF et diffLQ de signe différent

        if ((diffLQ<=0) & (((diffLQ*lambda)*(diffLQ*lambda)<0.25) | (IT>ITquotaExch))) GoOn = false; //on ne s'arrête que si diffLQ<=0 (Quota respecté)

        IT++;

        DIFF = diffLQ;

//        Rprintf("\n");
        Rprintf("diffLQ %f \n",diffLQ);
//
        Rprintf("PxQ_1 %f ",PxQ);

        if (diffLQ<0) {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmin2(lambda*diffLQ,-0.1),pxQuMax));
        } else {
        PxQ = fmax2(pxQuMin,fmin2(PxQ + fmax2(lambda*diffLQ,0.1),pxQuMax));
        }

        Rprintf("PxQ_2 %f \n",PxQ);
        Rprintf("lambda %f \n",lambda);
        //... et on recommence

        UNPROTECT(3);

    }


    double *pxQuot = REAL(VECTOR_ELT(out_PQuot_et,spQ));
    pxQuot[ind_t] = PxQ;

	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]-nbdsFM_G[ind_f+nbF*1]);
        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];

    }

    free_vector(multF,1,nbF);
    free_vector(z,1,2);
    free_vector(x,1,1);

    }

	return 0;
}

}



//
////Version customReport
//
//extern "C" {
//
//int BioEcoPar::QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t)
//{
//
//    if (ind_t<delay) {
//
//    } else {
//
//    FOTHoptim_use = true;
//
//    IND_T = ind_t;
//    PxQ = pxQuIni;
//    spQ = spp;
//
//	int ITquotaExch = 20;
//
//    Rprintf("TIME %i \n",IND_T);
//
//    double *multF = NRvector(1,nbF);
//    double *z = NRvector(1,2);
//    double *x = NRvector(1,1);
//
//    bool GoOn = true;
//
//	BEfn1_F foo2 = &BioEcoPar::fxMaxProf_FT_customReport; //fonction calculant le profit d'une flottille donnée à un instant donné en fonction du prix du quota d'une espèce donnée (dépend de IND_F)
//
//    double DIFF = 0.0;
//    int IT = 0;
//
//    while (GoOn) {
//
//    //for (int IT = 0 ; IT < ITquotaExch ; IT++){
//
//        //1ère étape : maximisation du profit par l'effort par flottille en fonction de PxQ et spQ
//    //Rprintf("IT %i \n",IT);
//
//        Zoptim_use = true;
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        IND_F = ind_f;
//
//        x[1]=1.0;
//        z[1]=(this->*foo2)(x);Rprintf("z1 %f \n",z[1]);
//        x[1]=350.0;
//        z[2]=(this->*foo2)(x);Rprintf("z2 %f \n",z[2]);
//
//        if (z[1]>z[2]) multF[ind_f+1] = 350.0; else multF[ind_f+1] = 0.0;
//
//        //if (IT==0)  Rprintf("Mult %f \n",multF[ind_f+1]);
//
//        }
//
//        //2ème étape : redéfinition de PxQ
//
//      Zoptim_use = false;
//        SEXP listTemp, nDimF;
//
//    PROTECT(listTemp = duplicate(list));
//    PROTECT(eVarCopy = duplicate(eVar));
//
//    double *g_nbdsFM = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f_m"));
//    double *g_nbdsF = REAL(getListElement(getListElement(listTemp, "Fleet"), "nbds_f"));
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        if (var==1) {
//            g_nbdsFM[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]);
//            g_nbdsFM[ind_f+nbF*1] = effortIni[ind_f] - g_nbdsFM[ind_f+nbF*0];
//        }
//        if (var==1) g_nbdsF[ind_f] = g_nbdsFM[ind_f+nbF*0] + g_nbdsFM[ind_f+nbF*1];
//
//    }
//
//
//        Mortalite(listTemp, IND_T, eVarCopy,0);
//        DynamicPop(listTemp, IND_T, eVarCopy,0);
//        CatchDL(listTemp, IND_T, eVarCopy,0);
//        Mortalite(listTemp, IND_T, eVarCopy,1);
//        DynamicPop(listTemp, IND_T, eVarCopy,1);
//        CatchDL(listTemp, IND_T, eVarCopy,1);
//        Mortalite(listTemp, IND_T, eVarCopy,2);
//        DynamicPop(listTemp, IND_T, eVarCopy,2);
//        CatchDL(listTemp, IND_T, eVarCopy,2);
//        Mortalite(listTemp, IND_T, eVarCopy,3);
//        DynamicPop(listTemp, IND_T, eVarCopy,3);
//        CatchDL(listTemp, IND_T, eVarCopy,3);
//
//        PROTECT(nDimF = allocVector(INTSXP,4));
//        int *ndF = INTEGER(nDimF); ndF[0] = nbF; ndF[1] = 0; ndF[2] = 0; ndF[3] = nbT;
//        double *totF = REAL(aggregObj(VECTOR_ELT(out_Y_efmit, spQ),nDimF));
//
//            //calcul de diffLQ
//        double diffLQ = 0.0;
//
//        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//         diffLQ = diffLQ + (totF[ind_f + nbF*IND_T] - TAC_byFleet[ind_f + (nbF+1)*IND_T]);
//
//        }
//
//        if ((DIFF*diffLQ)<0) lambda = lambda/3; //ie si DIFF et diffLQ de signe différent
//
//        if ((diffLQ<=0) & (((diffLQ*lambda)*(diffLQ*lambda)<0.25) | (IT>ITquotaExch))) GoOn = false; //on ne s'arrête que si diffLQ<=0 (Quota respecté)
//
//        IT++;
//
//        DIFF = diffLQ;
//
////        Rprintf("\n");
//        Rprintf("diffLQ %f \n",diffLQ);
////
//        Rprintf("PxQ_1 %f ",PxQ);
//
//        if (diffLQ<0) {
//            PxQ = fmax2(pxQuMin,fmin2(PxQ + fmin2(lambda*diffLQ,-0.1),pxQuMax));
//        } else {
//            PxQ = fmax2(pxQuMin,fmin2(PxQ + fmax2(lambda*diffLQ,0.1),pxQuMax));
//        }
//
//
////        PxQ = fmax2(pxQuMin,fmin2(PxQ + lambda*diffLQ,pxQuMax));
//
//        Rprintf("PxQ_2 %f \n",PxQ);
//        Rprintf("lambda %f \n",lambda);
//        //... et on recommence
//
//        UNPROTECT(3);
//
//    }
//
//
//    double *pxQuot = REAL(VECTOR_ELT(out_PQuot_et,spQ));
//    pxQuot[ind_t] = PxQ;
//
//	double *nbdsFM_G = REAL(getListElement(FList, "nbds_f_m"));
//    double *nbdsF_G = REAL(getListElement(FList, "nbds_f"));
//
//    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//
//        if (var==1) nbdsFM_G[ind_f+nbF*0] = fmin2(fmax2(multF[ind_f+1],0.0),effortIni[ind_f]-nbdsFM_G[ind_f+nbF*1]);
//        if (var==1) nbdsF_G[ind_f] = nbdsFM_G[ind_f+nbF*0] + nbdsFM_G[ind_f+nbF*1];
//
//    }
//
//    free_vector(multF,1,nbF);
//    free_vector(z,1,2);
//    free_vector(x,1,1);
//
//    }
//
//	return 0;
//}
//
//}
//



//------------------------------------------
// Module 'Economie' DCF
//------------------------------------------

extern "C" {

void BioEcoPar::EcoDCF(SEXP list, int ind_t, int adj, int lev, int ue_choice, int oths, int othsFM, int perscCalc, int report, double dr)
{

    SEXP Flist;
    PROTECT(Flist = getListElement(list, "Fleet"));

    PROTECT(out_EcoDCF);

//2 protect
    SEXP dimCstF, DimF, dimnamesF, dimCstFM, dimCstFini, dimCstFMini, DimFM, DimFMini, dimnamesFM, dimnamesFMini; //formatage des objets résultats

    SEXP eFACTf, eFACTfm, elmt;

    SEXP    Lref_f, GVLref_f, GVLref_f_m, nbv_f, nbv_f_m, lc_f, nbh_f, ue_f, ue_f_m, fc_f, fc_f_m, vf_f, vf_f_m, ovcDCF_f, ovcDCF_f_m,
            cnb_f, cshr_f, eec_f, mwh_f, rep_f, dep_f, ic_f, K_f, persc_f, fixc_f, mwhg_f;

    SEXP    dc_Lref_f, dc_GVLref_f, dc_GVLref_f_m, dc_nbv_f, dc_nbv_f_m, dc_lc_f, dc_nbh_f, dc_ue_f, dc_ue_f_m,
            dc_fc_f, dc_fc_f_m, dc_vf_f, dc_vf_f_m, dc_ovcDCF_f, dc_ovcDCF_f_m, dc_cnb_f, dc_cshr_f, dc_eec_f, dc_mwh_f, dc_rep_f,
            dc_dep_f, dc_ic_f, dc_K_f, dc_persc_f, dc_fixc_f, dc_mwhg_f;

    int *dCF,*dCFM,*dCFini,*dCFMini,*DF,*DFM, *DFMini;

    int     *dim_Lref_f, *dim_GVLref_f, *dim_GVLref_f_m, *dim_nbv_f, *dim_nbv_f_m, *dim_lc_f, *dim_nbh_f, *dim_ue_f,
            *dim_ue_f_m, *dim_fc_f, *dim_fc_f_m, *dim_vf_f, *dim_vf_f_m, *dim_ovcDCF_f, *dim_ovcDCF_f_m, *dim_cnb_f, *dim_cshr_f,
            *dim_eec_f, *dim_mwh_f, *dim_rep_f, *dim_dep_f, *dim_ic_f, *dim_K_f, *dim_persc_f, *dim_fixc_f, *dim_mwhg_f;

    double  *r_Lref_f, *r_GVLref_f, *r_GVLref_f_m, *r_nbv_f, *r_nbv_f_m, *r_lc_f, *r_nbh_f, *r_ue_f, *r_ue_f_m,
            *r_fc_f, *r_fc_f_m, *r_vf_f, *r_vf_f_m, *r_ovcDCF_f, *r_ovcDCF_f_m, *r_cnb_f, *r_cshr_f, *r_eec_f, *r_mwh_f,
            *r_rep_f, *r_dep_f, *r_ic_f, *r_K_f, *r_persc_f, *r_fixc_f, *r_mwhg_f;

    double  *r_GVLtot_f_m_out, *r_GVLav_f_m_out, *r_GVLtot_f_out, *r_GVLav_f_out, *r_NGVLav_f_m_out, *r_NGVLav_f_out,
            *r_rtbs_f_out, *r_rtbs_f_m_out, *r_cshrT_f_out, *r_sshr_f_out, *r_ncshr_f_out, *r_oclg_f_out, *r_ocl_f_out,
            *r_csg_f_out, *r_cs_f_out, *r_gva_f_out, *r_ccw_f_out, *r_ccwCr_f_out,
            *r_wageg_f_out, *r_wagen_f_out, *r_gcf_f_out, *r_ngcf_f_out, *r_gp_f_out, *r_ps_f_out, *r_sts_f_out,
            *r_ratio_gva_GVL_f_out, *r_ratio_gcf_GVL_f_out, *r_ratio_fc_GVL_f_out,
            *r_ratio_rep_GVL_f_out, *r_ratio_fvol_GVL_f_out, *r_ratio_fvol_gva_f_out, *r_ratio_gcf_gva_f_out,
            *r_ratio_K_cnb_f_out, *r_ratio_GVL_K_f_out, *r_ratio_gcf_K_f_out, *r_ratio_ngcf_K_f_out, *r_ratio_gp_K_f_out,
            *r_ratio_GVL_cnb_ue_f_out,
            *r_rtbsAct_f_out, *r_csAct_f_out, *r_gvaAct_f_out, *r_gcfAct_f_out, *r_psAct_f_out, *r_stsAct_f_out;

//définition des dimensions
    PROTECT(dimCstF = allocVector(INTSXP, 4));
    PROTECT(DimF = allocVector(INTSXP, 2));
    PROTECT(dimnamesF = allocVector(VECSXP,2));
    PROTECT(dimCstFM = allocVector(INTSXP, 4));
    PROTECT(DimFM = allocVector(INTSXP, 3));
    PROTECT(dimnamesFM = allocVector(VECSXP,3));

    PROTECT(dimCstFini = allocVector(INTSXP, 4));
    PROTECT(dimCstFMini = allocVector(INTSXP, 4));
    PROTECT(DimFMini = allocVector(INTSXP, 2));
    PROTECT(dimnamesFMini = allocVector(VECSXP,2));

    SET_VECTOR_ELT(dimnamesF, 0, fleetList); SET_VECTOR_ELT(dimnamesF, 1, times);
    SET_VECTOR_ELT(dimnamesFM, 0, fleetList); SET_VECTOR_ELT(dimnamesFM, 1, metierListEco); SET_VECTOR_ELT(dimnamesFM, 2, times);
    SET_VECTOR_ELT(dimnamesFMini, 0, fleetList); SET_VECTOR_ELT(dimnamesFMini, 1, metierListEco);

    dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
    dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;
    dCFini = INTEGER(dimCstFini) ; dCFini[0] = nbF; dCFini[1] = 0; dCFini[2] = 0; dCFini[3] = 0;
    dCFMini = INTEGER(dimCstFMini) ; dCFMini[0] = nbF; dCFMini[1] = nbMe; dCFMini[2] = 0; dCFMini[3] = 0;

    DF = INTEGER(DimF) ; DF[0] = nbF; DF[1] = nbT;
    DFM = INTEGER(DimFM) ; DFM[0] = nbF; DFM[1] = nbMe; DFM[2] = nbT;
    DFMini = INTEGER(DimFMini) ; DFMini[0] = nbF; DFMini[1] = nbMe;

    //facteurs des indices génériques F/FM
    PROTECT(eFACTf = iDim(dCF));
    PROTECT(eFACTfm = iDim(dCFM));

//12 protect -> 14

    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    PROTECT(Lref_f = getListElement(Flist, "Lref_f"));      PROTECT(dc_Lref_f = iDim(INTEGER(getAttrib(Lref_f, install("DimCst")))));
    PROTECT(GVLref_f = getListElement(Flist, "GVLref_f"));    PROTECT(dc_GVLref_f = iDim(INTEGER(getAttrib(GVLref_f, install("DimCst")))));
    PROTECT(GVLref_f_m = getListElement(Flist, "GVLref_f_m")); PROTECT(dc_GVLref_f_m = iDim(INTEGER(getAttrib(GVLref_f_m, install("DimCst")))));
    PROTECT(nbv_f = getListElement(Flist, "nbv_f"));        PROTECT(dc_nbv_f = iDim(INTEGER(getAttrib(nbv_f, install("DimCst")))));
    PROTECT(nbv_f_m = getListElement(Flist, "nbv_f_m"));    PROTECT(dc_nbv_f_m = iDim(INTEGER(getAttrib(nbv_f_m, install("DimCst")))));
    PROTECT(lc_f = getListElement(Flist, "lc_f"));          PROTECT(dc_lc_f = iDim(INTEGER(getAttrib(lc_f, install("DimCst")))));
    PROTECT(nbh_f = getListElement(Flist, "nbh_f"));        PROTECT(dc_nbh_f = iDim(INTEGER(getAttrib(nbh_f, install("DimCst")))));
    PROTECT(fc_f = getListElement(Flist, "fc_f"));          PROTECT(dc_fc_f = iDim(INTEGER(getAttrib(fc_f, install("DimCst")))));
    PROTECT(fc_f_m = getListElement(Flist, "fc_f_m"));      PROTECT(dc_fc_f_m = iDim(INTEGER(getAttrib(fc_f_m, install("DimCst")))));
    PROTECT(vf_f = getListElement(Flist, "vf_f"));          PROTECT(dc_vf_f = iDim(INTEGER(getAttrib(vf_f, install("DimCst")))));
    PROTECT(vf_f_m = getListElement(Flist, "vf_f_m"));      PROTECT(dc_vf_f_m = iDim(INTEGER(getAttrib(vf_f_m, install("DimCst")))));
    PROTECT(ovcDCF_f = getListElement(Flist, "ovcDCF_f"));  PROTECT(dc_ovcDCF_f = iDim(INTEGER(getAttrib(ovcDCF_f, install("DimCst")))));
    PROTECT(ovcDCF_f_m = getListElement(Flist, "ovcDCF_f_m")); PROTECT(dc_ovcDCF_f_m = iDim(INTEGER(getAttrib(ovcDCF_f_m, install("DimCst")))));
    PROTECT(cnb_f = getListElement(Flist, "cnb_f"));        PROTECT(dc_cnb_f = iDim(INTEGER(getAttrib(cnb_f, install("DimCst")))));
    PROTECT(cshr_f = getListElement(Flist, "cshr_f"));      PROTECT(dc_cshr_f = iDim(INTEGER(getAttrib(cshr_f, install("DimCst")))));
    PROTECT(eec_f = getListElement(Flist, "eec_f"));        PROTECT(dc_eec_f = iDim(INTEGER(getAttrib(eec_f, install("DimCst")))));
    PROTECT(mwh_f = getListElement(Flist, "mwh_f"));        PROTECT(dc_mwh_f = iDim(INTEGER(getAttrib(mwh_f, install("DimCst")))));
    PROTECT(rep_f = getListElement(Flist, "rep_f"));        PROTECT(dc_rep_f = iDim(INTEGER(getAttrib(rep_f, install("DimCst")))));
    PROTECT(dep_f = getListElement(Flist, "dep_f"));        PROTECT(dc_dep_f = iDim(INTEGER(getAttrib(dep_f, install("DimCst")))));
    PROTECT(ic_f = getListElement(Flist, "ic_f"));          PROTECT(dc_ic_f = iDim(INTEGER(getAttrib(ic_f, install("DimCst")))));
    PROTECT(K_f = getListElement(Flist, "K_f"));            PROTECT(dc_K_f = iDim(INTEGER(getAttrib(K_f, install("DimCst")))));
    PROTECT(persc_f = getListElement(Flist, "persc_f"));    PROTECT(dc_persc_f = iDim(INTEGER(getAttrib(persc_f, install("DimCst")))));
    PROTECT(fixc_f = getListElement(Flist, "fixc_f"));      PROTECT(dc_fixc_f = iDim(INTEGER(getAttrib(fixc_f, install("DimCst")))));
    PROTECT(mwhg_f = getListElement(Flist, "mwhg_f"));      PROTECT(dc_mwhg_f = iDim(INTEGER(getAttrib(mwhg_f, install("DimCst")))));
//42 protect  ->56

if (ue_choice == 1) {

    PROTECT(ue_f = getListElement(Flist, "nbds_f"));
    PROTECT(ue_f_m = getListElement(Flist, "nbds_f_m"));

} else {

    if (ue_choice == 2) {

        PROTECT(ue_f = getListElement(Flist, "nbh_f"));
        PROTECT(ue_f_m = getListElement(Flist, "nbh_f_m"));

    } else {

        PROTECT(ue_f = getListElement(Flist, "nbtrip_f"));
        PROTECT(ue_f_m = getListElement(Flist, "nbtrip_f_m"));

    }
}

    PROTECT(dc_ue_f = iDim(INTEGER(getAttrib(ue_f, install("DimCst")))));
    PROTECT(dc_ue_f_m = iDim(INTEGER(getAttrib(ue_f_m, install("DimCst")))));
//4 protect -> 60

    dim_Lref_f = INTEGER(dc_Lref_f);                        r_Lref_f = REAL(Lref_f);
    dim_GVLref_f = INTEGER(dc_GVLref_f);                    r_GVLref_f = REAL(GVLref_f);
    dim_GVLref_f_m = INTEGER(dc_GVLref_f_m);                r_GVLref_f_m = REAL(GVLref_f_m);
    dim_nbv_f = INTEGER(dc_nbv_f);                          r_nbv_f = REAL(nbv_f);
    dim_nbv_f_m = INTEGER(dc_nbv_f_m);                      r_nbv_f_m = REAL(nbv_f_m);
    dim_lc_f = INTEGER(dc_lc_f);                            r_lc_f = REAL(lc_f);
    dim_nbh_f = INTEGER(dc_nbh_f);                          r_nbh_f = REAL(nbh_f);
    dim_ue_f = INTEGER(dc_ue_f);                            r_ue_f = REAL(ue_f);
    dim_ue_f_m = INTEGER(dc_ue_f_m);                        r_ue_f_m = REAL(ue_f_m);
    dim_fc_f = INTEGER(dc_fc_f);                            r_fc_f = REAL(fc_f);
    dim_fc_f_m = INTEGER(dc_fc_f_m);                        r_fc_f_m = REAL(fc_f_m);
    dim_vf_f = INTEGER(dc_vf_f);                            r_vf_f = REAL(vf_f);
    dim_vf_f_m = INTEGER(dc_vf_f_m);                        r_vf_f_m = REAL(vf_f_m);
    dim_ovcDCF_f = INTEGER(dc_ovcDCF_f);                    r_ovcDCF_f = REAL(ovcDCF_f);
    dim_ovcDCF_f_m = INTEGER(dc_ovcDCF_f_m);                r_ovcDCF_f_m = REAL(ovcDCF_f_m);
    dim_cnb_f = INTEGER(dc_cnb_f);                          r_cnb_f = REAL(cnb_f);
    dim_cshr_f = INTEGER(dc_cshr_f);                        r_cshr_f = REAL(cshr_f);
    dim_eec_f = INTEGER(dc_eec_f);                          r_eec_f = REAL(eec_f);
    dim_mwh_f = INTEGER(dc_mwh_f);                          r_mwh_f = REAL(mwh_f);
    dim_rep_f = INTEGER(dc_rep_f);                          r_rep_f = REAL(rep_f);
    dim_dep_f = INTEGER(dc_dep_f);                          r_dep_f = REAL(dep_f);
    dim_ic_f = INTEGER(dc_ic_f);                            r_ic_f = REAL(ic_f);
    dim_K_f = INTEGER(dc_K_f);                              r_K_f = REAL(K_f);
    dim_persc_f = INTEGER(dc_persc_f);                      r_persc_f = REAL(persc_f);
    dim_fixc_f = INTEGER(dc_fixc_f);                        r_fixc_f = REAL(fixc_f);
    dim_mwhg_f = INTEGER(dc_mwhg_f);                        r_mwhg_f = REAL(mwhg_f);

    int nbI, nbC;

if (ind_t==0) {

    SEXP  GVL_f_m_e_out, GVLtot_f_m_out, GVLav_f_m_out, GVLtot_f_out, GVLav_f_out, NGVLav_f_m_out, NGVLav_f_out,
          rtbs_f_out, rtbs_f_m_out, cshrT_f_out, sshr_f_out, ncshr_f_out, oclg_f_out, ocl_f_out, csg_f_out, cs_f_out,
          gva_f_out, ccw_f_out, ccwCr_f_out, wageg_f_out, wagen_f_out, gcf_f_out,
          ngcf_f_out, gp_f_out, ps_f_out, sts_f_out, ratio_gva_GVL_f_out, ratio_gcf_GVL_f_out, ratio_fc_GVL_f_out,
          ratio_rep_GVL_f_out, ratio_fvol_GVL_f_out, ratio_fvol_gva_f_out,
          ratio_gcf_gva_f_out, ratio_K_cnb_f_out, ratio_GVL_K_f_out, ratio_gcf_K_f_out, ratio_ngcf_K_f_out,
          ratio_gp_K_f_out, ratio_GVL_cnb_ue_f_out,
          rtbsAct_f_out, csAct_f_out, gvaAct_f_out, gcfAct_f_out, psAct_f_out, stsAct_f_out;

    SEXP  GVLtot_f_m_e, GVLreftot_f_m_e, GVLreftot_f_m, GVLoths_f_m, GVLothsref_f_m, GVLothsue_f_m, GVLothsrefue_f_m,
          GVLtot_f_e, GVLreftot_f, GVLoths_f, GVLothsue_f, GVLothmet_f, GVLothmetue_f,
          fvolue_f, fvolue_f_m, ovcDCFue_f, ovcDCFue_f_m, GVLav_f, rtbs_f, ccwr_f, opersc_f, eco_names;

    double  *r_GVLreftot_f_m, *r_GVLoths_f_m, *r_GVLothsref_f_m, *r_GVLothsue_f_m,
            *r_GVLothsrefue_f_m, *r_GVLreftot_f, *r_GVLoths_f, *r_GVLothsue_f, *r_GVLothmet_f,
            *r_GVLothmetue_f, *r_fvolue_f, *r_fvolue_f_m, *r_ovcDCFue_f, *r_ovcDCFue_f_m, *r_GVLav_f, *r_rtbs_f, *r_ccwr_f, *r_opersc_f;



//-------------------------
// Stade préliminaire (temps initial)
//-------------------------

    PROTECT(GVLreftot_f_m = NEW_NUMERIC(nbF*nbMe));          r_GVLreftot_f_m = REAL(GVLreftot_f_m);
    PROTECT(GVLoths_f_m = NEW_NUMERIC(nbF*nbMe));            r_GVLoths_f_m = REAL(GVLoths_f_m);
    PROTECT(GVLothsref_f_m = NEW_NUMERIC(nbF*nbMe));         r_GVLothsref_f_m = REAL(GVLothsref_f_m);
    PROTECT(GVLothsue_f_m = NEW_NUMERIC(nbF*nbMe));          r_GVLothsue_f_m = REAL(GVLothsue_f_m);
    PROTECT(GVLothsrefue_f_m = NEW_NUMERIC(nbF*nbMe));       r_GVLothsrefue_f_m = REAL(GVLothsrefue_f_m);
    PROTECT(GVLreftot_f = NEW_NUMERIC(nbF));                 r_GVLreftot_f = REAL(GVLreftot_f);
    PROTECT(GVLoths_f = NEW_NUMERIC(nbF));                   r_GVLoths_f = REAL(GVLoths_f);
    PROTECT(GVLothsue_f = NEW_NUMERIC(nbF));                 r_GVLothsue_f = REAL(GVLothsue_f);
    PROTECT(GVLothmet_f = NEW_NUMERIC(nbF));                 r_GVLothmet_f = REAL(GVLothmet_f);
    PROTECT(GVLothmetue_f = NEW_NUMERIC(nbF));               r_GVLothmetue_f = REAL(GVLothmetue_f);
    PROTECT(fvolue_f = NEW_NUMERIC(nbF));                   r_fvolue_f = REAL(fvolue_f);
    PROTECT(fvolue_f_m = NEW_NUMERIC(nbF*nbMe));             r_fvolue_f_m = REAL(fvolue_f_m);
    PROTECT(ovcDCFue_f = NEW_NUMERIC(nbF));                  r_ovcDCFue_f = REAL(ovcDCFue_f);
    PROTECT(ovcDCFue_f_m = NEW_NUMERIC(nbF*nbMe));           r_ovcDCFue_f_m = REAL(ovcDCFue_f_m);
    PROTECT(GVLav_f = NEW_NUMERIC(nbF));                     r_GVLav_f = REAL(GVLav_f);
    PROTECT(rtbs_f = NEW_NUMERIC(nbF));                     r_rtbs_f = REAL(rtbs_f);
    PROTECT(ccwr_f = NEW_NUMERIC(nbF));                     r_ccwr_f = REAL(ccwr_f);
    PROTECT(opersc_f = NEW_NUMERIC(nbF));                   r_opersc_f = REAL(opersc_f);
//16 protect  -> 16 (t0)


    for (int e = 0 ; e < nbE ; e++) {

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        nbC = length(getListElement(elmt, "modC"));

        PROTECT(GVLtot_f_m_e = NEW_NUMERIC(nbF*nbMe*nbT));
        PROTECT(GVLreftot_f_m_e = NEW_NUMERIC(nbF*nbMe));
        PROTECT(GVLtot_f_e = NEW_NUMERIC(nbF));

        double *r_GVLtot_f_m_e = REAL(GVLtot_f_m_e);
        double *r_GVLreftot_f_m_e = REAL(GVLreftot_f_m_e);
        double *r_GVLtot_f_e = REAL(GVLtot_f_e);
        double *r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct2, e));
        double *r_GVLref_f_m_e = REAL(getListElement(elmt, "GVLref_f_m_e"));
        double *r_GVLref_f_e = REAL(getListElement(elmt, "GVLref_f_e"));
        double *r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
        int *dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct2, e), install("DimCst")))));
        int *dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));


        //------------------------------
        //équations de la table "p"
        //------------------------------


        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){   //on rappelle ici que ind_t est en fait égal à 0

    double countGVLtotf = 0.0; //pour sommer GVLtot_f_m_e sur les métiers

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

         //-- 3. GVLtot_f_m_e

             double count = 0.0;

             for (int ind_c = 0 ; ind_c < nbC ; ind_c++){

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]))

                count = count +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] ;

             }

            countGVLtotf = countGVLtotf + count;
            r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = count;


         //-- 4. GVLreftot_f_m_e

                r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLref_f_m_e[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                  r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]; //rappel : ind_t = 0 ici




        //-- 5. GVLreftot_f_m

            if (e==0) {

                r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLref_f_m[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                  r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]; //rappel : ind_t = 0 ici

            }

         //-- 6. GVLoths_f_m

            if (e==0) {

                if (!ISNA(r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }

            } else {

                if (!ISNA(r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }



        //-- 7. GVLothsref_f_m

            if (e==0) {

                if (!ISNA(r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                } else {

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }

            } else {

                if (!ISNA(r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                    r_GVLreftot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }

            }

    //-- 10. GVLtot_f_e

    if (lev==2 & adj==1) {

            r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = countGVLtotf;

    } else {

            r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLref_f_e[ind_f*dim_GVLref_f[0] + 0*dim_GVLref_f[1] + 0*dim_GVLref_f[2] + ind_t*dim_GVLref_f[3]] *
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];
    }


    //-- 11. GVLreftot_f

    if (e==0) {

            r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLref_f[ind_f*dim_GVLref_f[0] + 0*dim_GVLref_f[1] + 0*dim_GVLref_f[2] + ind_t*dim_GVLref_f[3]] *
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

    }

    //-- 14. GVLoths_f

            if (e==0) {

                if (!ISNA(r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]])) {

                r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                } else {

                r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                }

            } else {

                if (!ISNA(r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]))

                r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_GVLtot_f_e[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            }



        }

        //on formatte le(s) résultat(s) et on les intègre à 'eVar'

        setAttrib(GVLtot_f_m_e, R_DimSymbol, DimFM);
        setAttrib(GVLtot_f_m_e, R_DimNamesSymbol, dimnamesFM);
        setAttrib(GVLtot_f_m_e, install("DimCst"), dimCstFM);
        SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 41, GVLtot_f_m_e);

        UNPROTECT(4);

}

    // à ce stade, plus de considération d'espèce pour les variables à initialiser

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

 double countEff = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){


        //-- 14. GVLothmet_f

            if (lev==2) {

                r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

            } else {

                if (ind_m==0) {

                    if (!ISNA(r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                    } else {

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                    }

                } else {

                    if (!ISNA(r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]))

                    r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        r_GVLreftot_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                }
            }

           //-- terme : SUM_m (ue_f_m * nbv_f_m)

            countEff = countEff +
                r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

            }

        //-- 13. GVLothsue_f

        if (report==0) {

                r_GVLothsue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) );


        } else {

                r_GVLothsue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLoths_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    ((r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) - countEff) );

        }


        //-- 15. GVLothmetue_f

        if (lev==2) {

                r_GVLothmetue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

        } else {

                r_GVLothmetue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_GVLothmet_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    ((r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]]) - countEff) );

        }


   for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 8. GVLothsue_f_m

                r_GVLothsue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_GVLoths_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    (r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) );


        //-- 9. GVLothsrefue_f_m

                r_GVLothsrefue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_GVLothsref_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    (r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]) );




      //-- 22. fvolue_f_m

               r_fvolue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]] /
                    (r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]] *
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]]) );

      //-- 23. ovcDCFue_f_m

                r_ovcDCFue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_ovcDCF_f_m[ind_f*dim_ovcDCF_f_m[0] + ind_m*dim_ovcDCF_f_m[1] + 0*dim_ovcDCF_f_m[2] + ind_t*dim_ovcDCF_f_m[3]] /
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] );

   }

        //-- 16. fvolue_f

                r_fvolue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_fc_f[ind_f*dim_fc_f[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]] /
                    (r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] *
                     r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) );


        //-- 17. GVLav_f

            r_GVLav_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLreftot_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


        //-- 18. rtbs_f

            r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_GVLav_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                r_fc_f[ind_f*dim_fc_f[0] + 0*dim_fc_f[1] + 0*dim_fc_f[2] + ind_t*dim_fc_f[3]] -
                r_ovcDCF_f[ind_f*dim_ovcDCF_f[0] + 0*dim_ovcDCF_f[1] + 0*dim_ovcDCF_f[2] + ind_t*dim_ovcDCF_f[3]];


        //-- 19. ccwr_f

            r_ccwr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] /
                r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 20. opersc_f

            r_opersc_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] -
                (0.01 * r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                r_rtbs_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]);

        //-- 21. ovcDCFue_f

                r_ovcDCFue_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    finite( r_ovcDCF_f[ind_f*dim_ovcDCF_f[0] + 0*dim_ovcDCF_f[1] + 0*dim_ovcDCF_f[2] + ind_t*dim_ovcDCF_f[3]] /
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] );


}



//on formatte le(s) résultat(s) et on intègre à fVar

        setAttrib(GVLoths_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLoths_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLoths_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 1, GVLoths_f_m);

        setAttrib(GVLothsref_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsref_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsref_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 2, GVLothsref_f_m);

        setAttrib(GVLothsue_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 23, GVLothsue_f_m);

        setAttrib(GVLothsrefue_f_m, R_DimSymbol, DimFMini);
        setAttrib(GVLothsrefue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(GVLothsrefue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 24, GVLothsrefue_f_m);

        setAttrib(GVLothmet_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothmet_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 3, GVLothmet_f);

        setAttrib(GVLothmetue_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothmetue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 25, GVLothmetue_f);

        setAttrib(GVLothsue_f, R_NamesSymbol, fleetList);
        setAttrib(GVLothsue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 26, GVLothsue_f);

        setAttrib(fvolue_f_m, R_DimSymbol, DimFMini);
        setAttrib(fvolue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(fvolue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 4, fvolue_f_m);

        setAttrib(ovcDCFue_f_m, R_DimSymbol, DimFMini);
        setAttrib(ovcDCFue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(ovcDCFue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 10, ovcDCFue_f_m);

        setAttrib(fvolue_f, R_NamesSymbol, fleetList);
        setAttrib(fvolue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 17, fvolue_f);

        setAttrib(ovcDCFue_f, R_NamesSymbol, fleetList);
        setAttrib(ovcDCFue_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 18, ovcDCFue_f);

        setAttrib(ccwr_f, R_NamesSymbol, fleetList);
        setAttrib(ccwr_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 27, ccwr_f);

        setAttrib(opersc_f, R_NamesSymbol, fleetList);
        setAttrib(opersc_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 28, opersc_f);

        setAttrib(GVLoths_f, R_NamesSymbol, fleetList);
        setAttrib(GVLoths_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 29, GVLoths_f);

        SET_VECTOR_ELT(fVar, 31, rtbs_f);

//enfin, on initialise l'output

    PROTECT(GVL_f_m_e_out = allocVector(VECSXP, nbE));
    setAttrib(GVL_f_m_e_out, R_NamesSymbol, sppList);
    SET_VECTOR_ELT(out_EcoDCF, 0, GVL_f_m_e_out);

    PROTECT(GVLtot_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLtot_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLtot_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLtot_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 1, GVLtot_f_m_out);

    PROTECT(GVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 2, GVLav_f_m_out);

    PROTECT(GVLtot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLtot_f_out, R_DimSymbol, DimF);
    setAttrib(GVLtot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLtot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 3, GVLtot_f_out);

    PROTECT(GVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLav_f_out, R_DimSymbol, DimF);
    setAttrib(GVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 4, GVLav_f_out);

    PROTECT(NGVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(NGVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(NGVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(NGVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 5, NGVLav_f_m_out);

    PROTECT(NGVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(NGVLav_f_out, R_DimSymbol, DimF);
    setAttrib(NGVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(NGVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 6, NGVLav_f_out);

    PROTECT(rtbs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbs_f_out, R_DimSymbol, DimF);
    setAttrib(rtbs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 7, rtbs_f_out);

    PROTECT(cshrT_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cshrT_f_out, R_DimSymbol, DimF);
    setAttrib(cshrT_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cshrT_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 8, cshrT_f_out);

    PROTECT(sshr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(sshr_f_out, R_DimSymbol, DimF);
    setAttrib(sshr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(sshr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 9, sshr_f_out);

    PROTECT(ncshr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ncshr_f_out, R_DimSymbol, DimF);
    setAttrib(ncshr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ncshr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 10, ncshr_f_out);

    PROTECT(oclg_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(oclg_f_out, R_DimSymbol, DimF);
    setAttrib(oclg_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(oclg_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 11, oclg_f_out);

    PROTECT(ocl_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ocl_f_out, R_DimSymbol, DimF);
    setAttrib(ocl_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ocl_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 12, ocl_f_out);

    PROTECT(csg_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csg_f_out, R_DimSymbol, DimF);
    setAttrib(csg_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csg_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 13, csg_f_out);

    PROTECT(cs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cs_f_out, R_DimSymbol, DimF);
    setAttrib(cs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 14, cs_f_out);

    PROTECT(gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gva_f_out, R_DimSymbol, DimF);
    setAttrib(gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 15, gva_f_out);

    PROTECT(ccw_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccw_f_out, R_DimSymbol, DimF);
    setAttrib(ccw_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccw_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 16, ccw_f_out);

    PROTECT(ccwCr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccwCr_f_out, R_DimSymbol, DimF);
    setAttrib(ccwCr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccwCr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 17, ccwCr_f_out);

    PROTECT(wageg_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wageg_f_out, R_DimSymbol, DimF);
    setAttrib(wageg_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wageg_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 18, wageg_f_out);

    PROTECT(wagen_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wagen_f_out, R_DimSymbol, DimF);
    setAttrib(wagen_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wagen_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 19, wagen_f_out);

    PROTECT(gcf_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gcf_f_out, R_DimSymbol, DimF);
    setAttrib(gcf_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gcf_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 20, gcf_f_out);

    PROTECT(ngcf_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ngcf_f_out, R_DimSymbol, DimF);
    setAttrib(ngcf_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ngcf_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 21, ngcf_f_out);

    PROTECT(gp_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gp_f_out, R_DimSymbol, DimF);
    setAttrib(gp_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gp_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 22, gp_f_out);

    PROTECT(ps_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ps_f_out, R_DimSymbol, DimF);
    setAttrib(ps_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ps_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 23, ps_f_out);

    PROTECT(sts_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(sts_f_out, R_DimSymbol, DimF);
    setAttrib(sts_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(sts_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 24, sts_f_out);

    PROTECT(ratio_gva_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gva_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gva_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gva_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 25, ratio_gva_GVL_f_out);

    PROTECT(ratio_gcf_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 26, ratio_gcf_GVL_f_out);

    PROTECT(ratio_fc_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fc_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fc_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fc_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 27, ratio_fc_GVL_f_out);

    PROTECT(ratio_rep_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_rep_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_rep_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_rep_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 28, ratio_rep_GVL_f_out);

    PROTECT(ratio_fvol_GVL_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_GVL_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_GVL_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_GVL_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 29, ratio_fvol_GVL_f_out);

    PROTECT(ratio_fvol_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 30, ratio_fvol_gva_f_out);

    PROTECT(ratio_gcf_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 31, ratio_gcf_gva_f_out);

    PROTECT(ratio_K_cnb_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_K_cnb_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_K_cnb_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_K_cnb_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 32, ratio_K_cnb_f_out);

    PROTECT(ratio_GVL_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 33, ratio_GVL_K_f_out);

    PROTECT(ratio_gcf_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gcf_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gcf_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gcf_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 34, ratio_gcf_K_f_out);

    PROTECT(ratio_ngcf_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_ngcf_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_ngcf_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_ngcf_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 35, ratio_ngcf_K_f_out);

    PROTECT(ratio_gp_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gp_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gp_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gp_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 36, ratio_gp_K_f_out);

    PROTECT(ratio_GVL_cnb_ue_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_cnb_ue_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 37, ratio_GVL_cnb_ue_f_out);


    PROTECT(rtbsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbsAct_f_out, R_DimSymbol, DimF);
    setAttrib(rtbsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 38, rtbsAct_f_out);

    PROTECT(csAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csAct_f_out, R_DimSymbol, DimF);
    setAttrib(csAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 39, csAct_f_out);

    PROTECT(gvaAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gvaAct_f_out, R_DimSymbol, DimF);
    setAttrib(gvaAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gvaAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 40, gvaAct_f_out);

    PROTECT(gcfAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gcfAct_f_out, R_DimSymbol, DimF);
    setAttrib(gcfAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gcfAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 41, gcfAct_f_out);

    PROTECT(psAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(psAct_f_out, R_DimSymbol, DimF);
    setAttrib(psAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(psAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 42, psAct_f_out);

    PROTECT(stsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(stsAct_f_out, R_DimSymbol, DimF);
    setAttrib(stsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(stsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 43, stsAct_f_out);

    PROTECT(rtbs_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(rtbs_f_m_out, R_DimSymbol, DimFM);
    setAttrib(rtbs_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(rtbs_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 44, rtbs_f_m_out);


    //on nomme les éléments de out_EcoDCF
    const char *namesEco[45] = {"GVL_f_m_e","GVLtot_f_m","GVLav_f_m","GVLtot_f","GVLav_f","NGVLav_f_m","NGVLav_f",
                          "rtbs_f","cshrT_f","sshr_f","ncshr_f","oclg_f","ocl_f","csg_f","cs_f",
                          "gva_f","ccw_f","ccwCr_f","wageg_f","wagen_f","gcf_f","ngcf_f","gp_f",
                          "ps_f","sts_f","ratio_gva_GVL_f","ratio_gcf_GVL_f","ratio_fc_GVL_f",
                          "ratio_rep_GVL_f","ratio_fvol_GVL_f","ratio_fvol_gva_f","ratio_gcf_gva_f",
                          "ratio_K_cnb_f","ratio_GVL_K_f","ratio_gcf_K_f","ratio_ngcf_K_f",
                          "ratio_gp_K_f","ratio_GVL_cnb_ue_f","rtbsAct_f","csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","rtbs_f_m"};

    PROTECT(eco_names = allocVector(STRSXP, 45));

    for(int ct = 0; ct < 45; ct++) SET_STRING_ELT(eco_names, ct, mkChar(namesEco[ct]));

    setAttrib(out_EcoDCF, R_NamesSymbol, eco_names);

//39 protect    --> 55
}

//on importe les outputs afin de les mettre à jour à l'instant ind_t

    r_GVLtot_f_m_out = REAL(VECTOR_ELT(out_EcoDCF,1));
    r_GVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF,2));
    r_GVLtot_f_out = REAL(VECTOR_ELT(out_EcoDCF,3));
    r_GVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF,4));
    r_NGVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF,5));
    r_NGVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF,6));
    r_rtbs_f_out = REAL(VECTOR_ELT(out_EcoDCF,7));
    r_cshrT_f_out = REAL(VECTOR_ELT(out_EcoDCF,8));
    r_sshr_f_out = REAL(VECTOR_ELT(out_EcoDCF,9));
    r_ncshr_f_out = REAL(VECTOR_ELT(out_EcoDCF,10));
    r_oclg_f_out = REAL(VECTOR_ELT(out_EcoDCF,11));
    r_ocl_f_out = REAL(VECTOR_ELT(out_EcoDCF,12));
    r_csg_f_out = REAL(VECTOR_ELT(out_EcoDCF,13));
    r_cs_f_out = REAL(VECTOR_ELT(out_EcoDCF,14));
    r_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF,15));
    r_ccw_f_out = REAL(VECTOR_ELT(out_EcoDCF,16));
    r_ccwCr_f_out = REAL(VECTOR_ELT(out_EcoDCF,17));
    r_wageg_f_out = REAL(VECTOR_ELT(out_EcoDCF,18));
    r_wagen_f_out = REAL(VECTOR_ELT(out_EcoDCF,19));
    r_gcf_f_out = REAL(VECTOR_ELT(out_EcoDCF,20));
    r_ngcf_f_out = REAL(VECTOR_ELT(out_EcoDCF,21));
    r_gp_f_out = REAL(VECTOR_ELT(out_EcoDCF,22));
    r_ps_f_out = REAL(VECTOR_ELT(out_EcoDCF,23));
    r_sts_f_out = REAL(VECTOR_ELT(out_EcoDCF,24));
    r_ratio_gva_GVL_f_out = REAL(VECTOR_ELT(out_EcoDCF,25));
    r_ratio_gcf_GVL_f_out = REAL(VECTOR_ELT(out_EcoDCF,26));
    r_ratio_fc_GVL_f_out = REAL(VECTOR_ELT(out_EcoDCF,27));
    r_ratio_rep_GVL_f_out = REAL(VECTOR_ELT(out_EcoDCF,28));
    r_ratio_fvol_GVL_f_out = REAL(VECTOR_ELT(out_EcoDCF,29));
    r_ratio_fvol_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF,30));
    r_ratio_gcf_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF,31));
    r_ratio_K_cnb_f_out = REAL(VECTOR_ELT(out_EcoDCF,32));
    r_ratio_GVL_K_f_out = REAL(VECTOR_ELT(out_EcoDCF,33));
    r_ratio_gcf_K_f_out = REAL(VECTOR_ELT(out_EcoDCF,34));
    r_ratio_ngcf_K_f_out = REAL(VECTOR_ELT(out_EcoDCF,35));
    r_ratio_gp_K_f_out = REAL(VECTOR_ELT(out_EcoDCF,36));
    r_ratio_GVL_cnb_ue_f_out = REAL(VECTOR_ELT(out_EcoDCF,37));
    r_rtbsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,38));
    r_csAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,39));
    r_gvaAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,40));
    r_gcfAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,41));
    r_psAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,42));
    r_stsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF,43));
    r_rtbs_f_m_out = REAL(VECTOR_ELT(out_EcoDCF,44));

    double *r_GVLoths_f_m2 = REAL(VECTOR_ELT(fVar,1));
    double *r_GVLothsref_f_m2 = REAL(VECTOR_ELT(fVar,2));
    double *r_GVLothsue_f_m2 = REAL(VECTOR_ELT(fVar,23));
    double *r_GVLothsrefue_f_m2 = REAL(VECTOR_ELT(fVar,24));
    double *r_GVLothmet_f2 = REAL(VECTOR_ELT(fVar,3));
    double *r_GVLothmetue_f2 = REAL(VECTOR_ELT(fVar,25));
    double *r_GVLothsue_f2 = REAL(VECTOR_ELT(fVar,26));
    double *r_fvolue_f2 = REAL(VECTOR_ELT(fVar,17));
    double *r_fvolue_f_m2 = REAL(VECTOR_ELT(fVar,4));
    double *r_ovcDCFue_f2 = REAL(VECTOR_ELT(fVar,18));
    double *r_ovcDCFue_f_m2 = REAL(VECTOR_ELT(fVar,10));
    double *r_ccwr_f2 = REAL(VECTOR_ELT(fVar,27));
    //double *r_opersc_f2 = REAL(VECTOR_ELT(fVar,28));
    double *r_GVLoths_f2 = REAL(VECTOR_ELT(fVar,29));
    double *r_rtbs_f2 = REAL(VECTOR_ELT(fVar,31));

//for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
//for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){
//Rprintf("%f -- ",r_ovcDCFue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
//Rprintf("%f -- ",r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
////Rprintf("%f -- ",r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]);
////Rprintf("%f -- ",r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]]);
//Rprintf("AA\n ");
//}
//PrintValue(VECTOR_ELT(fVar,10));



   SEXP countGVLf;
   PROTECT(countGVLf = NEW_NUMERIC(nbF)); // --> 61
   double *r_countGVLf = REAL(countGVLf);


   for (int e = 0 ; e < nbE ; e++) {

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        nbC = length(getListElement(elmt, "modC"));

        double *r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
        double *r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct2, e));
        double *r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
        int *dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct2, e), install("DimCst")))));
        int *dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));


        //---------------------
        //équations de la table "t"
        //---------------------

  for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        r_countGVLf[ind_f] = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 1. GVL_f_m_e

             double count = 0.0;

             for (int ind_c = 0 ; ind_c < nbC ; ind_c++){

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]))

                count = count +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

             }

            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = count / pow(1+0.0,ind_t) ;
            r_countGVLf[ind_f] = r_countGVLf[ind_f] + (count / pow(1+0.0,ind_t));


        //-- 2. GVLtot_f_m


            if (e==0) {

                if (othsFM==1) {

                    if (adj==2) {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLothsrefue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                            r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                            r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]
                             / pow(1+0.0,ind_t)) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] ;

                    } else {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLothsue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                            r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                            r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]]
                            / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                    }

                } else {

                    if (adj==2) {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            ( r_GVLothsref_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] ;

                    } else {

                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                            (r_GVLoths_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] / pow(1+0.0,ind_t) ) +
                            r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                    }

                }

            } else {

                r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
                  r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

            }
        }
   }

  SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,0), e, VECTOR_ELT(VECTOR_ELT(eVar, e),41));

  UNPROTECT(1);

}


    // à ce stade, plus de considération d'espèce pour les indicateurs

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            double countEff = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++)

                countEff = countEff +
                    r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] *
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];


            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            //-- 3. GVLav_f_m

                r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + 0*dim_nbv_f_m[3]];

             //-- 4. GVLtot_f

                if (ind_m==0) {

                    if (lev==1) {

                        if (report==0) {

                            if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmet_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) ) +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                            } else {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmet_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t) );
                            }

                        } else {

                            if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmetue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                    ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff) +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                            } else {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    (r_GVLothmetue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                    ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff);

                            }
                        }

                    } else { //lev = 2

                        if (oths==0) {

                            r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                r_countGVLf[ind_f] +
                                (r_GVLoths_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));

                        } else {

                           if (report==0) {

                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                        r_countGVLf[ind_f] +
                                        (r_GVLothsue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                        ((r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]) - countEff);

                            } else {

                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                        r_countGVLf[ind_f] +
                                        ((r_GVLothsue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t)) *
                                        r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]] *
                                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]);

                            }
                        }
                    }

                } else { //ind_m>0

                    if (lev==1) {

                        if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
                        }
                    }
                }



            //-- 6. NGVLav_f_m

                r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                    (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]]);


            //-- 8bis. rtbs_f_m

                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                        ((r_ovcDCFue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] +
                        r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] *
                        r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]) *
                        r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] / pow(1+0.0,ind_t));
//Rprintf("%f -- ",r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
//Rprintf("%f -- ",r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
//Rprintf("%f -- ",r_ovcDCFue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
//Rprintf("%f -- ",r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);
//Rprintf("%f -- ",r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]);
//Rprintf("%f -- ",r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]]);
//Rprintf("BB\n ");
            } //on sort de la boucle sur les niveaux métiers
//Rprintf("\n ");
//Rprintf("\n ");
//if (ind_f<5 & ind_t<5) PrintValue(VECTOR_ELT(fVar,4));
//if (ind_f<5 & ind_t<5) PrintValue(VECTOR_ELT(fVar,10));

            //-- 5. GVLav_f

                r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + 0*dim_nbv_f[3]];


            //-- 7. NGVLav_f

                r_NGVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    (1 - 0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]]);


            //-- 8. rtbs_f

                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                        ((r_ovcDCFue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] +
                        r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                        r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]]) *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t));

                //version actualisée
            r_rtbsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);




            //-- 9. cshrT_f

            if (perscCalc==0) {

               r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];

            } else {

                 r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];
            }

            //-- 10. sshr_f


                 r_sshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 11. ncshr_f

                r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_eec_f[ind_f*dim_eec_f[0] + 0*dim_eec_f[1] + 0*dim_eec_f[2] + ind_t*dim_eec_f[3]] / pow(1+0.0,ind_t) );

             //-- 12. oclg_f

                r_oclg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_mwhg_f[ind_f*dim_mwhg_f[0] + 0*dim_mwhg_f[1] + 0*dim_mwhg_f[2] + 0*dim_mwhg_f[3]] *
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                    r_nbh_f[ind_f*dim_nbh_f[0] + 0*dim_nbh_f[1] + 0*dim_nbh_f[2] + ind_t*dim_nbh_f[3]] / pow(1+0.0,ind_t);


             //-- 13. ocl_f

                r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_mwh_f[ind_f*dim_mwh_f[0] + 0*dim_mwh_f[1] + 0*dim_mwh_f[2] + 0*dim_mwh_f[3]] *
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                    r_nbh_f[ind_f*dim_nbh_f[0] + 0*dim_nbh_f[1] + 0*dim_nbh_f[2] + ind_t*dim_nbh_f[3]] / pow(1+0.0,ind_t);

             //-- 14. csg_f

                r_csg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_oclg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 15. cs_f

                r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                //version actualisée
            r_csAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



             //-- 16. gva_f

                r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] +
                    r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t*dim_fixc_f[3]]) / pow(1+0.0,ind_t) ;


            //version actualisée
            r_gvaAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 17. ccw_f

            if (perscCalc==1) {

                r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ; //+
                    //(r_opersc_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));

            } else {

                if (perscCalc==2) {

                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ccwr_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                } else {

                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] ;//+
                    //(r_opersc_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] / pow(1+0.0,ind_t));
                }
            }


             //-- 18. ccwCr_f

             r_ccwCr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];


            //-- 19. wageg_f

             r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];

            //-- 20. wagen_f

             r_wagen_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
             r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
             r_cnb_f[ind_f* dim_cnb_f[0] + 0* dim_cnb_f[1] + 0* dim_cnb_f[2] + ind_t* dim_cnb_f[3]];


             //-- 21. gcf_f

                r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


            //version actualisée
            r_gcfAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



             //-- 22. ngcf_f

                r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]] / pow(1+0.0,ind_t) );


             //-- 23. gp_f

                r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]] / pow(1+0.0,ind_t) );

             //-- 24. ps_f

                r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                     r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]) *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


            //version actualisée
            r_psAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 25. sts_f

                r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    0.01*r_lc_f[ind_f*dim_lc_f[0] + 0*dim_lc_f[1] + 0*dim_lc_f[2] + ind_t*dim_lc_f[3]] *
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


            //version actualisée
            r_stsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



             //-- 26. ratio_gva_GVL_f

                r_ratio_gva_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 27. ratio_gcf_GVL_f

                r_ratio_gcf_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 28. ratio_fc_GVL_f

                r_ratio_fc_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                        r_vf_f[ind_f*dim_vf_f[0] + 0*dim_vf_f[1] + 0*dim_vf_f[2] + ind_t*dim_vf_f[3]] *
                        r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] / pow(1+0.0,ind_t)) /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 29. ratio_rep_GVL_f

                r_ratio_rep_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] / pow(1+0.0,ind_t) )/
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 30. ratio_fvol_GVL_f

                r_ratio_fvol_GVL_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] /
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 31. ratio_fvol_gva_f

                r_ratio_fvol_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_fvolue_f2[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]] /
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 32. ratio_gcf_gva_f

                r_ratio_gcf_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 33. ratio_K_cnb_f

                r_ratio_K_cnb_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) ) /
                    r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]];

            //-- 34. ratio_GVL_K_f

                r_ratio_GVL_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

           //-- 35. ratio_gcf_K_f

                r_ratio_gcf_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

           //-- 36. ratio_ngcf_K_f

                r_ratio_ngcf_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ngcf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

          //-- 37. ratio_gp_K_f

                r_ratio_gp_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

          //-- 38. ratio_GVL_cnb_ue_f

                r_ratio_GVL_cnb_ue_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_cnb_f[ind_f*dim_cnb_f[0] + 0*dim_cnb_f[1] + 0*dim_cnb_f[2] + ind_t*dim_cnb_f[3]] *
                     r_ue_f[ind_f*dim_ue_f[0] + 0*dim_ue_f[1] + 0*dim_ue_f[2] + ind_t*dim_ue_f[3]]);



        }

if (ind_t==0) UNPROTECT(64);
UNPROTECT(67);
//if (ind_t>10) PrintValue(VECTOR_ELT(out_EcoDCF, 44));
//if (ind_t>10) PrintValue(VECTOR_ELT(out_EcoDCF, 2));
//if (ind_t>10) PrintValue(VECTOR_ELT(fVar, 10));
//if (ind_t>10) PrintValue(VECTOR_ELT(fVar, 4));
//if (ind_t>10) PrintValue(vf_f_m);
//if (ind_t>10) PrintValue(ue_f_m);
//if (ind_t==(nbT-1)) {PrintValue(VECTOR_ELT(fVar,4));PrintValue(VECTOR_ELT(fVar,10));}

}}


//détection d'un caractère donné dans un objet SEXP de type SXPSTR

extern "C" {

bool isCharIn(SEXP names, const char *str)
{

    int i;
    bool test = false;

    for (i = 0; i < length(names); i++)
        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
            test = true;
            break;
        }

    return test;
}

}


extern "C" {
SEXP IAM(SEXP listInput, SEXP listSpec, SEXP listStochastic, SEXP listScen,
            SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, SEXP Bootstrp, SEXP nbBoot,
            SEXP GestInd, SEXP mF, SEXP mOth, SEXP bounds, SEXP TAC, SEXP FBAR, SEXP GestParam, SEXP EcoDcf,
            SEXP EcoInd, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyF, SEXP parBHV, SEXP parQEX, SEXP bootVar = R_NilValue)
{
//Rprintf("OO");
    if (INTEGER(Bootstrp)[0]==0) {

        BioEcoPar *object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                            RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                            GestInd, mF, mOth, bounds, TAC, FBAR, GestParam, EcoDcf,
                                            EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, parBHV, parQEX);

        SEXP output, out_names, out_Foth;
        PROTECT(output = allocVector(VECSXP, 23));
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
        SET_VECTOR_ELT(output, 13, object->out_L_efmct2);
        SET_VECTOR_ELT(output, 14, object->out_P_t);
        if (INTEGER(EcoDcf)[0]==0) SET_VECTOR_ELT(output, 15, object->out_Eco); else SET_VECTOR_ELT(output, 15, object->out_EcoDCF);
        PROTECT(out_Foth = allocVector(VECSXP, object->nbE));
        setAttrib(out_Foth, R_NamesSymbol, object->sppList);
        for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 44));
        SET_VECTOR_ELT(output, 16, out_Foth);
        SET_VECTOR_ELT(output, 17, object->mu_nbds);
        SET_VECTOR_ELT(output, 18, object->mu_nbv);
        SET_VECTOR_ELT(output, 19, object->out_effort);
        SET_VECTOR_ELT(output, 20, object->out_Fr_fmi);
        SET_VECTOR_ELT(output, 21, VECTOR_ELT(object->fVar, 29));
        SET_VECTOR_ELT(output, 22, object->out_PQuot_et);

        //on nomme les éléments de output
        const char *namesOut[23] = {"F","Z","Fbar","N","B","SSB","C","Ctot","Y","Ytot","D","Li","Lc","Lcm","P","E","Fothi","mu_nbds","mu_nbv","Eff","Fr","GVLoths_f","PQuot"};

        PROTECT(out_names = allocVector(STRSXP, 23));

        for(int ct = 0; ct < 23; ct++) SET_STRING_ELT(out_names, ct, mkChar(namesOut[ct]));

        setAttrib(output, R_NamesSymbol, out_names);

        UNPROTECT(3);
        return(output);
        delete object;


    } else {

        //on n'oublie pas d'activer les parties stochastiques pour que ça ait un sens

        //on commence par créer l'objet qui va accueillir la donnée  (3 outputs pour l'instant : biomasse, SSB, captures --> à développer selon les besoins)
        SEXP output, out_names, out_Foth, emptyObj;
        PROTECT(output = allocVector(VECSXP, 38)); //36
        SEXP eBoot;

        for (int ind = 0 ; ind < 38 ; ind++) { //36

            PROTECT(eBoot = allocVector(VECSXP, INTEGER(nbBoot)[0]));
            SET_VECTOR_ELT(output, ind, eBoot);

        }

        //on commence le bootstrap

        BioEcoPar *object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                    RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                    GestInd, mF, mOth, bounds, TAC, FBAR, GestParam, EcoDcf,
                                    EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, parBHV, parQEX);
        for (int it = 0 ; it < INTEGER(nbBoot)[0] ; it++) {

            if (it>0) object = new BioEcoPar(listInput, listSpec, listStochastic, listScen,
                                    RecType1, RecType2, RecType3, Scenarii, Bootstrp, nbBoot,
                                    GestInd, mF, mOth, bounds, TAC, FBAR, GestParam, EcoDcf,
                                    EcoInd, dr, SRind, listSR, TypeSR, mFM, TACbyF, parBHV, parQEX);

            //objet vide pour garder la structuration malgré la non-sélection de la variable en question
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
            for (int i = 0; i < object->nbE; i++) SET_VECTOR_ELT(out_Foth, i, VECTOR_ELT(VECTOR_ELT(object->eVar, i), 44));
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


            UNPROTECT(2);

        }

        //on nomme les éléments de output
        const char *namesOut[38] = {"B","SSB","Ctot","Ytot","Yfmi","Ffmi","Zeit","Fbar","Foth","mu_nbds","mu_nbv","N","Eff",
                                    "GVL_fme","GVLtot_fm","GVLav_f","rtbs_f","gp_f","ps_f","gcf_f","gva_f","cs_f","sts_f","rtbsAct_f",
                                    "csAct_f","gvaAct_f","gcfAct_f","psAct_f","stsAct_f","ccwCr_f","GVLtot_f","wagen_f","L_efmit","D_efmit","Fr_fmi","C_efmit","vcst_f","vcst_fm"};

        PROTECT(out_names = allocVector(STRSXP, 38));

        for(int ct = 0; ct < 38; ct++) SET_STRING_ELT(out_names, ct, mkChar(namesOut[ct]));

        setAttrib(output, R_NamesSymbol, out_names);

        UNPROTECT(2+38);
        return(output);
        delete object;

    }



}
}



////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
////------------------------------------------------------------------------------------------------
//
//
//using namespace boost;
//
//
////------------------------------------------
//// accesseur à un élément d'une liste donnée (list = liste en question , str {caractère} = intitulé de l'élément de la liste)
////------------------------------------------
//extern "C" {
//
//SEXP getListElt(SEXP list, const char *str)
//{
//    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
//    int i;
//
//    for (i = 0; i < length(list); i++)
//        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
//            elmt = VECTOR_ELT(list, i);
//            break;
//        }
//
//    return elmt;
//}
//
//}
//
//
////------------------------------------------
//// Transcripteur des sorties .txt de la fonction de conversion du fichier IN_IAM.r
//// Renvoie les paramètres au format IAM sous R (listes imbriquées)
//// Format standard des fichiers .txt input (séparateur "\t") :
//
//    // list	character	NA	NA	NA	NA	NA	NA	Langoustine	Merlu_commun	Sole_commune	Fleet
//    // list	character	NA	NA	NA	NA	NA	NA	modI	modL	modC	icat	alk	fm	mm	M_i	mat_i	wStock_i	wL_i	wD_i	N_it0	N_i0t	F_fmi	B_i	Y_mi	C_mi	Y_i	C_i	d_i	doth_i	sr	SelRef	P_fmce	alpha_fmce	beta_fmce	gamma_fmce	TAC	Fbar	FmaxTarget	Lref_f_e	Lref_f_m_e	GVLref_f_e	GVLref_f_m_e
//    // var	character	NA	NA	NA	NA	NA	NA	1	2	3	4	5	6	7	8	+gp
//    // var	double	0	0	0	0	NA	NA	NA
//    // var	character	NA	NA	NA	NA	NA	NA	10	20	30	40
//    // var	double	NA	NA	NA	NA	9	4	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	1	0	0	0	1	1	1	1	0	0	0	1	1	0	0	0	0	0	0	0
//    // var	double	0	0	0	0	NA	NA	NA
//
//// Descriptif :
//
//    // variable ou nouvelle liste imbriquée -- type de variable -- attributs DimCst (dim n°1 -- dim n°2 -- dim n°3 -- dim n°4) -- Codage Métier -- Dim Age ou Catégorie -- Elements énumérés de la variable ou des noms des éléments de la liste -- -- -- -- ...
//
////------------------------------------------
//
//
//
//extern "C" {
//SEXP Fun(SEXP File, SEXP Specific = R_NilValue) //file : character décrivant le path du fichier à transcrire  -----  specific : optionnel, la sortie de Fun appliquée au fichier 'specific.txt'
// {
//
//   SEXP imbricOBJ, file, specific;
//   PROTECT(file = duplicate(File));
//   PROTECT(specific = duplicate(Specific));
//   PROTECT(imbricOBJ = allocVector(VECSXP,5));  //on donne jusqu'à 5 niveaux d'imbrication pour l'objet
//   int IMAX[5]; //nbre d'élément à traiter
//   int I[5]; //élément traité par niveau
//   int n = 0; //niveau en cours de traitement
//
//   string ligne;
//   ifstream fichier(CHAR(STRING_ELT(file,0)));
//   if (!fichier) error("Can't read file !!\n");
//
//    vector< string > vec;
//    SEXP out, names, var, dimcst, Dim, DimNam;
//
//
//    //étape d'initialisation --> 1ère ligne
//
//    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
//    boost::char_separator<char> sep("\t" );
//
//    getline (fichier, ligne);
//    tokenizer tokens(ligne, sep);
//
//    vec.clear();
//    vec.assign(tokens.begin(),tokens.end());
//    std::vector<char const*> v( vec.size() );
//    for( size_t i = 8; i < v.size(); ++i ) v[i] = vec[i].c_str();
//
//    IMAX[n] = vec.size()-8;
//    I[n] = 0;
//    PROTECT(out = allocVector(VECSXP,IMAX[n]));
//    //on ajoute les noms des éléments
//    PROTECT(names = allocVector(STRSXP,IMAX[n]));
//    int count=0;
//    for( size_t i = 8; i < v.size(); ++i ) {SET_STRING_ELT(names,count,mkChar(v[i]));count++;}
//    setAttrib(out, R_NamesSymbol, names);
//    SET_VECTOR_ELT(imbricOBJ,0,out);
//
//
//
//   while (getline (fichier, ligne))
//   {
//
//     // tokenisation de la ligne courante
//     tokenizer tokens(ligne, sep);
//
//        vec.clear();
//        vec.assign(tokens.begin(),tokens.end());
//
//        std::vector<char const*> v( vec.size() );
//
//        for( size_t i = 8; i < v.size(); ++i ) v[i] = vec[i].c_str();
//
//        //1ère colonne qui détermine si on doit formater ou remplir
//
//        if (strcmp(vec[0].c_str(), "list") == 0) {
//
//            //on monte d'un niveau
//            n++;
//            //on fixe le nouveau imax pour le niveau en cours
//            IMAX[n] = vec.size()-8; //le premier élément de la ligne est seulement descriptif
//            //on initialise le compteur pour ce niveau
//            I[n] = 0;
//            //on crée une copie de l'objet à modifier après l'avoir intégré à l'output
//            SET_VECTOR_ELT(VECTOR_ELT(imbricOBJ,n-1),I[n-1],allocVector(VECSXP,IMAX[n]));
//            SET_VECTOR_ELT(imbricOBJ,n,VECTOR_ELT(VECTOR_ELT(imbricOBJ,n-1),I[n-1]));
//
//            PROTECT(names = allocVector(STRSXP,IMAX[n]));
//            int count=0;
//            for(size_t i = 8; i < v.size(); ++i ) {SET_STRING_ELT(names,count,mkChar(v[i]));count++;}
//            setAttrib(VECTOR_ELT(imbricOBJ,n), R_NamesSymbol, names);
//            UNPROTECT(1);
//
//            //on pointe sur cet élément dorénavant
//            I[n-1]++;
//
//        } else {
//
//            if (strcmp(vec[0].c_str(), "var") == 0) {
//
//            //on acccède à l'objet à modifier et on assigne un vecteur 'character' de taille donnée
//            PROTECT(var = allocVector(STRSXP,vec.size()-8));
//            //et on remplit
//            int count=0;
//            for(size_t i = 8; i < v.size(); ++i ) {SET_STRING_ELT(var,count,mkChar(v[i]));count++;}
//
//            //on prévoit d'ores et déjà le reformatage --> attribut 'DimCst'
//            if (strcmp(vec[2].c_str(), "NA") != 0) {
//
//               if (specific == NULL) error("missing 'specific' input object!!");
//
//               PROTECT(dimcst = allocVector(STRSXP,4)) ;
//               int count=0;
//               for(size_t i = 2; i < 6; ++i ) {SET_STRING_ELT(dimcst,count,mkChar(vec[i].c_str()));count++;}
//               int *dimC = INTEGER(AS_INTEGER(dimcst));
//               setAttrib(var, install("DimCst"),AS_INTEGER(dimcst));
//
//                //dimensions
//                int nbDim = 0; for (int j=0; j<4; j++) if (dimC[j]>0) nbDim++;
//
//                //si nbDim<1, pas de formatage
//
//                //si nbDim>1 formatage matriciel
//                if (nbDim>1) {
//
//                    PROTECT(Dim = allocVector(INTSXP,nbDim));
//                    int *dim = INTEGER(Dim); int index=0;
//                    for (int j=0; j<4; j++) if (dimC[j]>0) {dim[index] = dimC[j];index++;}
//                    setAttrib(var,R_DimSymbol,Dim);
//
//                    //noms de dimensions : si nbDim<1, pas de noms
//                    PROTECT(DimNam = allocVector(VECSXP,nbDim));
//                    //et on remplit en fonction des colonnes Dim et de l'objet Specific en input
//                    int rank=0;
//
//                    if (strcmp(vec[2].c_str(), "0") != 0) {
//                        SET_VECTOR_ELT(DimNam,rank,getListElt(specific, "Fleet"));rank++;
//                    }
//
//                    if (strcmp(vec[3].c_str(), "0") != 0) {
//                        if (strcmp(vec[6].c_str(), "B") == 0) {
//                        SET_VECTOR_ELT(DimNam,rank,getListElt(specific, "Metier"));rank++;
//                        } else {
//                        SET_VECTOR_ELT(DimNam,rank,getListElt(specific, "MetierEco"));rank++;
//                        }
//                    }
//
//                    if (strcmp(vec[4].c_str(), "0") != 0) { //il y a alors une dimension Espèce
//                        //on part du niveau 'n' et on remonte dans imbricObj pour les éléments I pour retrouver une modalité "espèce"
//                        int indSp = length(getListElt(specific, "Ages"));
//                        int flag = -1;
//
//                        for (int ii = n; ii>=0 ; ii--) {
//                            for (int J = 0; J<indSp; J++) {
//                                int val = I[ii]-1;
//                                if (ii==n) val++;
//
//                                if (strcmp( CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),J)) ,
//                                            CHAR(STRING_ELT(getAttrib(VECTOR_ELT(imbricOBJ,ii), R_NamesSymbol),val))) == 0) {
//                                    flag = J;
//                                    break;
//                                }
//                            }
//                            if (flag>=0) break;
//                        }
//
//                        if (strcmp(vec[7].c_str(), "A") == 0) {
//                            SET_VECTOR_ELT(DimNam,rank,
//                                            getListElt(getListElt(specific, "Ages"),
//                                                CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),flag))));
//                            rank++;
//                        } else {
//                            SET_VECTOR_ELT(DimNam,rank,
//                                            getListElt(getListElt(specific, "Cat"),
//                                                CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),flag))));
//                            rank++;
//                        }
//                    }
//
//                    if (strcmp(vec[5].c_str(), "0") != 0) {
//                        SET_VECTOR_ELT(DimNam,rank,AS_CHARACTER(getListElt(specific, "times")));rank++;
//                    }
//
//                    //et on assigne le résultat
//                    setAttrib(var,R_DimNamesSymbol,DimNam);
//
//                    UNPROTECT(2);
//                }
//
//                //si nbDim==1 formatage vectoriel
//                if (nbDim==1) {
//
//                    if (strcmp(vec[2].c_str(), "0") != 0) {
//                        setAttrib(var,R_NamesSymbol,getListElt(specific, "Fleet"));
//                    }
//
//                    if (strcmp(vec[3].c_str(), "0") != 0) {
//                        if (strcmp(vec[6].c_str(), "B") == 0) {
//                            setAttrib(var,R_NamesSymbol,getListElt(specific, "Metier"));
//                        } else {
//                            setAttrib(var,R_NamesSymbol,getListElt(specific, "MetierEco"));
//                        }
//                    }
//
//                    if (strcmp(vec[4].c_str(), "0") != 0) {
//                        //on part du niveau 'n' et on remonte dans imbricObj pour les éléments I pour retrouver une modalité "espèce"
//                        int indSp = length(getListElt(specific, "Ages"));
//                        int flag = -1;
//                        for (int ii = n; ii>=0 ; ii--) {
//                            for (int J = 0; J<indSp; J++) {
//
//                                if (strcmp( CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),J)) ,
//                                            CHAR(STRING_ELT(getAttrib(VECTOR_ELT(imbricOBJ,ii), R_NamesSymbol),I[ii]-1))) == 0) {
//                                    flag = J;
//                                    break;
//                                }
//                            }
//                            if (flag>=0) break;
//                        }
//
//                        if (strcmp(vec[7].c_str(), "A") == 0) {
//                            setAttrib(var,R_NamesSymbol,
//                                            getListElt(getListElt(specific, "Ages"),
//                                                CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),flag))));
//                        } else {
//                            setAttrib(var,R_NamesSymbol,
//                                            getListElt(getListElt(specific, "Cat"),
//                                                CHAR(STRING_ELT(getAttrib(getListElt(specific, "Ages"),R_NamesSymbol),flag))));
//                        }
//                    }
//
//                    if (strcmp(vec[5].c_str(), "0") != 0) {
//                        setAttrib(var,R_NamesSymbol,AS_CHARACTER(getListElt(specific, "times")));
//                    }
//                }
//
//                UNPROTECT(1);
//
//            } else { //matrice icat ou alk
//
//                PROTECT(Dim = allocVector(INTSXP,2));
//                int *dim = INTEGER(Dim);
//                PROTECT(DimNam = allocVector(VECSXP,2));
//
//                if (strcmp( CHAR(STRING_ELT(getAttrib(VECTOR_ELT(imbricOBJ,n), R_NamesSymbol),I[n])), "icat") == 0) {
//
//                    dim[0] = length(getListElt(VECTOR_ELT(imbricOBJ,n), "modI"));
//                    dim[1] = length(getListElt(VECTOR_ELT(imbricOBJ,n), "modC"));
//                    SET_VECTOR_ELT(DimNam,0,getListElt(VECTOR_ELT(imbricOBJ,n), "modI"));
//                    SET_VECTOR_ELT(DimNam,1,getListElt(VECTOR_ELT(imbricOBJ,n), "modC"));
//
//                    setAttrib(var,R_DimSymbol,Dim);
//                    setAttrib(var,R_DimNamesSymbol,DimNam);
//                }
//
//                if (strcmp( CHAR(STRING_ELT(getAttrib(VECTOR_ELT(imbricOBJ,n), R_NamesSymbol),I[n])), "alk") == 0) {
//
//                    dim[0] = length(getListElt(VECTOR_ELT(imbricOBJ,n), "modL"));
//                    dim[1] = length(getListElt(VECTOR_ELT(imbricOBJ,n), "modI"));
//                    SET_VECTOR_ELT(DimNam,0,getListElt(VECTOR_ELT(imbricOBJ,n), "modL"));
//                    SET_VECTOR_ELT(DimNam,1,getListElt(VECTOR_ELT(imbricOBJ,n), "modI"));
//
//                    setAttrib(var,R_DimSymbol,Dim);
//                    setAttrib(var,R_DimNamesSymbol,DimNam);
//                }
//
//                UNPROTECT(2);
//
//            }
//
//            //on insère le résultat dans 'out' en fonction du type de sortie
//            if (strcmp(vec[1].c_str(), "double") == 0) {
//
//                SET_VECTOR_ELT(VECTOR_ELT(imbricOBJ,n),I[n],AS_NUMERIC(var));
//
//            } else {
//
//                if (strcmp(vec[1].c_str(), "integer") == 0) {
//
//                    SET_VECTOR_ELT(VECTOR_ELT(imbricOBJ,n),I[n],AS_INTEGER(var));
//
//                } else {
//
//                    SET_VECTOR_ELT(VECTOR_ELT(imbricOBJ,n),I[n],var);
//
//                }
//            }
//
//            UNPROTECT(1);
//            I[n]++;
//            while (I[n]>=IMAX[n]) n--;
//
//        }
//        }
//
//   }
//    fichier.close();
//    UNPROTECT(5);
//    return out;
//}
//}
//
//
//
////------------------------------------------
//// fonction similaire à IAM, mais appelant des fichiers .txt de paramètres créés à partir des fonctions 'unl'
////------------------------------------------
//
//extern "C" {
//SEXP IAM_txt(SEXP fileParam, SEXP fileSpec, SEXP fileInput, SEXP fileScen /*= R_NilValue*/, SEXP fileStoch /*= R_NilValue*/) {
//
//
//SEXP outp, listInput_txt, listSpec_txt, listStoch_txt, listScen_txt, listParam_txt, inScen;
//
//
//PROTECT(listSpec_txt = Fun(fileSpec));
//PROTECT(listInput_txt = Fun(fileInput,listSpec_txt));
//PROTECT(listParam_txt = Fun(fileParam,listSpec_txt));
//if (fileScen == NULL) {
//        PROTECT(listScen_txt = R_NilValue);
//} else {
//        PROTECT(listScen_txt = Fun(fileScen,listSpec_txt));
//}
//
//if (fileStoch == NULL) {
//        PROTECT(listStoch_txt = R_NilValue);
//} else {
//        PROTECT(listStoch_txt = Fun(fileStoch,listSpec_txt));
//}
//
//
//if (length(getListElt(listParam_txt,"scenario"))==0) {
//    PROTECT(inScen = R_NilValue);
//} else {
//    PROTECT(inScen = getListElt(listScen_txt, CHAR(STRING_ELT(getListElt(listParam_txt,"scenario"),0))));
//}
//
//PROTECT(outp = IAM(listInput_txt,
//                  listSpec_txt,
//                  listStoch_txt,
//                  inScen,
//                  getListElt(listParam_txt,"RecType1"),
//                  getListElt(listParam_txt,"RecType2"),
//                  getListElt(listParam_txt,"RecType3"),
//                  getListElt(listParam_txt,"Scenarii"),
//                  getListElt(listParam_txt,"Bootstrp"),
//                  getListElt(listParam_txt,"nbBoot"),
//                  getListElt(listParam_txt,"GestInd"),
//                  getListElt(listParam_txt,"mF"),
//                  getListElt(listParam_txt,"mOth"),
//                  getListElt(listParam_txt,"bounds"),
//                  getListElt(listParam_txt,"TAC"),
//                  getListElt(listParam_txt,"FBAR"),
//                  getListElt(listParam_txt,"GestParam"),
//                  getListElt(listParam_txt,"EcoDcf"),
//                  getListElt(listParam_txt,"EcoInd"),
//                  getListElt(listParam_txt,"dr"),
//                  getListElt(listParam_txt,"SRind"),
//                  getListElt(listParam_txt,"listSR"),
//                  getListElt(listParam_txt,"TypeSR"),
//                  getListElt(listParam_txt,"mFM"),
//                  getListElt(listParam_txt,"bootVar")
//                  )
//        );
//
//
//UNPROTECT(7);
//
//return(outp);
//
//}
//}
//
//
////------------------------------------------
//// fonction d'exportation d'une variable de sortie au format data.frame, dans un fichier .txt avec séparateurs '\t'
////------------------------------------------
//
//SEXP IDim(int *dimInput) {
//
//    SEXP Tab;
//    PROTECT(Tab = allocVector(INTSXP,4));
//    int *tab = INTEGER(Tab);
//
//    tab[0] = (dimInput[0]>0);
//    tab[1] = (dimInput[1]>0)*(1 + (dimInput[0]-1)*(dimInput[0]>0));
//    tab[2] = (dimInput[2]>0)*(1 + (dimInput[1]-1)*(dimInput[1]>0))*(1 + (dimInput[0]-1)*(dimInput[0]>0));
//    tab[3] = (dimInput[3]>0)*(1 + (dimInput[2]-1)*(dimInput[2]>0))*(1 + (dimInput[1]-1)*(dimInput[1]>0))*(1 + (dimInput[0]-1)*(dimInput[0]>0));
//
//    UNPROTECT(1);
//    return(Tab);
//
//}
//
//
//extern "C" {
//SEXP IAM_export(SEXP vrbl, SEXP fileExp, SEXP replic, SEXP species) { // vrbl : variable de sortie
//                                                             // fileExp : nom du fichier en sortie
//                                                             // rep : 1/0 itérations ou non
//                                                             // spp : 1/0 par espèce ou non
//
//    SEXP vrblType, val = R_NilValue, dimnam = R_NilValue, namSpp = R_NilValue;
//    int nbIter = 0, nbSpp = 0, rep = INTEGER(replic)[0], spp = INTEGER(species)[0];
//    int *dimCst, *ind;
//    int index[4];
//    double *values;
//
//
//    string const fichier(CHAR(STRING_ELT(fileExp,0)));
//    ofstream flux(fichier.c_str());
//    if (!flux) error("Can't open file!!\n");
//
//    if (rep) {
//        nbIter = length(vrbl);
//        if (spp) {
//            nbSpp = length(VECTOR_ELT(vrbl,0));//Rprintf("A1 ");
//            PROTECT(namSpp=getAttrib(VECTOR_ELT(vrbl,0),R_NamesSymbol));//Rprintf("A2 ");
//            PROTECT(vrblType=VECTOR_ELT(VECTOR_ELT(vrbl,0),0));//Rprintf("A3 ");
//        } else {
//            PROTECT(vrblType=VECTOR_ELT(vrbl,0));//Rprintf("A4 ");
//        }
//    } else {
//        if (spp) {
//            nbSpp = length(vrbl);//Rprintf("A5 ");
//            PROTECT(namSpp=getAttrib(vrbl,R_NamesSymbol));//Rprintf("A6 ");
//            PROTECT(vrblType=VECTOR_ELT(vrbl,0));//Rprintf("A7 ");
//        } else {
//            PROTECT(vrblType=vrbl);//Rprintf("A8 \n");
//        }
//    }
////Rprintf("AA\n");
//    dimCst = INTEGER(getAttrib(vrblType, install("DimCst")));
////Rprintf("BB\n");
//
//
////Rprintf("CC\n");
//
//    if (flux)
//    {   //en-têtes
//        if (rep) flux << "iter" << '\t';
//        if (spp) flux << "spp" << '\t';
//        if (dimCst[0]>0) flux << "fleet" << '\t';
//        if (dimCst[1]>0) flux << "metier" << '\t';
//        if (dimCst[2]>0) flux << "age" << '\t';
//        if (dimCst[3]>0) flux << "year" << '\t';
//        flux << "value" << endl;
//
//        for (int it = 0; it < imax2(nbIter,1); it++)
//            for (int sp = 0; sp < imax2(nbSpp,1); sp++) {
//
//                if (rep & spp) PROTECT(val=VECTOR_ELT(VECTOR_ELT(vrbl,it),sp));//Rprintf("DD1 ");
//
//                if (rep & !spp) PROTECT(val=VECTOR_ELT(vrbl,it));//Rprintf("DD2 ");
//
//                if (!rep & spp) PROTECT(val=VECTOR_ELT(vrbl,sp));//Rprintf("DD3 ");
//
//                if (!rep & !spp) PROTECT(val=vrbl);//Rprintf("DD4 \n");
//
//                dimCst = INTEGER(getAttrib(val, install("DimCst")));
//                ind = INTEGER(IDim(dimCst));
//
//                //on veut les indices des dimensions selon l'attribut dimCst
//                int nb = 0;
//                for (int i=0; i<4; i++)  {
//
//                    if (dimCst[i]>0) {
//
//                        index[i] = nb; nb++;
//
//                    } else {
//
//                        index[i] = -1;
//                    }
//                }
//                    //on en déduit les intitulés de la variable
//                if (nb>1) {
//                    PROTECT(dimnam = getAttrib(val,R_DimNamesSymbol));
//                } else {
//                    PROTECT(dimnam = getAttrib(val,R_NamesSymbol));
//                }
//
//                values = REAL(val);
//
//                for (int fl = 0; fl < imax2(dimCst[0],1); fl++)
//                    for (int met = 0; met < imax2(dimCst[1],1); met++)
//                        for (int ag = 0; ag < imax2(dimCst[2],1); ag++)
//                            for (int yr = 0; yr < imax2(dimCst[3],1); yr++) {
//
//                                if (rep & spp) flux << it+1 << '\t' << CHAR(STRING_ELT(namSpp,sp)) << '\t';
//
//                                if (rep & !spp) flux << it+1 << '\t';
//
//                                if (!rep & spp) flux << CHAR(STRING_ELT(namSpp,sp)) << '\t';
//
//                                if (dimCst[0]>0) {
//                                    if (nb>1) {
//                                       flux << CHAR(STRING_ELT(VECTOR_ELT(dimnam,index[0]),fl)) << '\t';
//                                    } else {
//                                       flux << CHAR(STRING_ELT(dimnam,fl)) << '\t';
//                                    }
//                                }
//
//                                if (dimCst[1]>0) {
//                                    if (nb>1) {
//                                       flux << CHAR(STRING_ELT(VECTOR_ELT(dimnam,index[1]),met)) << '\t';
//                                    } else {
//                                       flux << CHAR(STRING_ELT(dimnam,met)) << '\t';
//                                    }
//                                }
//
//                                if (dimCst[2]>0) {
//                                    if (nb>1) {
//                                       flux << CHAR(STRING_ELT(VECTOR_ELT(dimnam,index[2]),ag)) << '\t';
//                                    } else {
//                                       flux << CHAR(STRING_ELT(dimnam,ag)) << '\t';
//                                    }
//                                }
//
//                                if (dimCst[3]>0) {
//                                    if (nb>1) {
//                                       flux << CHAR(STRING_ELT(VECTOR_ELT(dimnam,index[3]),yr)) << '\t';
//                                    } else {
//                                       flux << CHAR(STRING_ELT(dimnam,yr)) << '\t';
//                                    }
//                                }
//
//                                flux << values[fl*ind[0] + met*ind[1] + ag*ind[2] + yr*ind[3]] << endl;
//
//                                }
//
//                     UNPROTECT(2);
//
//                }
//    }
//
//    UNPROTECT(1);
//    if (spp) UNPROTECT(1);
//
//    return(fileExp);
//
//}
//}
//
////------------------------------------------------------------------------------
////-------------------------------------------------------------------------------
////------------------------------------------------------------------------------
////-------------------------------------------------------------------------------
////------------------------------------------------------------------------------
////-------------------------------------------------------------------------------
//
//
//extern "C" {
//
//SEXP IAM_txtIN_txtOUT(){
//
//
//SEXP nmsIN, nmsOUT, result, reP, spP, out;
//SEXP tempIN, tempOUT, spec, args;
//int *sppInt, *repInt;
//bool is_scen, is_sto;
//
//
////const char *inputFiles[5] = {"C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Input\\argsCPP.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Input\\specific.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Input\\input.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Input\\scenario.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Input\\stochastic.txt"};
//
//const char *inputFiles[5] = {"C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Input\\argsCPP.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Input\\specific.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Input\\input.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Input\\scenario.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Input\\stochastic.txt"};
//
//PROTECT(nmsIN = allocVector(VECSXP, 5));
//
//for(int k = 0; k < 5; k++) {
//    PROTECT(tempIN = allocVector(STRSXP,1));
//    SET_STRING_ELT(tempIN, 0, mkChar(inputFiles[k]));
//    SET_VECTOR_ELT(nmsIN,k,tempIN);
//    UNPROTECT(1);
//}
//
//
////const char *outputFiles[11] = {"C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\SSB.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\Fbar.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\Ctot.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\Ytot.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\Y.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\L.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\D.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\GVL.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\GVA.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\GCF.txt",
////                            "C:\\Documents and Settings\\mmerzere\\Bureau\\COST_R\\IAMwdSIAD\\Output\\PS.txt"};
//
//const char *outputFiles[11] = {"C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\SSB.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\Fbar.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\Ctot.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\Ytot.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\Y.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\L.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\D.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\GVL.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\GVA.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\GCF.txt",
//                            "C:\\Documents and Settings\\mmerzere\\Bureau\\SiteSIAD_Report\\Output\\PS.txt"};
//
//
//PROTECT(nmsOUT = allocVector(VECSXP, 11));
//
//for(int k = 0; k < 11; k++) {
//    PROTECT(tempOUT = allocVector(STRSXP,1));
//    SET_STRING_ELT(tempOUT, 0, mkChar(outputFiles[k]));
//    SET_VECTOR_ELT(nmsOUT,k,tempOUT);
//    UNPROTECT(1);
//}
//
//ifstream is_scenar(inputFiles[3]);//est-ce que le fichier 'scenario.txt' existe?
//ifstream is_stochastic(inputFiles[4]);//est-ce que le fichier 'stochastic.txt' existe?
//
//if (is_scenar) is_scen=true; else is_scen=false; is_scenar.close();
//if (is_stochastic) is_sto=true; else is_sto=false; is_stochastic.close();
//
////on doit determiner à partir du fichier 'arguments' si on est en présence de réplicats ou non
//PROTECT(spec = Fun(VECTOR_ELT(nmsIN,1)));
//PROTECT(args = Fun(VECTOR_ELT(nmsIN,0),spec));
//PROTECT(reP = getListElt(args,"Bootstrp")); repInt = INTEGER(reP);
//
//if (is_scen) {
//    if (is_sto) {
//        PROTECT(result = IAM_txt(VECTOR_ELT(nmsIN,0), VECTOR_ELT(nmsIN,1), VECTOR_ELT(nmsIN,2), VECTOR_ELT(nmsIN,3), VECTOR_ELT(nmsIN,4)));
//    } else {
//        PROTECT(result = IAM_txt(VECTOR_ELT(nmsIN,0), VECTOR_ELT(nmsIN,1), VECTOR_ELT(nmsIN,2), VECTOR_ELT(nmsIN,3), NULL));
//    }
//} else {
//    if (is_sto) {
//        PROTECT(result = IAM_txt(VECTOR_ELT(nmsIN,0), VECTOR_ELT(nmsIN,1), VECTOR_ELT(nmsIN,2), NULL, VECTOR_ELT(nmsIN,4)));
//    } else {
//        PROTECT(result = IAM_txt(VECTOR_ELT(nmsIN,0), VECTOR_ELT(nmsIN,1), VECTOR_ELT(nmsIN,2), NULL, NULL));
//    }
//}
//
////on exporte maintenant les résultats
//
//PROTECT(spP = allocVector(INTSXP,1)); sppInt = INTEGER(spP); //à remettre à jour pour chacun des cas (chaque variable)
//
////SSB
//sppInt[0] = 1; //variable Sp -> double contrôle à opérer sur la disponibilité de la donnée
//if (repInt[0]>0) {
//    if (length(getListElt(result,"SSB"))>0 & length(VECTOR_ELT(getListElt(result,"SSB"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"SSB"),0),0))>0)
//        out = IAM_export(getListElt(result,"SSB"),  VECTOR_ELT(nmsOUT,0), reP, spP);
//} else {
//    if (length(getListElt(result,"SSB"))>0 & length(VECTOR_ELT(getListElt(result,"SSB"),0))>0)
//        out = IAM_export(getListElt(result,"SSB"),  VECTOR_ELT(nmsOUT,0), reP, spP);
//}
//
////Fbar
//sppInt[0] = 1; // <- intutile mais à but illustratif
//if (repInt[0]>0) {
//    if (length(getListElt(result,"Fbar"))>0 & length(VECTOR_ELT(getListElt(result,"Fbar"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"Fbar"),0),0))>0)
//        out = IAM_export(getListElt(result,"Fbar"),  VECTOR_ELT(nmsOUT,1), reP, spP);
//} else {
//    if (length(getListElt(result,"Fbar"))>0 & length(VECTOR_ELT(getListElt(result,"Fbar"),0))>0)
//        out = IAM_export(getListElt(result,"Fbar"),  VECTOR_ELT(nmsOUT,1), reP, spP);
//}
//
////Ctot
//sppInt[0] = 1;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"Ctot"))>0 & length(VECTOR_ELT(getListElt(result,"Ctot"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"Ctot"),0),0))>0)
//        out = IAM_export(getListElt(result,"Ctot"),  VECTOR_ELT(nmsOUT,2), reP, spP);
//} else {
//    if (length(getListElt(result,"Ctot"))>0 & length(VECTOR_ELT(getListElt(result,"Ctot"),0))>0)
//        out = IAM_export(getListElt(result,"Ctot"),  VECTOR_ELT(nmsOUT,2), reP, spP);
//}
//
////Ytot
//sppInt[0] = 1;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"Ytot"))>0 & length(VECTOR_ELT(getListElt(result,"Ytot"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"Ytot"),0),0))>0)
//        out = IAM_export(getListElt(result,"Ytot"),  VECTOR_ELT(nmsOUT,3), reP, spP);
//} else {
//    if (length(getListElt(result,"Ytot"))>0 & length(VECTOR_ELT(getListElt(result,"Ytot"),0))>0)
//        out = IAM_export(getListElt(result,"Ytot"),  VECTOR_ELT(nmsOUT,3), reP, spP);
//}
//
////Y
//sppInt[0] = 1;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"Yfmi"))>0 & length(VECTOR_ELT(getListElt(result,"Yfmi"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"Yfmi"),0),0))>0)
//        out = IAM_export(getListElt(result,"Yfmi"),  VECTOR_ELT(nmsOUT,4), reP, spP);
//} else {
//    if (length(getListElt(result,"Y"))>0 & length(VECTOR_ELT(getListElt(result,"Y"),0))>0)
//        out = IAM_export(getListElt(result,"Y"),  VECTOR_ELT(nmsOUT,4), reP, spP);
//}
//
//
////L
//sppInt[0] = 1;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"L_efmit"))>0 & length(VECTOR_ELT(getListElt(result,"L_efmit"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"L_efmit"),0),0))>0)
//        out = IAM_export(getListElt(result,"L_efmit"),  VECTOR_ELT(nmsOUT,5), reP, spP);
//} else {
//    if (length(getListElt(result,"Li"))>0 & length(VECTOR_ELT(getListElt(result,"Li"),0))>0)
//        out = IAM_export(getListElt(result,"Li"),  VECTOR_ELT(nmsOUT,5), reP, spP);
//}
//
//
////D
//sppInt[0] = 1;
//if (repInt[0]>0) {
//     if (length(getListElt(result,"D_efmit"))>0 & length(VECTOR_ELT(getListElt(result,"D_efmit"),0))>0 & length(VECTOR_ELT(VECTOR_ELT(getListElt(result,"D_efmit"),0),0))>0)
//        out = IAM_export(getListElt(result,"D_efmit"),  VECTOR_ELT(nmsOUT,6), reP, spP);
//} else {
//    if (length(getListElt(result,"D"))>0 & length(VECTOR_ELT(getListElt(result,"D"),0))>0)
//        out = IAM_export(getListElt(result,"D"),  VECTOR_ELT(nmsOUT,6), reP, spP);
//}
//
//
////GVLav
//sppInt[0] = 0;
//if (repInt[0]>0) {
//     if (length(getListElt(result,"GVLav_f"))>0 & length(VECTOR_ELT(getListElt(result,"GVLav_f"),0))>0)
//        out = IAM_export(getListElt(result,"GVLav_f"),  VECTOR_ELT(nmsOUT,7), reP, spP);
//} else {
//    if (length(getListElt(getListElt(result,"E"),"GVLav_f"))>0)
//        out = IAM_export(getListElt(getListElt(result,"E"),"GVLav_f"),  VECTOR_ELT(nmsOUT,7), reP, spP);
//}
//
//
////GVA
//sppInt[0] = 0;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"gva_f"))>0 & length(VECTOR_ELT(getListElt(result,"gva_f"),0))>0)
//        out = IAM_export(getListElt(result,"gva_f"),  VECTOR_ELT(nmsOUT,8), reP, spP);
//} else {
//    if (length(getListElt(getListElt(result,"E"),"gva_f"))>0)
//        out = IAM_export(getListElt(getListElt(result,"E"),"gva_f"),  VECTOR_ELT(nmsOUT,8), reP, spP);
//}
//
//
////GCF
//sppInt[0] = 0;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"gcf_f"))>0 & length(VECTOR_ELT(getListElt(result,"gcf_f"),0))>0)
//        out = IAM_export(getListElt(result,"gcf_f"),  VECTOR_ELT(nmsOUT,9), reP, spP);
//} else {
//    if (length(getListElt(getListElt(result,"E"),"gcf_f"))>0)
//        out = IAM_export(getListElt(getListElt(result,"E"),"gcf_f"),  VECTOR_ELT(nmsOUT,9), reP, spP);
//}
//
////PS
//sppInt[0] = 0;
//if (repInt[0]>0) {
//    if (length(getListElt(result,"ps_f"))>0 & length(VECTOR_ELT(getListElt(result,"ps_f"),0))>0)
//        out = IAM_export(getListElt(result,"ps_f"),  VECTOR_ELT(nmsOUT,10), reP, spP);
//} else {
//    if (length(getListElt(getListElt(result,"E"),"ps_f"))>0)
//        out = IAM_export(getListElt(getListElt(result,"E"),"ps_f"),  VECTOR_ELT(nmsOUT,10), reP, spP);
//}
//
//
//UNPROTECT(7);
//
////return(R_NilValue);
//return(result);
//
//}
//}
//

//------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------------------------------------------------------------


int main()
{

    return 0;
}

