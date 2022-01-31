#ifndef PARAM_H_INCLUDED
#define PARAM_H_INCLUDED

//classe regroupant les parametres a integrer dans le modele, les variables intermediaires (discretisation du processus),
//ainsi que les sorties du modele
class BioEcoPar
{
public: //normalement, selon les conventions,  les attributs doivent �tre "private"

typedef  double (BioEcoPar::*BEfn1)(double mult);
typedef  double (BioEcoPar::*BEfn1_F)(double *x);

//   INPUTS   ------------------

SEXP    list;       //liste d'objets R constituant la donn�e entr�e du mod�le. Certaines variables seront remises � jour dans le cadre
                    //de mesures de gestion ou de mod�lisation de comportement.


//   OUTPUTS  ------------------

//outputs des diff�rents modules (� initialiser pour t=0)
SEXP    out_F_fmi,  //mortalit� "captures" par p�che (par esp�ce)
        out_Fr_fmi,  //mortalit� totale (corrig�e de la survie) par p�che (par esp�ce)
        out_Z_eit,  //coefficient de mortalit� totale
        out_Fbar_et,  //Fbar par esp�ce (t)
        out_N_eit,  //effectifs en nombre
        out_B_et,   //biomasse (t)
        out_SSB_et, //biomasse de reproducteurs (t)
        out_C_efmit,//captures en nombres
        out_C_eit,  //captures totales en nombres
        out_Y_efmit,//captures en poids (t)
        out_Y_eit,  //captures totales en poids (t)
        out_D_efmit,//rejets en poids (t)
        out_L_efmit,//d�barquements en poids aux �ges(t)
        out_L_efmct,//d�barquements en poids par cat�gories(t)
        out_L_eit,//d�barquements en poids par cat�gories(t) pour le codage m�tier eco
        out_L_et , //debarquements totaux par espece
        out_L_pt , // debarquements totaux par produit marche

        out_F_fmi_G1,  //mortalit� "captures" par p�che (par esp�ce)
        out_F_fmi_G2,
        out_Fr_fmi_G1,  //mortalit� totale (corrig�e de la survie) par p�che (par esp�ce)
        out_Fr_fmi_G2,
        out_Z_eit_G1,  //coefficient de mortalit� totale
        out_Z_eit_G2,
        out_N_eit_G1,  //effectifs en nombre
        out_N_eit_G2,
        out_C_efmit_G1,
        out_C_efmit_G2,
        out_C_eit_G1,
        out_C_eit_G2,
        out_D_efmit_G1,
        out_D_efmit_G2,

        out_oqD_eft,//rejets over-quotas par flottille (esp�ces dynamiques)
        out_oqD_et,//rejets over-quotas total (esp�ces dynamiques)

        out_Ystat,  //captures totales en poids (t) pour les esp�ces statiques
        out_Lstat,  //d�barquements totaux en poids (t) pour les esp�ces statiques
        out_Dstat,  //rejets totaux en poids (t) pour les esp�ces statiques

        out_oqDstat, //rejets over-quotas par flottille esp�ces statiques

        out_P_t,    //prix moyen (en euros) (niveau m�tier �co si dispo)
        out_Pstat,    //prix moyen (en euros) (niveau m�tier �co si dispo) pour les esp�ces statiques
        out_CA_eft, //chiffre d'affaires par esp�ce dynamique (en euros)
        out_CAstat, //chiffre d'affaires par esp�ce statique (en euros)
        out_CAT_ft, //chiffre d'affaires total (en euros)
        out_CA_ft,  //chiffre d'affaires moyen par navire (en euros)
        out_RAP_ft, //reste � partager (en euros)
        out_EBE_ft, //exc�dent brut d'exploitation (en euros)
        out_ENE_ft, //exc�dent net d'exploitation (en euros)
        out_SA_ft,  //surplus du capital (en euros)
        out_PS_t,   //surplus producteur (en euros)
        out_ES_t,   //surplus de l'Etat (en euros)
        out_Eco,
        out_EcoDCF,
        out_effort, //variables d'effort utilis�es lors de la simulation
        out_allocEff_fm,
        out_SRmod,
        out_PQuot_et, out_QuotaTrade_fe, out_diffLQ, out_PQuot_temp,
        out_typeGest,
        out_F_fmi_S1M1, out_F_fmi_S1M2, out_F_fmi_S1M3, out_F_fmi_S1M4, out_F_fmi_S2M1, out_F_fmi_S2M2, out_F_fmi_S2M3, out_F_fmi_S2M4,
        out_F_fmi_S3M1, out_F_fmi_S3M2, out_F_fmi_S3M3, out_F_fmi_S3M4, out_F_fmi_S4M1, out_F_fmi_S4M2, out_F_fmi_S4M3, out_F_fmi_S4M4,  //mortalit� "captures" par p�che (par esp�ce)
        out_Fr_fmi_S1M1, out_Fr_fmi_S1M2, out_Fr_fmi_S1M3, out_Fr_fmi_S1M4, out_Fr_fmi_S2M1, out_Fr_fmi_S2M2, out_Fr_fmi_S2M3, out_Fr_fmi_S2M4,
        out_Fr_fmi_S3M1, out_Fr_fmi_S3M2, out_Fr_fmi_S3M3, out_Fr_fmi_S3M4, out_Fr_fmi_S4M1, out_Fr_fmi_S4M2, out_Fr_fmi_S4M3, out_Fr_fmi_S4M4,  //mortalit� totale (corrig�e de la survie) par p�che (par esp�ce)
        out_FRWT_fmi_S1M1, out_FRWT_fmi_S1M2, out_FRWT_fmi_S1M3, out_FRWT_fmi_S1M4, out_FRWT_fmi_S2M1, out_FRWT_fmi_S2M2, out_FRWT_fmi_S2M3, out_FRWT_fmi_S2M4,
        out_FRWT_fmi_S3M1, out_FRWT_fmi_S3M2, out_FRWT_fmi_S3M3, out_FRWT_fmi_S3M4, out_FRWT_fmi_S4M1, out_FRWT_fmi_S4M2, out_FRWT_fmi_S4M3, out_FRWT_fmi_S4M4,
        out_FDWT_fmi_S1M1, out_FDWT_fmi_S1M2, out_FDWT_fmi_S1M3, out_FDWT_fmi_S1M4, out_FDWT_fmi_S2M1, out_FDWT_fmi_S2M2, out_FDWT_fmi_S2M3, out_FDWT_fmi_S2M4,
        out_FDWT_fmi_S3M1, out_FDWT_fmi_S3M2, out_FDWT_fmi_S3M3, out_FDWT_fmi_S3M4, out_FDWT_fmi_S4M1, out_FDWT_fmi_S4M2, out_FDWT_fmi_S4M3, out_FDWT_fmi_S4M4,
        out_Z_eit_S1M1, out_Z_eit_S1M2, out_Z_eit_S1M3, out_Z_eit_S1M4, out_Z_eit_S2M1, out_Z_eit_S2M2, out_Z_eit_S2M3, out_Z_eit_S2M4,
        out_Z_eit_S3M1, out_Z_eit_S3M2, out_Z_eit_S3M3, out_Z_eit_S3M4, out_Z_eit_S4M1, out_Z_eit_S4M2, out_Z_eit_S4M3, out_Z_eit_S4M4,  //coefficient de mortalit� totale
        out_N_eit_S1M1, out_N_eit_S1M2, out_N_eit_S1M3, out_N_eit_S1M4, out_N_eit_S2M1, out_N_eit_S2M2, out_N_eit_S2M3, out_N_eit_S2M4,
        out_N_eit_S3M1, out_N_eit_S3M2, out_N_eit_S3M3, out_N_eit_S3M4, out_N_eit_S4M1, out_N_eit_S4M2, out_N_eit_S4M3, out_N_eit_S4M4,  //effectifs en nombre;

        out_Ytot_fm, out_DD_efmi, out_DD_efmc, out_LD_efmi, out_LD_efmc, out_statDD_efm, out_statLD_efm,
        out_DD_efmi_G1, out_DD_efmi_G2, out_LD_efmi_G1, out_LD_efmi_G2,
        out_statLDst_efm, out_statLDor_efm,
        intermBIOMspict, //effort1_fm et effort2_fm sont d�sormais inclus dans out_effort
        intermGoFish,
        multPrice;
//    VARIABLES  ---------------

//parties des inputs
SEXP    FList, sppList, sppListStat, sppListQ,sppListQM,sppListQM_dyn, pList, sppListAll, fleetList, metierList, metierListEco, namDC, t_init, times, Q, S, itListQ,
        NBVF, NBVFM, NBDSF, NBDSFM, EFF2F, EFF2FM, dnmsF, dnmsFM,dnmsIter, nmsEF, mu_nbds, mu_nbv, //mulitplicateurs d'effort
        m_f, m_fm, m_oth, eVar, eVarCopy, eStatVar, //variables interm�diaires par esp�ces
        fVar /*variable interm�diaire flottilles*/, list_copy, FList_copy, eVar_copy, fVar_copy, othSpSupList, effSupMat, listQR, listQR_f, TACbyF, TAC, reconcilSPP, reconcilSPP_copy, recList, recParamList, ParamSPMlist;

SEXP inpFtarg, inpW_Ftarg, inpMeanRec_Ftarg;
SEXP Qholdings;

int     nbT, nbF, nbM, nbMe, nbE, nbEstat, nbP,nbEall,nbEQuota,nbEQuotaMarket, nbEQuotaMarket_dyn,//dimensions
        curQ, spQ, scen, //application du sc�nario??
        bhv_active /*application du module report d'effort*/, type, /*boot, nbBoot,*/ ecodcf, typeGest, //special request ICES 2013 : pistage des r�gles de sc�nario int�gr� dans la variable out_typeGest
        var, trgt, delay, upd, gestInd, gestyp/*Module de gestion*/, activeQR,
        IND_T, IND_F, eTemp, fTemp /*indicateurs de temps, d'esp�ces et de flottilles consid�r�s*/, corVarTACnby_CPP, Blim_trigger, maxIter, t_stop,
        *SRInd, /**EcoIndCopy,*/ *Qvec, *recType1, *recType2, *recType3, *Svec; //indicateur conditionnant l'utilisation d'un recrutement al�atoire d�fini par la m�thode impl�ment�e RecAlea


double  PxQ, expEff, X1, X2, drCopy, tolVarTACinf_CPP, tolVarTACsup_CPP, corVarTACval_CPP, Blim_CPP, Bmax_CPP, //module de traitement stochastique de mod�le de prix
        *TAC_byFleet, *TAC_glob, *Fbar_trgt, diffZmax, lambda,//param�tres de contr�le du TAC pour l'analyse des mesures de gestion pour la requ�te CIEM 2013 sur la SOLE GG
        *effortIni, *effort1Ini, *Zoptim, *FOTHoptim, *Ztemp, *Etemp, *Einterm_fm, *EffsupTMP_fm, *Einterm_fm_copy, *multFOTHinterm_e;//Z fix� pour r�soudre l'ajustement par flottille (dimension �ge)


bool    Zoptim_use, FOTHoptim_use, boolQ, ZoptSS3, //indicateur de pr�sence de donn�es d'effort disponible (calcul capturabilit�,...)
        constMM, //indicateur qui d�termine si les niveaux m�tiers de la partie bio et de la partie �co sont les m�mes (utilisation de l'effort par flottille-m�tier pour calculer la capturabilit� dans le module 'Mortalit�')
        fUpdate, dUpdate, cUpdate, pflex, eUpdate; //indicateur de mise � jour des variables de calcul



int     *SPPstatOPT, *SPPspictOPT, *SPPdynOPT, N_SPPstatOPT, N_SPPspictOPT, N_SPPdynOPT;
SEXP    ZtempList;

//    METHODES  ---------------

	//constructeur
    BioEcoPar(SEXP list, SEXP listSpec, SEXP listStochastic, SEXP listScen,
                SEXP RecType1, SEXP RecType2, SEXP RecType3, SEXP Scenarii, /*SEXP Bootstrp, SEXP nbBootstrp , */ // TODO : remove unused arg
                SEXP GestInd, SEXP mOth, SEXP bounds, SEXP TAC, SEXP FBAR, SEXP othSpSup, SEXP effSup, SEXP GestParam, SEXP EcoDcf,
                SEXP persCalc, SEXP dr, SEXP SRind, SEXP listSR, SEXP TypeSR, SEXP mFM, SEXP TACbyFL, SEXP Ftarg, SEXP W_Ftarg, SEXP MeanRec_Ftarg,
                SEXP parBHV, SEXP parQEX,
                SEXP tacCTRL, SEXP stochPrice, SEXP updateE, SEXP parOQD, int VERBOSE, int force_T);

	//destructeur
    ~BioEcoPar();

    //accesseur d'un �l�ment de l'input
    SEXP getListElement(SEXP list, const char *str);

    //indice d'un �l�ment dans une liste ou un vecteur R
    int getListIndex(SEXP list, const char *str);

    //indice d'un �l�ment dans un vecteur non nomme R
    int getVectorIndex(SEXP vect, const char *str);

    //analyse des NAs dans un objet SEXP
    int all_is_na(SEXP object);

    double finite(double value);

    //indices multipliateurs pour les concordances de dimensions
    SEXP iDim(int *dimInput);

    //fontion d'agr�gation d'un objet R accompagn� de son attribut 'DimCst'
    SEXP aggregObj(SEXP object, SEXP newDim);

    // fonction de ventilation de la mortalit� en fonction d'une matrice de donn�es "capture"
    SEXP allocMortality(SEXP mortality, SEXP capture, SEXP captureTot);

    // fonction de calcul de l'indice de capturabilit� en fonction de la mortalit� par p�che et d'une variable d'effort quelconque
    SEXP calcCapturabilite(SEXP adjustedMortal, SEXP effortIni);

    void RecAlea(SEXP list, SEXP listSto, int ind_t, int type, int *recTyp);

    void SRmod(SEXP list, SEXP listSR, int ind_t, SEXP TypeSR, int *srind);

    void Scenario(SEXP list, SEXP listScen, int ind_t);

    // MODULES :
    //----------

    // Module 'Mortalit� par p�che et survie des rejets'
    void Mortalite(SEXP list, int ind_t, SEXP EVAR, int VERBOSE = 0);

    // Module 'Dynamique de population'
    void DynamicPop(SEXP list, int ind_t, SEXP EVAR, bool Reality, int VERBOSE = 0);

    // Module 'Captures, rejets et d�barquements'
    void CatchDL(SEXP list, int ind_t, SEXP EVAR, int VERBOSE = 0);

    // Module 'March�' : 'modCatch'
    void Marche(SEXP list, int ind_t);

    // Module 'Economie' DCF
    void EcoDCF(SEXP list, int ind_t,int persCalc,double dr);

    // Module Gestion
    double fxTAC_glob(double mult);
    void Gestion(SEXP list, int ind_t, int VERBOSE = 0);
    double zbrent(BEfn1 fx, double x1, double x2, double tol);
    void zbrak(BEfn1 fx, double x1, double x2, int n, double xb1[], double xb2[], int *nb);

    double func(double *x);
    //double fxTAC_F(double *x);

    int MinimizeF(double **p, double y[], int ndim, double ftol);
    void amoeba(BEfn1_F funk, double **p, double y[], int ndim, double ftol, int *nfunk);
    double amotry(BEfn1_F funk, double **p, double y[], double psum[], int ndim, int ihi, double fac);

    double *NRvector(long nl, long nh);
    double **NRmatrix(long nrl, long nrh, long ncl, long nch);
    void free_vector(double *v, long nl, long nh);
    void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);

    //int GestionF(double **p, double y[], int ndim, double ftol, int ind_t);

    void FleetBehav(SEXP list, int ind_t, SEXP paramBehav);

    void QuotaMarket(SEXP list, SEXP pQuotaIni, SEXP pQuotaMin, SEXP pQuotaMax, double lambdaQ, double sdmax, double ftol, int itmax, SEXP paramBehav, int ind_t, int persCalc );

    //int QuotaExch(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    //double fxMaxProf_FT(double *x);

    //double fxTAC_F_customCst(double *x);

    //double fxTAC_F_customReport(double *x);

    //double fxMaxProf_FT_customCst(double *x);

    //double fxMaxProf_FT_customReport(double *x);

    int EstimationTACfromF(int ind_t);

    // Module GestionF2
    void GestionF2(int ind_t);
    void abv_GestionF2(int ind_t);

    //double fxTAC_F_customCst2(double *x);

    //int GestionF2report(int spp, int ind_t);

    //double fxTAC_F_customReport2(double *x);

    //int QuotaExchV2(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    //double fxMaxProf_FT_customCstV2(double *x);

    //int QuotaExchV2Report(double pxQuIni, double pxQuMin, double pxQuMax, double lambda, int spp, double ftol, int ind_t);

    //double fxMaxProf_FT_customReportV2(double *x);

    //void PriceAlea(SEXP list, SEXP stPrice, int ind_t);


};

#endif // PARAM_H_INCLUDED

