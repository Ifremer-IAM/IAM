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


//extern "C" : pour �viter le "name mangling process" qui renomme les fonctions export�es dans les dll.


//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// UTILITIES
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//------------------------------------------
// accesseur � un �l�ment d'une liste donn�e (list = liste en question , str {caract�re} = intitul� de l'�l�ment de la liste)
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

extern "C" {

int BioEcoPar::getListIndex(SEXP list, const char *str) //fonctionne aussi pour les vecteurs nomm�s
{
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) break;

    return i;
}

}

extern "C" {

int BioEcoPar::getVectorIndex(SEXP vect, const char *str) //fonctionne aussi pour les vecteurs nomm�s
{
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(vect,i)), str) == 0) break;

    return i;
}

}


//------------------------------------------
// fonction all.is.na (teste si tous les �l�ments d'un objet sont � NA ou non)
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
// (permet la g�n�ricit� des �quations en assurant la compatibilit� des variables en pr�sence,
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
// fonction d'agr�gation d'un objet attribu� type ('object'), en fonction d'un nouveau vecteur dimension DimCst ('newDim')
// NB : toutes les valeurs de 'newDim' doivent �tre au plus �gales aux dimensions correspondantes de l'objet pour que la fonction s'applique
//------------------------------------------

extern "C" {

SEXP BioEcoPar::aggregObj(SEXP object, SEXP newDim)
{
    PROTECT(object=object);
    PROTECT(newDim=newDim);

    SEXP ans, dimObj, dimnames, Dim;

    int *dim, *ndim, *rdim;
    double *rans = &NA_REAL, *robj = &NA_REAL;

    PROTECT(dimObj = getAttrib(object, install("DimCst")));

    dim = INTEGER(dimObj); ndim = INTEGER(newDim); //Rprintf("in aggegObj:") ;PrintValue(dimObj); PrintValue(newDim);

    //tests sur les dimensions
    if ((dim[0]==0) & (dim[1]==0) & (dim[2]==0) & (dim[3]==0)) {  //c'est termin�, rien � agr�ger
        UNPROTECT(3);
        return(object);

    } else {
        //Rprintf("Dans aggregobj dim[0] %i ndim[0] %i test %d ;  dim[1] %i ndim[1] %i test %d \n", dim[0],ndim[0],dim[0]<ndim[0],dim[1],ndim[1],dim[1]<ndim[1])
        if ((dim[0]<ndim[0]) | (dim[1]<ndim[1]) | (dim[2]<ndim[2]) | (dim[3]<ndim[3]))
        {
            error("Check input dimensions in 'aggregObj'!!\n");
        }

        //on calcule le nombre de cellules � remplir et le nombre de dimensions nulles
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

            //en-t�tes
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

        //il ne reste plus qu'� effectuer l'agr�gation
        for (int ind_f = 0 ; ind_f < (1 + (dim[0] - 1)*(dim[0]>0)) ; ind_f++)
        for (int ind_m = 0 ; ind_m < (1 + (dim[1] - 1)*(dim[1]>0)) ; ind_m++)
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
// fonction de calcul de l'indice de capturabilit� en fonction de la mortalit� par p�che et d'une variable d'effort donn�e : � op�rer � t=0
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
    if (((dimE[0]!=0) & (dimM[0]!=0) & (dimE[0]!=dimM[0])) | ((dimE[1]!=0) & (dimM[1]!=0) & (dimE[1]!=dimM[1])) |
        ((dimE[2]!=0) & (dimM[2]!=0) & (dimE[2]!=dimM[2])) | ((dimE[3]!=0) & (dimM[3]!=0) & (dimE[3]!=dimM[3])))
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
    dimEffort[3] = dimE[3]; //on n'agr�ge pas sur le temps puisque on ne consid�re ensuite que l'instant initial

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
                finite(rMort[ind_f*fact1[0] + ind_m*fact1[1] + ind_i*fact1[2]] /
                rEff[ind_f*fact2[0] + ind_m*fact2[1] + ind_i*fact2[2]]) ;

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




//
////------------------------------------------
//// Module 'March�'
////------------------------------------------
//
extern "C" {

void BioEcoPar::Marche(SEXP list, int ind_t)
{

//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.Marche.txt", ios::out | ios::trunc);

    SEXP    elmt, intC, v_P_fmce, v_icat, v_L_efmit, dimCst_P_fmce, dimCst_L_efmit, dimCst_L_efmct, Dim_L_efmct, //Dim_P_fmce,
            ans_L_efmct = R_NilValue, ans_P_fmce = R_NilValue, dimnames_Lc = R_NilValue, dimnames_P = R_NilValue,
            rnames_Esp, cFACTc, cFACTi, cFACTp, cFACTpini, cFACTpStat, cFACTpStatini,
            ans_DD_efmc = R_NilValue, ans_LD_efmc = R_NilValue, v_DD_efmit, v_LD_efmit;
    SEXP    v_P_eStat, dimCst_P_eStat, dimCst_P_eStatR, ans_P_eStat = R_NilValue, dimnames_Pstat = R_NilValue;

    int *dim_P_fmce, *dim_P_fmcet, *dim_L_efmit, *dim_icat, *dim_L_efmct, *dim_P_eStat, *dim_P_eStat_t, *dim_P_eStat_tR, *dimLc;//, *dimP;

    int nbI, nbC;

    double *rans_L_efmct, *r_L_efmit, *r_P_fmce, *r_icat, *r_Pstat, *r_P_fmceIni, *r_PstatIni,
            *rans_DD_efmc, *rans_LD_efmc, *r_DD_efmit, *r_LD_efmit;
    double mult_p; //multiplicateur de prix (1 par defaut sinon fonction de relation prix quantite
    SEXP MarketList= R_NilValue, v_ep= R_NilValue, v_beta_pp= R_NilValue;
    double *r_beta_pp = &NA_REAL;
    int *r_ep = &NA_INTEGER;
    int ind_p;
    SEXP ans_MultPrice;


////Rprintf("CCC1");

//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 5));////Rprintf("Mort20.2\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 6));////Rprintf("Mort20.3\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 7));////Rprintf("Mort20.4\n");
//PrintValue(VECTOR_ELT(VECTOR_ELT(eVar, 0), 61));////Rprintf("Mort20.5\n");



    if (ind_t==0) {

        PROTECT(rnames_Esp = getAttrib(out_L_efmit,R_NamesSymbol));
        setAttrib(out_L_efmct, R_NamesSymbol, rnames_Esp);
        setAttrib(out_DD_efmc, R_NamesSymbol, rnames_Esp);
        setAttrib(out_LD_efmc, R_NamesSymbol, rnames_Esp);
        setAttrib(out_P_t, R_NamesSymbol, rnames_Esp);
        if (nbEstat>0) setAttrib(out_Pstat, R_NamesSymbol, sppListStat);

    }



if(pUpdate){


    PROTECT(MarketList = getListElement(list,"Market"));
    PROTECT(v_ep = getListElement(MarketList,"ep")); //PrintValue(v_ep);
    PROTECT(v_beta_pp = getListElement(MarketList,"beta_pp")); //PrintValue(v_beta_pp);

    r_beta_pp = REAL(v_beta_pp);
    r_ep = INTEGER(v_ep);

    //---------
    // Calcul debarquements produits
    //---------

    SEXP rnames_p = R_NilValue, ans_L_pt = R_NilValue, dimCstL_pt, dimNamL_pt;
    int *int_dimCstL_pt;
    double L_p;
    double *rans_L_pt;

    if (ind_t==0) {
        PROTECT(rnames_p = allocVector(STRSXP, nbP));
        setAttrib(out_L_pt, R_NamesSymbol, rnames_p);
    }

    for (int p = 0; p < nbP; p++) {

            if (ind_t==0) {
            PROTECT(ans_L_pt = NEW_NUMERIC(nbT));

            PROTECT(dimCstL_pt = allocVector(INTSXP, 1));
            int_dimCstL_pt = INTEGER(dimCstL_pt);  int_dimCstL_pt[0] = nbT;
            setAttrib(ans_L_pt, R_DimSymbol, dimCstL_pt);

            PROTECT( dimNamL_pt = allocVector(VECSXP,1));
            SET_VECTOR_ELT(dimNamL_pt, 0, times);
            setAttrib(ans_L_pt, R_DimNamesSymbol, dimNamL_pt);

            rans_L_pt = REAL(ans_L_pt);
            } else {
                rans_L_pt = REAL(VECTOR_ELT(out_L_pt, p));
                    }

            //Colonne de la matrice ep : 1 si espece correspond au meme produit, 0 sinon
            L_p = 0.0;
            for (int k = 0; k < nbEall; k++) {
                    L_p = L_p + r_ep[p*nbEall + k] * REAL(VECTOR_ELT(out_L_et, k))[ind_t];
                    //fichier << "p = " << CHAR(STRING_ELT(pList,p)) << ", k = " << CHAR(STRING_ELT(sppListAll,k)) << ", r_ep = " << r_ep[p*nbEall + k] << ", L_p = " << L_p << endl;
            }

            rans_L_pt [ind_t] = L_p; //fichier << "p = " << CHAR(STRING_ELT(pList,p)) << ", rans_L_pt = " << rans_L_pt [ind_t] << endl;

            if (ind_t==0) {
                SET_VECTOR_ELT(out_L_pt, p, ans_L_pt);
                SET_STRING_ELT(rnames_p, p, STRING_ELT(pList,p));
                UNPROTECT(3);
                }

        }
    }

if (nbE>0) {

    for (int e = 0 ; e < nbE ; e++) {
//Rprintf("M2\n");fichier << "M2" << endl;
//fichier << "e = " << CHAR(STRING_ELT(sppList,e)) << endl;

        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));

        nbI = length(getListElement(elmt, "modI"));
        intC = getListElement(elmt, "modC");
        nbC = length(intC);

        PROTECT(v_P_fmce = getListElement(elmt, "P_fmce"));
        PROTECT(v_icat = getListElement(elmt, "icat"));
        PROTECT(v_L_efmit = getListElement(out_L_efmit, CHAR(STRING_ELT(sppList,e))));
        PROTECT(v_DD_efmit = getListElement(out_DD_efmi, CHAR(STRING_ELT(sppList,e))));
        PROTECT(v_LD_efmit = getListElement(out_LD_efmi, CHAR(STRING_ELT(sppList,e))));

        PROTECT(dimCst_P_fmce = getAttrib(v_P_fmce, install("DimCst")));
        PROTECT(dimCst_L_efmit = getAttrib(v_L_efmit, install("DimCst")));
//Rprintf("M3\n");fichier << "M3" << endl;

        //tests sur les dimensions :
        dim_P_fmce = INTEGER(dimCst_P_fmce);////Rprintf("AAA1");
        if (((dim_P_fmce[0]!=0) & (dim_P_fmce[0]!=nbF)) | ((dim_P_fmce[1]!=0) & (dim_P_fmce[1]!=nbMe)) |
            ((dim_P_fmce[2]!=0) & (dim_P_fmce[2]!=nbC)) | ((dim_P_fmce[3]!=0) & (dim_P_fmce[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in P_fmce element. Check .ini biological parameters files !!\n");
        }

        dim_L_efmit = INTEGER(dimCst_L_efmit);////Rprintf("AAA2");
        if (((dim_L_efmit[0]!=0) & (dim_L_efmit[0]!=nbF)) | ((dim_L_efmit[1]!=0) & (dim_L_efmit[1]!=nbM)) |
            ((dim_L_efmit[2]!=0) & (dim_L_efmit[2]!=nbI)) | ((dim_L_efmit[3]!=0) & (dim_L_efmit[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in L_efmit element. Check .ini biological parameters files !!\n");
        }

        dim_icat = INTEGER(getAttrib(v_icat, R_DimSymbol));////Rprintf("AAA3");
        if ((dim_icat[0]!=nbI) & (dim_icat[1]!=nbC))
        {
            error("Non_homogeneous dimensions in icat element. Check .ini biological parameters files !!\n");
        }
//Rprintf("M4\n");fichier << "M4" << endl;

//        ---------
//         Ajustement prix avec relation prix/quantite
//        ---------
        mult_p = 1.0;
        if(pUpdate){
//Rprintf("M5\n");fichier << "M5" << endl;
                int ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e))); //indice dans matrice ep
                //fichier << "e = " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

                double *r_L_pt;

                //Trouver l'indice du produit correspondant
                for (ind_p = 0 ; ind_p < nbP ; ind_p++) {if (r_ep[ind_p*nbEall + ind_e] == 1) break;}

                //fichier << "p = " << CHAR(STRING_ELT(pList,ind_p)) << endl;

                for ( int k = 0; k < nbP; k++) {

                    r_L_pt = REAL(VECTOR_ELT(out_L_pt, k));
                    if (r_L_pt[0]!=0) mult_p = mult_p * pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) ;
                    //fichier << "k = " << CHAR(STRING_ELT(pList,k)) << ", ratio L = " << r_L_pt[ind_t] / r_L_pt[0] << ", beta = " << r_beta_pp[k*nbP + ind_p] << ", pow = "<<  pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) << ", mult = " << mult_p << endl;

                }

                PROTECT(ans_MultPrice = VECTOR_ELT(multPrice, ind_e));
                REAL(ans_MultPrice)[ind_t] = mult_p;
                UNPROTECT(1);

        }




        //---------
        // calcul de L_efmct
        //---------

        PROTECT(dimCst_L_efmct = allocVector(INTSXP, 4));
        dim_L_efmct = INTEGER(dimCst_L_efmct);
        dim_L_efmct[0] = dim_L_efmit[0] ; dim_L_efmct[1] = dim_L_efmit[1] ; dim_L_efmct[2] = nbC; dim_L_efmct[3] = dim_L_efmit[3];

        int count = 0, prod = 1, count2 = 0, count3 = 0;
        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                count++;
                prod = prod * dim_L_efmct[k];
            }

        }

        PROTECT(Dim_L_efmct = allocVector(INTSXP, count));
        dimLc = INTEGER(Dim_L_efmct);

//Rprintf("M5\n");fichier << "M5" << endl;
        for (int k = 0 ; k < 4 ; k++) {

            if (dim_L_efmct[k]>0) {
                dimLc[count2] = dim_L_efmct[k];
                count2++;
            }

        }

//--------------------------

        PROTECT(dimCst_P_fmce = allocVector(INTSXP, 4));
        dim_P_fmcet = INTEGER(dimCst_P_fmce);
        dim_P_fmcet[0] = nbF ; dim_P_fmcet[1] = nbMe ; dim_P_fmcet[2] = nbC; dim_P_fmcet[3] = nbT;



if (ind_t==0){
//Rprintf("M6\n");fichier << "M6" << endl;
        //on cr�e le tableau r�sultat pour l'esp�ce en question
        PROTECT(ans_L_efmct = NEW_NUMERIC(prod));
        setAttrib(ans_L_efmct, R_DimSymbol, Dim_L_efmct);

        PROTECT(ans_DD_efmc = NEW_NUMERIC(prod));
        setAttrib(ans_DD_efmc, R_DimSymbol, Dim_L_efmct);

        PROTECT(ans_LD_efmc = NEW_NUMERIC(prod));
        setAttrib(ans_LD_efmc, R_DimSymbol, Dim_L_efmct);

        PROTECT(dimnames_Lc = allocVector(VECSXP,count));
        if (dim_L_efmct[0]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, fleetList) ; count3++;}
        if (dim_L_efmct[1]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, metierList) ; count3++;}
        if (dim_L_efmct[2]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, intC) ; count3++;}
        if (dim_L_efmct[3]>0) {SET_VECTOR_ELT(dimnames_Lc, count3, times) ; count3++;}

        rans_L_efmct = REAL(ans_L_efmct);
        rans_DD_efmc = REAL(ans_DD_efmc);
        rans_LD_efmc = REAL(ans_LD_efmc);

        PROTECT(ans_P_fmce = NEW_NUMERIC(nbF*nbMe*nbC*nbT));
        setAttrib(ans_P_fmce, R_DimSymbol, dimCst_P_fmce);

        PROTECT(dimnames_P = allocVector(VECSXP,4));
        SET_VECTOR_ELT(dimnames_P, 0, fleetList);
        SET_VECTOR_ELT(dimnames_P, 1, metierListEco);
        SET_VECTOR_ELT(dimnames_P, 2, intC);
        SET_VECTOR_ELT(dimnames_P, 3, times);

        r_P_fmce = REAL(ans_P_fmce);

//Rprintf("M7\n");fichier << "M7" << endl;
} else {

        rans_L_efmct = REAL(VECTOR_ELT(out_L_efmct, e));
        rans_DD_efmc = REAL(VECTOR_ELT(out_DD_efmc, e));
        rans_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
        r_P_fmce = REAL(VECTOR_ELT(out_P_t, e));

//Rprintf("M8\n");fichier << "M8" << endl;
}

        r_L_efmit = REAL(v_L_efmit);
        r_DD_efmit = REAL(v_DD_efmit);
        r_LD_efmit = REAL(v_LD_efmit);
        r_icat = REAL(v_icat);

//Rprintf("M9\n");fichier << "M9" << endl;
        //facteurs des indices
        PROTECT(cFACTc = iDim(dim_L_efmct));
        PROTECT(cFACTi = iDim(dim_L_efmit));
        PROTECT(cFACTp = iDim(dim_P_fmcet));
        PROTECT(cFACTpini = iDim(dim_P_fmce));

        int *fact_Cc = INTEGER(cFACTc);
        int *fact_Ci = INTEGER(cFACTi);
        int *fact_P = INTEGER(cFACTp);
        int *fact_Pini = INTEGER(cFACTpini);

        r_P_fmceIni = REAL(v_P_fmce);

        //�quation n�1 : conversion �ge/catg�gorie

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_m = 0 ; ind_m < nbM ; ind_m++)
        for (int ind_c = 0 ; ind_c < nbC ; ind_c++) {

            for (int ind_i = 0 ; ind_i < nbI ; ind_i++) {

                if (ind_i ==0) {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

           rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_LD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

           rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                r_DD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                } else {

            rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_L_efmct[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_L_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

            rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_LD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_LD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

            rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] =
                rans_DD_efmc[ind_f*fact_Cc[0] + ind_m*fact_Cc[1] + ind_c*fact_Cc[2] + ind_t*fact_Cc[3]] +
                r_DD_efmit[ind_f*fact_Ci[0] + ind_m*fact_Ci[1] + ind_i*fact_Ci[2] + ind_t*fact_Ci[3]] * r_icat[ind_i + nbI*ind_c];

                }
            }
        // valable car nbM = nbMe (M=Me)
            //r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] =
            //    r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + ind_t*fact_Pini[3]];

            r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] =
                r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + 0*fact_Pini[3]] * mult_p; //Florence 05/2019

                //if ((ind_f==0) & (ind_m==0)) fichier << "r_P_fmceIni = " << r_P_fmceIni[ind_f*fact_Pini[0] + ind_m*fact_Pini[1] + ind_c*fact_Pini[2] + 0*fact_Pini[3]] << ", mult_p = " << mult_p << ", r_P_fmce = " << r_P_fmce[ind_f*fact_P[0] + ind_m*fact_P[1] + ind_c*fact_P[2] + ind_t*fact_P[3]] << endl;

        }


//Rprintf("M10\n");fichier << "M10" << endl;

if (ind_t==0) {

        setAttrib(ans_L_efmct, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_L_efmct, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_L_efmct, e, ans_L_efmct);

        setAttrib(ans_DD_efmc, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_DD_efmc, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_DD_efmc, e, ans_DD_efmc);

        setAttrib(ans_LD_efmc, R_DimNamesSymbol, dimnames_Lc);
        setAttrib(ans_LD_efmc, install("DimCst"), dimCst_L_efmct);

        SET_VECTOR_ELT(out_LD_efmc, e, ans_LD_efmc);

      //  SET_VECTOR_ELT(out_L_efmct2, e, ans_L_efmct);

        setAttrib(ans_P_fmce, R_DimNamesSymbol, dimnames_P);
        setAttrib(ans_P_fmce, install("DimCst"), dimCst_P_fmce);

        SET_VECTOR_ELT(out_P_t, e, ans_P_fmce);

}

//Rprintf("M11\n");fichier << "M11" << endl;
if (ind_t==0) UNPROTECT(6);
UNPROTECT(15);

}
}


if (nbEstat>0) {
    for (int e = 0 ; e < nbEstat ; e++) {
//Rprintf("M2\n");fichier << "M2" << endl;
        PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e))));

//        ---------
//         Ajustement prix avec relation prix/quantite
//        ---------
        mult_p = 1.0;
        if(pUpdate){
//Rprintf("M5\n");fichier << "M5" << endl;
                int ind_e = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e))); //indice dans matrice ep
                //fichier << "e = " << CHAR(STRING_ELT(sppListAll,ind_e)) << endl;

                double *r_L_pt;

                //Trouver l'indice du produit correspondant
                for (ind_p = 0 ; ind_p < nbP ; ind_p++) {if (r_ep[ind_p*nbEall + ind_e] == 1) break;}

                //fichier << "p = " << CHAR(STRING_ELT(pList,ind_p)) << endl;

                for ( int k = 0; k < nbP; k++) {

                    r_L_pt = REAL(VECTOR_ELT(out_L_pt, k));
                    if (r_L_pt[0]!=0) mult_p = mult_p * pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) ;
                    //fichier << "k = " << CHAR(STRING_ELT(pList,k)) << ", ratio L = " << r_L_pt[ind_t] / r_L_pt[0] << ", beta = " << r_beta_pp[k*nbP + ind_p] << ", pow = "<<  pow(r_L_pt[ind_t] / r_L_pt[0], r_beta_pp[k*nbP + ind_p]) << ", mult = " << mult_p << endl;

                }
            PROTECT(ans_MultPrice = VECTOR_ELT(multPrice, ind_e));
            REAL(ans_MultPrice)[ind_t] = mult_p;
            UNPROTECT(1);


        }


        PROTECT(v_P_eStat = getListElement(elmt, "P_fme"));
        PROTECT(dimCst_P_eStat = getAttrib(v_P_eStat, install("DimCst")));

        //tests sur les dimensions :
        dim_P_eStat = INTEGER(dimCst_P_eStat);////Rprintf("AAA1");
        if (((dim_P_eStat[0]!=0) & (dim_P_eStat[0]!=nbF)) | ((dim_P_eStat[1]!=0) & (dim_P_eStat[1]!=nbMe)) |
            (dim_P_eStat[2]!=0) | ((dim_P_eStat[3]!=0) & (dim_P_eStat[3]!=nbT)))
        {
            error("Non_homogeneous dimensions in P_fme element. Check .ini biological parameters files !!\n");
        }



        PROTECT(dimCst_P_eStat = allocVector(INTSXP, 4));
        dim_P_eStat_t = INTEGER(dimCst_P_eStat);
        dim_P_eStat_t[0] = nbF ; dim_P_eStat_t[1] = nbMe ; dim_P_eStat_t[2] = 0; dim_P_eStat_t[3] = nbT;
        PROTECT(dimCst_P_eStatR = allocVector(INTSXP, 3));
        dim_P_eStat_tR = INTEGER(dimCst_P_eStatR);
        dim_P_eStat_tR[0] = nbF ; dim_P_eStat_tR[1] = nbMe ; dim_P_eStat_tR[2] = nbT;


if (ind_t==0){
//Rprintf("M6\n");fichier << "M6" << endl;

        PROTECT(ans_P_eStat = NEW_NUMERIC(nbF*nbMe*nbT));
        setAttrib(ans_P_eStat, R_DimSymbol, dimCst_P_eStatR);

        PROTECT(dimnames_Pstat = allocVector(VECSXP,3));
        SET_VECTOR_ELT(dimnames_Pstat, 0, fleetList);
        SET_VECTOR_ELT(dimnames_Pstat, 1, metierListEco);
        SET_VECTOR_ELT(dimnames_Pstat, 2, times);

        r_Pstat = REAL(ans_P_eStat);

//Rprintf("M7\n");fichier << "M7" << endl;
} else {

        r_Pstat = REAL(VECTOR_ELT(out_Pstat, e));

//Rprintf("M8\n");fichier << "M8" << endl;
}
//Rprintf("M9\n");fichier << "M9" << endl;
        //facteurs des indices
        PROTECT(cFACTpStat = iDim(dim_P_eStat_t));
        PROTECT(cFACTpStatini = iDim(dim_P_eStat));

        int *fact_Pstat = INTEGER(cFACTpStat);
        int *fact_Pstatini = INTEGER(cFACTpStatini);

        r_PstatIni = REAL(v_P_eStat);

        //�quation n�1 : conversion �ge/catg�gorie

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_m = 0 ; ind_m < nbM ; ind_m++) {

        // valable car nbM = nbMe (M=Me)
           // r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] =
             //   r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + ind_t*fact_Pstatini[3]];

            r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] =
                r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + 0*fact_Pstatini[3]] * mult_p;// Florence 05/2019


                //if ((ind_f==0) & (ind_m==0)) fichier << "r_PstatIni = " << r_PstatIni[ind_f*fact_Pstatini[0] + ind_m*fact_Pstatini[1] + 0*fact_Pstatini[2] + 0*fact_Pstatini[3]] << ", mult_p = " << mult_p << ", r_Pstat = " << r_Pstat[ind_f*fact_Pstat[0] + ind_m*fact_Pstat[1] + 0*fact_Pstat[2] + ind_t*fact_Pstat[3]] << endl;


        }

//Rprintf("M10\n");fichier << "M10" << endl;

if (ind_t==0) {

        setAttrib(ans_P_eStat, R_DimNamesSymbol, dimnames_Pstat);
        setAttrib(ans_P_eStat, install("DimCst"), dimCst_P_eStat);

        SET_VECTOR_ELT(out_Pstat, e, ans_P_eStat);

}

//Rprintf("M11\n");fichier << "M11" << endl;
if (ind_t==0) UNPROTECT(2);
UNPROTECT(7);

}

}

if (pUpdate){
        if(ind_t == 0) UNPROTECT(1);
    UNPROTECT(3);
}
if (ind_t==0) UNPROTECT(1);


//fichier.close();

}}




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





//---------------------------------
//
// Module de mod�lisation de relations S/R
//
//---------------------------------



extern "C" {

void BioEcoPar::SRmod(SEXP list, SEXP listSR, int ind_t, SEXP TypeSR, int *srind)
        //list : liste des param�tres d'entr�e;
        //listSR : liste des param�tres a,b&c du mod�le SR + e.t bruit normal ou lognormal + type de bruit (1=normal, 2=lognormal) (un vecteur de longueur 5 par esp�ce mod�lis�e contenant des "doubles")
        //type : type de relation Stock-Recrutement : (liste de longueur "nb d'esp�ces mod�lis�es" contenant des entiers)
        //                                              1 -> recrutement constant moyen (rec~a)
        //                                              2 -> Hockey stick (rec ~ (si (ssb<=b) a*ssb sinon a*b))
        //                                              3 -> Beverton & Holt (rec ~ a*ssb/(b+ssb))
        //                                              4 -> Ricker (rec ~ a*ssb*exp(-b*ssb))
        //                                              5 -> Shepherd (rec ~ a*ssb/(1+ (ssb/b)^c))
        //                                              6 -> Hockey Stick Quadratic (rec ~ (si (ssb<=b*(1-c)) a*ssb ; si (b*(1-c)<ssb<b*(1+c)) a*(ssb-((ssb-b*(1-c))^2)/(4*b*c)) ; sinon a*b))
        //                                              7 -> Hockey Stick Smooth (rec ~ a*(ssb+sqrt(b^2+g)-sqrt((ssb-b)^2+g)), avec g=0.001 )
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

    if (ind_t==0) { //deuxi�me �tape d'initialisation (niveau esp�ce)

        PROTECT(ans = NEW_NUMERIC(nbT));
        setAttrib(ans, R_NamesSymbol, times);
        SET_VECTOR_ELT(out_SRmod, e, ans);
        SET_STRING_ELT(rnames, e, STRING_ELT(sppList,e));
    }

    rans = REAL(VECTOR_ELT(out_SRmod, e));
    paramet = REAL(VECTOR_ELT(listSR, e));
    typeSR = INTEGER(VECTOR_ELT(TypeSR, e))[0];
    if (ind_t>0) ssb = REAL(VECTOR_ELT(out_SSB_et, e)); else ssb = &NA_REAL;

    //il nous faut aussi le d�calage temporel d� au premier �ge mod�lis� -> un SSB(t) g�n�rera un R(t+age0) #correction de R(t+age0+1)
    fstAge = CHAR(STRING_ELT(VECTOR_ELT(namDC, e),0))[0] - '0'; //++fstAge;

    //on en profite pour initialiser l'objet pour les premi�res ann�es pour lesquelles on devra aller chercher l'info dans Ni0
    if (((typeSR!=1) & (ind_t<fstAge)) | (ind_t==0)) {  //deuxi�me condition : si t initial et recrutement d�duit de la ssb de la m�me ann�e, on part des param�tres initiaux et non de la relation SR

        rans[ind_t] = NA_REAL;

    } else {

    switch (typeSR) {

        case 1 :

        rans[ind_t] = paramet[0*nbT + ind_t]; break;  //modif MM 27/08/2013 pour donner la possibilit� de d�finir plusieurs param�trages pour la relation SR au cours de la simu
                                                        //indices : 0 --> 0*nbT + ind_t
        case 2 :

        if (ssb[ind_t-fstAge] <= paramet[1*nbT + ind_t])
            rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge];
        else
            rans[ind_t] = paramet[0*nbT + ind_t] * paramet[1*nbT + ind_t]; break;

        case 3 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] / (paramet[1*nbT + ind_t] + ssb[ind_t-fstAge]); break;

        case 4 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] * exp(-1.0 * paramet[1*nbT + ind_t] * ssb[ind_t-fstAge]); break;

        case 5 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] / (1 + pow(ssb[ind_t-fstAge] / paramet[1*nbT + ind_t] , paramet[2*nbT + ind_t])); break;

        case 6 :

        if (ssb[ind_t-fstAge] <= (paramet[1*nbT + ind_t]*(1-paramet[2*nbT + ind_t])))

            rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge];

        else {

            if (ssb[ind_t-fstAge] >= (paramet[1*nbT + ind_t]*(1+paramet[2*nbT + ind_t])))

                rans[ind_t] = paramet[0*nbT + ind_t] * paramet[1*nbT + ind_t];

            else

                rans[ind_t] = paramet[0*nbT + ind_t]*(ssb[ind_t-fstAge] - (pow(ssb[ind_t-fstAge] - paramet[1*nbT + ind_t]*(1-paramet[2*nbT + ind_t]),2.0)/(4*paramet[1*nbT + ind_t]*paramet[2*nbT + ind_t])));

            } break;

        case 7 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ( ssb[ind_t-fstAge] + sqrt( pow(paramet[1*nbT + ind_t],2) + 0.001) - sqrt( pow(ssb[ind_t-fstAge] - paramet[1*nbT + ind_t],2) + 0.001)); break;

        default :

        rans[ind_t] = NA_REAL;

    }
    //il ne reste plus qu'� ajouter le bruit blanc issue de N(0,sigma) avec sigma = paramet[3]
    double v_alea = 0.0;
GetRNGstate();
        if (!ISNA(paramet[3*nbT + ind_t])) v_alea = rnorm(0.0,paramet[3*nbT + ind_t]); //////Rprintf("%f ",v_alea);//Rprintf("%f ",rnorm(0.0,0.157));
PutRNGstate();
    if (paramet[4*nbT + ind_t]==1)  //bruit de loi normale

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
// Module de gestion des sc�narios
//
//---------------------------------



extern "C" {

void BioEcoPar::Scenario(SEXP list, SEXP listScen, int ind_t) //list : liste des param�tres d'entr�e ; listScen : liste des multiplicateurs pour un sc�nario donn�
{

//1er niveau de la liste de multiplicateurs : Fleet ou Species --> on cible la partie de "list" correspondante

SEXP mult_lvl_1, target_lvl_1, mult_lvl_2, target_lvl_2=R_NilValue, namVar, namElt, dimMult, dimTarget, fTarg, fMult;

int nbElt = length(listScen);

for (int elt = 0 ; elt < nbElt ; elt++) {

    PROTECT(namElt = STRING_ELT(getAttrib(listScen, R_NamesSymbol), elt));

    PROTECT(mult_lvl_1 = getListElement(listScen, CHAR(namElt))); //Rprintf("%i \n",elt); PrintValue(namElt);

    if (mult_lvl_1 != NULL) {

        PROTECT(target_lvl_1 = getListElement(list, CHAR(namElt)));
        int nbVar = length(mult_lvl_1);

        for (int i = 0 ; i < nbVar ; i++) {

            PROTECT(namVar = STRING_ELT(getAttrib(mult_lvl_1, R_NamesSymbol), i)); //Rprintf("%i \n",i); PrintValue(namVar);
            PROTECT(mult_lvl_2 = getListElement(mult_lvl_1, CHAR(namVar)));

            //ici, selon que la variable consid�r�e est un input ou une variable interne (ex : Foth_i), on agit diff�remment
int indic = 0;

if (strcmp(CHAR(namVar), "Ffmi_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 100)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 101)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 102)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 103)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 104)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 105)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 106)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 107)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 108)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 109)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 110)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 111)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 112)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 113)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 114)); indic=1;}
if (strcmp(CHAR(namVar), "Ffmi_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 115)); indic=1;}

if (strcmp(CHAR(namVar), "Foth_i_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 116)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 117)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 118)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 119)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 120)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 121)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 122)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 123)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 124)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 125)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 126)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 127)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 128)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 129)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 130)); indic=1;}
if (strcmp(CHAR(namVar), "Foth_i_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 131)); indic=1;}

if (strcmp(CHAR(namVar), "FLWfmi_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 160)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 161)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 162)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 163)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 164)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 165)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 166)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 167)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 168)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 169)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 170)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 171)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 172)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 173)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 174)); indic=1;}
if (strcmp(CHAR(namVar), "FLWfmi_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 175)); indic=1;}

if (strcmp(CHAR(namVar), "FRWToth_i_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 176)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 177)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 178)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 179)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 180)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 181)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 182)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 183)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 184)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 185)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 186)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 187)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 188)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 189)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 190)); indic=1;}
if (strcmp(CHAR(namVar), "FRWToth_i_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 191)); indic=1;}

if (strcmp(CHAR(namVar), "FDWfmi_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 192)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 193)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 194)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 195)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 196)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 197)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 198)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 199)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 200)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 201)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 202)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 203)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 204)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 205)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 206)); indic=1;}
if (strcmp(CHAR(namVar), "FDWfmi_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 207)); indic=1;}

if (strcmp(CHAR(namVar), "FDWToth_i_S1M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 208)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S1M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 209)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S1M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 210)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S1M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 211)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S2M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 212)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S2M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 213)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S2M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 214)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S2M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 215)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S3M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 216)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S3M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 217)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S3M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 218)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S3M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 219)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S4M1") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 220)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S4M2") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 221)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S4M3") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 222)); indic=1;}
if (strcmp(CHAR(namVar), "FDWToth_i_S4M4") == 0) {PROTECT(target_lvl_2 = VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 223)); indic=1;}

            if (strcmp(CHAR(namVar), "Foth_i") == 0) {

                ////Rprintf("%i",IS_NUMERIC(VECTOR_ELT(getListElement(eVar, CHAR(namElt)), 44)));
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

                                            if (indic==0) PROTECT(target_lvl_2 = getListElement(target_lvl_1, CHAR(namVar)));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

                int *dimM, *dimT;
                int typeSc = 0;

                PROTECT(dimMult = getAttrib(mult_lvl_2, install("DimCst"))); //Rprintf("gg\n");
                if (getAttrib(mult_lvl_2, install("type")) != NULL) typeSc = INTEGER(getAttrib(mult_lvl_2, install("type")))[0];
                //si 'target_lvl_2' est un �l�ment de eVar, s'assurer au pr�alable de l'existence de l'attribut DimCst
                PROTECT(dimTarget = getAttrib(target_lvl_2, install("DimCst")));

                dimM = INTEGER(dimMult); dimT = INTEGER(dimTarget);//Rprintf("hh");

            //tests sur les dimensions
                if ((dimM[0]>dimT[0]) | (dimM[1]>dimT[1]) | (dimM[2]>dimT[2])) error("Wrong dimensions specification in 'Scenario' !!\n");

                PROTECT(fTarg = iDim(dimT));
                PROTECT(fMult = iDim(dimM));

                int *ftarg = INTEGER(fTarg); //Rprintf("ii");
                int *fmult = INTEGER(fMult); //Rprintf("jj");

                double *target = REAL(target_lvl_2), *mult = REAL(mult_lvl_2);

            //et on applique la mise � jour selon typeSc

            if ((typeSc==0) | (typeSc==1)) {

                for (int ind_f = 0 ; ind_f < imax2(1,dimT[0]) ; ind_f++)
                for (int ind_m = 0 ; ind_m < imax2(1,dimT[1]) ; ind_m++)
                for (int ind_i = 0 ; ind_i < imax2(1,dimT[2]) ; ind_i++) {
                    if (!ISNA(mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]])) {

                        target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] =
                        target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] *
                        mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]];

                    }
                }
            }

           if (typeSc==2) {

                for (int ind_f = 0 ; ind_f < imax2(1,dimT[0]) ; ind_f++)
                for (int ind_m = 0 ; ind_m < imax2(1,dimT[1]) ; ind_m++)
                for (int ind_i = 0 ; ind_i < imax2(1,dimT[2]) ; ind_i++) {
                    if (!ISNA(mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]])) {

                        target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] =
                        target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] +
                        mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]];

                    }
                }
            }

           if (typeSc==3) {

                for (int ind_f = 0 ; ind_f < imax2(1,dimT[0]) ; ind_f++)
                for (int ind_m = 0 ; ind_m < imax2(1,dimT[1]) ; ind_m++)
                for (int ind_i = 0 ; ind_i < imax2(1,dimT[2]) ; ind_i++) {
                    if (!ISNA(mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]])) {

                        target[ind_f*ftarg[0] + ind_m*ftarg[1] + ind_i*ftarg[2] + ind_t*ftarg[3]] =
                        mult[ind_f*fmult[0] + ind_m*fmult[1] + ind_i*fmult[2] + ind_t*fmult[3]];

                    }
                }
            }

           UNPROTECT(7);

        }

        UNPROTECT(1);
    }

    UNPROTECT(2);
}

}}

//---------------------------------
//
// Module de marche de quotas
//
//---------------------------------

extern "C" {

void BioEcoPar::QuotaMarket(SEXP list, SEXP pQuotaIni, SEXP pQuotaMin, SEXP pQuotaMax, double lambdaQ, double sdmax, double ftol, int itmax, SEXP paramBehav, int ind_t, int persCalc ) //ind_t>0
{

// ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.QuotaMarket.txt", ios::out | ios::trunc);
//ofstream fichier;
//if (ind_t ==1) fichier.open ("C:\\Users\\fbriton\\Dropbox\\These\\IAM_Dvt\\test.QuotaMarket.txt");
//fichier << "Start" << endl;
 //time_t my_time;

    SEXP listTemp, eVarCopy, alphaBhv, pQuota ,nam_eQuota,nam_eQuota_dyn,
        dimCstF,dimCstFM, nDimT,nDimT2, eFACTf, eFACTfm;

    PROTECT(listTemp = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));

    double *g_effSup = REAL(effSupMat); // specified in objArgs@arguments$Gestion$effSup = relates to effort1*effort2
//fichier << "QM0.1"  << endl;
    double *r_out_effort1_f = REAL(NBDSF); // efforts from output for past values
    double *r_out_effort2_f = REAL(EFF2F);
    double *r_out_effort1_f_m = REAL(NBDSFM);
    double *r_out_effort2_f_m = REAL(EFF2FM);
//fichier << "QM0.2"  << endl;

    double *g_effort1FM = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f_m")); // efforts from list to update to then send to other modules
    double *g_effort1F = REAL(getListElement(getListElement(listTemp, "Fleet"), "effort1_f"));
//fichier << "QM0.3"  << endl;
    double* r_L_f_m_e;

    SEXP rtbs_f_m_out, rtbs_f_out, ccw_f_out, rep_f, gc_f, fixc_f, dep_f, ic_f, GVLtot_f_m_out, cshrT_f_m_out,GVLtot_f_out;
    double *r_rtbs_f_m_out, *r_rtbs_f_out, *r_ccw_f_out, *r_rep_f, *r_gc_f, *r_fixc_f, *r_dep_f, *r_ic_f, *r_GVLtot_f_m_out, *r_GVLtot_f_m_e_ref, *r_P_f_m_e,*r_cshrT_f_m_out,*r_GVLtot_f_out;
    PROTECT(rep_f = getListElement(getListElement(listTemp, "Fleet"), "rep_f"));
    r_rep_f = REAL(rep_f);
    PROTECT(gc_f = getListElement(getListElement(listTemp, "Fleet"), "gc_f"));
    r_gc_f = REAL(gc_f);
    PROTECT(fixc_f = getListElement(getListElement(listTemp, "Fleet"), "fixc_f"));
    r_fixc_f = REAL(fixc_f);
    PROTECT(dep_f = getListElement(getListElement(listTemp, "Fleet"), "dep_f"));
    r_dep_f = REAL(dep_f);
    PROTECT(ic_f = getListElement(getListElement(listTemp, "Fleet"), "ic_f"));
    r_ic_f = REAL(ic_f);
    PROTECT(GVLtot_f_m_out = VECTOR_ELT(out_EcoDCF, 6));
    r_GVLtot_f_m_out = REAL(GVLtot_f_m_out);
    PROTECT(rtbs_f_m_out = VECTOR_ELT(out_EcoDCF, 15));
    r_rtbs_f_m_out = REAL(rtbs_f_m_out);
    PROTECT(ccw_f_out = VECTOR_ELT(out_EcoDCF, 29));
    r_ccw_f_out = REAL(ccw_f_out);
    PROTECT(cshrT_f_m_out = VECTOR_ELT(out_EcoDCF, 18));
    r_cshrT_f_m_out = REAL(cshrT_f_m_out);
    PROTECT(rtbs_f_out = VECTOR_ELT(out_EcoDCF, 16));
    r_rtbs_f_out = REAL(rtbs_f_out);

    SEXP dc_rep_f ,dc_gc_f,dc_fixc_f,dc_dep_f;
    int *dim_rep_f,*dim_gc_f,*dim_fixc_f,*dim_dep_f;
    PROTECT(dc_rep_f = iDim(INTEGER(getAttrib(rep_f, install("DimCst")))));
    PROTECT(dc_gc_f = iDim(INTEGER(getAttrib(gc_f, install("DimCst")))));
    PROTECT(dc_fixc_f = iDim(INTEGER(getAttrib(fixc_f, install("DimCst")))));
    PROTECT(dc_dep_f = iDim(INTEGER(getAttrib(dep_f, install("DimCst")))));
    dim_rep_f = INTEGER(dc_rep_f);
    dim_gc_f = INTEGER(dc_gc_f);
    dim_fixc_f = INTEGER(dc_fixc_f);
    dim_dep_f = INTEGER(dc_dep_f);

    SEXP PQuot_temp, diffLQ;
    double *r_PQuot_temp, *r_diffLQ;


//fichier << "QM0.4"  << endl;

    SEXP ProfUE_f_m, ProfUE_f_m_ctr, ProfUE_f, ExpProfUE_f, Profmin_f;
    PROTECT(ProfUE_f_m = NEW_NUMERIC(nbF*nbMe));
    double *r_ProfUE_f_m = REAL(ProfUE_f_m);
    PROTECT(ProfUE_f_m_ctr = NEW_NUMERIC(nbF*nbMe));
    double *r_ProfUE_f_m_ctr = REAL(ProfUE_f_m_ctr);
    PROTECT(ProfUE_f = NEW_NUMERIC(nbF));
    double *r_ProfUE_f = REAL(ProfUE_f);
    double *r_allocEff_f_m = REAL(out_allocEff_fm);
    PROTECT(ExpProfUE_f = NEW_NUMERIC(nbF));
    double *r_ExpProfUE_f = REAL(ExpProfUE_f);
//fichier << "QM0.4"  << endl;
    double *r_GoFish_f = REAL(intermGoFish);
//    PrintValue(intermGoFish);
    PROTECT(Profmin_f = NEW_NUMERIC(nbF));
    double *r_Profmin_f = REAL(Profmin_f);
    double ratio_p;
    double *rans_multPrice;
    double GVLtot_temp_f_m;
//fichier << "QM0.5"  << endl;
    int ind_t_last, ind_t_price;

    PROTECT(alphaBhv = getListElement(paramBehav, "ALPHA"));
    double *r_alphaBhv = REAL(alphaBhv);
//fichier << "QM0.6"  << endl;

    PROTECT(dimCstF = allocVector(INTSXP, 4));
    PROTECT(dimCstFM = allocVector(INTSXP, 4));
    int *dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
    int *dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;
    PROTECT(eFACTf = iDim(dCF));
    PROTECT(eFACTfm = iDim(dCFM));
    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    //fichier << "dCFM[0]: " << dCFM[0] << "; dCFM[1]: " << dCFM[1] << "; dCFM[2]: " << dCFM[2] << "; dCFM[3]: " << dCFM[3] << endl;
    //fichier << "eF_fm[0]: " << eF_fm[0] << "; eF_fm[1]: " << eF_fm[1] << "; eF_fm[2]: " << eF_fm[2] << "; eF_fm[3]: " << eF_fm[3] << endl;

    PROTECT(nDimT = allocVector(INTSXP,4));
    int *ndT = INTEGER(nDimT); ndT[0] = 0; ndT[1] = 0; ndT[2] = 0; ndT[3] = nbT;
    PROTECT(nDimT2 = allocVector(INTSXP,2));
    int *ndT2 = INTEGER(nDimT2); ndT2[0] = 0; ndT2[1] = nbT;

    double *r_pQuota;
    double* r_pQuotaIni;
    SEXP ans_multPrice;
    int ind_ePrice;
    SEXP Lmodel,TACmodel;
    double *r_Lmodel, *r_TACmodel  ;
    SEXP reducelambda;
    PROTECT(reducelambda = NEW_LOGICAL(nbEQuotaMarket_dyn));
    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) LOGICAL(reducelambda)[int_eQuota_dyn] = false;

    double lambdaQ_iter;

    //---------------------------------------------
    // 1 - Tatonnement of quota market
    //---------------------------------------------

    // Initialize quota price
//    fichier << "QM0.7"  << endl;
    //fichier << "nbEQuotaMarket:" << nbEQuotaMarket  << endl;

    // year used to calculate ProfUE_fm and base fishing decisions:
    ind_t_last=0 ; // initial year
    //if you want it to be the last year the vessel has been active, uncomment lines 14774-14783

    // year used to retrieve fish prices and calculate GVL (mostly to account for dynamics in fish price while using ind_t_last=0
    if (ind_t == 1){ind_t_price = 0;} else {ind_t_price = 1;} // to capture market dynamics
    //ind_t_price = ind_t-1; // previous year (less stable than previous option)
//    fichier << "ind_t_price: " << ind_t_price << endl;

    for (int int_eQuota = 0 ; int_eQuota  < nbEQuotaMarket ; int_eQuota++) { // initialize quota prices
            ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListQM,int_eQuota))); // index of species in MultPrice
            PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
            PROTECT(pQuota = VECTOR_ELT(out_PQuot_et,int_eQuota));
            r_pQuota = REAL(pQuota);
            r_pQuotaIni = REAL(VECTOR_ELT(pQuotaIni,int_eQuota));
            r_pQuota[ind_t] = r_pQuotaIni[ind_t] * REAL(ans_multPrice)[ind_t_price]; // adjust quota price for change in fish price in previous year (only value when entering the loop, but to keep ratio quota price/fish price stable)
            UNPROTECT(2);
    }
//                PrintValue(sppListAll);
//                PrintValue(multPrice);
//                PrintValue(sppListQM);
//                PrintValue(out_PQuot_et);

//    fichier << "QM0.8"  << endl;

    bool toinit;
    bool keep;
    bool GoOn = true;

    int itQ=0;
    double min_sum_diffLQ, sum_diffLQ, nb_overTAC, max_diffLQ, min_max_diffLQ, max_nb_overTAC;
    int min_itQ;
    //fichier << "itmax" << itmax  << endl;

    //Calculate L_fme for quoted species
    SEXP L_f_m_e, L_f_m_e_i, ans_L_f_m_e;
    PROTECT(L_f_m_e = allocVector(VECSXP, nbEQuotaMarket));
    setAttrib(L_f_m_e, R_NamesSymbol, sppListQM);

     for (int ind_eQ = 0 ; ind_eQ< nbEQuotaMarket ; ind_eQ++) {
        PROTECT(nam_eQuota = STRING_ELT(sppListQM,ind_eQ));

        if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))){ // espece dyn
                PROTECT(L_f_m_e_i = getListElement(out_L_efmit, CHAR(nam_eQuota)));
                PROTECT(ans_L_f_m_e = aggregObj(L_f_m_e_i,dimCstFM));
        } else{ // espece stat
                PROTECT(ans_L_f_m_e = getListElement(out_Lstat, CHAR(nam_eQuota)));
        }

        SET_VECTOR_ELT(L_f_m_e, ind_eQ, ans_L_f_m_e);

        UNPROTECT(2);
        if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))) UNPROTECT(1);
     }
     //PrintValue(L_f_m_e);


    while (GoOn & (itQ<itmax)){
//            fichier << "itQ: " << itQ << endl;
            //Rprintf("itQ = %i \n", itQ);

//            fichier << "QM1.1"  << endl;
            // Initialization ProfUE_f_m
            for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = 0.0;

                }
                r_ProfUE_f[ind_f] = 0.0;
                r_ExpProfUE_f[ind_f] = 0.0;
            }
            //fichier << "QM1.2"  << endl;
           //Rprintf("QM1.2\n");


        // Calculate metier profitability for each fleet and metier
       // my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                //if (ind_f==5) fichier << "f: " << CHAR(STRING_ELT(fleetList,ind_f)) << endl;
                toinit = true;

                // retrieve last ind_t when active
//                keep = true;
//                ind_t_last = ind_t-1;
//                while(keep){
//                    if (r_GoFish_f[ind_f+nbF*ind_t_last]>0){
//                            keep = false;
//                    } else{
//                    ind_t_last = ind_t_last -1;}
//
//                }

                //fichier << "QM1.3"  << endl;
                //Rprintf("QM1.3\n");
//                if ((itQ==0) & (ind_t<5)) Rprintf("f: %d, ind_t_last: %d \n", ind_f, ind_t_last);

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
//                    if (ind_f==6) fichier << "ind_m: " << ind_m << endl;

                     // Adjust GVL by species to account for market (fish price) dynamics
                        for (int e = 0 ; e < nbE+nbEstat ; e++) {
                                if ((nbE>0) & (e<nbE)) {
                                    r_GVLtot_f_m_e_ref = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
                                    ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppList,e)));
                                    //Rprintf("e = %s, ind_ePrice = %d \n", CHAR(STRING_ELT(sppList,e)),ind_ePrice);
                                    PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
                                    rans_multPrice = REAL(ans_multPrice);

                                } if ((nbEstat>0) & (e>=nbE)) {
                                    r_GVLtot_f_m_e_ref = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));
                                    ind_ePrice = getVectorIndex(sppListAll,CHAR(STRING_ELT(sppListStat,e-nbE)));
                                    //Rprintf("e = %s, ind_ePrice = %d \n", CHAR(STRING_ELT(sppListStat,e-nbE)),ind_ePrice);
                                    PROTECT(ans_multPrice = VECTOR_ELT(multPrice,ind_ePrice));
                                    rans_multPrice = REAL(ans_multPrice);
                                    }
//                                    fichier << "QM1.4"  << endl;

                        ratio_p = rans_multPrice[ind_t_price] / rans_multPrice[ind_t_last];

                        if (!ISNA(r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]* ratio_p)){
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] +
                                r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * ratio_p;} //adjustment for variations in fish price

                        UNPROTECT(1);

//                        if (ind_f==6) fichier << "e = " <<  e <<
//                            "; GVL init = " << r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]  <<
//                            "; ratio = " << ratio_p <<
//                            "; GVL actual = " << r_GVLtot_f_m_e_ref[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * ratio_p <<
//                            "; ProfUE_f_m = " << r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
                        }
                        GVLtot_temp_f_m = r_ProfUE_f_m[ind_f + nbF*ind_m]; // temporary save to use to calculate crew share if persCalc=5

                        // Deduce variable costs : RTBS
                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                            (r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] - r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]);
//                       if (ind_f==6) fichier << "Deduce var costs : r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;

                        //Deduce costs
                         // 1- crew costs
                        if ((persCalc == 1) | (persCalc == 2) | (persCalc == 3) | (persCalc == 4)) { // crew share RTBS
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] *
                                        (1 - r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]);
                                        }else if (persCalc == 5){ // crew share GVL
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                        (r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] /
                                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *
                                        GVLtot_temp_f_m);
                        } else { // fixed wages
                            r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                                          r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]  *
                                                          (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) /
                                                          (r_out_effort1_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] * r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]);
                        }

//                        if (ind_f==6) {
//                                if ((persCalc == 1) | (persCalc == 2) | (persCalc == 3) | (persCalc == 4)) {
//                                fichier << " cshrT = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " ; rtbs = " <<  r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " crew share = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl; } else if (persCalc == 5){
//                                    fichier << " cshrT ref = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " ; GVL ref = " <<  r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                 " crew share = " << r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] / r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] <<
//                                "; GVL  =  " <<  GVLtot_temp_f_m <<
//                                    "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                                    } else{
//                                fichier << " crew costs_f =  " <<  r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] <<
//                                "; ratio Effort" << (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) /
//                                                          (r_out_effort1_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]] * r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t_last*eF_f[3]]) <<
//                                "; Deduce crew costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                                }}

                        // 2- Quota costs

                        for (int ind_eQ = 0 ; ind_eQ< nbEQuotaMarket ; ind_eQ++) { //deduce quota expenses for dynamic markets
                                PROTECT(nam_eQuota = STRING_ELT(sppListQM,ind_eQ));

                                r_L_f_m_e = REAL(getListElement(L_f_m_e, CHAR(nam_eQuota)));

//                                if (ind_f==6) fichier << "nam_eQuota: " << CHAR(nam_eQuota) << endl;


                                if (!ISNA(r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]])){

                                    PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota)));
                                    r_pQuota = REAL(pQuota);

                                    r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m]  -
                                                                r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *1000 * r_pQuota[ind_t];
//                                    if (ind_f==6){
//                                            fichier << "index: " << ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3] << endl;
//                                            fichier << "r_L_f_m_e: " << r_L_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] *1000 << endl;
//                                            fichier << "r_pQuota: " << r_pQuota[ind_t] << endl;
//                                            fichier << "Deduce quota costs: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;}


                                    UNPROTECT(1);
                               }
                               UNPROTECT(1);

                        }


                        // Divide by effort
                        if(r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]>0){
                                r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] /
                                (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]]) ;
                        }

//                        if (ind_f==6){
//                            fichier << "Effort_f_m = " << r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t_last*eF_fm[3]] << endl;
//                            fichier << "Divide by effort: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;
//                       }



                        // 3- fixed costs
//                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
//                                                          (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
//                                                           r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
//                                                           r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]) /
//                                                          g_effSup[ind_f + nbF*ind_t];

                        r_ProfUE_f_m[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m] -
                                                          (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
                                                           r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
                                                           r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]+
                                                           r_dep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]]+
                                                           r_ic_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]]) /
                                                          g_effSup[ind_f + nbF*ind_t];

//                        if (ind_f==6) fichier << " fixed costs_f =  " <<  r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] +
//                                                                            r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t_last*dim_fixc_f[3]] +
//                                                                            r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t_last*dim_gc_f[3]]+
//                                                                            r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t_last*dim_dep_f[3]]+
//                                                                            r_ic_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t_last*dim_rep_f[3]] <<
//                                "; Eff sup = " << g_effSup[ind_f + nbF*ind_t] <<
//                                "; Deduce fixed costs per UE: r_ProfUE_f_m =  " <<  r_ProfUE_f_m[ind_f + nbF*ind_m] << endl;




                        if (!ISNA(r_ProfUE_f_m[ind_f + nbF*ind_m])){
                            if (toinit){
                                r_Profmin_f[ind_f] = r_ProfUE_f_m[ind_f + nbF*ind_m];
                                toinit=false;
                            }
                            if(r_ProfUE_f_m[ind_f + nbF*ind_m] < r_Profmin_f[ind_f]) r_Profmin_f[ind_f] =  r_ProfUE_f_m[ind_f + nbF*ind_m];
                        }
                }
        }
//        fichier << "QM1.4"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
        //Rprintf("QM1.4\n");

        // Center Profitabilities from the minimal value so that they are all positive
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

            for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] = r_ProfUE_f_m[ind_f + nbF*ind_m]  - r_Profmin_f[ind_f];
                       if(!ISNA(r_ProfUE_f_m_ctr[ind_f + nbF*ind_m])) r_ProfUE_f[ind_f] = r_ProfUE_f[ind_f] + r_ProfUE_f_m_ctr[ind_f + nbF*ind_m];
            }
            if (r_ProfUE_f[ind_f] ==0.0) { //in case all metiers have a profitability equal to the minimum one
                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       if(r_ProfUE_f_m_ctr[ind_f + nbF*ind_m]==0.0) {
                            r_ProfUE_f_m_ctr[ind_f + nbF*ind_m]=1.0;
                            r_ProfUE_f[ind_f] = r_ProfUE_f[ind_f] + r_ProfUE_f_m_ctr[ind_f + nbF*ind_m];
                       }
                }
            }
        }
//        fichier << "QM2"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
        //Rprintf("QM2\n");

        // Update effort allocation: (weighted average between profitability and habit)
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {
                       r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_alphaBhv[ind_f + nbF*ind_t] * r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] / r_ProfUE_f[ind_f] +
                                                           (1-r_alphaBhv[ind_f + nbF*ind_t]) * (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) /
                                                           (r_out_effort1_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]  * r_out_effort2_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]);
                        if (ISNA(r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t])) r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.0;
//                        if(ind_f ==6) fichier << "m= " << ind_m << "; attract = " << r_ProfUE_f_m_ctr[ind_f + nbF*ind_m] / r_ProfUE_f[ind_f] << "; trad = " << (r_out_effort1_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]] * r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) /
//                                                          (r_out_effort1_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]  * r_out_effort2_f[ind_f*eF_f[0] + ind_m*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]]) << "; allocEff = " << r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] << endl;
                        if (!ISNA(r_ProfUE_f_m[ind_f + nbF*ind_m])) r_ExpProfUE_f[ind_f] = r_ExpProfUE_f[ind_f] + r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * r_ProfUE_f_m[ind_f + nbF*ind_m];
//                        if (ind_f ==6){ fichier << "r_allocEff_f_m = " << r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] << "; r_ProfUE_f_m = " << r_ProfUE_f_m[ind_f + nbF*ind_m] << "; r_ExpProfUE_f = " << r_ExpProfUE_f[ind_f] << endl;}

                }

          //Decision to go fishing or not based on expected profitability
                if(r_ExpProfUE_f[ind_f]>=0.0) {
                        r_GoFish_f[ind_f+nbF*ind_t] = 1.0;
                } else {
                    r_GoFish_f[ind_f+nbF*ind_t] = 0.0;
                    }
        }
//        fichier << "QM3"  << endl;
        //my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;
         //Rprintf("QM3\n");

        // Update effort
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

                for (int ind_m = 0 ; ind_m< nbMe ; ind_m++) {

                       g_effort1FM[ind_f + nbF*ind_m] = r_GoFish_f[ind_f+nbF*ind_t] * g_effSup[ind_f + nbF*ind_t] * r_allocEff_f_m[ind_f + nbF*ind_m + nbF*nbMe*ind_t] / r_out_effort2_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]];
//                      if(ind_f ==6) fichier << "m= " << ind_m << "; Eff = " << g_effort1FM[ind_f + nbF*ind_m] << endl;
                }
                g_effort1F[ind_f] = r_GoFish_f[ind_f+nbF*ind_t] * g_effSup[ind_f + nbF*ind_t] / r_out_effort2_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]];
//                if(ind_f ==6) fichier << "; Eff_f = " << g_effort1F[ind_f] << endl;

        }
        //PrintValue(getListElement(getListElement(listTemp, "Fleet"), "effort1_f"));
        //PrintValue(getListElement(getListElement(listTemp, "Fleet"), "effort1_f_m"));


//       fichier << "QM4.1"  << endl;
        // Rprintf("QM4.1\n");

        Mortalite(listTemp, ind_t, eVarCopy);
//        fichier << "QM4.2"  << endl;
        // Rprintf("QM4.2\n");


        DynamicPop(listTemp, ind_t, eVarCopy,false);
//        fichier << "QM4.3"  << endl;
        // Rprintf("QM4.3\n");

        CatchDL(listTemp, ind_t, eVarCopy);
//        fichier << "QM4.4"  << endl;
        // Rprintf("QM4.4\n");
       // my_time = time(NULL);
        //fichier << "time: " << ctime(&my_time) << endl;

        // Adjust quota prices for all species
        GoOn = false;
        //lambdaQ_iter = lambdaQ * ( 1 - ((double)rand() / (double)RAND_MAX) * pow(((double)itQ/(double)itmax),2.0) * sdmax);
        lambdaQ_iter = lambdaQ * ( 1 - pow(((double)itQ/(double)itmax),2.0) * sdmax); //increase lambda with itQ to avoid oscillating between the same values towards the end of convergence (optionnal, you can also keep it constant when commenting this line)

        //fichier << "LambdaQ = "<< lambdaQ <<  "; itQ = "  << itQ  <<  "; Lambda iter = "<< lambdaQ_iter <<endl;
        sum_diffLQ = 0.0;
        max_diffLQ=0.0;
        nb_overTAC = 0.0;

        for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {

                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                //fichier << "nam_eQuota_dyn: " << CHAR(nam_eQuota_dyn) << endl;
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));

                if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(Lmodel = aggregObj(getListElement(out_L_eit, CHAR(nam_eQuota_dyn)),nDimT));
                                    //PROTECT(Lmodel = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),nDimT));
                                } else{ // espece stat
                                    PROTECT(Lmodel = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),nDimT));
                                }

                //PrintValue(Lmodel);
                r_Lmodel = REAL(Lmodel);
//                fichier << "QM5.1"  << endl;
                PROTECT(TACmodel = getListElement(TAC, CHAR(nam_eQuota_dyn)));
                //PROTECT(TACmodel = aggregObj(getListElement(TACbyF, CHAR(nam_eQuota_dyn)),nDimT)); //ne fonctionne plus si quota pas detenu par navires modelises
                r_TACmodel = REAL(TACmodel);
                r_pQuota = REAL(pQuota);

//                fichier << "QM5.2"  << endl;

                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                r_PQuot_temp = REAL(PQuot_temp);
                PROTECT(diffLQ = getListElement(out_diffLQ, CHAR(nam_eQuota_dyn)));
                r_diffLQ = REAL(diffLQ);
//                fichier << "QM5.3"  << endl;

                r_diffLQ[itQ + itmax*ind_t] =  (r_Lmodel[ind_t] - r_TACmodel[ind_t]) / r_TACmodel[ind_t]; //save value for current iteration to check algorithm convergence
                sum_diffLQ = sum_diffLQ + fabs(r_diffLQ[itQ + itmax*ind_t]);

//                fichier << "L_eit = " << r_Lmodel[ind_t] << "; L_efmit = " <<  REAL(aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),nDimT))[ind_t] <<
//                 "; TAC = " << r_TACmodel[ind_t] << "; diff = " << r_diffLQ[itQ + itmax*ind_t] << endl;

                if (int_eQuota_dyn ==0){
                      max_diffLQ = fabs(r_diffLQ[itQ + itmax*ind_t]);
                } else if (fabs(r_diffLQ[itQ + itmax*ind_t]) > max_diffLQ) max_diffLQ = fabs(r_diffLQ[itQ + itmax*ind_t]);

                if (r_diffLQ[itQ + itmax*ind_t]>0) nb_overTAC = nb_overTAC+1;
                r_PQuot_temp[itQ + itmax*ind_t] = r_pQuota[ind_t]; //save value for current iteration to check algorithm convergence

//                fichier << "QM5"  << endl;
                    if ((!LOGICAL(reducelambda)[int_eQuota_dyn]) & (r_diffLQ[itQ + itmax*ind_t]*r_diffLQ[itQ-1 + itmax*ind_t] < 0)) LOGICAL(reducelambda)[int_eQuota_dyn] = true;
                    //Rprintf("Reduce lambda at iter %i: %d \n",itQ,LOGICAL(reducelambda)[int_eQuota_dyn]);


                    if(!LOGICAL(reducelambda)[int_eQuota_dyn]){
                        if ((lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t]) > 0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 1.2);
                        } else if ((lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t]) < -0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 0.8);
                        } else {r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * (1 + lambdaQ_iter * r_diffLQ[itQ + itmax*ind_t])); }
                    } else {
                        if ((lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t]) > 0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 1.2);
                        } else if ((lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t]) < -0.2){
                            r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * 0.8);
                        } else {r_pQuota[ind_t] = fmax2(0, r_pQuota[ind_t] * (1 + lambdaQ_iter/5 * r_diffLQ[itQ + itmax*ind_t])); }}


                    if((fabs(r_diffLQ[itQ + itmax*ind_t]) > ftol) | (itQ < floor(itmax/2))) GoOn=true;

                    //fichier << "itQ: " << itQ << ", diffLQ = " << r_diffLQ[itQ + itmax*ind_t] << "Quota price = "<< r_pQuota[ind_t] << "GoOn =" << GoOn << endl;
                UNPROTECT(6);
        }
//        fichier << "QM5"  << endl;
//        if (itQ == 0){ //floor(itmax/2)
//            min_sum_diffLQ = sum_diffLQ;
//            min_itQ = itQ;
//            max_nb_overTAC = nb_overTAC;
//        } else if (itQ > 0){
//            if(nb_overTAC > max_nb_overTAC) {
//                min_sum_diffLQ = sum_diffLQ;
//                min_itQ = itQ;
//                max_nb_overTAC = nb_overTAC;
//            } else if ((nb_overTAC == max_nb_overTAC) & (sum_diffLQ < min_sum_diffLQ)){
//                min_sum_diffLQ = sum_diffLQ;
//                min_itQ = itQ;
//            }

//fichier << "itQ = "  << itQ << ", max_diffLQ = " << max_diffLQ << ", min_max_diffLQ = " << min_max_diffLQ << endl;
//Save best iteration (the one with the minimum max_diffLQ (=max difference between landing and TAC across stocks)):
        if (itQ == 0){ //floor(itmax/2)
            min_itQ = itQ;
            min_max_diffLQ = fabs(max_diffLQ);
            min_sum_diffLQ = sum_diffLQ;
            } else if (fabs(max_diffLQ) < min_max_diffLQ) {
                //Rprintf("Cas 1 , old diffLQ = %f; new diffLQ = %f; min itQ = %d \n",min_max_diffLQ,fabs(max_diffLQ) ,itQ);
                min_itQ = itQ;
                min_max_diffLQ = fabs(max_diffLQ);
                min_sum_diffLQ = sum_diffLQ;
            } else if ((fabs(max_diffLQ) == min_max_diffLQ) & (sum_diffLQ < min_sum_diffLQ)){
                //Rprintf("Cas 2 , diffLQ = %f; old sum_diffLQ = %f; new sum_diffLQ = %f; min itQ = %d \n",fabs(max_diffLQ),min_sum_diffLQ,sum_diffLQ ,itQ);
                min_itQ = itQ;
                min_max_diffLQ = fabs(max_diffLQ);
                min_sum_diffLQ = sum_diffLQ;
            }
//fichier << "min_max_diffLQ = " << min_max_diffLQ << endl;

         //rerun the best iteration which will be saved at index itmax-1 (last iteration)
        if (((GoOn == false) & (itQ<itmax-1)) | (itQ == itmax-2)){
            GoOn = true;
            itQ = itmax-2;

            for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                REAL(pQuota)[ind_t] = REAL(PQuot_temp)[min_itQ + itmax*ind_t]; // retrieve quota price for the best iteration min_itQ

                UNPROTECT(3);
            }
        }

        // for the last iteration, save the final price
        if (itQ == itmax-1){
            for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
                PROTECT(PQuot_temp = getListElement(out_PQuot_temp, CHAR(nam_eQuota_dyn)));
                REAL(pQuota)[ind_t] = REAL(PQuot_temp)[itQ + itmax*ind_t];

                UNPROTECT(3);
            }
        }

        itQ = itQ+1;


    }

    //---------------------------------------------
    // 2 - Final adjustment quota prices to account for current stock status and market
    //---------------------------------------------
    // Run Market and EcoDCF modules
    double max_ratio_GOS_QuotaExp;
    SEXP GOS_f,QuotaExp_f;
    double *r_GOS_f, *r_QuotaExp_f;

    Marche(listTemp, ind_t);

    GoOn=TRUE;
    itQ=0;

    while(GoOn && (itQ<10)){
//            fichier << "itQ = " << itQ << endl;

        EcoDCF(listTemp, ind_t, EcoIndCopy[4], drCopy);

        PROTECT(GOS_f = VECTOR_ELT(out_EcoDCF, 35));
        r_GOS_f = REAL(GOS_f);
        PROTECT(QuotaExp_f = VECTOR_ELT(out_EcoDCF, 59));
        r_QuotaExp_f = REAL(QuotaExp_f);

        max_ratio_GOS_QuotaExp = 1.0;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                if ((r_GoFish_f[ind_f+nbF*ind_t]==1) && (r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] > max_ratio_GOS_QuotaExp)){
                    max_ratio_GOS_QuotaExp = r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]  ;
                }
//                fichier << "GoFish = " << r_GoFish_f[ind_f+nbF*ind_t] <<
//                "; GOS = " <<  r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]<<
//                "; QuotaExp = " << r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] <<
//                "; ratio = " << r_QuotaExp_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GOS_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]<<
//                "; max ratio = " << max_ratio_GOS_QuotaExp << endl;
        }

        if (max_ratio_GOS_QuotaExp == 1.0){ //OK
                GoOn=FALSE;
        } else{// Adjust Quota price by the ratio
                for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
                        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
                        PROTECT(pQuota = getListElement(out_PQuot_et,CHAR(nam_eQuota_dyn)));
//                        fichier << "Pquot before = " << REAL(pQuota)[ind_t] << endl;
                        REAL(pQuota)[ind_t] = REAL(pQuota)[ind_t]/(1.05 * max_ratio_GOS_QuotaExp);
//                        fichier << " PQuota after = " << REAL(pQuota)[ind_t] << endl;

                        UNPROTECT(2);
                    }
        }
//        fichier << "GoOn = " << GoOn << endl;
        itQ ++ ;
        UNPROTECT(2);
    }



    //---------------------------------------------
    // 3 - Trade quotas
    //---------------------------------------------

    // Rank fleet in order of profitability
    SEXP rank_Prof_f, ExpProfUE_f_copy;
    PROTECT(rank_Prof_f = NEW_INTEGER(nbF));
    int* r_rank_Prof_f = INTEGER(rank_Prof_f);
    ExpProfUE_f_copy = duplicate(ExpProfUE_f);
    double* r_ExpProfUE_f_copy = REAL(ExpProfUE_f_copy);

    double maxTemp = r_ExpProfUE_f_copy[0];
    int ind_maxTemp ;

    for (int ind_rank = 0 ; ind_rank < nbF ; ind_rank++){
        toinit=true;
        ind_maxTemp = 0;
        while ((ind_maxTemp <nbF) & toinit){
            if (!ISNA(r_ExpProfUE_f_copy[ind_maxTemp])){
                toinit=false;
                maxTemp = r_ExpProfUE_f_copy[ind_maxTemp];
            } else{ind_maxTemp++;}
        }

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            if ((!ISNA(r_ExpProfUE_f_copy[ind_f])) & (r_ExpProfUE_f_copy[ind_f] >= maxTemp)){
                maxTemp = r_ExpProfUE_f_copy[ind_f] ;
                ind_maxTemp = ind_f;
            }
        }
        r_ExpProfUE_f_copy[ind_maxTemp] = NA_REAL;
        r_rank_Prof_f[ind_rank] = ind_maxTemp; //index of the fleet of rank ind_rank
    }
    //PrintValue(ExpProfUE_f);
    //PrintValue(rank_Prof_f);

//    fichier << "QM6"  << endl;


    if (false){ //Old version

    SEXP demandQ_f_e, offerQ_f_e, L_ef,TACbyF_e,TACbyF_ex, QuotaTrade_fe ;
    PROTECT(demandQ_f_e = NEW_NUMERIC(nbF));
    PROTECT(offerQ_f_e = NEW_NUMERIC(nbF));
    double *r_L_ef;
    double *r_demandQ_f_e = REAL(demandQ_f_e);
    double *r_offerQ_f_e = REAL(offerQ_f_e);
    double *r_TACbyF_e;
    double *r_TACbyF_ex;
    double *r_QuotaTrade_fe;
    int ind_buyer, ind_seller;

    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
        //fichier << "nam_eQuota_dyn: " << CHAR(nam_eQuota_dyn) << endl;
        if (!isNull(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(L_ef = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),dimCstF));
                } else{ // espece stat
                                    PROTECT(L_ef = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),dimCstF));
                                }
        r_L_ef = REAL(L_ef);
        PROTECT(TACbyF_e = getListElement(TACbyF, CHAR(nam_eQuota_dyn)));
        r_TACbyF_e = REAL(TACbyF_e);
        PROTECT(TACbyF_ex = duplicate(TACbyF_e));
        r_TACbyF_ex = REAL(TACbyF_ex);

        //fichier << "QM7.1"  << endl;
//        if (ind_t==4){
//                Rprintf("sp = %s Expected catches: \n",CHAR(nam_eQuota_dyn));
//                PrintValue(L_ef);
//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Expected catches:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_L_ef[ind_f+nbF*ind_t]  << endl;
//                }
//        }

        // Calculate net demand for quota by fleet = expected catches - holdings
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            r_demandQ_f_e[ind_f] = r_L_ef[ind_f+nbF*ind_t] - r_TACbyF_e[ind_f+nbF*ind_t] ;
            r_offerQ_f_e[ind_f] = - r_demandQ_f_e[ind_f] ;
        }
        //fichier << "QM7.2"  << endl;
        //Rprintf("TAC before trade\n");
        //PrintValue(TACbyF_e);
        //Rprintf("Expected catches\n");
        //PrintValue(L_ef);
        //Rprintf("Demand before trade\n");
        //PrintValue(demandQ_f_e);

        //Trade quotas
        int rank_buyer = 0;
        ind_buyer = r_rank_Prof_f[rank_buyer];
        int rank_seller = nbF-1;
        ind_seller = r_rank_Prof_f[rank_seller];
        double trade;



                while ((rank_buyer < nbF) &&  (rank_seller >=0)){
                        //fichier << "rank buyer = " << rank_buyer  << "; ind buyer = " << ind_buyer << endl;
                    while ((rank_seller >=0) && (r_demandQ_f_e[ind_buyer] > 0)){
                        //fichier << "rank seller = " << rank_seller << "; ind seller = " << ind_seller << endl;

                        if(r_offerQ_f_e[ind_seller] > 0){
                            trade = fmin2(r_demandQ_f_e[ind_buyer],r_offerQ_f_e[ind_seller]);
                            //fichier << "demand = " << r_demandQ_f_e[ind_buyer] << ", offer = " << r_offerQ_f_e[ind_seller] << ", trade = " << trade << endl;
                            r_demandQ_f_e[ind_buyer] = r_demandQ_f_e[ind_buyer] - trade;
                            r_offerQ_f_e[ind_seller] = r_offerQ_f_e[ind_seller] - trade;

                            r_TACbyF_e[ind_buyer+nbF*ind_t] = r_TACbyF_e[ind_buyer+nbF*ind_t] + trade;
                            r_TACbyF_e[ind_seller+nbF*ind_t] = r_TACbyF_e[ind_seller+nbF*ind_t] - trade;
                        } else {
                            rank_seller --;
                            ind_seller = r_rank_Prof_f[rank_seller];
                        }
                    }
                    rank_buyer ++;
                    ind_buyer = r_rank_Prof_f[rank_buyer];

                }



//        fichier << "QM7.4"  << endl;
        //Rprintf("Demand after trade\n");
        //PrintValue(demandQ_f_e);
        //Rprintf("Offer after trade\n");
        //PrintValue(offerQ_f_e);
        //Rprintf("TAC after trade\n");
        //PrintValue(TACbyF_e);
//
//        if (ind_t==4){
//                Rprintf("sp = %s , Quota after trade: \n",CHAR(nam_eQuota_dyn));
//                PrintValue(TACbyF_e);
//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Quota after trade:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_TACbyF_e[ind_f+nbF*ind_t]  << endl;
//                }
//        }

        // REcord quota trades
         PROTECT(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota_dyn)));
         r_QuotaTrade_fe= REAL(QuotaTrade_fe);

         for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                r_QuotaTrade_fe[ind_f+ nbF*ind_t] = r_TACbyF_e[ind_f+ nbF*ind_t] - r_TACbyF_ex[ind_f+ nbF*ind_t];
         }
         //PrintValue(QuotaTrade_fe);

        UNPROTECT(5);
    }
    UNPROTECT(2);
    } else{ // New

    SEXP demandQ_f_e, offerQ_f_e, L_ef,TACbyF_e, QuotaTrade_fe,Qholdings_e ;
    PROTECT(demandQ_f_e = NEW_NUMERIC(nbF+1));
    PROTECT(offerQ_f_e = NEW_NUMERIC(nbF+1));
    double *r_L_ef;
    double *r_demandQ_f_e = REAL(demandQ_f_e);
    double *r_offerQ_f_e = REAL(offerQ_f_e);
    double *r_TACbyF_e, *r_Qholdings;
    double *r_QuotaTrade_fe;
    int ind_buyer, ind_seller;

    for (int int_eQuota_dyn = 0 ; int_eQuota_dyn  < nbEQuotaMarket_dyn ; int_eQuota_dyn++) {
        PROTECT(nam_eQuota_dyn=STRING_ELT(sppListQM_dyn,int_eQuota_dyn));
        if (!isNull(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)))){ // espece dyn
                                    PROTECT(L_ef = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota_dyn)),dimCstF));
                } else{ // espece stat
                                    PROTECT(L_ef = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota_dyn)),dimCstF));
                                }
        r_L_ef = REAL(L_ef);
        PROTECT(TACbyF_e = getListElement(TACbyF, CHAR(nam_eQuota_dyn)));
        r_TACbyF_e = REAL(TACbyF_e);
        PROTECT(Qholdings_e = getListElement(Qholdings, CHAR(nam_eQuota_dyn)));
        r_Qholdings = REAL(Qholdings_e);

//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Expected catches:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_L_ef[ind_f+nbF*ind_t]  << endl;
//                }


        // Calculate net demand for quota by fleet = expected catches - holdings
        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
            r_demandQ_f_e[ind_f] = r_L_ef[ind_f+nbF*ind_t] - r_Qholdings[ind_f+(nbF+1)*ind_t] ;
            r_offerQ_f_e[ind_f] = - r_demandQ_f_e[ind_f] ;
        }

        // And for external quota holders
        r_demandQ_f_e[nbF] = 0.0;
        r_offerQ_f_e[nbF] = r_Qholdings[nbF+(nbF+1)*ind_t];

        //Trade quotas
        int rank_buyer = 0;
        ind_buyer = r_rank_Prof_f[rank_buyer];
        int rank_seller = nbF; // start with external holders
        ind_seller = nbF;
        double trade;



                while ((rank_buyer < nbF) &&  (rank_seller >=0)){
                        //fichier << "rank buyer = " << rank_buyer  << "; ind buyer = " << ind_buyer << endl;
                    while ((rank_seller >=0) && (r_demandQ_f_e[ind_buyer] > 0)){
                        //fichier << "rank seller = " << rank_seller << "; ind seller = " << ind_seller << endl;

                        if(r_offerQ_f_e[ind_seller] > 0){
                            trade = fmin2(r_demandQ_f_e[ind_buyer],r_offerQ_f_e[ind_seller]);
                            //fichier << "demand = " << r_demandQ_f_e[ind_buyer] << ", offer = " << r_offerQ_f_e[ind_seller] << ", trade = " << trade << endl;
                            r_demandQ_f_e[ind_buyer] = r_demandQ_f_e[ind_buyer] - trade;
                            r_offerQ_f_e[ind_seller] = r_offerQ_f_e[ind_seller] - trade;

                            r_TACbyF_e[ind_buyer+nbF*ind_t] = r_TACbyF_e[ind_buyer+nbF*ind_t] + trade;

                            if(ind_seller<nbF)// Only for vessels, not external holders
                            r_TACbyF_e[ind_seller+nbF*ind_t] = r_TACbyF_e[ind_seller+nbF*ind_t] - trade;
                        } else {
                            rank_seller --;
                            ind_seller = r_rank_Prof_f[rank_seller];
                        }
                    }
                    rank_buyer ++;
                    ind_buyer = r_rank_Prof_f[rank_buyer];

                }

//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Quota after trade:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_TACbyF_e[ind_f+nbF*ind_t]  << endl;
//                }


        // REcord quota trades
         PROTECT(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota_dyn)));
         r_QuotaTrade_fe= REAL(QuotaTrade_fe);

         for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
                r_QuotaTrade_fe[ind_f+ nbF*ind_t] = r_TACbyF_e[ind_f+ nbF*ind_t] - r_Qholdings[ind_f+ (nbF+1)*ind_t];
         }


//                fichier << "sp = " << CHAR(nam_eQuota_dyn) <<", Traded quota:" << endl;
//                for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
//                    fichier << ind_f << ";" << r_QuotaTrade_fe[ind_f+nbF*ind_t]  << endl;
//                }


        UNPROTECT(5);
    }
    UNPROTECT(2);
    }
//    fichier << "QM7.5"  << endl;



    UNPROTECT(31);
//  fichier.close();
}
}


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
	if (!m) //Rprintf("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) //Rprintf("allocation failure 2 in matrix()");
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
            } else if ((y[i] > y[inhi]) && (i != ihi)) inhi=i;
        }
        rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);//Rprintf("ihi = %i ilo = %i rtol = %f\n",ihi,ilo,rtol);
        //Rprintf("rtol %f \n",rtol);
        if (rtol < ftol) {
            SWAP(y[1],y[ilo])
            for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
            break;
        }
        if (*nfunk >= NMAX) //Rprintf("NMAX exceeded : rtol = %f\n",rtol); //break;}//'//Rprintf' remplace 'nrerror'
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




extern "C" {

int BioEcoPar::EstimationTACfromF(int ind_t)
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
//ofstream fichier("C:\\Users\\BRI281\\Dropbox\\These\\IAM_Dvt\\test.EstimationTAC.txt", ios::out | ios::trunc);
//fichier << "D�but" << endl;


    if (ind_t<delay) {

    } else {
//Rprintf("A1\n");
    IND_T = ind_t;
//spQ = spp;

    //double *totFM, *totFM2, *totF, *totF2, *totFF, *totFF2, *tot, *totMod, *totMod2;

    SEXP listTempP, nDim;

    PROTECT(nDim = allocVector(INTSXP,4));
    int *nd = INTEGER(nDim); nd[0] = 0;  nd[1] = 0; nd[2] = 0; nd[3] = nbT;

    PROTECT(listTempP = duplicate(list));
    PROTECT(eVarCopy = duplicate(eVar));


    double *g_effort1FM = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f_m"));
    double *g_effort1F = REAL(getListElement(getListElement(listTempP, "Fleet"), "effort1_f"));
    double *g_nbTripFM = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f_m"));
    double *g_nbTripF = REAL(getListElement(getListElement(listTempP, "Fleet"), "nbTrip_f"));

//Rprintf("A2\n");
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



//for (int ind_f = 0 ; ind_f<nbF ; ind_f++)
//for (int ind_m = 0 ; ind_m<nbMe ; ind_m++) {
//    if ((ind_t==1) & (ind_f==0)) {
//
//        std::stringstream ggg2;
//        ggg2 << g_effort1FM[ind_f + nbF*ind_m];
//
//        fichier << "effort_step2T1" << ggg2.str() << endl;
//
//    }
//}

//if (false) PrintValue(inpFtarg);
//if (false) {


int nbEtarg = length(getAttrib(inpFtarg, R_NamesSymbol));

int denom=0, denom2=0;

double *r_N_eit_S1M1=&NA_REAL, *r_N_eit_S2M2=&NA_REAL, *r_N_eit_S3M3=&NA_REAL, *r_N_eit_S4M4=&NA_REAL, *rans_N_eit=&NA_REAL, *r_N_eit_G1=&NA_REAL, *r_N_eit_G2=&NA_REAL,
       *r_N_e0t_S1M1=&NA_REAL, *r_N_e0t_S2M2=&NA_REAL, *r_N_e0t_S3M3=&NA_REAL, *r_N_e0t_S4M4=&NA_REAL, *r_N_e0t=&NA_REAL, *r_N_e0t_G1=&NA_REAL, *r_N_e0t_G2=&NA_REAL, *recValues=&NA_REAL, *LTOT=&NA_REAL,
       *TAC_byFleet=&NA_REAL, *TAC_glob=&NA_REAL, *r_W_Ftarg=&NA_REAL, *r_Qholdings;

double *Fothi2=&NA_REAL, *Fothi2_G1=&NA_REAL, *Fothi2_G2=&NA_REAL;

double *Fothi2_S1M1=&NA_REAL,*Fothi2_S1M2=&NA_REAL,*Fothi2_S1M3=&NA_REAL,*Fothi2_S1M4=&NA_REAL,
       *Fothi2_S2M1=&NA_REAL,*Fothi2_S2M2=&NA_REAL,*Fothi2_S2M3=&NA_REAL,*Fothi2_S2M4=&NA_REAL,
       *Fothi2_S3M1=&NA_REAL,*Fothi2_S3M2=&NA_REAL,*Fothi2_S3M3=&NA_REAL,*Fothi2_S3M4=&NA_REAL,
       *Fothi2_S4M1=&NA_REAL,*Fothi2_S4M2=&NA_REAL,*Fothi2_S4M3=&NA_REAL,*Fothi2_S4M4=&NA_REAL,
       *FRWTothi2_S1M1=&NA_REAL,*FRWTothi2_S1M2=&NA_REAL,*FRWTothi2_S1M3=&NA_REAL,*FRWTothi2_S1M4=&NA_REAL,
       *FRWTothi2_S2M1=&NA_REAL,*FRWTothi2_S2M2=&NA_REAL,*FRWTothi2_S2M3=&NA_REAL,*FRWTothi2_S2M4=&NA_REAL,
       *FRWTothi2_S3M1=&NA_REAL,*FRWTothi2_S3M2=&NA_REAL,*FRWTothi2_S3M3=&NA_REAL,*FRWTothi2_S3M4=&NA_REAL,
       *FRWTothi2_S4M1=&NA_REAL,*FRWTothi2_S4M2=&NA_REAL,*FRWTothi2_S4M3=&NA_REAL,*FRWTothi2_S4M4=&NA_REAL,
       *FDWTothi2_S1M1=&NA_REAL,*FDWTothi2_S1M2=&NA_REAL,*FDWTothi2_S1M3=&NA_REAL,*FDWTothi2_S1M4=&NA_REAL,
       *FDWTothi2_S2M1=&NA_REAL,*FDWTothi2_S2M2=&NA_REAL,*FDWTothi2_S2M3=&NA_REAL,*FDWTothi2_S2M4=&NA_REAL,
       *FDWTothi2_S3M1=&NA_REAL,*FDWTothi2_S3M2=&NA_REAL,*FDWTothi2_S3M3=&NA_REAL,*FDWTothi2_S3M4=&NA_REAL,
       *FDWTothi2_S4M1=&NA_REAL,*FDWTothi2_S4M2=&NA_REAL,*FDWTothi2_S4M3=&NA_REAL,*FDWTothi2_S4M4=&NA_REAL;

double newRec=0.0, newRec_Q1=0.0, newRec_Q2=0.0, newRec_Q3=0.0, newRec_Q4=0.0, newRec_G1=0.0, newRec_G2=0.0;
//Rprintf("%i",nbEtarg);


for (int intEspTarg = 0 ; intEspTarg < nbEtarg ; intEspTarg++) {
//Rprintf("A3\n");
    SEXP namVarTarg, elmt, v_N_e0t, v_N_e0t_S1M1, v_N_e0t_S2M2, v_N_e0t_S3M3, v_N_e0t_S4M4, v_N_e0t_G1, v_N_e0t_G2, v_MeanRec_Ftarg, v_W_Ftarg, v_out_L_eit;
    PROTECT(namVarTarg=STRING_ELT(getAttrib(inpFtarg, R_NamesSymbol),intEspTarg));

    //calcul du ratio Ftarg/Fbar
    double r_Ftarg = REAL(getListElement(inpFtarg, CHAR(namVarTarg)))[IND_T];
    double r_Fbar_prev = REAL(getListElement(out_Fbar_et, CHAR(namVarTarg)))[IND_T-1];
    double r_Fbar_init = REAL(getListElement(out_Fbar_et, CHAR(namVarTarg)))[0];
    //Rprintf("Ftarg: %f, Fbar_init: %f, Fbar_prev: %f \n",r_Ftarg, r_Fbar_init,r_Fbar_prev); fichier << "Ftarg:" << r_Ftarg << ", Fbar_prev:" << r_Fbar_prev << ", Fbar_init:" << r_Fbar_init << endl;

    PROTECT(elmt = getListElement(listTempP, CHAR(namVarTarg)));
    int nbI = length(getListElement(elmt, "modI"));

//Rprintf("A4\n");
    //correction des efforts par le ratio pr�c�dent
    for (int indF = 0 ; indF < nbF ; indF++) {

        for (int indM = 0 ; indM<nbMe ; indM++) {

            g_effort1FM[indF + nbF*indM] = g_effort1FM[indF + nbF*indM] * r_Ftarg / r_Fbar_prev;
            g_nbTripFM[indF + nbF*indM] = g_nbTripFM[indF + nbF*indM] * r_Ftarg / r_Fbar_prev;

        }

        //fichier << "Avant T" << ind_t << ": Effort1 Fleet " << indF << " = " << g_effort1F[indF] << "; Nbtrip = " << g_nbTripF[indF]<< endl;

        g_effort1F[indF] = g_effort1F[indF] * r_Ftarg / r_Fbar_prev;
        g_nbTripF[indF] = g_nbTripF[indF] * r_Ftarg / r_Fbar_prev;

        //fichier << "Apres T" << ind_t << ": Effort1 Fleet " << indF << " = " << g_effort1F[indF] << "; Nbtrip = " << g_nbTripF[indF]<< endl;

    }


    //et correction des mortalit�s autres pour les esp�ces dynamiques XSA, Spict et SS3


            int nbi = length(getListElement(getListElement(list, CHAR(namVarTarg)), "modI"));

            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global

                    // Dans eVarCopy
                    Fothi2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44)); //Rprintf("Dans EVARcopy (l.14478), Fothi2 = "); PrintValue(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi age " << ag <<"=" << Fothi2[ag + IND_T*nbi] << endl;
                            Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi age" << ag << "=" << Fothi2[ag + IND_T*nbi] << endl;
                    }

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 44)); //Rprintf("Dans EVARcopy (l.14478), Fothi2 = "); PrintValue(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 44));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi age " << ag <<"=" << Fothi2[ag + IND_T*nbi] << endl;
                            Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi age" << ag << "=" << Fothi2[ag + IND_T*nbi] << endl;
                    }

            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)) {//Age and sex-based

                    // Dans EvarCopy
                    Fothi2_G1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 224));
                    Fothi2_G2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 225));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi_G1 age " << ag <<"=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                            Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi_G1 age" << ag << "=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                    }

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2_G1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 224));
                    Fothi2_G2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 225));
                    for (int ag = 0; ag < nbi; ag++) {
                            //fichier << "Avant T" << ind_t << "; Fothi_G1 age " << ag <<"=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                            Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                            //fichier << "Apres T" << ind_t << "; Fothi_G1 age" << ag << "=" << Fothi2_G1[ag + IND_T*nbi] << "/ Fothi_G2 age " << ag <<"=" << Fothi2_G2[ag + IND_T*nbi] << endl;
                    }
            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly

                    // Dans eVarCopy
                    Fothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 116));
                    Fothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 117));
                    Fothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 118));
                    Fothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 119));
                    Fothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 120));
                    Fothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 121));
                    Fothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 122));
                    Fothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 123));
                    Fothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 124));
                    Fothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 125));
                    Fothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 126));
                    Fothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 127));
                    Fothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 128));
                    Fothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 129));
                    Fothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 130));
                    Fothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 131));

                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FRWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 176));
                    FRWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 177));
                    FRWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 178));
                    FRWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 179));
                    FRWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 180));
                    FRWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 181));
                    FRWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 182));
                    FRWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 183));
                    FRWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 184));
                    FRWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 185));
                    FRWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 186));
                    FRWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 187));
                    FRWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 188));
                    FRWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 189));
                    FRWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 190));
                    FRWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 191));

                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FDWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 208));
                    FDWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 209));
                    FDWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 210));
                    FDWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 211));
                    FDWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 212));
                    FDWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 213));
                    FDWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 214));
                    FDWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 215));
                    FDWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 216));
                    FDWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 217));
                    FDWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 218));
                    FDWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 219));
                    FDWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 220));
                    FDWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 221));
                    FDWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 222));
                    FDWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVarCopy, CHAR(namVarTarg)), 223));

                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    // Dans eVar pour usage hors de cette fonction
                    Fothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 116));
                    Fothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 117));
                    Fothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 118));
                    Fothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 119));
                    Fothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 120));
                    Fothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 121));
                    Fothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 122));
                    Fothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 123));
                    Fothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 124));
                    Fothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 125));
                    Fothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 126));
                    Fothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 127));
                    Fothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 128));
                    Fothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 129));
                    Fothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 130));
                    Fothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 131));

                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FRWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 176));
                    FRWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 177));
                    FRWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 178));
                    FRWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 179));
                    FRWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 180));
                    FRWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 181));
                    FRWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 182));
                    FRWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 183));
                    FRWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 184));
                    FRWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 185));
                    FRWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 186));
                    FRWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 187));
                    FRWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 188));
                    FRWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 189));
                    FRWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 190));
                    FRWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 191));

                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;

                    FDWTothi2_S1M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 208));
                    FDWTothi2_S1M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 209));
                    FDWTothi2_S1M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 210));
                    FDWTothi2_S1M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 211));
                    FDWTothi2_S2M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 212));
                    FDWTothi2_S2M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 213));
                    FDWTothi2_S2M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 214));
                    FDWTothi2_S2M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 215));
                    FDWTothi2_S3M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 216));
                    FDWTothi2_S3M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 217));
                    FDWTothi2_S3M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 218));
                    FDWTothi2_S3M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 219));
                    FDWTothi2_S4M1 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 220));
                    FDWTothi2_S4M2 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 221));
                    FDWTothi2_S4M3 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 222));
                    FDWTothi2_S4M4 = REAL(VECTOR_ELT(getListElement(eVar, CHAR(namVarTarg)), 223));

                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Ftarg / r_Fbar_init;
                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Ftarg / r_Fbar_init;


            }

    // envoi du module Mortalit�
    //Rprintf("call.Mortalite\n");
   // fichier << "call.Mortalite" << endl;
    Mortalite(listTempP, IND_T, eVarCopy);
    //Rprintf("End.Mortalite\n");
    //fichier << "End.Mortalite" << endl;

    // calcul recrutement
    if ((nbI>1) && !isNull(inpMeanRec_Ftarg) && (!isNull(getListElement(inpMeanRec_Ftarg,CHAR(namVarTarg))))){ // si MeanRec_Ftarg renseigne: forcage selon type 1 (Moyenne sur X dernieres annees) ou 2 (For�age avec valeurs renseignees)

        //fichier << "Recruitment in HCR = from MeanRecFtarg" << endl;

            PROTECT(v_MeanRec_Ftarg = getListElement(inpMeanRec_Ftarg, CHAR(namVarTarg)));
            //PrintValue(v_MeanRec_Ftarg);
        //Rprintf("A5\n");
            if (length(v_MeanRec_Ftarg)==1) {

                denom = INTEGER(v_MeanRec_Ftarg)[0];
                denom2 = denom;

                if (!ISNA(denom)){
                    if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global
        //Rprintf("A6\n");
                      newRec = 0.0;
                      rans_N_eit = REAL(getListElement(out_N_eit,CHAR(namVarTarg)));
                      for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec = newRec + rans_N_eit[(IND_T-index)*nbI] ;
                                                                 }}
                      newRec = newRec / denom2;
                      PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                      r_N_e0t = REAL(v_N_e0t);
                      r_N_e0t[IND_T] = newRec;
                      rans_N_eit[IND_T*nbI] = newRec;
                      UNPROTECT(1);
        //Rprintf("A7\n");
        //fichier << "Type 1 (moyenne) MeanRec: " << rans_N_eit[IND_T*nbI] << endl;

                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)) {//Ageand sex-based
                        r_N_eit_G1 = REAL(getListElement(out_N_eit_G1,CHAR(namVarTarg)));
                        r_N_eit_G2 = REAL(getListElement(out_N_eit_G2,CHAR(namVarTarg)));
                        newRec_G1 = 0.0; newRec_G2 = 0.0;

                        for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec_G1 = newRec_G1 + r_N_eit_G1[(IND_T-index)*nbI] ;
                                                                   newRec_G2 = newRec_G2 + r_N_eit_G2[(IND_T-index)*nbI] ;
                                                                 }}
                        newRec_G1 = newRec_G1 / denom2 ;
                        newRec_G2 = newRec_G2 / denom2 ;

                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2"));

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = newRec_G1;
                        r_N_e0t_G2[IND_T] = newRec_G2;

                        UNPROTECT(2);
       // fichier << "Type 1 (moyenne) MeanRecG1: " << r_N_e0t_G1[IND_T] << "/ MeanRecG2: " << r_N_e0t_G2[IND_T]<< endl;
                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly
       // Rprintf("A8\n");
                      //Quarter 1
                        r_N_eit_S1M1 = REAL(getListElement(out_N_eit_S1M1,CHAR(namVarTarg)));
                        r_N_eit_S2M2 = REAL(getListElement(out_N_eit_S2M2,CHAR(namVarTarg)));
                        r_N_eit_S3M3 = REAL(getListElement(out_N_eit_S3M3,CHAR(namVarTarg)));
                        r_N_eit_S4M4 = REAL(getListElement(out_N_eit_S4M4,CHAR(namVarTarg)));
                        newRec_Q1 = 0.0; newRec_Q2 = 0.0; newRec_Q3 = 0.0; newRec_Q4 = 0.0;
        //Rprintf("A9\n");
                        for (int index=1; index<=denom; index++) {if ((IND_T-index)<0) {
                                                                   denom2 = denom2 - 1 ;
                                                                 } else {
                                                                   newRec_Q1 = newRec_Q1 + r_N_eit_S1M1[(IND_T-index)*nbI] ;
                                                                   newRec_Q2 = newRec_Q2 + r_N_eit_S2M2[(IND_T-index)*nbI] ;
                                                                   newRec_Q3 = newRec_Q3 + r_N_eit_S3M3[(IND_T-index)*nbI] ;
                                                                   newRec_Q4 = newRec_Q4 + r_N_eit_S4M4[(IND_T-index)*nbI] ;
                                                                 }}
                        newRec_Q1 = newRec_Q1 / denom2 ;
                        newRec_Q2 = newRec_Q2 / denom2 ;
                        newRec_Q3 = newRec_Q3 / denom2 ;
                        newRec_Q4 = newRec_Q4 / denom2 ;
        //Rprintf("A10\n");
                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = newRec_Q1;
                        r_N_e0t_S2M2[0] = newRec_Q2;
                        r_N_e0t_S3M3[0] = newRec_Q3;
                        r_N_e0t_S4M4[0] = newRec_Q4;

                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...
                        UNPROTECT(4);

                    }

                }
            } else {

              if (length(v_MeanRec_Ftarg)>1) {  //historique XSA ou SS3
        //Rprintf("A11\n");
                recValues = REAL(v_MeanRec_Ftarg);
       // Rprintf("A110\n");
                if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {//Age-based + global
        //Rprintf("A12\n");
                      PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                      r_N_e0t = REAL(v_N_e0t);
                      r_N_e0t[IND_T] = recValues[IND_T];
                      UNPROTECT(1);
        //Rprintf("A13\n");
        //fichier << "Type 2 (historique) MeanRec: " << rans_N_eit[IND_T*nbI] << endl;

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){//age and sex-based
        //Rprintf("A14\n");
                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1")); //PrintValue(v_N_e0t_G1);
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2")); //PrintValue(v_N_e0t_G2);

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = recValues[2*IND_T];
                        r_N_e0t_G2[IND_T] = recValues[2*IND_T + 1];
       // Rprintf("A15\n");
       // fichier << "Type 2 (historique) MeanRecG1: " << r_N_e0t_G1[IND_T] << "/ MeanRecG2: " << r_N_e0t_G2[IND_T]<< endl;

                       UNPROTECT(2);

                    } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){//Quarterly
     //   Rprintf("A14\n");
                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = recValues[4*IND_T];
                        r_N_e0t_S2M2[0] = recValues[4*IND_T + 1];
                        r_N_e0t_S3M3[0] = recValues[4*IND_T + 2];
                        r_N_e0t_S4M4[0] = recValues[4*IND_T + 3];
        //Rprintf("A15\n");
                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...

                       UNPROTECT(4);

                    }

              }

            }

            UNPROTECT(1);

    } else if ((nbI>1) && (!isNull(getListElement(recParamList,CHAR(namVarTarg)))) && (IND_T>0)) { // recrutement = recrutement attendu avec la relation SR (sans incertitude)
//fichier << "Recruitment in HCR = expected from SR" << endl;

            double  *rans_SSB_et = REAL(getListElement(out_SSB_et,CHAR(namVarTarg)));

            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
                double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type")); //Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));

                if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                    double recr = 0.0;

                    if (typeSR[IND_T]==1){ // Hockey Stick
                        if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T - 1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                            recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                        } else{
                            recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                        }
                    } else if (typeSR[IND_T]==2){ // Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                        recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T -1 - del]) /
                        (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T -1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                        param[IND_T + 3*nbT] ;
                    }

                    PROTECT(v_N_e0t = getListElement(elmt, "N_i0t"));
                    r_N_e0t = REAL(v_N_e0t);
                    r_N_e0t[IND_T] = recr; //fichier << "e = " << CHAR(namVarTarg) << " recr = " << recr << ", N_i0t = " << r_N_e0t[IND_T] << endl;
                    UNPROTECT(1);
                }

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
                    double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                    int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));// Rprintf("type = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                    int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));
                    double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil")); //Rprintf("ventil = "); //PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil"));

                    if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                        double recr = 0.0;

                        if (typeSR[IND_T]==1){ // Hockey Stick
                            if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T - 1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                                recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                            } else{
                                recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                            }
                        } else if (typeSR[IND_T]==2){ //Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                            recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T -1 - del]) /
                            (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T -1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                            param[IND_T + 3*nbT];
                        }

                        PROTECT(v_N_e0t_S1M1 = getListElement(elmt, "Ni0_S1M1"));
                        PROTECT(v_N_e0t_S2M2 = getListElement(elmt, "Ni0_S2M2"));
                        PROTECT(v_N_e0t_S3M3 = getListElement(elmt, "Ni0_S3M3"));
                        PROTECT(v_N_e0t_S4M4 = getListElement(elmt, "Ni0_S4M4"));

                        r_N_e0t_S1M1 = REAL(v_N_e0t_S1M1);
                        r_N_e0t_S2M2 = REAL(v_N_e0t_S2M2);
                        r_N_e0t_S3M3 = REAL(v_N_e0t_S3M3);
                        r_N_e0t_S4M4 = REAL(v_N_e0t_S4M4);

                        r_N_e0t_S1M1[0] = recr*ventil[0];
                        r_N_e0t_S2M2[0] =  recr*ventil[1];
                        r_N_e0t_S3M3[0] =  recr*ventil[2];
                        r_N_e0t_S4M4[0] =  recr*ventil[3];
                       // ce serait bien de mettre aussi � jour "N0t_S1M1[0]",...

                       UNPROTECT(4);

                    }

                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){
                    double *param = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param")); //Rprintf("param = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"param"));
                    int *typeSR = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));// Rprintf("type = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"type"));
                    int del = INTEGER(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"))[0]; //Rprintf("delay = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"delay"));
                    double *ventil = REAL(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil")); //Rprintf("ventil = "); PrintValue(getListElement(getListElement(recParamList,CHAR(namVarTarg)),"ventil"));

                    if ((!ISNA(param[IND_T])) & (IND_T>=del)) {
                        double recr = 0.0;
                        //fichier << "e = " << CHAR(namVarTarg) << " SSB = " << rans_SSB_et[IND_T - 1 - del] << endl;

                        if (typeSR[IND_T]==1){ // Hockey Stick
                            if ((1/param[IND_T + 1*nbT])>rans_SSB_et[IND_T -1 - del]) { // on prend SSB [t-1] comme approximation de SSB [t] pour calcul recrutement
                                recr = param[IND_T + 0*nbT] * rans_SSB_et[IND_T - 1 - del] * param[IND_T + 2*nbT];
                            } else{
                                recr = param[IND_T + 0*nbT] * param[IND_T + 2*nbT] / param[IND_T + 1*nbT];
                            }
                        } else if (typeSR[IND_T]==2){ //Beverton-Holt (param = [h,R0,B0,exp(RecDev in log scale),delta = exp(norm(0,sigma^2))])
                            //fichier << " Beverton Holt" << endl;
                            recr = (4*param[IND_T + 0*nbT] * param[IND_T + 1*nbT] * rans_SSB_et[IND_T - 1 - del]) /
                            (param[IND_T + 2*nbT]*(1-param[IND_T + 0*nbT]) + rans_SSB_et[IND_T - 1 - del]*(5*param[IND_T + 0*nbT]-1)) *
                            param[IND_T + 3*nbT];
                        }

                        PROTECT(v_N_e0t_G1 = getListElement(elmt, "N_i0t_G1"));  //PrintValue(v_N_e0t_G1);
                        PROTECT(v_N_e0t_G2 = getListElement(elmt, "N_i0t_G2")); //PrintValue(v_N_e0t_G2);

                        r_N_e0t_G1 = REAL(v_N_e0t_G1);
                        r_N_e0t_G2 = REAL(v_N_e0t_G2);

                        r_N_e0t_G1[IND_T] = recr*ventil[0]; //fichier << "e = " << CHAR(namVarTarg) << " recr = " << recr << ", N_i0t_G1 = " << r_N_e0t_G1[IND_T] << endl;
                        r_N_e0t_G2[IND_T] = recr*ventil[1];

                         UNPROTECT(2);

                    }
                }

        }



//Rprintf("A16\n");
 //   Rprintf("call.DynamicPop\n");
 //   fichier << "call.DynamicPop" << endl;
    DynamicPop(listTempP, IND_T, eVarCopy, false);
 //   Rprintf("end.DynamicPop\n");
 //   fichier << "end.DynamicPop" << endl;
//Rprintf("A17\n");
  //  Rprintf("call.CatchDL\n");
  //  fichier << "call.CatchDL" << endl;
    CatchDL(listTempP, IND_T, eVarCopy);
 //   Rprintf("end.CatchDL\n");
 //   fichier << "end.CatchDL" << endl;
//Rprintf("A18\n");
    //on peut d�sormais d�duire des d�barquements mod�lis�s les TAC par flottille et totaux

    PROTECT(v_W_Ftarg = getListElement(inpW_Ftarg, CHAR(namVarTarg)));
//Rprintf("A19\n");
    PROTECT(v_out_L_eit = getListElement(out_L_eit, CHAR(namVarTarg)));
//Rprintf("A20\n");
    LTOT = REAL(aggregObj(v_out_L_eit,nDim));
//Rprintf("A21\n");
    //if (IND_T==2) {Rprintf("AA\n"); PrintValue(TACbyF);PrintValue(TAC);}

    TAC_byFleet = REAL(getListElement(TACbyF, CHAR(namVarTarg)));
    TAC_glob = REAL(getListElement(TAC, CHAR(namVarTarg)));
    r_W_Ftarg = REAL(v_W_Ftarg);
    r_Qholdings = REAL(getListElement(Qholdings, CHAR(namVarTarg)));
//Rprintf("A22\n");

    TAC_glob[IND_T] = LTOT[IND_T];
    //fichier << "Ltot t-1: " << LTOT[IND_T-1] << endl;
    //fichier << "Tactot t: " << TAC_glob[IND_T] << endl;

    for (int indF = 0 ; indF < nbF ; indF++) TAC_byFleet[indF + nbF*IND_T] = r_W_Ftarg[indF + (nbF+1)*IND_T] * LTOT[IND_T]; // for use in Gestion F2: only modelled fleets

    for (int indF = 0 ; indF <= nbF ; indF++) r_Qholdings[indF + (nbF+1)*IND_T] = r_W_Ftarg[indF + (nbF+1)*IND_T] * LTOT[IND_T]; // for use in quota trading, contains also external investors
    //PrintValue(getListElement(Qholdings, CHAR(namVarTarg)));
        //re-correction des efforts par l'inverse du ratio pr�c�dent
    for (int indF = 0 ; indF < nbF ; indF++) {

        for (int indM = 0 ; indM<nbMe ; indM++) {

            g_effort1FM[indF + nbF*indM] = g_effort1FM[indF + nbF*indM] * r_Fbar_prev/ r_Ftarg;
            g_nbTripFM[indF + nbF*indM] = g_nbTripFM[indF + nbF*indM] * r_Fbar_prev / r_Ftarg;

        }

        g_effort1F[indF] = g_effort1F[indF] * r_Fbar_prev / r_Ftarg;
        g_nbTripF[indF] = g_nbTripF[indF] * r_Fbar_prev / r_Ftarg;

    }

        //et re-correction des mortalit�s autres pour les esp�ces dynamiques XSA, Spict et SS3

//            if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)) {
//
//                    for (int ag = 0; ag < nbi; ag++) Fothi2[ag + IND_T*nbi] = Fothi2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//            } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==0) & (Svec[getListIndex(S, CHAR(namVarTarg))]==1)){
//                for (int ag = 0; ag < nbi; ag++) {
//                        Fothi2_G1[ag + IND_T*nbi] = Fothi2_G1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                        Fothi2_G2[ag + IND_T*nbi] = Fothi2_G2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                }
//
//                } else if ((Qvec[getListIndex(Q, CHAR(namVarTarg))]==1) & (Svec[getListIndex(S, CHAR(namVarTarg))]==0)){
//
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M1[ag + IND_T*nbi] = Fothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M2[ag + IND_T*nbi] = Fothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M3[ag + IND_T*nbi] = Fothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S1M4[ag + IND_T*nbi] = Fothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M1[ag + IND_T*nbi] = Fothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M2[ag + IND_T*nbi] = Fothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M3[ag + IND_T*nbi] = Fothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S2M4[ag + IND_T*nbi] = Fothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M1[ag + IND_T*nbi] = Fothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M2[ag + IND_T*nbi] = Fothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M3[ag + IND_T*nbi] = Fothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S3M4[ag + IND_T*nbi] = Fothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M1[ag + IND_T*nbi] = Fothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M2[ag + IND_T*nbi] = Fothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M3[ag + IND_T*nbi] = Fothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) Fothi2_S4M4[ag + IND_T*nbi] = Fothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M1[ag + IND_T*nbi] = FRWTothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M2[ag + IND_T*nbi] = FRWTothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M3[ag + IND_T*nbi] = FRWTothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S1M4[ag + IND_T*nbi] = FRWTothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M1[ag + IND_T*nbi] = FRWTothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M2[ag + IND_T*nbi] = FRWTothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M3[ag + IND_T*nbi] = FRWTothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S2M4[ag + IND_T*nbi] = FRWTothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M1[ag + IND_T*nbi] = FRWTothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M2[ag + IND_T*nbi] = FRWTothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M3[ag + IND_T*nbi] = FRWTothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S3M4[ag + IND_T*nbi] = FRWTothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M1[ag + IND_T*nbi] = FRWTothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M2[ag + IND_T*nbi] = FRWTothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M3[ag + IND_T*nbi] = FRWTothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FRWTothi2_S4M4[ag + IND_T*nbi] = FRWTothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M1[ag + IND_T*nbi] = FDWTothi2_S1M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M2[ag + IND_T*nbi] = FDWTothi2_S1M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M3[ag + IND_T*nbi] = FDWTothi2_S1M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S1M4[ag + IND_T*nbi] = FDWTothi2_S1M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M1[ag + IND_T*nbi] = FDWTothi2_S2M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M2[ag + IND_T*nbi] = FDWTothi2_S2M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M3[ag + IND_T*nbi] = FDWTothi2_S2M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S2M4[ag + IND_T*nbi] = FDWTothi2_S2M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M1[ag + IND_T*nbi] = FDWTothi2_S3M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M2[ag + IND_T*nbi] = FDWTothi2_S3M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M3[ag + IND_T*nbi] = FDWTothi2_S3M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S3M4[ag + IND_T*nbi] = FDWTothi2_S3M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M1[ag + IND_T*nbi] = FDWTothi2_S4M1[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M2[ag + IND_T*nbi] = FDWTothi2_S4M2[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M3[ag + IND_T*nbi] = FDWTothi2_S4M3[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//                    for (int ag = 0; ag < nbi; ag++) FDWTothi2_S4M4[ag + IND_T*nbi] = FDWTothi2_S4M4[ag + 0*nbi] * r_Fbar_init / r_Ftarg;
//
//            }

     UNPROTECT(4);


}
//UNPROTECT(2);
 //}
     UNPROTECT(3);


    }

   return(0);
   //fichier.close();


  }
}




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

    if (ind_t<delay) {

    } else {

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
    CatchDL(listTempP, IND_T, eVar);
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
        CatchDL(listTempP, IND_T, eVar); //hors boucle esp�ce � optimiser
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

    }

    ZoptSS3 = false;
	return 0;

//fichier.close();
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
	//Rprintf("\nNumber of function evaluations: %3d\n",nfunc);
	//Rprintf("Vertices of final 3-d simplex and\n");
	//Rprintf("function values at the vertices:\n\n");
	//Rprintf("%3s %10s %12s %12s %14s\n\n","i","x[i]","y[i]","z[i]","function");
	//for (i=1;i<=MP;i++) {
		//Rprintf("%3d ",i);
	//	for (j=1;j<=NP;j++) Rprintf("%12.6f ",p[i][j]);
		//Rprintf("%12.6f\n",y[i]);
	//}
	//Rprintf("\nTrue minimum is at (0.5,0.6,0.7)\n");
	//free_matrix(p,1,MP,1,NP);
	//free_vector(y,1,MP);
	free_vector(x,1,NP);
	return 0;
}

}



extern "C" {

void BioEcoPar::EcoDCF(SEXP list, int ind_t, int perscCalc, double dr)
{

//ofstream fichier;
//if (ind_t ==4) fichier.open ("C:\\Users\\fbriton\\Dropbox\\These\\IAM_Dvt\\EcoDCF.txt", ios::out | ios::trunc);

//extern "C" {
//void BioEcoPar::EcoDCF(SEXP list, int ind_t, int perscCalc, double dr)
//{
//Rprintf("\nJ1\n");fichier << "J1" << endl;

    SEXP Flist;
    PROTECT(Flist = getListElement(list, "Fleet"));

    PROTECT(out_EcoDCF);

    SEXP dimCstF, DimF, dimnamesF, dimCstFM, dimCstFini, dimCstFMini, DimFM, DimFMini, dimnamesFM, dimnamesFMini; //formatage des objets r�sultats

    SEXP eFACTf, eFACTfm, elmt;

    SEXP    nbv_f, nbv_f_m, lc_f_m, lcd_f_m, tripLgth_f, tripLgth_f_m, nbTrip_f, nbTrip_f_m, nbds_f, nbds_f_m,
            effort1_f, effort1_f_m, effort2_f, effort2_f_m, Lref_f_m, cnb_f_m, ovcDCF_f_m, fc_f_m, vf_f_m, cshr_f_m, cshr_f, cnb_f, persc_f,
            eec_f, mwh_f, rep_f, gc_f, fixc_f, FTE_f, dep_f, ic_f, K_f, inv_f, FTE_f_m, GVLref_f_m, ue_f, ue_f_m;

    SEXP    dc_nbv_f, dc_nbv_f_m, dc_lc_f_m, dc_lcd_f_m, dc_tripLgth_f, dc_tripLgth_f_m, dc_nbTrip_f, dc_nbTrip_f_m, dc_nbds_f, dc_nbds_f_m,
            dc_effort1_f, dc_effort1_f_m, dc_effort2_f, dc_effort2_f_m, dc_Lref_f_m, dc_cnb_f_m, dc_ovcDCF_f_m, dc_fc_f_m, dc_vf_f_m, dc_cshr_f_m, dc_cshr_f, dc_cnb_f, dc_persc_f,
            dc_eec_f, dc_mwh_f, dc_rep_f, dc_gc_f, dc_fixc_f, dc_FTE_f, dc_dep_f, dc_ic_f, dc_K_f, dc_inv_f, dc_FTE_f_m, dc_GVLref_f_m, dc_ue_f, dc_ue_f_m;

    int *dCF,*dCFM,*dCFini,*dCFMini,*DF,*DFM, *DFMini;

    int     *dim_nbv_f, *dim_nbv_f_m, *dim_lc_f_m, *dim_tripLgth_f_m, *dim_nbTrip_f_m,
            *dim_Lref_f_m, *dim_cnb_f_m, *dim_ovcDCF_f_m, *dim_fc_f_m, *dim_vf_f_m, *dim_cshr_f_m, *dim_cshr_f, *dim_persc_f,
            *dim_eec_f, *dim_mwh_f, *dim_rep_f, *dim_gc_f, *dim_fixc_f, *dim_FTE_f, *dim_dep_f, *dim_ic_f, *dim_K_f, *dim_inv_f, *dim_GVLref_f_m,
            *dim_ue_f_m;

    double  *r_nbv_f, *r_nbv_f_m, *r_lc_f_m, *r_lcd_f_m, *r_tripLgth_f_m, *r_nbTrip_f_m,
            *r_Lref_f_m, *r_cnb_f_m, *r_ovcDCF_f_m, *r_fc_f_m, *r_vf_f_m, *r_cshr_f_m, *r_cshr_f, *r_persc_f,
            *r_eec_f, *r_mwh_f, *r_rep_f, *r_gc_f, *r_fixc_f, *r_FTE_f, *r_dep_f, *r_ic_f, *r_K_f, *r_inv_f, *r_GVLref_f_m;

    double  *r_ET_f_m_out,
            *r_GVLcom_f_m_e_out,*r_GVLst_f_m_e_out, *r_GVLtot_f_m_out, *r_GVLav_f_m_out, *r_GVLtot_f_out,
            *r_GVLav_f_out, *r_NGVLav_f_m_out, *r_NGVLav_f_out, *r_cnb_f_m_out, *r_cnb_f_out,
            *r_rtbs_f_m_out, *r_rtbs_f_out, *r_cshrT_f_m_out, *r_cshrT_f_out, *r_ncshr_f_out, *r_ocl_f_out, *r_cs_f_out, *r_csTot_f_out, *r_gva_f_out, *r_gvamargin_f_out,
            *r_gva_FTE_f_out, *r_ccw_f_out, *r_ccwCr_f_out, *r_wageg_f_out, *r_wagen_f_out, *r_wageg_FTE_f_out, *r_wageg_h_f_out, *r_gp_f_out, *r_gpmargin_f_out,
            *r_ncf_f_out, *r_np_f_out, *r_npmargin_f_out, *r_prof_f_out, *r_npmargin_trend_f_out, *r_ssTot_f_out, *r_ps_f_out, *r_sts_f_out, *r_BER_f_out, *r_CR_BER_f_out,
            *r_fuelEff_f_out, *r_ratio_fvol_gva_f_out, *r_ratio_gp_gva_f_out, *r_ratio_GVL_K_f_out, *r_ratio_gp_K_f_out, *r_RoFTA_f_out, *r_ROI_f_out,
            *r_ratio_np_K_f_out, *r_ratio_GVL_cnb_ue_f_out,
            *r_rtbsAct_f_out, *r_csAct_f_out, *r_gvaAct_f_out, *r_gpAct_f_out, *r_psAct_f_out, *r_stsAct_f_out, *r_QuotaExp_f_out;


//Rprintf("Eco 2");fichier << "Eco2" << endl;
//d�finition des dimensions


    PROTECT(dimnamesF = allocVector(VECSXP,2));
    PROTECT(dimnamesFM = allocVector(VECSXP,3));
    PROTECT(dimnamesFMini = allocVector(VECSXP,2));

    SET_VECTOR_ELT(dimnamesF, 0, fleetList); SET_VECTOR_ELT(dimnamesF, 1, times);
    SET_VECTOR_ELT(dimnamesFM, 0, fleetList); SET_VECTOR_ELT(dimnamesFM, 1, metierListEco); SET_VECTOR_ELT(dimnamesFM, 2, times);
    SET_VECTOR_ELT(dimnamesFMini, 0, fleetList); SET_VECTOR_ELT(dimnamesFMini, 1, metierListEco);


    PROTECT(dimCstF = allocVector(INTSXP, 4));
    PROTECT(dimCstFini = allocVector(INTSXP, 4));
    PROTECT(dimCstFM = allocVector(INTSXP, 4));
    PROTECT(dimCstFMini = allocVector(INTSXP, 4));


    dCF = INTEGER(dimCstF) ; dCF[0] = nbF; dCF[1] = 0; dCF[2] = 0; dCF[3] = nbT;
    dCFM = INTEGER(dimCstFM) ; dCFM[0] = nbF; dCFM[1] = nbMe; dCFM[2] = 0; dCFM[3] = nbT;
    dCFini = INTEGER(dimCstFini) ; dCFini[0] = nbF; dCFini[1] = 0; dCFini[2] = 0; dCFini[3] = 0;
    dCFMini = INTEGER(dimCstFMini) ; dCFMini[0] = nbF; dCFMini[1] = nbMe; dCFMini[2] = 0; dCFMini[3] = 0;


    PROTECT(DimF = allocVector(INTSXP, 2));
    PROTECT(DimFM = allocVector(INTSXP, 3));
    PROTECT(DimFMini = allocVector(INTSXP, 2));

    DF = INTEGER(DimF) ; DF[0] = nbF; DF[1] = nbT;
    DFM = INTEGER(DimFM) ; DFM[0] = nbF; DFM[1] = nbMe; DFM[2] = nbT;
    DFMini = INTEGER(DimFMini) ; DFMini[0] = nbF; DFMini[1] = nbMe;

    // facteurs des indices g�n�riques F/FM

    PROTECT(eFACTf = iDim(dCF));
    PROTECT(eFACTfm = iDim(dCFM));

    //Rprintf("Eco 3");fichier << "Eco3" << endl;
    // protect.root -> 14

// ---> P = 14

    int *eF_f = INTEGER(eFACTf);
    int *eF_fm = INTEGER(eFACTfm);

    PROTECT(nbv_f = getListElement(Flist, "nbv_f"));                PROTECT(dc_nbv_f = iDim(INTEGER(getAttrib(nbv_f, install("DimCst"))))); //Rprintf("Eco 31");
    PROTECT(nbv_f_m = getListElement(Flist, "nbv_f_m"));            PROTECT(dc_nbv_f_m = iDim(INTEGER(getAttrib(nbv_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(lc_f_m = getListElement(Flist, "lc_f_m"));              PROTECT(dc_lc_f_m = iDim(INTEGER(getAttrib(lc_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(lcd_f_m = getListElement(Flist, "lcd_f_m"));            PROTECT(dc_lcd_f_m = iDim(INTEGER(getAttrib(lcd_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(tripLgth_f = getListElement(Flist, "tripLgth_f"));      PROTECT(dc_tripLgth_f = iDim(INTEGER(getAttrib(tripLgth_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(tripLgth_f_m = getListElement(Flist, "tripLgth_f_m"));  PROTECT(dc_tripLgth_f_m = iDim(INTEGER(getAttrib(tripLgth_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(nbTrip_f = getListElement(Flist, "nbTrip_f"));          PROTECT(dc_nbTrip_f = iDim(INTEGER(getAttrib(nbTrip_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(nbTrip_f_m = getListElement(Flist, "nbTrip_f_m"));      PROTECT(dc_nbTrip_f_m = iDim(INTEGER(getAttrib(nbTrip_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(nbds_f = getListElement(Flist, "nbds_f"));              PROTECT(dc_nbds_f = iDim(INTEGER(getAttrib(nbds_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(nbds_f_m = getListElement(Flist, "nbds_f_m"));          PROTECT(dc_nbds_f_m = iDim(INTEGER(getAttrib(nbds_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(effort1_f = getListElement(Flist, "effort1_f"));        PROTECT(dc_effort1_f = iDim(INTEGER(getAttrib(effort1_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(effort1_f_m = getListElement(Flist, "effort1_f_m"));    PROTECT(dc_effort1_f_m = iDim(INTEGER(getAttrib(effort1_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(effort2_f = getListElement(Flist, "effort2_f"));        PROTECT(dc_effort2_f = iDim(INTEGER(getAttrib(effort2_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(effort2_f_m = getListElement(Flist, "effort2_f_m"));    PROTECT(dc_effort2_f_m = iDim(INTEGER(getAttrib(effort2_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(Lref_f_m = getListElement(Flist, "Lref_f_m"));          PROTECT(dc_Lref_f_m = iDim(INTEGER(getAttrib(Lref_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(cnb_f_m = getListElement(Flist, "cnb_f_m"));            PROTECT(dc_cnb_f_m = iDim(INTEGER(getAttrib(cnb_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(ovcDCF_f_m = getListElement(Flist, "ovcDCF_f_m"));      PROTECT(dc_ovcDCF_f_m = iDim(INTEGER(getAttrib(ovcDCF_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(fc_f_m = getListElement(Flist, "fc_f_m"));              PROTECT(dc_fc_f_m = iDim(INTEGER(getAttrib(fc_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(vf_f_m = getListElement(Flist, "vf_f_m"));              PROTECT(dc_vf_f_m = iDim(INTEGER(getAttrib(vf_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(cshr_f_m = getListElement(Flist, "cshr_f_m"));          PROTECT(dc_cshr_f_m = iDim(INTEGER(getAttrib(cshr_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(cshr_f = getListElement(Flist, "cshr_f"));              PROTECT(dc_cshr_f = iDim(INTEGER(getAttrib(cshr_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(cnb_f = getListElement(Flist, "cnb_f"));                PROTECT(dc_cnb_f = iDim(INTEGER(getAttrib(cnb_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(persc_f = getListElement(Flist, "persc_f"));            PROTECT(dc_persc_f = iDim(INTEGER(getAttrib(persc_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(eec_f = getListElement(Flist, "eec_f"));                PROTECT(dc_eec_f = iDim(INTEGER(getAttrib(eec_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(mwh_f = getListElement(Flist, "mwh_f"));                PROTECT(dc_mwh_f = iDim(INTEGER(getAttrib(mwh_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(rep_f = getListElement(Flist, "rep_f"));                PROTECT(dc_rep_f = iDim(INTEGER(getAttrib(rep_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(gc_f = getListElement(Flist, "gc_f"));                  PROTECT(dc_gc_f = iDim(INTEGER(getAttrib(gc_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(fixc_f = getListElement(Flist, "fixc_f"));              PROTECT(dc_fixc_f = iDim(INTEGER(getAttrib(fixc_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(FTE_f = getListElement(Flist, "FTE_f"));                PROTECT(dc_FTE_f = iDim(INTEGER(getAttrib(FTE_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(dep_f = getListElement(Flist, "dep_f"));                PROTECT(dc_dep_f = iDim(INTEGER(getAttrib(dep_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(ic_f = getListElement(Flist, "ic_f"));                  PROTECT(dc_ic_f = iDim(INTEGER(getAttrib(ic_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(K_f = getListElement(Flist, "K_f"));                    PROTECT(dc_K_f = iDim(INTEGER(getAttrib(K_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(inv_f = getListElement(Flist, "inv_f"));                PROTECT(dc_inv_f = iDim(INTEGER(getAttrib(inv_f, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(FTE_f_m = getListElement(Flist, "FTE_f_m"));            PROTECT(dc_FTE_f_m = iDim(INTEGER(getAttrib(FTE_f_m, install("DimCst")))));//Rprintf("Eco 31");
    PROTECT(GVLref_f_m = getListElement(Flist, "GVLref_f_m"));      PROTECT(dc_GVLref_f_m = iDim(INTEGER(getAttrib(GVLref_f_m, install("DimCst")))));//Rprintf("Eco 31");

// ---> P = 14 + 35*2 = 84
//Rprintf("Eco 4");fichier << "Eco4" << endl;

    PROTECT(ue_f = NEW_NUMERIC(nbF));
    setAttrib(ue_f, R_DimSymbol, getAttrib(getListElement(Flist, "effort1_f"), R_DimSymbol));
    setAttrib(ue_f, R_DimNamesSymbol, getAttrib(getListElement(Flist, "effort1_f"), R_DimNamesSymbol));
    setAttrib(ue_f, install("DimCst"), getAttrib(getListElement(Flist, "effort1_f"), install("DimCst")));

    PROTECT(ue_f_m = NEW_NUMERIC(nbF*nbMe));
    setAttrib(ue_f_m, R_DimSymbol, getAttrib(getListElement(Flist, "effort1_f_m"), R_DimSymbol));
    setAttrib(ue_f_m, R_DimNamesSymbol, getAttrib(getListElement(Flist, "effort1_f_m"), R_DimNamesSymbol));
    setAttrib(ue_f_m, install("DimCst"), getAttrib(getListElement(Flist, "effort1_f_m"), install("DimCst")));

    double *r_ue_f = REAL(ue_f); double *reff1_f = REAL(getListElement(Flist, "effort1_f")) ; double *reff2_f = REAL(getListElement(Flist, "effort2_f"));
    double *r_ue_f_m = REAL(ue_f_m); double *reff1 = REAL(getListElement(Flist, "effort1_f_m")) ; double *reff2 = REAL(getListElement(Flist, "effort2_f_m"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        r_ue_f[ind_f] = reff1_f[ind_f]*reff2_f[ind_f];
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) r_ue_f_m[ind_f + nbF*ind_m] = reff1[ind_f + nbF*ind_m]*reff2[ind_f + nbF*ind_m];
    }

//Rprintf("Eco 5");fichier << "Eco5" << endl;
    PROTECT(dc_ue_f = iDim(INTEGER(getAttrib(ue_f, install("DimCst")))));
    PROTECT(dc_ue_f_m = iDim(INTEGER(getAttrib(ue_f_m, install("DimCst")))));
    dim_ue_f_m = INTEGER(dc_ue_f_m);
// ---> P = 84 + 4 = 88

    dim_nbv_f = INTEGER(dc_nbv_f);                          r_nbv_f = REAL(nbv_f);//Rprintf("Eco 51\n");
    dim_nbv_f_m = INTEGER(dc_nbv_f_m);                      r_nbv_f_m = REAL(nbv_f_m);//Rprintf("Eco 51\n");
    dim_lc_f_m = INTEGER(dc_lc_f_m);                        r_lc_f_m = REAL(lc_f_m);//Rprintf("Eco 51\n");
    /*dim_lcd_f_m = INTEGER(dc_lcd_f_m); */                 r_lcd_f_m = REAL(lcd_f_m);//Rprintf("Eco 51\n");
    //dim_tripLgth_f = INTEGER(dc_tripLgth_f);                r_tripLgth_f = REAL(tripLgth_f);//Rprintf("Eco 51\n");
    dim_tripLgth_f_m = INTEGER(dc_tripLgth_f_m);            r_tripLgth_f_m = REAL(tripLgth_f_m);//Rprintf("Eco 51\n");
    //dim_nbTrip_f = INTEGER(dc_nbTrip_f);                    r_nbTrip_f = REAL(nbTrip_f);//Rprintf("Eco 51\n");
    dim_nbTrip_f_m = INTEGER(dc_nbTrip_f_m);                r_nbTrip_f_m = REAL(nbTrip_f_m);//Rprintf("Eco 51\n");
    //dim_nbds_f = INTEGER(dc_nbds_f);                        r_nbds_f = REAL(nbds_f);//Rprintf("Eco 51\n");
    //dim_nbds_f_m = INTEGER(dc_nbds_f_m);                    r_nbds_f_m = REAL(nbds_f_m);//Rprintf("Eco 51\n");
    //dim_effort1_f = INTEGER(dc_effort1_f);                  r_effort1_f = REAL(effort1_f);//Rprintf("Eco 51\n");
    //dim_effort1_f_m = INTEGER(dc_effort1_f_m);              r_effort1_f_m = REAL(effort1_f_m);//Rprintf("Eco 51\n");
    //dim_effort2_f = INTEGER(dc_effort2_f);                  r_effort2_f = REAL(effort2_f);//Rprintf("Eco 51\n");
    //dim_effort2_f_m = INTEGER(dc_effort2_f_m);              r_effort2_f_m = REAL(effort2_f_m);//Rprintf("Eco 51\n");
    dim_Lref_f_m = INTEGER(dc_Lref_f_m);                    r_Lref_f_m = REAL(Lref_f_m);//Rprintf("Eco 51\n");
    dim_cnb_f_m = INTEGER(dc_cnb_f_m);                      r_cnb_f_m = REAL(cnb_f_m);//Rprintf("Eco 51\n");
    dim_ovcDCF_f_m = INTEGER(dc_ovcDCF_f_m);                r_ovcDCF_f_m = REAL(ovcDCF_f_m);//Rprintf("Eco 51\n");
    dim_fc_f_m = INTEGER(dc_fc_f_m);                        r_fc_f_m = REAL(fc_f_m);//Rprintf("Eco 51\n");
    dim_vf_f_m = INTEGER(dc_vf_f_m);                        r_vf_f_m = REAL(vf_f_m);//Rprintf("Eco 51\n");
    dim_cshr_f_m = INTEGER(dc_cshr_f_m);                    r_cshr_f_m = REAL(cshr_f_m);//Rprintf("Eco 51\n");
    dim_cshr_f = INTEGER(dc_cshr_f);                        r_cshr_f = REAL(cshr_f);//Rprintf("Eco 51\n");
    //dim_cnb_f = INTEGER(dc_cnb_f);                          r_cnb_f = REAL(cnb_f);//Rprintf("Eco 51\n");
    dim_persc_f = INTEGER(dc_persc_f);                      r_persc_f = REAL(persc_f);//Rprintf("Eco 51\n");
    dim_eec_f = INTEGER(dc_eec_f);                          r_eec_f = REAL(eec_f);//Rprintf("Eco 51\n");
    dim_mwh_f = INTEGER(dc_mwh_f);                          r_mwh_f = REAL(mwh_f);//Rprintf("Eco 51\n");
    dim_rep_f = INTEGER(dc_rep_f);                          r_rep_f = REAL(rep_f);//Rprintf("Eco 51\n");
    dim_gc_f = INTEGER(dc_gc_f);                            r_gc_f = REAL(gc_f);//Rprintf("Eco 51\n");
    dim_fixc_f = INTEGER(dc_fixc_f);                        r_fixc_f = REAL(fixc_f);//Rprintf("Eco 51\n");
    dim_FTE_f = INTEGER(dc_FTE_f);                          r_FTE_f = REAL(FTE_f);//Rprintf("Eco 51\n");
    dim_dep_f = INTEGER(dc_dep_f);                          r_dep_f = REAL(dep_f);//Rprintf("Eco 51\n");
    dim_ic_f = INTEGER(dc_ic_f);                            r_ic_f = REAL(ic_f);//Rprintf("Eco 51\n");
    dim_K_f = INTEGER(dc_K_f);                              r_K_f = REAL(K_f);//Rprintf("Eco 51\n");
    dim_inv_f = INTEGER(dc_inv_f);                          r_inv_f = REAL(inv_f);//Rprintf("Eco 51\n");
    //dim_FTE_f_m = INTEGER(dc_FTE_f_m);                    r_FTE_f_m = REAL(FTE_f_m);//Rprintf("Eco 51\n");
    dim_GVLref_f_m = INTEGER(dc_GVLref_f_m);                r_GVLref_f_m = REAL(GVLref_f_m);//Rprintf("Eco 51\n");


    int nbC=0;
    //int nbI=0;

//Rprintf("Eco 6");fichier << "Eco6" << endl;

if (ind_t==0) {

    SEXP ETini_f_m, fvolue_f_m, ovcDCFue_f_m, rtbsIni_f, ccwr_f, opersc_f, eco_names,
         GVLcom_f_m_e_out, GVLcom_f_m_eStat_out, GVLcom_f_m_e, GVLst_f_m_e_out, GVLst_f_m_eStat_out, GVLst_f_m_e, GVL_f_m_e_out, GVL_f_m_eStat_out, GVLtot_f_m_e,
         GVLtot_f_m_out, GVLav_f_m_out, GVLtot_f_out, GVLav_f_out, NGVLav_f_m_out, NGVLav_f_out, ET_f_m_out,
         cnb_f_m_out, cnb_f_out, rtbs_f_m_out, rtbs_f_out, rtbsAct_f_out, cshrT_f_m_out, cshrT_f_out, ncshr_f_out, ocl_f_out, cs_f_out, csAct_f_out, csTot_f_out,
         gva_f_out, gvaAct_f_out, gvamargin_f_out, gva_FTE_f_out, ccw_f_out, ccwCr_f_out, wageg_f_out, wagen_f_out, wageg_FTE_f_out, wageg_h_f_out,
         gp_f_out, gpAct_f_out, gpmargin_f_out, ncf_f_out, np_f_out, npmargin_f_out, prof_f_out, npmargin_trend_f_out,
         ssTot_f_out, ps_f_out, psAct_f_out, sts_f_out, stsAct_f_out, BER_f_out, CR_BER_f_out, fuelEff_f_out,
         ratio_fvol_gva_f_out, ratio_gp_gva_f_out, ratio_GVL_K_f_out, ratio_gp_K_f_out, RoFTA_f_out, ROI_f_out, ratio_np_K_f_out, ratio_GVL_cnb_ue_f_out, QuotaExp_f_out;

    double  *r_ETini_f_m, *r_fvolue_f_m, *r_ovcDCFue_f_m, *r_rtbsIni_f, *r_ccwr_f, *r_opersc_f;

//-------------------------
// Stade pr�liminaire (temps initial)
//-------------------------

    PROTECT(ETini_f_m = NEW_NUMERIC(nbF*nbMe));                 r_ETini_f_m = REAL(ETini_f_m);
    PROTECT(fvolue_f_m = NEW_NUMERIC(nbF*nbMe));                r_fvolue_f_m = REAL(fvolue_f_m);
    PROTECT(ovcDCFue_f_m = NEW_NUMERIC(nbF*nbMe));              r_ovcDCFue_f_m = REAL(ovcDCFue_f_m);
    PROTECT(rtbsIni_f = NEW_NUMERIC(nbF));                      r_rtbsIni_f = REAL(rtbsIni_f);
    PROTECT(ccwr_f = NEW_NUMERIC(nbF));                         r_ccwr_f = REAL(ccwr_f);
    PROTECT(opersc_f = NEW_NUMERIC(nbF));                       r_opersc_f = REAL(opersc_f);
// ---> P(t0) = 6
//Rprintf("Eco 7");fichier << "Eco7" << endl;

// on cr�e ETini
//    double *rnbTrip = REAL(getListElement(Flist, "nbTrip_f_m"));
//    double *rtripLgth = REAL(getListElement(Flist, "tripLgth_f_m"));
//    double *rnbv = REAL(getListElement(Flist, "nbv_f_m"));
//    double *rcnb = REAL(getListElement(Flist, "cnb_f_m"));

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++) {
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) {
            r_ETini_f_m[ind_f + nbF*ind_m] =
               finite(
               r_Lref_f_m[ind_f*dim_Lref_f_m[0] + ind_m*dim_Lref_f_m[1] + 0*dim_Lref_f_m[2] + ind_t*dim_Lref_f_m[3]] /
               (r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]] *
                r_nbTrip_f_m[ind_f*dim_nbTrip_f_m[0] + ind_m*dim_nbTrip_f_m[1] + 0*dim_nbTrip_f_m[2] + ind_t*dim_nbTrip_f_m[3]] *
                r_tripLgth_f_m[ind_f*dim_tripLgth_f_m[0] + ind_m*dim_tripLgth_f_m[1] + 0*dim_tripLgth_f_m[2] + ind_t*dim_tripLgth_f_m[3]] *
                r_cnb_f_m[ind_f*dim_cnb_f_m[0] + ind_m*dim_cnb_f_m[1] + 0*dim_cnb_f_m[2] + ind_t*dim_cnb_f_m[3]]));
            //if (ISNA(r_ETini_f_m[ind_f + nbF*ind_m])) r_ETini_f_m[ind_f + nbF*ind_m] = 0.0;
        }
    }


//Rprintf("Eco 8");fichier << "Eco8" << endl;
 for (int e = 0 ; e < nbE+nbEstat ; e++) {

        if (e<nbE) {
         PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e))));
        } else {
         PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e-nbE))));
        }

        //if (e<nbE) nbI = length(getListElement(elmt, "modI"));
        if (e<nbE) nbC = length(getListElement(elmt, "modC"));

        PROTECT(GVLtot_f_m_e = NEW_NUMERIC(nbF*nbMe*nbT));
        PROTECT(GVLcom_f_m_e = NEW_NUMERIC(nbF*nbMe*nbT));
        PROTECT(GVLst_f_m_e = NEW_NUMERIC(nbF*nbMe*nbT));

        double *r_GVLtot_f_m_e = REAL(GVLtot_f_m_e);
        double *r_GVLcom_f_m_e = REAL(GVLcom_f_m_e);
        double *r_GVLst_f_m_e = REAL(GVLst_f_m_e);

        double *r_Lbio_f_m_e ,  *r_P_f_m_e, r_Pst_e=NA_REAL, *r_LD_efmc=&NA_REAL, *r_statLDor_efm=&NA_REAL, *r_statLDst_efm=&NA_REAL, r_theta_e;
        int *dim_Lbio_e, *dim_P_e;
//Rprintf("Eco 9");fichier << "Eco9" << endl;
        if (e<nbE) {
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct, e));
            //r_Lbio_f_sum_e = REAL(aggregObj(VECTOR_ELT(out_L_efmct, e),dimCstF));
            r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
            r_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct, e), install("DimCst")))));
            dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));
        } else {
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_Lstat, e-nbE));
            //r_Lbio_f_sum_e = REAL(aggregObj(VECTOR_ELT(out_Lstat, e-nbE),dimCstF));
            r_P_f_m_e = REAL(VECTOR_ELT(out_Pstat, e-nbE));
            r_statLDor_efm = REAL(VECTOR_ELT(out_statLDor_efm, e-nbE));
            r_statLDst_efm = REAL(VECTOR_ELT(out_statLDst_efm, e-nbE));
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_Lstat, e-nbE), install("DimCst")))));
            dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_Pstat, e-nbE), install("DimCst")))));
            r_Pst_e = REAL(getListElement(elmt, "Pst_e"))[0];
        }

        //------------------------------
        //�quations de la table "p"
        //------------------------------
//Rprintf("Eco 10");fichier << "Eco10" << endl;

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){   //on rappelle ici que ind_t est en fait �gal � 0

    //double countGVLtotf = 0.0; //pour sommer GVLtot_f_m_e sur les m�tiers

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

         //-- 3. GVLtot_f_m_e

    double countCom = 0.0;

    if (e<nbE) {

             if (ISNA(r_theta_e)) r_theta_e = 1.0;

             for (int ind_c = 0 ; ind_c < (nbC-1) ; ind_c++){ //sur les classes non sous-tailles

                if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]]))
                    r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

                countCom = countCom +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                  r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

             }

             }

             if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]]))
                    r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;

             //if (ISNA(r_theta_e)) r_theta_e = 1.0;

             if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

                  r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                   r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                   r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                   r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                   r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

             } else {

                  r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;

             }

               //if (e==1 & ind_f==0 & ind_m==4) PrintValue(ETini_f_m_out);

    } else {

        if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]]))
                    r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;

        if (ISNA(r_theta_e)) r_theta_e = 1.0;

        if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

            countCom = r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                  r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                  finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]) ;

        }

        if (!ISNA(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

        if (ISNA(r_Pst_e)) r_Pst_e = 0.0;

            r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_Pst_e * 1000 * r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]];

        }

    }

            r_GVLcom_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = countCom;
            r_GVLtot_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
            r_GVLcom_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
            r_GVLst_f_m_e[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

//Rprintf("Eco 11");fichier << "Eco11" << endl;


        }
        }

        //on formatte le(s) r�sultat(s) et on les int�gre � 'eVar'
//Rprintf("Eco 13");fichier << "Eco13" << endl;
        setAttrib(GVLtot_f_m_e, R_DimSymbol, DimFM);
        setAttrib(GVLtot_f_m_e, R_DimNamesSymbol, dimnamesFM);
        setAttrib(GVLtot_f_m_e, install("DimCst"), dimCstFM);
        if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 41, GVLtot_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 1, GVLtot_f_m_e);

        setAttrib(GVLcom_f_m_e, R_DimSymbol, DimFM);
        setAttrib(GVLcom_f_m_e, R_DimNamesSymbol, dimnamesFM);
        setAttrib(GVLcom_f_m_e, install("DimCst"), dimCstFM);
        if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 246, GVLcom_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 8, GVLcom_f_m_e);

        setAttrib(GVLst_f_m_e, R_DimSymbol, DimFM);
        setAttrib(GVLst_f_m_e, R_DimNamesSymbol, dimnamesFM);
        setAttrib(GVLst_f_m_e, install("DimCst"), dimCstFM);
        if (e<nbE) SET_VECTOR_ELT(VECTOR_ELT(eVar, e), 247, GVLst_f_m_e); else SET_VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE), 9, GVLst_f_m_e);

//Rprintf("Eco 14");fichier << "Eco14" << endl;
        UNPROTECT(4);

}


//Rprintf("Eco 8");fichier << "Eco8" << endl;

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

       double countRTBSnum = 0.0;

       for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        countRTBSnum = countRTBSnum + finite(r_GVLref_f_m[ind_f*dim_GVLref_f_m[0] + ind_m*dim_GVLref_f_m[1] + 0*dim_GVLref_f_m[2] + ind_t*dim_GVLref_f_m[3]] *
                                      r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]] *
                                      (1 - finite(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]))) -
                                      finite(r_ovcDCF_f_m[ind_f*dim_ovcDCF_f_m[0] + ind_m*dim_ovcDCF_f_m[1] + 0*dim_ovcDCF_f_m[2] + ind_t*dim_ovcDCF_f_m[3]]) -
                                      finite(r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]]);

      //-- 4. fvolue_f_m

        r_fvolue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_fc_f_m[ind_f*dim_fc_f_m[0] + ind_m*dim_fc_f_m[1] + 0*dim_fc_f_m[2] + ind_t*dim_fc_f_m[3]] /
                    (r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]] *
                    r_ue_f_m[ind_f + ind_m*nbF]) );

      //-- 5. ovcDCFue_f_m

        r_ovcDCFue_f_m[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    finite( r_ovcDCF_f_m[ind_f*dim_ovcDCF_f_m[0] + ind_m*dim_ovcDCF_f_m[1] + 0*dim_ovcDCF_f_m[2] + ind_t*dim_ovcDCF_f_m[3]] /
                    r_ue_f_m[ind_f + ind_m*nbF] );

        }

//Rprintf("Eco 18");fichier << "Eco18" << endl;

        //-- 6. rtbsIni_f

            r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                  countRTBSnum / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


        //-- 7. ccwr_f

            r_ccwr_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] /
                r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        //-- 8. opersc_f

            r_opersc_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_persc_f[ind_f*dim_persc_f[0] + 0*dim_persc_f[1] + 0*dim_persc_f[2] + ind_t*dim_persc_f[3]] -
                (0.01 * r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] * //cshr_f en %
                r_rtbsIni_f[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]);

    }

//Rprintf("Eco 19\n");fichier << "Eco19" << endl;

//on formatte le(s) r�sultat(s) et on int�gre � fVar

        setAttrib(fvolue_f_m, R_DimSymbol, DimFMini);
        setAttrib(fvolue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(fvolue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 4, fvolue_f_m);

        setAttrib(ovcDCFue_f_m, R_DimSymbol, DimFMini);
        setAttrib(ovcDCFue_f_m, R_DimNamesSymbol, dimnamesFMini);
        setAttrib(ovcDCFue_f_m, install("DimCst"), dimCstFMini);
        SET_VECTOR_ELT(fVar, 10, ovcDCFue_f_m);

        setAttrib(ccwr_f, R_NamesSymbol, fleetList);
        setAttrib(ccwr_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 27, ccwr_f);

        setAttrib(opersc_f, R_NamesSymbol,  fleetList);
        setAttrib(opersc_f, install("DimCst"), dimCstFini);
        SET_VECTOR_ELT(fVar, 28, opersc_f);

        SET_VECTOR_ELT(fVar, 31, rtbsIni_f);

        SET_VECTOR_ELT(fVar, 33, ETini_f_m);
//Rprintf("Eco 20\n");fichier << "Eco20" << endl;
//enfin, on initialise l'output


    PROTECT(GVLcom_f_m_e_out = allocVector(VECSXP, nbE));
    setAttrib(GVLcom_f_m_e_out, R_NamesSymbol, sppList);
    SET_VECTOR_ELT(out_EcoDCF, 0, GVLcom_f_m_e_out);//Rprintf("Eco 20\n");

    PROTECT(GVLcom_f_m_eStat_out = allocVector(VECSXP, nbEstat));
    setAttrib(GVLcom_f_m_eStat_out, R_NamesSymbol, sppListStat);
    SET_VECTOR_ELT(out_EcoDCF, 1, GVLcom_f_m_eStat_out);//Rprintf("Eco 20\n");

    PROTECT(GVLst_f_m_e_out = allocVector(VECSXP, nbE));
    setAttrib(GVLst_f_m_e_out, R_NamesSymbol, sppList);
    SET_VECTOR_ELT(out_EcoDCF, 2, GVLst_f_m_e_out);//Rprintf("Eco 20\n");

    PROTECT(GVLst_f_m_eStat_out = allocVector(VECSXP, nbEstat));
    setAttrib(GVLst_f_m_eStat_out, R_NamesSymbol, sppListStat);
    SET_VECTOR_ELT(out_EcoDCF, 3, GVLst_f_m_eStat_out);//Rprintf("Eco 20\n");

    PROTECT(GVL_f_m_e_out = allocVector(VECSXP, nbE));
    setAttrib(GVL_f_m_e_out, R_NamesSymbol, sppList);
    SET_VECTOR_ELT(out_EcoDCF, 4, GVL_f_m_e_out);//Rprintf("Eco 20\n");

    PROTECT(GVL_f_m_eStat_out = allocVector(VECSXP, nbEstat));
    setAttrib(GVL_f_m_eStat_out, R_NamesSymbol, sppListStat);
    SET_VECTOR_ELT(out_EcoDCF, 5, GVL_f_m_eStat_out);//Rprintf("Eco 20\n");

    PROTECT(GVLtot_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLtot_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLtot_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLtot_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 6, GVLtot_f_m_out);//Rprintf("Eco 20\n");

    PROTECT(GVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(GVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(GVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(GVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 7, GVLav_f_m_out);//Rprintf("Eco 20\n");

    PROTECT(GVLtot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLtot_f_out, R_DimSymbol, DimF);
    setAttrib(GVLtot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLtot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 8, GVLtot_f_out);//Rprintf("Eco 20\n");

    PROTECT(GVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(GVLav_f_out, R_DimSymbol, DimF);
    setAttrib(GVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(GVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 9, GVLav_f_out);//Rprintf("Eco 20\n");

    PROTECT(NGVLav_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(NGVLav_f_m_out, R_DimSymbol, DimFM);
    setAttrib(NGVLav_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(NGVLav_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 10, NGVLav_f_m_out);//Rprintf("Eco 20\n");

    PROTECT(NGVLav_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(NGVLav_f_out, R_DimSymbol, DimF);
    setAttrib(NGVLav_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(NGVLav_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 11, NGVLav_f_out);//Rprintf("Eco 20\n");

    PROTECT(ET_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(ET_f_m_out, R_DimSymbol, DimFM);
    setAttrib(ET_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(ET_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 12, ET_f_m_out);//Rprintf("Eco 20\n");

    PROTECT(cnb_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(cnb_f_m_out, R_DimSymbol, DimFM);
    setAttrib(cnb_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(cnb_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 13, cnb_f_m_out);//Rprintf("Eco 20\n");

    PROTECT(cnb_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cnb_f_out, R_DimSymbol, DimF);
    setAttrib(cnb_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cnb_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 14, cnb_f_out);//Rprintf("Eco 20\n");

    PROTECT(rtbs_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(rtbs_f_m_out, R_DimSymbol, DimFM);
    setAttrib(rtbs_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(rtbs_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 15, rtbs_f_m_out);

    PROTECT(rtbs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbs_f_out, R_DimSymbol, DimF);
    setAttrib(rtbs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 16, rtbs_f_out);//Rprintf("Eco 20\n");

    PROTECT(rtbsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(rtbsAct_f_out, R_DimSymbol, DimF);
    setAttrib(rtbsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(rtbsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 17, rtbsAct_f_out);

    PROTECT(cshrT_f_m_out = NEW_NUMERIC(nbF*nbMe*nbT));
    setAttrib(cshrT_f_m_out, R_DimSymbol, DimFM);
    setAttrib(cshrT_f_m_out, R_DimNamesSymbol, dimnamesFM);
    setAttrib(cshrT_f_m_out, install("DimCst"), dimCstFM);
    SET_VECTOR_ELT(out_EcoDCF, 18, cshrT_f_m_out);

    PROTECT(cshrT_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cshrT_f_out, R_DimSymbol, DimF);
    setAttrib(cshrT_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cshrT_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 19, cshrT_f_out);//Rprintf("Eco 20\n");

    PROTECT(ncshr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ncshr_f_out, R_DimSymbol, DimF);
    setAttrib(ncshr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ncshr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 20, ncshr_f_out);//Rprintf("Eco 20\n");

    PROTECT(ocl_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ocl_f_out, R_DimSymbol, DimF);
    setAttrib(ocl_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ocl_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 21, ocl_f_out);//Rprintf("Eco 20\n");

    PROTECT(cs_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(cs_f_out, R_DimSymbol, DimF);
    setAttrib(cs_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(cs_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 22, cs_f_out);//Rprintf("Eco 20\n");

    PROTECT(csAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csAct_f_out, R_DimSymbol, DimF);
    setAttrib(csAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 23, csAct_f_out);

    PROTECT(csTot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(csTot_f_out, R_DimSymbol, DimF);
    setAttrib(csTot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(csTot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 24, csTot_f_out);

    PROTECT(gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gva_f_out, R_DimSymbol, DimF);
    setAttrib(gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 25, gva_f_out);//Rprintf("Eco 20\n");

    PROTECT(gvaAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gvaAct_f_out, R_DimSymbol, DimF);
    setAttrib(gvaAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gvaAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 26, gvaAct_f_out);//Rprintf("Eco 20\n");

    PROTECT(gvamargin_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gvamargin_f_out, R_DimSymbol, DimF);
    setAttrib(gvamargin_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gvamargin_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 27, gvamargin_f_out);//Rprintf("Eco 20\n");

    PROTECT(gva_FTE_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gva_FTE_f_out, R_DimSymbol, DimF);
    setAttrib(gva_FTE_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gva_FTE_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 28, gva_FTE_f_out);//Rprintf("Eco 20\n");

    PROTECT(ccw_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccw_f_out, R_DimSymbol, DimF);
    setAttrib(ccw_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccw_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 29, ccw_f_out);//Rprintf("Eco 20\n");

    PROTECT(ccwCr_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ccwCr_f_out, R_DimSymbol, DimF);
    setAttrib(ccwCr_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ccwCr_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 30, ccwCr_f_out);//Rprintf("Eco 20\n");

    PROTECT(wageg_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wageg_f_out, R_DimSymbol, DimF);
    setAttrib(wageg_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wageg_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 31, wageg_f_out);//Rprintf("Eco 20\n");

    PROTECT(wagen_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wagen_f_out, R_DimSymbol, DimF);
    setAttrib(wagen_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wagen_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 32, wagen_f_out);//Rprintf("Eco 20\n");

    PROTECT(wageg_FTE_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wageg_FTE_f_out, R_DimSymbol, DimF);
    setAttrib(wageg_FTE_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wageg_FTE_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 33, wageg_FTE_f_out);//Rprintf("Eco 20\n");

    PROTECT(wageg_h_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(wageg_h_f_out, R_DimSymbol, DimF);
    setAttrib(wageg_h_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(wageg_h_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 34, wageg_h_f_out);//Rprintf("Eco 20\n");

    PROTECT(gp_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gp_f_out, R_DimSymbol, DimF);
    setAttrib(gp_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gp_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 35, gp_f_out);//Rprintf("Eco 20\n");

    PROTECT(gpAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gpAct_f_out, R_DimSymbol, DimF);
    setAttrib(gpAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gpAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 36, gpAct_f_out);//Rprintf("Eco 20\n");

    PROTECT(gpmargin_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(gpmargin_f_out, R_DimSymbol, DimF);
    setAttrib(gpmargin_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(gpmargin_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 37, gpmargin_f_out);//Rprintf("Eco 20\n");

    PROTECT(ncf_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ncf_f_out, R_DimSymbol, DimF);
    setAttrib(ncf_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ncf_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 38, ncf_f_out);//Rprintf("Eco 20\n");

    PROTECT(np_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(np_f_out, R_DimSymbol, DimF);
    setAttrib(np_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(np_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 39, np_f_out);//Rprintf("Eco 20\n");

    PROTECT(npmargin_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(npmargin_f_out, R_DimSymbol, DimF);
    setAttrib(npmargin_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(npmargin_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 40, npmargin_f_out);//Rprintf("Eco 20\n");

    PROTECT(prof_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(prof_f_out, R_DimSymbol, DimF);
    setAttrib(prof_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(prof_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 41, prof_f_out);//Rprintf("Eco 20\n");

    PROTECT(npmargin_trend_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(npmargin_trend_f_out, R_DimSymbol, DimF);
    setAttrib(npmargin_trend_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(npmargin_trend_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 42, npmargin_trend_f_out);//Rprintf("Eco 20\n");

    PROTECT(ssTot_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ssTot_f_out, R_DimSymbol, DimF);
    setAttrib(ssTot_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ssTot_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 43, ssTot_f_out);

    PROTECT(ps_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ps_f_out, R_DimSymbol, DimF);
    setAttrib(ps_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ps_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 44, ps_f_out);

    PROTECT(psAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(psAct_f_out, R_DimSymbol, DimF);
    setAttrib(psAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(psAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 45, psAct_f_out);

    PROTECT(sts_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(sts_f_out, R_DimSymbol, DimF);
    setAttrib(sts_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(sts_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 46, sts_f_out);

    PROTECT(stsAct_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(stsAct_f_out, R_DimSymbol, DimF);
    setAttrib(stsAct_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(stsAct_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 47, stsAct_f_out);

    PROTECT(BER_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(BER_f_out, R_DimSymbol, DimF);
    setAttrib(BER_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(BER_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 48, BER_f_out);

    PROTECT(CR_BER_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(CR_BER_f_out, R_DimSymbol, DimF);
    setAttrib(CR_BER_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(CR_BER_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 49, CR_BER_f_out);

    PROTECT(fuelEff_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(fuelEff_f_out, R_DimSymbol, DimF);
    setAttrib(fuelEff_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(fuelEff_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 50, fuelEff_f_out);

    PROTECT(ratio_fvol_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_fvol_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_fvol_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_fvol_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 51, ratio_fvol_gva_f_out);

    PROTECT(ratio_gp_gva_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gp_gva_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gp_gva_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gp_gva_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 52, ratio_gp_gva_f_out);

    PROTECT(ratio_GVL_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 53, ratio_GVL_K_f_out);

    PROTECT(ratio_gp_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_gp_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_gp_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_gp_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 54, ratio_gp_K_f_out);

    PROTECT(RoFTA_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(RoFTA_f_out, R_DimSymbol, DimF);
    setAttrib(RoFTA_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(RoFTA_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 55, RoFTA_f_out);

    PROTECT(ROI_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ROI_f_out, R_DimSymbol, DimF);
    setAttrib(ROI_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ROI_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 56, ROI_f_out);

    PROTECT(ratio_np_K_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_np_K_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_np_K_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_np_K_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 57, ratio_np_K_f_out);

    PROTECT(ratio_GVL_cnb_ue_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimSymbol, DimF);
    setAttrib(ratio_GVL_cnb_ue_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(ratio_GVL_cnb_ue_f_out, install("DimCst"), dimCstF);
    SET_VECTOR_ELT(out_EcoDCF, 58, ratio_GVL_cnb_ue_f_out);

    PROTECT(QuotaExp_f_out = NEW_NUMERIC(nbF*nbT));
    setAttrib(QuotaExp_f_out, R_DimSymbol, DimF);
    setAttrib(QuotaExp_f_out, R_DimNamesSymbol, dimnamesF);
    setAttrib(QuotaExp_f_out, install("DimCst"), dimCstF);
    for (int ind_f = 0 ; ind_f < nbF ; ind_f++)
        for (int ind_tt = 0 ; ind_tt < nbT ; ind_tt++) REAL(QuotaExp_f_out)[ind_f + nbF*ind_tt] = 0.0;
    SET_VECTOR_ELT(out_EcoDCF, 59, QuotaExp_f_out);

//Rprintf("Eco 20.8\n");fichier << "Eco 20.8" << endl;
    //on nomme les �l�ments de out_EcoDCF




    const char *namesEco[60] = {"GVLcom_f_m_e_out","GVLcom_f_m_eStat_out","GVLst_f_m_e_out","GVLst_f_m_eStat_out","GVL_f_m_e_out","GVL_f_m_eStat_out","GVLtot_f_m_out",
                                "GVLav_f_m_out","GVLtot_f_out","GVLav_f_out","NGVLav_f_m_out","NGVLav_f_out","ET_f_m_out","cnb_f_m_out","cnb_f_out","rtbs_f_m_out","rtbs_f_out",
                                "rtbsAct_f_out","cshrT_f_m_out","cshrT_f_out","ncshr_f_out","ocl_f_out","cs_f_out","csAct_f_out","csTot_f_out","gva_f_out","gvaAct_f_out",
                                "gvamargin_f_out","gva_FTE_f_out","ccw_f_out","ccwCr_f_out","wageg_f_out","wagen_f_out","wageg_FTE_f_out","wageg_h_f_out","gp_f_out",
                                "gpAct_f_out","gpmargin_f_out","ncf_f_out","np_f_out","npmargin_f_out","prof_f_out","npmargin_trend_f_out","ssTot_f_out","ps_f_out",
                                "psAct_f_out","sts_f_out","stsAct_f_out","BER_f_out","CR_BER_f_out","fuelEff_f_out","ratio_fvol_gva_f_out","ratio_gp_gva_f_out",
                                "ratio_GVL_K_f_out","ratio_gp_K_f_out","RoFTA_f_out","ROI_f_out","ratio_np_K_f_out","ratio_GVL_cnb_ue_f_out","QuotaExp_f_out"};

//Rprintf("Eco 20.9\n");fichier << "Eco 20.9" << endl;
    PROTECT(eco_names = allocVector(STRSXP, 60));

    for(int ct = 0; ct < 60; ct++) SET_STRING_ELT(eco_names, ct, mkChar(namesEco[ct])); //PrintValue(out_EcoDCF);
//Rprintf("Eco 20.95\n");fichier << "Eco20.95" << endl;
    setAttrib(out_EcoDCF, R_NamesSymbol, eco_names);
//Rprintf("Eco 21\n");fichier << "Eco21" << endl;

// ---> P(t0) = 6 + 60 = 66
}


//on importe les outputs afin de les mettre � jour � l'instant ind_t

//    r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 0));//Rprintf("Eco 20\n");
//    r_GVLcom_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 1));//Rprintf("Eco 20\n");
//    r_GVLst_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 2));//Rprintf("Eco 20\n");
//    r_GVLst_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 3));//Rprintf("Eco 20\n");
//    r_GVL_f_m_e_out = REAL(VECTOR_ELT(out_EcoDCF, 4));//Rprintf("Eco 20\n");
//    r_GVL_f_m_eStat_out = REAL(VECTOR_ELT(out_EcoDCF, 5));//Rprintf("Eco 20\n");
    r_GVLtot_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 6));//Rprintf("Eco 20\n");
    r_GVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 7));//Rprintf("Eco 20\n");
    r_GVLtot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 8));//Rprintf("Eco 20\n");
    r_GVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF, 9));//Rprintf("Eco 20\n");
    r_NGVLav_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 10));//Rprintf("Eco 20\n");
    r_NGVLav_f_out = REAL(VECTOR_ELT(out_EcoDCF, 11));//Rprintf("Eco 20\n");
    r_ET_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 12));//Rprintf("Eco 20\n");
    r_cnb_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 13));//Rprintf("Eco 20\n");
    r_cnb_f_out = REAL(VECTOR_ELT(out_EcoDCF, 14));//Rprintf("Eco 20\n");
    r_rtbs_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 15));
    r_rtbs_f_out = REAL(VECTOR_ELT(out_EcoDCF, 16));//Rprintf("Eco 20\n");
    r_rtbsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 17));
    r_cshrT_f_m_out = REAL(VECTOR_ELT(out_EcoDCF, 18));
    r_cshrT_f_out = REAL(VECTOR_ELT(out_EcoDCF, 19));//Rprintf("Eco 20\n");
    r_ncshr_f_out = REAL(VECTOR_ELT(out_EcoDCF, 20));//Rprintf("Eco 20\n");
    r_ocl_f_out = REAL(VECTOR_ELT(out_EcoDCF, 21));//Rprintf("Eco 20\n");
    r_cs_f_out = REAL(VECTOR_ELT(out_EcoDCF, 22));//Rprintf("Eco 20\n");
    r_csAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 23));
    r_csTot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 24));
    r_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 25));//Rprintf("Eco 20\n");
    r_gvaAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 26));//Rprintf("Eco 20\n");
    r_gvamargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 27));//Rprintf("Eco 20\n");
    r_gva_FTE_f_out = REAL(VECTOR_ELT(out_EcoDCF, 28));//Rprintf("Eco 20\n");
    r_ccw_f_out = REAL(VECTOR_ELT(out_EcoDCF, 29));//Rprintf("Eco 20\n");
    r_ccwCr_f_out = REAL(VECTOR_ELT(out_EcoDCF, 30));//Rprintf("Eco 20\n");
    r_wageg_f_out = REAL(VECTOR_ELT(out_EcoDCF, 31));//Rprintf("Eco 20\n");
    r_wagen_f_out = REAL(VECTOR_ELT(out_EcoDCF, 32));//Rprintf("Eco 20\n");
    r_wageg_FTE_f_out = REAL(VECTOR_ELT(out_EcoDCF, 33));//Rprintf("Eco 20\n");
    r_wageg_h_f_out = REAL(VECTOR_ELT(out_EcoDCF, 34));//Rprintf("Eco 20\n");
    r_gp_f_out = REAL(VECTOR_ELT(out_EcoDCF, 35));//Rprintf("Eco 20\n");
    r_gpAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 36));//Rprintf("Eco 20\n");
    r_gpmargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 37));//Rprintf("Eco 20\n");
    r_ncf_f_out = REAL(VECTOR_ELT(out_EcoDCF, 38));//Rprintf("Eco 20\n");
    r_np_f_out = REAL(VECTOR_ELT(out_EcoDCF, 39));//Rprintf("Eco 20\n");
    r_npmargin_f_out = REAL(VECTOR_ELT(out_EcoDCF, 40));//Rprintf("Eco 20\n");
    r_prof_f_out = REAL(VECTOR_ELT(out_EcoDCF, 41));//Rprintf("Eco 20\n");
    r_npmargin_trend_f_out = REAL(VECTOR_ELT(out_EcoDCF, 42));//Rprintf("Eco 20\n");
    r_ssTot_f_out = REAL(VECTOR_ELT(out_EcoDCF, 43));
    r_ps_f_out = REAL(VECTOR_ELT(out_EcoDCF, 44));
    r_psAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 45));
    r_sts_f_out = REAL(VECTOR_ELT(out_EcoDCF, 46));
    r_stsAct_f_out = REAL(VECTOR_ELT(out_EcoDCF, 47));
    r_BER_f_out = REAL(VECTOR_ELT(out_EcoDCF, 48));
    r_CR_BER_f_out = REAL(VECTOR_ELT(out_EcoDCF, 49));
    r_fuelEff_f_out = REAL(VECTOR_ELT(out_EcoDCF, 50));
    r_ratio_fvol_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 51));
    r_ratio_gp_gva_f_out = REAL(VECTOR_ELT(out_EcoDCF, 52));
    r_ratio_GVL_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 53));
    r_ratio_gp_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 54));
    r_RoFTA_f_out = REAL(VECTOR_ELT(out_EcoDCF, 55));
    r_ROI_f_out = REAL(VECTOR_ELT(out_EcoDCF, 56));
    r_ratio_np_K_f_out = REAL(VECTOR_ELT(out_EcoDCF, 57));
    r_ratio_GVL_cnb_ue_f_out = REAL(VECTOR_ELT(out_EcoDCF, 58));
    r_QuotaExp_f_out = REAL(VECTOR_ELT(out_EcoDCF, 59));

//Rprintf("Eco 22\n");fichier << "Eco22" << endl;

    double *r_fvolue_f_m2 = REAL(VECTOR_ELT(fVar,4));
    double *r_ovcDCFue_f_m2 = REAL(VECTOR_ELT(fVar,10));
    double *r_ccwr_f2 = REAL(VECTOR_ELT(fVar,27));
    double *r_opersc_f2 = REAL(VECTOR_ELT(fVar,28));

    SEXP countLf;
    PROTECT(countLf = NEW_NUMERIC(nbF)); // --> 67
    double *r_countLf = REAL(countLf);
    for (int INd_f = 0 ; INd_f < nbF ; INd_f++) r_countLf[INd_f] = 0.0; // pour le calcul de 'fuelEff'

// ---> P = 88 + 1 = 89



double *rnbv = REAL(getListElement(Flist, "nbv_f_m"));
double *rnbv_f = REAL(getListElement(Flist, "nbv_f"));
double *rnbTrip = REAL(getListElement(Flist, "nbTrip_f_m"));
double *rtripLgth = REAL(getListElement(Flist, "tripLgth_f_m"));
double *rnbTrip_f = REAL(getListElement(Flist, "nbTrip_f"));
double *rtripLgth_f = REAL(getListElement(Flist, "tripLgth_f"));


//Rprintf("Eco 23\n");fichier << "Eco23" << endl;


//initialisation de cnb, GVLtot et NGVLav, et remplissage de ET


    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        r_cnb_f_out[ind_f + nbF*ind_t] = 0.0;

        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            if (ind_t==0) {
              r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = REAL(VECTOR_ELT(fVar,33))[ind_f + nbF*ind_m];
            } else {
              r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*(ind_t-1)];
            }

            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = 0.001;

            r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;

            r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
        }
    }



// indicateurs esp�ces ---------------------------------------------------------

for (int e = 0 ; e < nbE+nbEstat ; e++) {//on assume qu'il y a au moins une esp�ce mod�lis�e, qu'elle soit dynamique ou non --> pas de if

        if (e<nbE) {
         PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppList,e)))); //esp�ce dynamique
        } else {
         PROTECT(elmt = getListElement(list, CHAR(STRING_ELT(sppListStat,e-nbE)))); //esp�ce statique
        }
//Rprintf("Eco X10\n");fichier << "EcoX10" << endl;
        //if (e<nbE) nbI = length(getListElement(elmt, "modI"));
        if (e<nbE) nbC = length(getListElement(elmt, "modC"));

        double *r_Lbio_f_sum_t_e=&NA_REAL, *r_GVLtot_f_m_e2=&NA_REAL, *r_Lbio_f_m_e=&NA_REAL, *r_P_f_m_e=&NA_REAL, *r_LD_efmc=&NA_REAL, r_theta_e=NA_REAL, *r_statLDor_efm=&NA_REAL,
               *r_statLDst_efm=&NA_REAL, r_Pst_e=NA_REAL, *r_LD_f_sum_t_e=&NA_REAL, *r_statLDor_f_sum_t_e=&NA_REAL, *r_statLDst_f_sum_t_e=&NA_REAL;
        int *dim_Lbio_e, *dim_P_e;

//Rprintf("Eco X11\n");fichier << "EcoX11" << endl;
       if ((nbE>0) & (e<nbE)) {
            SEXP gg1=R_NilValue, gg2=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue;
            PROTECT(Pgg1=VECTOR_ELT(out_L_efmct, e));
            PROTECT(Pgg2=VECTOR_ELT(out_LD_efmc, e));
            PROTECT(gg1=aggregObj(Pgg1,dimCstF));
            PROTECT(gg2=aggregObj(Pgg2,dimCstF));

            r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),41));
            r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),246));
            r_GVLst_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eVar, e),247));
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_L_efmct, e));
            r_Lbio_f_sum_t_e = REAL(gg1);
            r_LD_efmc = REAL(VECTOR_ELT(out_LD_efmc, e));
            r_LD_f_sum_t_e = REAL(gg2);
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            //r_Lref_f_e = REAL(getListElement(elmt, "Lref_f_e"));
            r_P_f_m_e = REAL(VECTOR_ELT(out_P_t, e));
            dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_L_efmct, e), install("DimCst")))));
            dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_P_t, e), install("DimCst")))));
//Rprintf("Eco X12\n");fichier << "EcoX12" << endl;
            UNPROTECT(4);
        }
        if ((nbEstat>0) & (e>=nbE)) {

            SEXP gg1=R_NilValue, gg2=R_NilValue, gg3=R_NilValue, Pgg1=R_NilValue, Pgg2=R_NilValue, Pgg3=R_NilValue;
            PROTECT(Pgg1=VECTOR_ELT(out_Lstat, e-nbE));
            PROTECT(Pgg2=VECTOR_ELT(out_statLDor_efm, e-nbE));
            PROTECT(Pgg3=VECTOR_ELT(out_statLDst_efm, e-nbE));
            PROTECT(gg1=aggregObj(Pgg1,dimCstF));
            PROTECT(gg2=aggregObj(Pgg2,dimCstF));
            PROTECT(gg3=aggregObj(Pgg3,dimCstF));

            r_GVLtot_f_m_e2 = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));
            r_GVLcom_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),8));
            r_GVLst_f_m_e_out = REAL(VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),9));
            r_Lbio_f_m_e = REAL(VECTOR_ELT(out_Lstat, e-nbE));
            r_Lbio_f_sum_t_e = REAL(gg1);
            r_statLDor_efm = REAL(VECTOR_ELT(out_statLDor_efm, e-nbE));
            r_statLDor_f_sum_t_e = REAL(gg2);
            r_statLDst_efm = REAL(VECTOR_ELT(out_statLDst_efm, e-nbE));
            r_statLDst_f_sum_t_e = REAL(gg3);
            r_theta_e = REAL(getListElement(elmt, "theta_e"))[0];
            //r_Lref_f_e = REAL(getListElement(elmt, "Lref_f_e"));
            r_P_f_m_e = REAL(VECTOR_ELT(out_Pstat, e-nbE));
            dim_Lbio_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_Lstat, e-nbE), install("DimCst")))));
            dim_P_e = INTEGER(iDim(INTEGER(getAttrib(VECTOR_ELT(out_Pstat, e-nbE), install("DimCst")))));
            r_Pst_e = REAL(getListElement(elmt, "Pst_e"))[0];
//Rprintf("Eco X13\n");fichier << "EcoX13" << endl;
            UNPROTECT(6);
        }



        //---------------------
        //�quations de la table "t"
        //---------------------

  for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

        if (e<nbE) {

          r_countLf[ind_f] = r_countLf[ind_f] + finite(r_Lbio_f_sum_t_e[ind_f]) + finite(r_LD_f_sum_t_e[ind_f]);

        } else {

          r_countLf[ind_f] = r_countLf[ind_f] + finite(r_Lbio_f_sum_t_e[ind_f]) + finite(r_statLDor_f_sum_t_e[ind_f]) + finite(r_statLDst_f_sum_t_e[ind_f]);

        }
////Rprintf("Eco X131\n");fichier << "EcoX131" << endl;
        for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

        //-- 1. GVL_f_m_e

        double countCom = 0.0;
////Rprintf("Eco X132\n");fichier << "EcoX132" << endl;

    if (e<nbE) {

             if (ISNA(r_theta_e)) r_theta_e = 1.0;

             for (int ind_c = 0 ; ind_c < (nbC-1) ; ind_c++){ //sur les classes non sous-tailles

                if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]]))
                        r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;

                if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {

               countCom = countCom +
                  r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                  r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + ind_c*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

               r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + ind_c*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                }

             }

////Rprintf("Eco X133\n");fichier << "EcoX133" << endl;
             if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]]))
                        r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;

////Rprintf("Eco X1331\n");fichier << "EcoX1331" << endl;
             if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])) {
////Rprintf("Eco X1332\n");fichier << "EcoX1332" << endl;
                  r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                   r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                   r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                   r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + (nbC-1)*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                   finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

                  r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                    r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                    finite(r_LD_efmc[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + (nbC-1)*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);
////Rprintf("Eco X1333\n");fichier << "EcoX1333" << endl;
             } else {
////Rprintf("Eco X1334\n");fichier << "EcoX1334" << endl;
                  r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;
////Rprintf("Eco X1335\n");fichier << "EcoX1335" << endl;
             }

    } else {

////Rprintf("Eco X134\n");fichier << "EcoX134" << endl;

       if (ISNA(r_theta_e)) r_theta_e = 1.0;

       if (ISNA(r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]]))
                        r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] = 0.0;


       if (!ISNA(r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])){

            countCom = r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 * //prix au kg
                  r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                  r_theta_e * r_P_f_m_e[ind_f*dim_P_e[0] + ind_m*dim_P_e[1] + 0*dim_P_e[2] + ind_t*dim_P_e[3]] * 1000 *
                  finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]) ;

            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                r_Lbio_f_m_e[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]] +
                finite(r_statLDor_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);
        }

        if (!ISNA(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]])){

            if (ISNA(r_Pst_e)) r_Pst_e = 0.0;

            r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                  r_Pst_e * 1000 * finite(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

            r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] +
                finite(r_statLDst_efm[ind_f*dim_Lbio_e[0] + ind_m*dim_Lbio_e[1] + 0*dim_Lbio_e[2] + ind_t*dim_Lbio_e[3]]);

        } else {

           r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = 0.0;

        }
////Rprintf("Eco X135\n");fichier << "EcoX135" << endl;

    }

////Rprintf("Eco X1338\n");fichier << "EcoX1338" << endl;
    r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = countCom;
////Rprintf("Eco X1339\n");fichier << "EcoX1339" << endl;
    r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
       r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] + r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];
////Rprintf("Eco X136\n");fichier << "EcoX136" << endl;

    //-- 2. GVLtot_f_m


    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
       r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
       finite(r_GVLtot_f_m_e2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]);

    double LC = 0.0, LCD = 0.0;
    if (!ISNA(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]))
       LC = r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]];
    if (!ISNA(r_lcd_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]))
       LCD = r_lcd_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]];
////Rprintf("Eco X137\n");fichier << "EcoX137" << endl;
    r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
       r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] +
       finite(r_GVLcom_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]) * (1 - 0.01*LC) +
       finite(r_GVLst_f_m_e_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]]) * (1 - 0.01*LCD);


   }
  }
  ////Rprintf("Eco X14\n");fichier << "EcoX14" << endl;

////Rprintf("aa1"); fichier << "aa1" << endl;
//if (ind_t==4 & e==26) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
////Rprintf("aa1");


  if (e<nbE) {
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,4), e, VECTOR_ELT(VECTOR_ELT(eVar, e),41));//Rprintf("Eco X1\n");
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,0), e, VECTOR_ELT(VECTOR_ELT(eVar, e),246));//Rprintf("Eco X2\n");
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,2), e, VECTOR_ELT(VECTOR_ELT(eVar, e),247));//Rprintf("Eco X3\n");
  } else {
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,5), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),1));//Rprintf("Eco X4\n");
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,1), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),8));//Rprintf("Eco X5\n");
    SET_VECTOR_ELT(VECTOR_ELT(out_EcoDCF,3), e-nbE, VECTOR_ELT(VECTOR_ELT(eStatVar, e-nbE),9));//Rprintf("Eco X6\n");
  }

  UNPROTECT(1);

  //Rprintf("\nJ3\n");fichier << "J3" << endl;

////Rprintf("aa1"); fichier << "aa1" << endl;
//if (ind_t==4 & e==26) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
////Rprintf("aa1");
//if (ind_t==4) PrintValue(VECTOR_ELT(out_EcoDCF, 13));


}

// Calcul quota costs
for (int ind_f = 0 ; ind_f < nbF ; ind_f++) //Reinitialisation si plusieurs appels au module
        r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

SEXP nam_eQuota, PQuot_et,QuotaTrade_fe;
double *r_PQuot_et, *r_QuotaTrade_fe;
for (int eQuota = 0 ; eQuota < nbEQuotaMarket ; eQuota++) {

    PROTECT(nam_eQuota = STRING_ELT(sppListQM,eQuota));
    PROTECT(PQuot_et = getListElement(out_PQuot_et, CHAR(nam_eQuota)));

//    if(!isNull(getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota)))){ //explicitely traded: amount traded = landings - holdings
//        PROTECT(QuotaTrade_fe = getListElement(out_QuotaTrade_fe, CHAR(nam_eQuota)));
//        r_QuotaTrade_fe = REAL(QuotaTrade_fe);
//    } else { //otherwise amount traded = landings
        if (!isNull (getListElement(out_L_efmit, CHAR(nam_eQuota)))){ // espece dyn
                                    PROTECT(QuotaTrade_fe = aggregObj(getListElement(out_L_efmit, CHAR(nam_eQuota)),dimCstF));
                                } else{ // espece stat
                                    PROTECT(QuotaTrade_fe = aggregObj(getListElement(out_Lstat, CHAR(nam_eQuota)),dimCstF));
                                }
        r_QuotaTrade_fe = REAL(QuotaTrade_fe);
//    }

    r_PQuot_et = REAL(PQuot_et);

    for (int ind_f = 0 ; ind_f < nbF ; ind_f++){
        r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                                                                                    r_PQuot_et[ind_t] * r_QuotaTrade_fe[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * 1000;
//        if (ind_f==6){
//            fichier << "e = " << CHAR(nam_eQuota) <<
//                        "; Pquot = " << r_PQuot_et[ind_t] <<
//                        "; Traded amount = " << r_QuotaTrade_fe[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]*1000 <<
//                        "; QuotaExp = " << r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] << endl;
//        }

    }
    UNPROTECT(3);
}


//Rprintf("aa");
//if (ind_t==4) PrintValue(VECTOR_ELT(out_EcoDCF, 13));
//Rprintf("bb");//Rprintf("%f\n",r_cnb_f_m_out[0 + nbF*7 + nbF*nbMe*4]);//Rprintf("%f\n",r_ET_f_m_out[0 + nbF*7 + nbF*nbMe*4]);
////Rprintf("%f\n",rnbv[0 + nbF*7]);//Rprintf("%f\n",rnbTrip[0 + nbF*7]);//Rprintf("%f\n",rtripLgth[0 + nbF*7]);
//if (ind_t==1) PrintValue(VECTOR_ELT(out_EcoDCF, 12));
//Rprintf("cc");

for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

    for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

////Rprintf("Eco X1336\n");fichier << "EcoX1336" << endl; //� ce moment, cnb contient les d�barquements totaux par flottille et m�tier
    r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] = finite(r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] /
        (r_ET_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * rnbv[ind_f + nbF*ind_m] *
         rnbTrip[ind_f + nbF*ind_m] * rtripLgth[ind_f + nbF*ind_m]));

////Rprintf("Eco X1337\n");fichier << "EcoX1337" << endl;//calcul du num�rateur de cnb_f

    r_cnb_f_out[ind_f + nbF*ind_t] = r_cnb_f_out[ind_f + nbF*ind_t] +
             finite(r_cnb_f_m_out[ind_f + nbF*ind_m + nbF*nbMe*ind_t] * rnbv[ind_f + nbF*ind_m] *
             rnbTrip[ind_f + nbF*ind_m] * rtripLgth[ind_f + nbF*ind_m]);


}}

//if (ind_t==1) PrintValue(VECTOR_ELT(out_EcoDCF, 13));

// --------------------------------------------------------------------------------------




    // � ce stade, plus de consid�ration d'esp�ce pour les indicateurs

        for (int ind_f = 0 ; ind_f < nbF ; ind_f++){

          double NGVLtot_f = 0.0, RTBStot_f = 0.0;

            for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

            //-- 5. GVLav_f_m

                r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                    r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

             //-- 6. GVLtot_f & NGVLav_f_m

                if (ind_m==0) {

                    if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                    } else {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
                    }

                } else {

                    if (!ISNA(r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                                r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                                    r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                    }
                }


                    if (!ISNA(r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                        NGVLtot_f = NGVLtot_f +
                         r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

                    }

                        r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] /
                        r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];



            //-- 11. rtbs_f_m

                    r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        r_NGVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] -
                        ((finite(r_ovcDCFue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) +
                        finite(r_fvolue_f_m2[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + 0*eF_fm[3]]) *
                        r_vf_f_m[ind_f*dim_vf_f_m[0] + ind_m*dim_vf_f_m[1] + 0*dim_vf_f_m[2] + ind_t*dim_vf_f_m[3]]) *
                        r_ue_f_m[ind_f*dim_ue_f_m[0] + ind_m*dim_ue_f_m[1] + 0*dim_ue_f_m[2] + ind_t*dim_ue_f_m[3]] / pow(1+0.0,ind_t));

                   if (!ISNA(r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]])) {

                        RTBStot_f = RTBStot_f +
                         r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] *
                         r_nbv_f_m[ind_f*dim_nbv_f_m[0] + ind_m*dim_nbv_f_m[1] + 0*dim_nbv_f_m[2] + ind_t*dim_nbv_f_m[3]];

                    }

                   if (perscCalc<2) {

                    r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        0.01 * r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                        r_rtbs_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                   } else if (perscCalc==5){

                       r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] =
                        0.01 * r_cshr_f_m[ind_f*dim_cshr_f_m[0] + ind_m*dim_cshr_f_m[1] + 0*dim_cshr_f_m[2] + ind_t*dim_cshr_f_m[3]] *
                        r_GVLtot_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]];

                       }else {

                    r_cshrT_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] = NA_REAL;

                   }

            } //on sort de la boucle sur les niveaux m�tiers


            //-- 7. GVLav_f

                r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];


            //-- 8.5. NGVLav_f

                r_NGVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    NGVLtot_f / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];

            //-- 10. cnb_f

                r_cnb_f_out[ind_f + nbF*ind_t] =
                    finite(r_cnb_f_out[ind_f + nbF*ind_t] / (rnbv_f[ind_f] * rtripLgth_f[ind_f] * rnbTrip_f[ind_f]));


            //-- 12. rtbs_f

                r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        RTBStot_f / r_nbv_f[ind_f*dim_nbv_f[0] + 0*dim_nbv_f[1] + 0*dim_nbv_f[2] + ind_t*dim_nbv_f[3]];



                //version actualis�e
                r_rtbsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 14. cshrT_f


        if (perscCalc==0) {  //salaires par marin fix�s

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    r_cnb_f_out[ind_f + nbF*ind_t] / r_cnb_f_out[ind_f + nbF*0];

        }

        if (perscCalc==1) {  //part �quipage constante (RAP)

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        }

        if (perscCalc==2) {  //part �quipage constante calcul�e (RAP) - ccwr

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_ccwr_f2[ind_f] *
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

        }


        if (perscCalc==3) {  //part �quipage constante (RAP) + salaire marin suppl�mentaire fix�

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                    (r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    (r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] - r_cnb_f_out[ind_f + nbF*0 + nbF*0]) /
                    r_cnb_f_out[ind_f + nbF*0 + nbF*0]);

        }


        if (perscCalc==4) {  //part �quipage constante calcul�e (RAP)- salaires marin suppl�mentaire fix�

            r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_ccwr_f2[ind_f] *
                    (r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + 0*eF_f[3]] *
                    (r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] - r_cnb_f_out[ind_f + nbF*0 + nbF*0]) /
                    r_cnb_f_out[ind_f + nbF*0 + nbF*0]);

        }

        if (perscCalc==5) {  //part �quipage constante (GVL)

                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                        0.01*r_cshr_f[ind_f*dim_cshr_f[0] + 0*dim_cshr_f[1] + 0*dim_cshr_f[2] + ind_t*dim_cshr_f[3]] *
                            r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ;

                }
             //-- 15. ncshr_f

                r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_eec_f[ind_f*dim_eec_f[0] + 0*dim_eec_f[1] + 0*dim_eec_f[2] + ind_t*dim_eec_f[3]];

             //-- 16. ocl_f

                r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_mwh_f[ind_f*dim_mwh_f[0] + 0*dim_mwh_f[1] + 0*dim_mwh_f[2] + 0*dim_mwh_f[3]] *
                    r_cnb_f_out[ind_f + nbF*0 + nbF*ind_t] * rtripLgth_f[ind_f] * rnbTrip_f[ind_f];


             //-- 17. cs_f

                r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    r_ocl_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

                //version actualis�e
                r_csAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 18. csTot_f

                r_csTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * rnbv_f[ind_f];


             //-- 19. gva_f

                r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] -
                    (r_rep_f[ind_f*dim_rep_f[0] + 0*dim_rep_f[1] + 0*dim_rep_f[2] + ind_t*dim_rep_f[3]] +
                    r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t*dim_fixc_f[3]] +
                    r_gc_f[ind_f*dim_gc_f[0] + 0*dim_gc_f[1] + 0*dim_gc_f[2] + ind_t*dim_gc_f[3]])/ pow(1+0.0,ind_t) ;//+
                    //r_QuotaExp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+0.0,ind_t) ;


                //version actualis�e
                r_gvaAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


            //-- 20. gvamargin_f

            r_gvamargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
               r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


            //-- 21. gva_FTE_f

            r_gva_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
               r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];


            //-- 22. ccw_f

            r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            if ( (perscCalc==0) | (perscCalc==1) | (perscCalc==3) ) {

            r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] + r_opersc_f2[ind_f];
            }


            //-- 23. ccwCr_f

             r_ccwCr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];


            //-- 24. wageg_f

             r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];


            //-- 25. wagen_f

             r_wagen_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_ncshr_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_cnb_f_out[ind_f + nbF*ind_t];

            //-- 26. wageg_FTE_f

             r_wageg_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
               r_cshrT_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];
               //r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_FTE_f[ind_f*dim_FTE_f[0] + 0*dim_FTE_f[1] + 0*dim_FTE_f[2] + ind_t*dim_FTE_f[3]];


            //-- 27. wageg_h_f

             r_wageg_h_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                r_wageg_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / (rtripLgth_f[ind_f] * rnbTrip_f[ind_f]);
               //r_wageg_FTE_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / (rtripLgth_f[ind_f] * rnbTrip_f[ind_f]);

            //-- 28. gp_f

                r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_ccw_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


              //version actualis�e
                r_gpAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 29. gpmargin_f

                r_gpmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

             //-- 30. ncf_f

                r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]];


             //-- 31. np_f

                r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]];

             //-- 32. npmargin_f

                r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 33. prof_f

                r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = -1.0;
                if (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]>=0) r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
                if (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]>0.1) r_prof_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 1.0;


             //-- 34. npmargin_trend_f

                r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = -1.0;
                if (ind_t>=5) {
                     double devTrend;
                     devTrend = r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                                (0.2 * (r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-5)*eF_f[3]] +
                                        r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-4)*eF_f[3]] +
                                        r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-3)*eF_f[3]] +
                                        r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-2)*eF_f[3]] +
                                        r_npmargin_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + (ind_t-1)*eF_f[3]]));
                     if (devTrend>(-0.05)) r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;
                     if (devTrend>0.05) r_npmargin_trend_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 1.0;
                }

             //-- 35. ssTot_f

                r_ssTot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * rnbv_f[ind_f];

             //-- 36. ps_f

                r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    rnbv_f[ind_f] * (r_cs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] + r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]]);

              //version actualis�e
                r_psAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_ps_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);


             //-- 37. sts_f

                r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = 0.0;

                for (int ind_m = 0 ; ind_m < nbMe ; ind_m++){

                  r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] +
                        (finite(r_lc_f_m[ind_f*dim_lc_f_m[0] + ind_m*dim_lc_f_m[1] + 0*dim_lc_f_m[2] + ind_t*dim_lc_f_m[3]]) *
                         r_GVLav_f_m_out[ind_f*eF_fm[0] + ind_m*eF_fm[1] + 0*eF_fm[2] + ind_t*eF_fm[3]] * rnbv[ind_f + ind_m*nbF]);

                }

               //version actualis�e
                 r_stsAct_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                       r_sts_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] / pow(1+dr,ind_t);



             //-- 38. ber_f

                r_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLtot_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] *
                    (r_fixc_f[ind_f*dim_fixc_f[0] + 0*dim_fixc_f[1] + 0*dim_fixc_f[2] + ind_t*dim_fixc_f[3]] +
                    r_dep_f[ind_f*dim_dep_f[0] + 0*dim_dep_f[1] + 0*dim_dep_f[2] + ind_t*dim_dep_f[3]] +
                    r_ic_f[ind_f*dim_ic_f[0] + 0*dim_ic_f[1] + 0*dim_ic_f[2] + ind_t*dim_ic_f[3]]) /
                    r_rtbs_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] ;

             //-- 39. CR_BER_f

                r_CR_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_BER_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];


             //-- 40. fuelEff_f
                double numFuelEff = 0.0;
                for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) numFuelEff = numFuelEff + finite(r_fvolue_f_m2[ind_f + ind_m*nbF] * r_ue_f_m[ind_f+ ind_m*nbF] * r_nbv_f_m[ind_f+ ind_m*nbF]);
                r_fuelEff_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = numFuelEff / r_countLf[ind_f];


            //-- 41. ratio_fvol_GVA_f
                double numFvolGVA = 0.0;
                for (int ind_m = 0 ; ind_m < nbMe ; ind_m++) numFvolGVA = numFvolGVA + finite(r_fvolue_f_m2[ind_f + ind_m*nbF] * r_ue_f_m[ind_f+ ind_m*nbF]);
                r_ratio_fvol_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = numFvolGVA / r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 42. ratio_gp_GVA_f

                r_ratio_gp_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_gva_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]];

            //-- 43. ratio_GVL_K_f

                r_ratio_GVL_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

            //-- 44. ratio_gp_K_f

                r_ratio_gp_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

            //-- 45. RoFTA_f

                r_RoFTA_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_ncf_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]] / pow(1+0.0,ind_t) );

            //-- 46. ROI_f
                r_ROI_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] = NA_REAL;
                if (finite(r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]])>0) {
                    r_ROI_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                      (r_gp_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] - r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]]) /
                      finite(r_inv_f[ind_f*dim_inv_f[0] + 0*dim_inv_f[1] + 0*dim_inv_f[2] + ind_t*dim_inv_f[3]]);
                }

            //-- 47. ratio_np_K_f

                r_ratio_np_K_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_np_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    r_K_f[ind_f*dim_K_f[0] + 0*dim_K_f[1] + 0*dim_K_f[2] + ind_t*dim_K_f[3]];



            //-- 48. ratio_GVL_cnb_ue_f

                r_ratio_GVL_cnb_ue_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] =
                    r_GVLav_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] /
                    (r_cnb_f_out[ind_f*eF_f[0] + 0*eF_f[1] + 0*eF_f[2] + ind_t*eF_f[3]] * r_ue_f[ind_f] );


        }

if (ind_t==0) UNPROTECT(67);
UNPROTECT(89);

//Rprintf("\nJ2\n");fichier << "J2" << endl;

//if (ind_t>10) //PrintValue(VECTOR_ELT(out_EcoDCF, 44));
//if (ind_t>10) //PrintValue(VECTOR_ELT(out_EcoDCF, 2));
//if (ind_t>10) //PrintValue(VECTOR_ELT(fVar, 10));
//if (ind_t>10) //PrintValue(VECTOR_ELT(fVar, 4));
//if (ind_t>10) //PrintValue(vf_f_m);
//if (ind_t>10) //PrintValue(ue_f_m);
//if (ind_t==(nbT-1)) {//PrintValue(VECTOR_ELT(fVar,4));//PrintValue(VECTOR_ELT(fVar,10));}

//fichier.close();

}}

