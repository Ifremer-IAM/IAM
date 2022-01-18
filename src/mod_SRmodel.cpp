// #include <stdlib.h>
// #include <stdio.h>
// #include <time.h>
// #include <vector>
// #include <math.h>
// #include <string>
// #include <sstream>
// #include <iostream>
// #include <fstream>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
// #include <Rcpp.h>

#include "BioEcoPar.h" // Class is defined in this file.
// #include "Modules.h" // Contient tout les modules

//using namespace Rcpp;
// using namespace std;

//---------------------------------
//
// Module de mod�lisation de relations S/R
//
//---------------------------------



extern "C" {

void BioEcoPar::SRmod(SEXP list, SEXP listSR, int ind_t, SEXP TypeSR, int *srind)
        //list : liste des parametres d'entree;
        //listSR : liste des parametres a,b&c du modele SR + e.t bruit normal ou lognormal + type de bruit (1=normal, 2=lognormal) (un vecteur de longueur 5 par espece modelisee contenant des "doubles")
        //type : type de relation Stock-Recrutement : (liste de longueur "nb d'especes modelisees" contenant des entiers)
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
    if (((typeSR!=1) & (ind_t<fstAge)) | (ind_t==0)) {  //deuxi�me condition : si t initial et recrutement deduit de la ssb de la meme annee, on part des param�tres initiaux et non de la relation SR

        rans[ind_t] = NA_REAL;

    } else {

    switch (typeSR) {

        case 1 :

        rans[ind_t] = paramet[0*nbT + ind_t]; break;  //modif MM 27/08/2013 pour donner la possibilit� de d�finir plusieurs param�trages pour la relation SR au cours de la simu
                                                        //indices : 0 --> 0*nbT + ind_t
        case 2 :

        if (ssb[ind_t-fstAge] <= paramet[1*nbT + ind_t]){
            rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge];
        }else{
            rans[ind_t] = paramet[0*nbT + ind_t] * paramet[1*nbT + ind_t]; 
        } break;

        case 3 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] / (paramet[1*nbT + ind_t] + ssb[ind_t-fstAge]); break;

        case 4 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] * exp(-1.0 * paramet[1*nbT + ind_t] * ssb[ind_t-fstAge]); break;

        case 5 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge] / (1 + pow(ssb[ind_t-fstAge] / paramet[1*nbT + ind_t] , paramet[2*nbT + ind_t])); break;

        case 6 :

        if (ssb[ind_t-fstAge] <= (paramet[1*nbT + ind_t]*(1-paramet[2*nbT + ind_t]))){

            rans[ind_t] = paramet[0*nbT + ind_t] * ssb[ind_t-fstAge];

        } else {

            if (ssb[ind_t-fstAge] >= (paramet[1*nbT + ind_t]*(1+paramet[2*nbT + ind_t]))){

                rans[ind_t] = paramet[0*nbT + ind_t] * paramet[1*nbT + ind_t];

            }else{

                rans[ind_t] = paramet[0*nbT + ind_t]*(ssb[ind_t-fstAge] - (pow(ssb[ind_t-fstAge] - paramet[1*nbT + ind_t]*(1-paramet[2*nbT + ind_t]),2.0)/(4*paramet[1*nbT + ind_t]*paramet[2*nbT + ind_t])));
            }
        } break;

        case 7 :

        rans[ind_t] = paramet[0*nbT + ind_t] * ( ssb[ind_t-fstAge] + sqrt( pow(paramet[1*nbT + ind_t],2) + 0.001) - sqrt( pow(ssb[ind_t-fstAge] - paramet[1*nbT + ind_t],2) + 0.001)); break;

        default :

        rans[ind_t] = NA_REAL;

    }
    //il ne reste plus qu'a ajouter le bruit blanc issue de N(0,sigma) avec sigma = paramet[3]
    double v_alea = 0.0;
GetRNGstate();
        if (!ISNA(paramet[3*nbT + ind_t])) v_alea = rnorm(0.0,paramet[3*nbT + ind_t]); //////Rprintf("%f ",v_alea);//Rprintf("%f ",rnorm(0.0,0.157));
PutRNGstate();
    if (paramet[4*nbT + ind_t]==1){  //bruit de loi normale
        rans[ind_t] = rans[ind_t] + v_alea;
    }else{                //bruit de loi lognormale
        rans[ind_t] = rans[ind_t] * exp(v_alea);
    }
    }
    
    if (ind_t==0) UNPROTECT(1);

    }}

    if (ind_t==0) UNPROTECT(1);

}}
