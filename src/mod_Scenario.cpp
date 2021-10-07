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
