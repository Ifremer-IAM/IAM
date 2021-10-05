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



// Numerical Recipes //----------------------------------------------------------------------------------------

// --------  d�termination racine (unidimensionnel)

// TODO beign remove
// TODO : remove this function
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

// TODO : remove this function
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

// TODO : end remove

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

// TODO : remove this function
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

