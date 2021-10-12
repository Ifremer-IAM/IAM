// #include <stdlib.h>
// #include <stdio.h>
// #include <time.h> //
// #include <vector>
// #include <math.h>
// #include <string> //
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

// function used in mod_Mortalite.cpp
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

// TODO begin remove
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


// TODO : if MinimizeF removed, remove this function
// TODO : numerical recipes, fonction recup autre part.
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

// TODO : if amoeba removec, remove this function
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


#define MP 4 // TODO : used in MinimzeF, could be removed.
#define NP 3 // TODO : idem
#define FTOL 0.0000001 // TODO : remove because unused

extern "C" {
// TODO : remove this function ?
double BioEcoPar::func(double *x)
{
	return ((x[1]-23.14)*(x[1]-23.14) + (x[2]-0.256)*(x[2]-0.256) + (x[3]+17.45)*(x[3]+17.45));
}

}



// TODO : used in constr, but don't understand why. 800l for nothing ?
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
    CatchDL(listTempP, IND_T, eVarCopy, 0);
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



// TODO : used in constr, but don't understand why. 1kl for nothing ?
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

    }

    ZoptSS3 = false;
	return 0;

//fichier.close();
}

}


// TODO : remove this function
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

