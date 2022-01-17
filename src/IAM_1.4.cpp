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


//extern "C" : pour eviter le "name mangling process" qui renomme les fonctions export�es dans les dll.


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

    for (i = 0; i < length(list); i++) // TODO : why is list here ?
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
// fonction d'agregation d'un objet attribue type ('object'), en fonction d'un nouveau vecteur dimension DimCst ('newDim')
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



// TODO : only call is commented
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

