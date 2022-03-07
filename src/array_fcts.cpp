#include <Rdefines.h>

/*
* Regroup mutliple function to work on 4 dimension array.
* These functions are not dependant on BioEcoPar class !
*/





/*
*------------------------------------------
* fonction de calcul de multiplicateurs d'indices selon les dimensions d'un objet 'array'
* (permet la genericite des equations en assurant la compatibilite des variables en presence,
*  quelles que soient leurs dimensions respectives)
* INPUT : attribut 'DimCst' de l'objet en question
*------------------------------------------
*/

SEXP iDim(int *dimInput) {

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



/*
*------------------------------------------
* fonction d'agregation d'un objet attribue type ('object'), en fonction d'un nouveau vecteur dimension DimCst ('newDim')
* NB : toutes les valeurs de 'newDim' doivent etre au plus egales aux dimensions correspondantes de l'objet pour que la fonction s'applique
*------------------------------------------
*/

SEXP aggregObj(SEXP object, SEXP newDim)
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
    } 

    //Rprintf("Dans aggregobj dim[0] %i ndim[0] %i test %d ;  dim[1] %i ndim[1] %i test %d \n", dim[0],ndim[0],dim[0]<ndim[0],dim[1],ndim[1],dim[1]<ndim[1])
    if ((dim[0]<ndim[0]) | (dim[1]<ndim[1]) | (dim[2]<ndim[2]) | (dim[3]<ndim[3])) {
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
    for (int i = 0 ; i < nbCell ; i++){
        rans[i] = 0.0;
    }

    if (nbDim>0) {
        //en-t�tes
        PROTECT(Dim = allocVector(INTSXP,nbDim));
        rdim = INTEGER(Dim);
        PROTECT(dimnames = allocVector(VECSXP,nbDim));
        for (int i = 0 ; i < 4 ; i++) {
            if (ndim[i]>0) {
                if (GET_DIMNAMES(object)!=R_NilValue){
                    SET_VECTOR_ELT(dimnames, incr, VECTOR_ELT(GET_DIMNAMES(object), incr2)) ;
                }
                rdim[incr] = ndim[i] ;
                incr++;}
            if (dim[i]>0){
                incr2++;
            }
        }

        setAttrib(ans, R_DimSymbol, Dim);
        if (GET_DIMNAMES(object)!=R_NilValue){
            setAttrib(ans, R_DimNamesSymbol, dimnames);
        }
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