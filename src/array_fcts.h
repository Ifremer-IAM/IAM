#ifndef ARRAYS_H_INCLUDED
#define ARRAYS_H_INCLUDED

//indices multipliateurs pour les concordances de dimensions
SEXP iDim(int *dimInput);

//fontion d'agregation d'un objet R accompagne de son attribut 'DimCst'
SEXP aggregObj(SEXP object, SEXP newDim);

#endif //ARRAYS_H_INCLUDED