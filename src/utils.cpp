#include <Rdefines.h>
#include <Rmath.h>

/*
*-------------------------------------------
* access one element of a given SEXP list 
* (list, str {char} = name of list element)
* return SEXP list element.
*-------------------------------------------
*/
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) {
            elmt = VECTOR_ELT(list, i);
            break;
        }

    return elmt;
}

int getListIndex(SEXP list, const char *str)
{
    SEXP names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
        if (strcmp(CHAR(STRING_ELT(names,i)), str) == 0) break;

    return i;
}

int getVectorIndex(SEXP vect, const char *str) // work on named vector
{
    int i;
    for (i = 0; i < length(vect); i++)
        if (strcmp(CHAR(STRING_ELT(vect,i)), str) == 0) break;

    return i;
}