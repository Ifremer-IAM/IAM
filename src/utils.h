#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

SEXP getListElement(SEXP list, const char *str);

int getListIndex(SEXP list, const char *str);

int getVectorIndex(SEXP vect, const char *str);

#endif // UTILS_H_INCLUDED