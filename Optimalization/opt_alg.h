//Ten plik nie powinien byæ edytowany

#pragma once

#include"solution.h"

//unsigned int counter = 0;

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN);

solution expansion(matrix(*ff)(matrix, matrix, matrix), long double x0, long double d, long double alpha, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution fib(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, long double epsilon, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution lag(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, long double epsilon, long double gamma, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double s, long double alpha, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, long double s, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, long double alpha, long double beta, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double c, long double dc, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double s, long double alpha, long double beta, long double gamma, long double delta, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix(*Hf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
solution golden(matrix(*ff)(matrix, matrix, matrix), matrix x0,matrix d,long double a, long double b, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, long double epsilon, int Nmax, matrix ud1 = NAN, matrix ud2 = NAN); // throw (string);
