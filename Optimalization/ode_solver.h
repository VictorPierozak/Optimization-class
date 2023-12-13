//Ten plik nie powinien byæ edytowany

#pragma once

#include"matrix.h"
#include"user_funs.h"

matrix* solve_ode(matrix(*)(long double, matrix, matrix, matrix), long double, long double, long double, matrix, matrix = NAN, matrix = NAN); // throw (string);
