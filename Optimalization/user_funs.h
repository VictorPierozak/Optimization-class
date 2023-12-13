#pragma once

#include"ode_solver.h"


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(long double, matrix, matrix = NAN, matrix = NAN);

matrix df1(long double, matrix, matrix, matrix);
matrix df2(long double t, matrix y, matrix ud1, matrix ud2);
matrix df3(long double t, matrix y, matrix ud1, matrix ud2);

matrix ff_lab1(matrix, matrix, matrix);
matrix ff_lab2(matrix, matrix, matrix);
matrix ff_lab3(matrix x, matrix ud1, matrix ud2);

matrix test_second(matrix, matrix = NAN, matrix = NAN);

matrix testowa_lab_2(matrix, matrix, matrix);
matrix testowa_lab_3_zew(matrix x, matrix ud1, matrix ud2);
matrix testowa_lab_3_wew(matrix x, matrix ud1, matrix ud2);
matrix testowa_lab_3_wew_clear(matrix y, matrix x, matrix ud1, matrix ud2);
matrix testowa_lab_3_clear(matrix x, matrix ud1, matrix ud2);