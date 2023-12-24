#include"user_funs.h"

#define PI 3.14159265359

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new long double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	long double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(long double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	long double m = 1, l = 0.5, b = 0.5, g = 9.81;
	long double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}

matrix df1(long double t, matrix y, matrix ud1, matrix ud2)
{
	long double Pa = 0.7;
	long double a = -0.98;
	long double b = 0.63;
	long double Db = 36.5665 / 10e4;
	long double Pb = 1;
	long double Fin = 0.01;
	long double Tin = 10;
	long double Ta = 90;

	matrix result(3, 1);
	long double Fa_out = y(0) > 0 ?  a * b * m2d(ud2) * sqrt(2 * 9.81 * y(0) / Pa) : 0;
	long double Fb_out = y(1) > 0 ?  a * b * Db * sqrt(2 * 9.81 * y(1) / Pb) : 0;
	result(0) = Fa_out;
	result(1) = Fb_out - Fa_out + Fin;
	result(2) = Fin / y(1) * (Tin - y(2)) - Fa_out / y(1) * (Ta - y(2));

	return result;
}

matrix df2(long double t, matrix y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	long double b = 0.5;
	long double l = 0.6;
	long double mc = 9.5;
	long double mr = 1;
	long double I =  l * l * ((1 / 3) * mr + mc);
	dY(0, 0) = y(1);
	dY(1, 0) = (ud2(0) * (ud1(0) - y(0)) + ud2(1) * (ud1(1) - y(1)) - 0.5 * y(1))/I;

	return dY;
}

matrix ff_lab1(matrix x, matrix ud1, matrix ud2)
{
	matrix y0 = matrix(3, new long double[3] {5, 1, 10});
	matrix* y = solve_ode(df1, 0, 1, 1000, y0, ud1, x);
	int n = get_len(y[0]);
	long double max = y[1](0, 2);
	for (int i = 0; i < 1000; i++)
	{
		if (max < y[1](i, 2))
			max = y[1](i, 2);
	}
	matrix res = abs(max - 50);
	return res;
}

matrix ff_lab2(matrix x, matrix ud1, matrix ud2)
{
	matrix y0 = ud1;
	long double dt = 0.01;

	matrix y_ref(2, 1);
	y_ref(0, 0) = 3.14;
	y_ref(1, 0) = 0.0;

	matrix y_0(2, 1);
	y_0(0, 0) = 0;
	y_0(1, 0) = 0;

	matrix* y_solver = solve_ode(df2, 0.0f, dt, 100, y_0, y_ref, x);
	matrix y;
	int n = get_len(y_solver[0]);
	for (int i = 0; i < n; i++)
	{
		y = y + 10 * pow(y_ref(0) - y_solver[1](i, 0), 2) + pow(y_ref(1) - y_solver[1](i, 1), 2) + pow(x(0) * (y_ref(0) - y_solver[1](i, 0)) + x(1) *
			(y_ref(1) - y_solver[1](i, 1)), 2) ;
	}
	delete[] y_solver;
	y = y * dt;
	return y;
}

matrix testowa_lab_2(matrix x, matrix ud1, matrix ud2)
{
	return x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) - cos(2.5 * 3.14 * x(0, 0)) - cos(2.5 * 3.14 * x(1, 0)) + 2;
}


matrix test_second(matrix x, matrix ud1, matrix ud2)
{
	return ( -1)*cos(0.1 * x(0, 0))* exp((-1)* (0.1 * x(0, 0) - 2 * 3.14) * (0.1 * x(0, 0) - 2 * 3.14)) + 0.002 * 0.1 * x(0, 0) * 0.1 * x(0, 0);
}

matrix testowa_lab_3_zew(matrix x, matrix ud1, matrix ud2)
{
	matrix y = sin(PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2))) / (PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2)));
	if (-x(0) + 1 > 0)
	{
		y = y + ud2 * pow(-x(0) + 1, 2);
	}
	if (-x(1) + 1 > 0)
	{
		y = y + ud2 * pow(-x(0) + 1, 2);
	}
	if (norm(x) - ud1 > 0)
	{
		y = y + ud2*pow(norm(x), 2);
	}
	return y;
}

matrix testowa_lab_3_clear( matrix x, matrix ud1, matrix ud2)
{
	matrix y = sin(PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2))) / (PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2)));

	return y;
}

matrix testowa_lab_3_wew(matrix x, matrix ud1, matrix ud2)
{
	matrix y = sin(PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2))) / (PI * sqrt(pow(x(0) / PI, 2) + pow(x(1) / PI, 2)));
	if (-x(0) + 1 > 0)
	{
		y = y + 10e4;
	}
	else y = y + ud2(0, 0) / (1 - x(0));

	if (-x(1) + 1 > 0)
	{
		y = y + 10e4;
	}
	else y = y + ud2(0, 0) / (1 - x(1));
	if (norm(x) - ud1 > 0)
	{
		y = y + 10e4;
	}
	else y = y + ud2(0, 0) / (norm(x) - ud1);
	return y;
}

matrix testowa_lab_3_wew_clear(matrix y, matrix x, matrix ud1, matrix ud2)
{
	if (-x(0) + 1 < 0)
	 y = y - ud2(0, 0) / (1 - x(0));

	if (-x(1) + 1 < 0)
	y = y - ud2(0, 0) / (1 - x(1));

	if (norm(x) - ud1 < 0)
	y = y - ud2(0, 0) / (norm(x) - ud1);

	return y;
}

matrix testowa_lab_4(matrix x, matrix ud1, matrix ud2)
{
	return pow(x(0, 0) + 2 * x(1, 0) - 7, 2) + pow(2 * x(0, 0) + x(1, 0) - 5, 2);
}

matrix grad_testowa_lab_4(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 1);
	y(0, 0) = 2.0 * (x(0, 0) + 2.0 * x(1, 0) - 7.0) + 4.0 * (2.0 * x(0, 0) + x(1, 0) - 5.0);
	y(1, 0) = 4.0 * (x(0, 0) + 2.0 * x(1, 0) - 7.0) + 2.0 * (2.0 * x(0, 0) + x(1, 0) - 5.0);
	//double normY = norm(y);
	//y = y / normY;
	
	return y;
}

matrix hess_testowa_lab_4(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 2);
	y(0, 0) = 10;	y(1, 0) = 8;
	y(0, 1) = 8;	y(1, 1) = 10;
	return y;
}

matrix df3(long double t, matrix y, matrix ud1, matrix ud2)
{
	matrix dy(4, 1);

	long double C = 0.47;
	long double ro = 1.2;
	long double r = 0.12;
	long double S = PI * r * r;

	long double Fmx = ro * y(3) * ud2(0) * r * r * r * PI;
	long double Dx = 0.5 * C * ro * S * y(1) * y(1);
	long double Fmy = ro * y(1) * ud2(0) * r * r * r * PI;
	long double Dy = 0.5 * C * ro * S * y(3) * y(3);
	long double m = 0.6;
	long double g = 9.81;

	dy(0) = y(1);
	dy(1) = (-1) * (Fmx + Dx) / m;
	dy(2) = y(3);
	dy(3) = (-1) * (Fmx + Dx + m * g) / m;

	return dy;
}

matrix ff_lab3(matrix x, matrix ud1, matrix ud2)
{
	matrix y0(4, 1);
	y0(0) = 0;
	y0(1) = x(0);
	y0(2) = 100;
	y0(3) = 0;

	matrix* y = solve_ode(df3, 0, 0.01, 7, y0, ud1, x(1));
	int i0 = 0;
	int i50 = 0;
	int n = get_len(y[0]);
	for (int i = 0; i < n; i++)
	{
		if (abs(y[1](i, 2)) < abs(y[1](i0, 2)))
		{
			i0 = i;
		}
		if (abs(y[1](i, 2) - 50) < abs(y[1](i50, 2) - 50))
		{
			i50 = i;
		}
	}
	matrix res = y[1](i0, 0);
	

	if (abs(x(0)) - 10 > 0)
		res = res - ud2 * pow(abs(x(0)) - 10, 2);
	if(abs(x(1)) - 23 > 0)
		res = res - ud2 * pow(abs(x(1)) - 23, 2);
	if( abs(y[1](i50,0) - 5 ) - 1 > 0)
		res = res - ud2 * pow(abs(y[1](i50, 0) - 5) - 1, 2);
	delete[] y;

	return res*(-1);
}

double h0(matrix x, matrix theta)
{
	return 1 / (1 + exp(-(trans(theta) * x)(0, 0)));
}


matrix ff_lab4(matrix theta, matrix y, matrix x)
{
	matrix J(1, 1, 0.0);

	int* sizeX = get_size(y);
	int m = sizeX[1];
	delete [] sizeX;

	matrix xi(3, 1);
	for (int i = 0; i < m; i++)
	{
		xi(0, 0) = x(0, i);
		xi(1, 0) = x(1, i);
		xi(2, 0) = x(2, i);

		J = J + y(0, i) * log(h0(xi, theta)) + (1 - y(0, i)) * log(1 - h0(xi, theta));
	}
	J = -J * 1.0/(double)m;

	return J;
}

matrix df_lab4(matrix theta, matrix y, matrix x)
{
	matrix dJ(3, 1, 0.0);

	int* sizeX = get_size(y);
	int m = sizeX[1];
	delete[] sizeX;

	matrix xi(3, 1);
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < m; i++)
		{
			xi(0, 0) = x(0, i);
			xi(1, 0) = x(1, i);
			xi(2, 0) = x(2, i);

			dJ(j, 0) += (h0(xi, theta) - y(0, i)) * x(j, i);
		}
	}

	dJ = dJ / (double)m;

	return dJ;
}