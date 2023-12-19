#include"opt_alg.h"
#include<vector>
solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

solution expansion(matrix(*ff)(matrix, matrix, matrix), long double x0, long double d, long double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution p;
		matrix mx0(1, 1, x0);
		matrix f_mx0(1, 1, ff(mx0, ud1, ud2)(0,0));

		matrix mx1(1, 1, x0 + d);
		matrix f_mx1(1, 1, ff(mx1, ud1, ud2)(0, 0));
		
		if (f_mx1 == f_mx0)
		{
			p.x = x0;
			p.y = x0 + d;
			return p;
		}
		else if (f_mx1 > f_mx0)
		{
			d = -d;
			mx1(0, 0) = x0 + d;
			f_mx1(0, 0) = ff(mx1, ud1, ud2)(0, 0);
			if (f_mx1 >= f_mx0)
			{
				p.x = mx1(0, 0);
				p.y = mx0(0, 0) - d;
				return p;
			}
		}
		matrix mx2(1, 1, x0 + 2 * d);
		matrix f_mx2(1, 1, ff(mx2, ud1, ud2)(0, 0));
		int i = 1;
		for (i = 1; i <= Nmax; i++)
		{
			if (f_mx1(0, 0) <= f_mx2(0, 0)) break;

			mx0(0,0) = mx1(0,0);
			f_mx0(0,0) = f_mx1(0,0);

			mx1(0,0) = mx2(0,0);
			f_mx1(0,0) = f_mx2(0,0);

			mx2(0, 0) = mx2(0, 0) + pow(alpha, i) * d;
			f_mx2(0, 0) = ff(mx2, ud1, ud2)(0, 0);
		}
		
		if (d > 0)
		{
			p.x = mx0(0, 0);
			p.y = mx2(0, 0);
			p.f_calls = i + 4;
		}
		else
		{
			p.x = mx2(0, 0);
			p.y = mx0(0, 0);
			p.f_calls = i + 4;
		}
		return p;
	}
	catch (string ex_info)
	{
		throw ("long double* expansion(...):\n" + ex_info);
	}
}

long double fib_series(int k)
{
	long double a = 1.0;
	long double b = 1.0;
	long double c = a;
	for (int i = 0; i < k; i++)
	{
		c = b;
		b = a + b;
		a = c;
	}
	return b;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, long double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		std::vector<long double> fib_num;
		int k = 1;
		for (k = 1; ; k++)
		{
			fib_num.push_back(fib_series(k));
			if (fib_num[k - 1] > ((b - a) / epsilon)) break;
		}
		long double c = b - (fib_num[k - 2] / fib_num[k - 1] )* (b - a);
		long double d = a + b - c;
		for (int i = 0; i < k - 2; i++)
		{
			//std::cout << b - a << std::endl;
			if (ff(c, ud1, ud2) < ff(d, ud1, ud2))
			{
				b = d;
			}
			else
			{
				a = c;
			}
			c = b - (fib_num[k - i - 3] /( fib_num[k - i - 2])) * (b - a);
			d = a + b - c;
		}
		Xopt.x = c;
		Xopt.y = ff(c, ud1, ud2);
		Xopt.f_calls = k;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), long double a, long double b, long double epsilon, long double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		
		long double c = a + (b - a) / 3;
		long double licznik, mianownik;
		long double d = b - (b - a) / 3;
		int i = 0;
		for (i = 0; i < Nmax ; i++)
		{
			//std::cout << b - a << std::endl;
			if (abs(b - a) < epsilon || abs(ff(b, ud1, ud2)(0, 0) - ff(a, ud1, ud2)(0, 0)) < gamma)
			{
				//std::cout << "b - a: " << abs(b - a) << " | f(b) - f(a) : " << abs(ff(b, ud1, ud2)(0, 0) - ff(a, ud1, ud2)(0, 0)) << std::endl;
				break;
			}
			//std::cout << "b - a: " << abs(b - a) << " | f(b) - f(a) : " << abs(ff(b, ud1, ud2)(0, 0) - ff(a, ud1, ud2)(0, 0)) << std::endl;
			long double f_a = ff(a, ud1, ud2)(0, 0);
			long double f_b = ff(b, ud1, ud2)(0, 0);
			long double f_c = ff(c, ud1, ud2)(0, 0);
			long double f_d = 0;
			licznik = f_a * (b * b - c * c) + f_b * (c * c - a * a) + f_c * (a * a - b * b);
			mianownik = f_a * (b - c) + f_b * (c - a) + f_c * (a - b);
			if (mianownik <= 0)
			{
				break;
			}
			d = 0.5 * licznik / mianownik;
			f_d = ff(d, ud1, ud2)(0,0);
			if (a < d && d < c)
			{
				if (f_d < f_c)
				{
		
					b = c;
					c = d;
				}
				else
				{
					a = d;
				}
			}
			else if (c < d && d < b)
			{
				if (f_d < f_c)
				{
					a = c;
					c = d;
				}
				else
				{
					b = d;
				}
			}
			
			

		}
		Xopt.x = d;
		Xopt.y = ff(d, ud1, ud2);
		Xopt.f_calls = i*8;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double s, long double alpha, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		matrix x_base = x0;
		matrix ff_base = ff(x0, ud1, ud2);
		matrix x_next(2, 1);
		matrix ff_next = ff_base;
		matrix x_temp(2,1);

		matrix direction(2, 1);
		bool explore_new_direction = true;
		int itr = 0;
		int additional_f_calls = 0;
		std::ofstream file;
		file.open("HJ.txt");
		while (s > epsilon && itr < Nmax)
		{
			file << x_base(0, 0) << '\t' << x_base(1, 0) << std::endl;
			itr++;
			for (int i = -1; i <= 1 && explore_new_direction; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					x_temp(0, 0) = x_base(0, 0) + s * j;
					x_temp(1, 0) = x_base(1, 0) + s * i;
					if (ff(x_temp, ud1, ud2) < ff_next)
					{
						x_next = x_temp;
						ff_next = ff(x_next, ud1, ud2);
						additional_f_calls++;
					}
					additional_f_calls++;
				}
			}
			if (explore_new_direction && ff_next < ff_base)
			{
				direction(0, 0) = x_next(0, 0) - x_base(0, 0);
				direction(1, 0) = x_next(1, 0) - x_base(1, 0);

				x_base = x_next;
				ff_base = ff_next;

				explore_new_direction = false;
				continue;
			}
			else if (explore_new_direction && ff_next >= ff_base)
			{
				s *= alpha;
				continue;
			}

			x_next = x_next + direction;
			ff_next = ff(x_next, ud1, ud2);

			if (ff(x_next, ud1, ud2) < ff_base)
			{
				x_base = x_next;
				ff_base = ff_next;
			}
			else
			{
				s *= alpha;
				explore_new_direction = true;
			}
		}
		file.close();
		Xopt.x = x_base;
		Xopt.y = ff_base;
		Xopt.f_calls = itr * 2 + additional_f_calls;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, long double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

matrix new_coordinate_system(matrix d, int lambda[])
{
	matrix Q(2, 2);
	Q(0, 0) = lambda[0];
	Q(1, 0) = lambda[1];
	Q(0, 1) = 0;
	Q(1, 1) = lambda[1];

	Q = d * Q;

	matrix v(2, 2);
	v(0, 0) = Q(0, 0);
	v(1, 0) = Q(1, 0);

	v(0, 1) = Q(0, 1) - (Q(0, 1) * v(0, 1) + Q(1, 1) * v(1, 1)) * d(0, 1);
	v(1, 1) = Q(1, 1) - (Q(0, 1) * v(0, 1) + Q(1, 1) * v(1, 1)) * d(1, 1);

	d(0, 0) = v(0, 0) / sqrt(v(0, 0) * v(0, 0) + v(1, 0) * v(1, 0));
	d(1, 0) = v(1, 0) / sqrt(v(1, 0) * v(1, 0) + v(0, 0) * v(0, 0));
	d(1, 1) = v(1, 1) / sqrt(v(1, 1) * v(1, 1) + v(0, 1) * v(0, 1));
	d(0, 1) = v(0, 1) / sqrt(v(1, 1) * v(1, 1) + v(0, 1) * v(0, 1));

	return d;
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, long double alpha, long double beta, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//int n = get_dim(x0);
		solution Xopt;
		Xopt.f_calls = 0;
		int i = 0;
		matrix xB(x0);
		matrix x_test(x0);
		//matrix(1,n);
		matrix d(2, 2);
		d(0, 0) = 1;
		d(1, 0) = 0;
		d(0, 1) = 0;
		d(1, 1) = 1;
		int lambda[2] = { 0,0 };
		int p[2] = { 0,0 };
		matrix s = s0;
		std::ofstream file;
		file.open("Rosen.txt");
		do
		{
			file << xB(0, 0) << '\t' << xB(1, 0) << std::endl;
			for (int j = 0; j < 2; j++)
			{
				if (ff(xB + s(j, 0) * d[j], ud1, ud2) < ff(xB, ud1, ud2))
				{
					Xopt.f_calls += 2;
					xB = xB + s(j, 0) * d[j];
					lambda[j] = lambda[j] + s(j, 0);
					s(j, 0) = alpha * s(j, 0);
				}
				else
				{
					s(j, 0) = -beta * s(j, 0);
					xB = xB + s(j, 0) * d[j];
					p[j] = p[j] + 1;
				}
			}
			i++;
			//xB(i, 0) = x;??? PROBLEMY ?????
			if (lambda[0] != 0 && lambda[1] != 0 && p[0] != 0 && p[1] != 0)
			{
				d = new_coordinate_system(d, lambda);
				lambda[0] = 0;
				lambda[1] = 0;
				p[0] = 0;
				p[1] = 0;
				s = s0;
			}
			if (i >= Nmax)break;
			
		} while (abs(s(0)) > epsilon || abs(s(1)) > epsilon);
		file.close();
		Xopt.x = xB;
		Xopt.y = ff(xB, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double c, long double dc, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

int max_vertex(solution* Xopt, int dim)
{
	int max = 0;
	for (int i = 1; i < dim; i++)
		if (Xopt[max].y < Xopt[i].y)
			max = i;
	return max;
}

int min_vertex(solution* Xopt, int dim)
{
	int min = 0;
	for (int i = 1; i < dim; i++)
		if (Xopt[min].y > Xopt[i].y)
			min = i;
	return min;
}

matrix mass_center(solution* Xopt, int dim, int max)
{
	matrix x(dim, 1);
	int mod = 0;
	for (int i = 0; i < dim + 1; i++)
	{
		if (i == max) continue;
		x = x + Xopt[i].x;
	}
	x = x / dim;
	return x;
}


solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double s, long double alpha, long double beta, long double gamma, long double delta, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int* size = get_size(x0);
		int dim = size[0];
		delete[] size;

		solution* Xopt = new solution[dim + 1];
		matrix e(dim, 1, 0.0f);

		for (int i = 0; i < dim; i++)
		{
			e(i, 0) = 1;
			Xopt[i].x = x0 + s*e;
			Xopt[i].y = ff(Xopt[i].x, ud1, ud2); solution::f_calls++;
			//e(i, 0) = 0;
		}
		Xopt[dim].x = x0;
		Xopt[dim].y = ff(x0, ud1, ud2);
		int min, max;
		matrix p_center;
		matrix p_odb;
		matrix p_e;
		matrix ff_min;
		matrix ff_max;
		matrix ff_odb;

		for (int itr = 0; itr < Nmax; itr++)
		{
			max = max_vertex(Xopt, dim);
			min = min_vertex(Xopt, dim);
			if (norm(Xopt[max].x - Xopt[min].x) < epsilon) break;
			//std::cout << itr << std::endl;
			p_center = mass_center(Xopt, dim, max);
			p_odb = p_center + alpha*(p_center - Xopt[max].x);

			ff_min = Xopt[min].fit_fun(ff, ud1, ud2);
			ff_max = Xopt[max].fit_fun(ff, ud1, ud2);
			ff_odb = ff(p_odb, ud1, ud2); solution::f_calls++;

			if (ff_odb < ff_min)
			{
				p_e = p_center + gamma * (p_odb - p_center);
				matrix ff_e = ff(p_e, ud1, ud2); solution::f_calls++;
				if (ff_e < ff_odb)
				{
					Xopt[max].x = p_e;
					Xopt[max].y = ff_e;
				}
				else
				{
					Xopt[max].x = p_odb;
					Xopt[max].y = ff_odb;
				}
			}
			else if (ff_min <= ff_odb && ff_odb < ff_max)
			{
				Xopt[max].x = p_odb;
				Xopt[max].y = ff_odb;
			}
			else
			{
				matrix p_z = p_center + beta * (Xopt[max].x - p_center);
				matrix ff_z = ff(p_z, ud1, ud2); solution::f_calls++;
				if (ff_z >= ff_max)
				{
					for (int i = 0; i < dim; i++)
					{
						if (i == min) continue;
						Xopt[i].x = delta * (Xopt[i].x + Xopt[min].x);
						Xopt[i].y = ff(Xopt[i].x, ud1, ud2); solution::f_calls++;
					}
				}
				else
				{
					Xopt[max].x = p_z;
					Xopt[max].y = ff_z;
				}
			}

		}
		min = min_vertex(Xopt, dim);
		return Xopt[min];
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution Xopt_n;
		Xopt.x = x0;
		matrix d(2,1);
		double h = h0;
		double a = 0;
		double b = 1;

		for (int i = 0; i < Nmax; i++)
		{
			d = -Xopt.grad(gf, ud1, ud2);
			if (h0 == -1)
			{
				//b /= (1.0 + (double)i / (double)Nmax);
				h = golden(ff, Xopt.x, d, a, b, epsilon, Nmax, ud1, ud2).x(0, 0);
			}
			Xopt_n.x = Xopt.x + h * d;
			if (norm(Xopt_n.x - Xopt.x) < epsilon) break;
			Xopt = Xopt_n;
		}
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution Xopt_n;
		Xopt.x = x0;
		matrix d = -gf(Xopt.x, ud1, ud2);

		double h = h0;

		double a = 0;
		double b = 1;
		if(h != -1)
			Xopt_n.x = x0 + d*h;
		else
			Xopt_n.x = x0 + d * 0.05;

		double beta = pow(norm(Xopt_n.grad(gf, ud1, ud2)) / norm(Xopt.grad(gf, ud1, ud2)), 2);
		Xopt.x = Xopt_n.x;
		for (int i = 0; i < Nmax; i++)
		{
			d = -Xopt.grad(gf, ud1, ud2) + beta * d;
			if (h0 == -1)
			{
				b /= (1.0 + (double)i / (double)Nmax);
				h = golden(ff, Xopt.x, d, a, b, epsilon, Nmax, ud1, ud2).x(0, 0);
			}
			Xopt_n.x = Xopt.x + h * d;

			if (norm(Xopt_n.x - Xopt.x) < epsilon) break;

			beta = pow(norm(Xopt_n.grad(gf, ud1, ud2)) / norm(Xopt.grad(gf, ud1, ud2)), 2);
			Xopt = Xopt_n;
		}

		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, long double h0, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution Xopt_n;
		Xopt.x = x0;
		matrix d = 0;
		double h = h0;

		double a = 0;
		double b = 1;

		for (int i = 0; i < Nmax; i++)
		{
			d = -det(inv(Xopt.hess(Hf, ud1, ud2))) * Xopt.grad(gf, ud1, ud2);
			if (h0 == -1)
			{
				//b /= (1.0 + (double)i / (double)Nmax);
				h = golden(ff, Xopt.x, d, a, b, epsilon, Nmax, ud1, ud2).x(0, 0);
			}
			Xopt_n.x = Xopt.x + h * d;
			if (norm(Xopt_n.x - Xopt.x) < epsilon) break;
			Xopt = Xopt_n;
		}

		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix grad, long double a, long double b, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		double alpha = (pow(5, 0.5) - 1) / 2;
		double c = b - alpha * (b - a);
		double d = a + alpha*(b-a);

		// Repeat the following steps until the convergence criterion is met
		while (abs(b - a) > epsilon) {
			if (ff(x0 + grad*c, ud1, ud2) < ff(x0 + grad*d, ud1, ud2)) 
			{
				c = b - alpha * (b - a);
				a = a;
				b = d;
				d = c;
			}
			else {
				a = c;
				b = b;
				c = d;
				d = a + alpha * (b - a);
			}

			i++;
			if (i > Nmax) break;
		}

		// Calculate the optimal solution
		Xopt.x = (a + b) / 2;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, long double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
