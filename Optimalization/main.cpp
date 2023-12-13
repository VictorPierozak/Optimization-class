/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	long double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	//a(0) = -1;
	//a(1) = 2;
	//opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	//cout << opt << endl << endl;
	//solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-3;
	lb = 0;
	ub = 5;
	long double teta_opt = 1;
	opt = fib(ff_lab1, 1, 100, epsilon, lb, ub );
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new long double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	std::ofstream file;

	solution opt;
	matrix ud1; // niewykorzystywane
	matrix ud2; //niewykorzystywane
	int max_iter = 1000;
	srand(time(NULL));

	//Testowa funkcja celu

	file.open("lab_1.txt");

	//for(long double alpha = 1.1; alpha <= 1.5; alpha+= 0.2)
	//{
	//	for (int i = 0; i < 100; i++)
	//	{
	//		cout << "Iteracja: " << i << endl;
	//		long double initialGuess = std::rand() % 200 - 100;
	//		file << initialGuess << "\t";
	//		solution range = expansion(test_second, initialGuess, 5, alpha, max_iter);
	//		file << range.x(0, 0) << "\t" << range.y(0, 0) << "\t" << range.f_calls << "\t";

	//		opt = lag(test_second, range.x(0,0), range.y(0, 0), 0.1, 0.0001, max_iter);
	//		file << opt.x << "\t" << opt.y << "\t" << opt.f_calls;
	//		if (opt.x < 50) file << "\tLOKALNE\t";
	//		else file << "\tGLOBALNE\t";

	//		opt = fib(test_second, range.x(0, 0), range.y(0, 0), 0.0001);
	//		file << opt.x << "\t" << opt.y << "\t" << opt.f_calls;
	//		if (opt.x < 50) file << "\tLOKALNE\t";
	//		else file << "\tGLOBALNE\t";
	//		file << endl;
	//	}
	//}
	file.close();
	
	// Problem rzeczywisty

	solution optFib = fib(ff_lab1,10e-4, 10e-2, 10e-8);
	cout << "Fibonacci: " << optFib << endl;
	solution optLag = lag(ff_lab1, 10e-4, 10e-2, 10e-8, 0.01, 1000);
	cout << "Lagrange: " << optLag << endl;

	// Symulacja Fibonacci

	file.open("symulacja_fib.txt");
	matrix x(optFib.x);
	matrix y0 = matrix(3, new long double[3] {5, 1, 10});
	matrix* y = solve_ode(df1, 0, 1, 1000, y0, ud1, x);
	int n = get_len(y[0]);
	long double max = y[1](0, 2);
	for (int i = 0; i < 1000; i++)
	{
		file << y[1](i, 0) << "\t" << y[1](i, 1) << "\t" << y[1](i, 2) << endl;
	}
	file.close();

	// Symulacja Lagrange

	file.open("symulacja_lag.txt");
	x = matrix(optLag.x);
	y = solve_ode(df1, 0, 1, 1000, y0, ud1, x);
	max = y[1](0, 2);
	for (int i = 0; i < 1000; i++)
	{
		file << y[1](i, 0) << "\t" << y[1](i, 1) << "\t" << y[1](i, 2) << endl;
	}
	file.close();

	solution::clear_calls();
}

void lab2()
{
	solution opt;
	matrix x0(2,1);
	long double div = 10e4;
	unsigned int range_size = 100;

	long double s = 0.05;
	matrix s0(2, 1);
	long double alpha = 0.75;
	long double epsilon = 0.000005;
	int Nmax = 1000;

	long double rosen_alpha = 2.5;
	long double rosen_beta = 0.75;
	
	srand(time(NULL));
	x0(0, 0) = static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / 2.0f)) - 1.0f;
	x0(1, 0) = static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX / 2.0f)) - 1.0f;


	opt = HJ(testowa_lab_2, x0, s, alpha, epsilon, Nmax);
	cout << opt << endl;
	//matrix s0(2, 1);
	s0(0, 0) = s;
	s0(1, 0) = s;
	opt = Rosen(testowa_lab_2, x0, s0, 2.0, 0.5, 0.0001, Nmax);
	
	std::ofstream file;
	/*file.open("lab_02.txt");
	for (s = 0.1; s >= 0.02; s /= 2)
	{
		s0(0,0) = s;
		s0(1,0) = s;
		for (int i = 0; i < 100; i++)
		{
			//cout << "Iteracja: " << i << endl;
			x0(0, 0) = static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX/2.0f)) - 1.0f;
			x0(1, 0) = static_cast <float> (std::rand()) / (static_cast <float> (RAND_MAX/2.0f)) - 1.0f;

			file << x0(0, 0) << '\t' << x0(1, 0) << '\t';
			opt = HJ(testowa_lab_2, x0, s, alpha, epsilon, Nmax);
			file << opt.x(0, 0) << "\t" << opt.x(1, 0) << "\t" << opt.y << "\t" << opt.f_calls;
			if ((abs(opt.x(0, 0)) > 10e-4) || (abs(opt.x(1, 0)) > 10e-4)) file << "\tNIE\t";
			else file << "\tTAK\t";

			solution::clear_calls();

			opt = Rosen(testowa_lab_2, x0, s0, rosen_alpha, rosen_beta, epsilon, Nmax);
			file << opt.x(0, 0) << "\t" << opt.x(1, 0) << "\t" << opt.y << "\t" << opt.f_calls;
			if ((abs(opt.x(0, 0)) > 10e-4) || (abs(opt.x(1, 0)) > 10e-4)) file << "\tNIE\t";
			else file << "\tTAK\t";

			file << endl;
		}
	}
	file.close();*/

	// Problem rzeczywisty
	s = 0.01;
	s0(0, 0) = s;
	s0(1, 0) = s;
	
	// HJ

	/*matrix k0(2, 1);
	k0(0, 0) = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 10.0f));
	k0(1, 0) = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 10.0f));
	opt = HJ(ff_lab2, k0, s, alpha, epsilon, 1000);
	long double dt = 0.1;
	
	matrix y_ref(2, 1);
	y_ref(0, 0) = 3.14;
	y_ref(1, 0) = 0.0; 
	matrix y_0(2, 1);
	y_0(0, 0) = 0;
	y_0(1, 0) = 0;
	matrix* y_solver = solve_ode(df2, 0.0f, dt, 100, y_0, y_ref, opt.x);
	file.open("lab_02_symulacja_HJ.txt");
	int n = get_len(y_solver[0]);
	for (int i = 0; i < n; i++)
	{
		file << y_solver[1](i, 0) << '\t' << y_solver[1](i, 1) << '\n';
	}
	delete[] y_solver;
	file.close();

	cout << opt << endl;

	solution::clear_calls();
	// Rosen

	opt = Rosen(ff_lab2, k0, s0, rosen_alpha, rosen_beta, epsilon, 1000);

	y_solver = solve_ode(df2, 0.0f, dt, 100, y_0, y_ref, opt.x);
	file.open("lab_02_symulacja_Rosen.txt");
	n = get_len(y_solver[0]);
	for (int i = 0; i < n; i++)
	{
		file << y_solver[1](i, 0) << '\t' << y_solver[1](i, 1) << '\n';
	}
	delete[] y_solver;
	file.close();*/
	
	cout << opt << endl;
	solution::clear_calls();
}

void lab3()
{
		// FUNCKJA TESTOWA //
	srand(time(NULL));
	solution opt;

		// Pocz�tkowy simpleks //
	matrix x0(2, 1);

		// Ograniczenia //
	double xLower = 1;
	double xUpper = 3;
	double yLower = 1;
	double yUpper = 3;

	x0(0, 0) = ((float)std::rand() / (float)RAND_MAX)* (xUpper - xLower) + xLower;
	x0(1, 0) = ((float)std::rand() / (float)RAND_MAX)* (yUpper - yLower) + yLower;
	std::cout << x0(0, 0) << ' ' << x0(1, 0) << std::endl;

		// Parametry algorytmu //
	long double s = 0.1;
	long double alpha = 1;
	long double beta = 0.8;
	long double gamma = 1.5;
	long double delta = 0.8;
	long double epislon = 10e-4;
	int Nmax = 1000;

		// Parametry funkcji testowe //
	matrix ud1(1, 1, 5);
	matrix ud2(1, 1, 10e4);
	float a[3] = { 4, 4.4934, 5 };
	std::ofstream file;
	
	file.open("testowa_lab4.txt");
	for (int a_idx = 0; a_idx< 3; a_idx++)
	{
		ud1(0, 0) = a[a_idx];
		
		for (int i = 0; i < 100; i++)
		{
				// Losowanie wierzcho�ka pocz�tkowego simpleksu //
			xUpper = a[a_idx];
			x0(0, 0) = ((float)std::rand() / (float)RAND_MAX) * (xUpper - xLower) + xLower;
			yUpper = sqrt(a[a_idx] * a[a_idx] - x0(0, 0) * x0(0, 0));
			x0(1, 0) = ((float)std::rand() / (float)RAND_MAX) * (yUpper - yLower) + yLower;
			file << x0(0, 0) << '\t' << x0(1, 0) << '\t';

			// Wywo�anie funkcji algorytmu //
			solution::clear_calls();
			ud2(0, 0) = 10e4;
			opt = sym_NM(testowa_lab_3_zew, x0, s, alpha, beta, gamma, delta, epislon, Nmax, ud1, ud2);
			file << opt.x(0, 0) << '\t' << opt.x(1, 0) << '\t' << norm(opt.x) << '\t' << 
				opt.y << '\t' << opt.f_calls << '\t';

			solution::clear_calls();
			ud2(0, 0) = 0.01;
			opt = sym_NM(testowa_lab_3_wew, x0, s, alpha, beta, gamma, delta, epislon, Nmax, ud1, ud2);
			file << opt.x(0, 0) << '\t' << opt.x(1, 0) << '\t' << norm(opt.x) << '\t' << 
				testowa_lab_3_wew_clear(opt.y, opt.x, ud1, ud2) << '\t' << opt.f_calls << '\t';
			file << std::endl;
		}
		
	}
	file.close();
	
	// PROBELM RZECZYWISTY //
	solution opt_min;
	opt_min.y = 10e5;

		// Wywo�anie funkcji algorytmu //
	x0(1, 0) = 8; ((float)std::rand() / (float)RAND_MAX) * 46 - 23;
	x0(0, 0) = -1;  ((float)std::rand() / (float)RAND_MAX) * 20 - 10;
		ud2(0) = 10e4;
		std::cout << x0(0, 0) << ' ' << x0(1, 0) << std::endl;
		solution::clear_calls();
		opt = sym_NM(ff_lab3, x0, s, alpha, beta, gamma, delta, epislon, Nmax, ud1, ud2);
		if (opt.y < opt_min.y)
			opt_min = opt;

		std::cout << opt << endl;


	matrix y0(4, 1);
	y0(0) = 0;
	y0(1) = opt.x(0,0);
	y0(2) = 100; 
	y0(3) = 0;

	matrix* y = solve_ode(df3, 0, 0.01, 7, y0, ud1, opt.x(1,0));
	int i0 = 0;
	int i50 = 0;
	int n = get_len(y[0]);
	file.open("simulation.txt");
	for (int i = 0; i < n; i++)
	{
		file << y[1](i, 0) << '\t' << y[1](i, 2) << std::endl;
		if (abs(y[1](i, 2)) < abs(y[1](i0, 2)))
		{
			i0 = i;
		}
		if (abs(y[1](i, 2) - 50) < abs(y[1](i50, 2) - 50))
		{
			i50 = i;
		}
	}
	file.close();
	matrix res = y[1](i0, 0);
	delete[] y;
		// Sprawdzenie wyniku //
	std::cout << opt << endl;
}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}