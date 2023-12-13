//Ten plik nie powinien by� edytowany

#pragma once

#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<random>
#include<chrono>
using namespace std;

#define SEP_SYMBOL ','

class matrix
{
	int n, m;
	long double** M;
	friend int* get_size(const matrix&);
	friend int get_len(const matrix&); // throw (string);
public:
	matrix(long double = 0.0);
	matrix(int, int, long double = 0.0); // throw (string);
	matrix(int, long double*); // throw (string);
	matrix(int, int, long double**); // throw (string);
	matrix(const matrix&);
	~matrix();
	matrix& operator=(const matrix&);
	matrix operator[](int) const; // throw (string);
	long double& operator()(int = 0, int = 0); // throw (string);
	long double operator()(int = 0, int = 0) const; // throw (string);
	void set_col(const matrix&, int); // throw (string);
	void set_row(const matrix&, int); // throw (string);
	void add_col(long double = 0.0);
	void add_row(long double = 0.0);
	void add_col(const matrix&); // throw (string);
	void add_row(const matrix&); // throw (string);
};

matrix operator+(const matrix&, const matrix&); // throw (string);
matrix operator-(const matrix&, const matrix&); // throw (string);
matrix operator*(const matrix&, const matrix&); // throw (string);
matrix operator/(const matrix&, const matrix&); // throw (string);
matrix operator-(const matrix&);
bool operator<(const matrix&, const matrix&); // throw (string);
bool operator>(const matrix&, const matrix&); // throw (string);
bool operator<=(const matrix&, const matrix&); // throw (string);
bool operator>=(const matrix&, const matrix&); // throw (string);
bool operator==(const matrix&, const matrix&); // throw (string);
bool operator!=(const matrix&, const matrix&); // throw (string);
matrix ident_mat(int = 1); // throw (string);
matrix rand_mat(int = 1, int = 1); // throw (string);
matrix randn_mat(int = 1, int = 1); // throw (string);
long double m2d(const matrix&); // throw (string);
long double det(const matrix&); // throw (string);
matrix inv(const matrix&); // throw (string);
matrix trans(const matrix&);
matrix pow(const matrix&, int = 2); // throw (string);
long double norm(const matrix&); // throw (string);
matrix hcat(const matrix&, const matrix&); // throw (string);
matrix vcat(const matrix&, const matrix&); // throw (string);
matrix get_col(const matrix&, int); // throw (string);
matrix get_row(const matrix&, int); // throw (string);
ostream& operator<<(ostream&, const matrix&);
istream& operator>>(istream&, matrix&); // throw (string);
