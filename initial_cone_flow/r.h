#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

#define PI 3.141592653589793
#define gamma 1.4

using namespace std;

vector<double> addVectors(const vector<double>& a, const vector<double>& b);

template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b);

template <typename T>
vector<T> operator*(const vector<T>& vec, T scalar);

template <typename T>
vector<T> operator*(T scalar, const vector<T>& vec);

double a_sqaured(vector<double> v, double C);

vector<double> func(double theta, vector<double> v, double C);

vector<double> RK(double t0, double T, double h, vector<double> y0, string filename, double C);

double secant_method(double p0, double rho0, double V1, double theta0,
	double beta0, double beta1,
	double h, double tol
	);

double newton_method(double p0, double rho0, double V1, double theta0,
	double beta0,
	double h, double tol
	);