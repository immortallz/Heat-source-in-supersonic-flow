#pragma once
#include <vector>
#include <string>

constexpr double PI = 3.141592653589793;
constexpr double GAMMA = 1.4;

std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b);

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b);

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, T scalar);

template <typename T>
std::vector<T> operator*(T scalar, const std::vector<T>& vec);

double a_squared(const std::vector<double> &v, double C);

std::vector<double> func(double theta, const std::vector<double> &v, double C);

std::vector<double> RK(double t0, double T, double h, std::vector<double> y0, std::string filename, double C);

double secant_method(double p0, double rho0, double V1, double theta0,
	double beta0, double beta1,
	double h, double tol
	);

double newton_method(double p0, double rho0, double V1, double theta0,
	double beta0,
	double h, double tol
	);
