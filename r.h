#pragma once
#include <array>
#include <initializer_list>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <ctime>

constexpr double Gamma = 1.4;
constexpr double Pi = 3.141592653589793;

enum class BodyType {
    Cylindrical,
    Cone,
	Parabolic,
	DoubleCone
};

struct Params {
    double Mach_inf;
    double p_inf;
    double rho_inf;
    double a_inf;
    double V_inf;
    bool is_adiabatic;
	BodyType bodyType;
};

extern Params params;

std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b);

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b);

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, T scalar);

template <typename T>
std::vector<T> operator*(T scalar, const std::vector<T>& vec);

class BaseArray {
	protected:
		static const int SIZE = 5;
	public:
		std::array<double, SIZE> data;
		
		BaseArray();
		BaseArray(const std::initializer_list<double>& list);
	
		BaseArray& operator=(const BaseArray& other);
		BaseArray operator+(const BaseArray& other) const;
		BaseArray operator-(const BaseArray& other) const;
		BaseArray operator*(double scalar) const;
		friend BaseArray operator*(double scalar, const BaseArray& arr);
	
		double& operator[](std::size_t index);
		const double& operator[](std::size_t index) const;
	
		virtual BaseArray process() const;
	
		virtual double get_rho(double r) const;
		virtual double get_p(double r) const;
		virtual double get_u() const;
		virtual double get_v() const;
		virtual double get_w() const;
	
		virtual ~BaseArray();
	
		void print() const;
	};

class E_array : public BaseArray {
public:
	E_array();
	E_array(const BaseArray& base);
	E_array(const std::initializer_list<double>& list);
	E_array& operator=(const BaseArray& other);
};

class F_array : public BaseArray {
public:
	F_array();
	F_array(const BaseArray& base);
	F_array(const std::initializer_list<double>& list);
	F_array& operator=(const BaseArray& other);
};

class G_array : public BaseArray {
public:
	G_array();
	G_array(const BaseArray& base);
	G_array(const std::initializer_list<double>& list);
	G_array& operator=(const BaseArray& other);
	double get_alpha() const;
	double get_rho(double r) const override;
	double get_p(double r) const override;
	double get_u() const override;
	double get_v() const override;
	double get_w() const override;
};

class R_array : public BaseArray {
public:
	R_array();
	R_array(const BaseArray& base);
	R_array(const std::initializer_list<double>& list);
	R_array& operator=(const BaseArray& other);
};

G_array predictor(E_array Em, E_array Ep, F_array Fm, F_array Fp, G_array G, R_array R, double dr, double dth, double dz);
G_array corrector(E_array Em, E_array Ep, F_array Fm, F_array Fp, G_array Gm, G_array Gp, R_array R, double dr, double dth, double dz);
E_array get_E(G_array G, double r);
F_array get_F(G_array G, double r);
R_array get_R(G_array G, double r, double q);

double r_from_xi(double xi, double r_s, double r_b);
double r_b(double z, BodyType bodyType);
double r_b_z(double z, BodyType bodyType);
double xi_r(double r_s, double r_b);
double xi_theta(double xi, double r_s, double r_b, double r_s_theta);
double xi_z(double xi, double r_s, double r_b, double r_s_z, double r_b_z);
double q(double r, double theta, double z, double x_q, double z_q, bool is_adiabatic);
double lambda_r(double rho, double p, double u, double v, double w);

std::vector<double> solver(double x_q, double z_q);
