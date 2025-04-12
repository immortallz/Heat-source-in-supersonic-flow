#include <array>
#include <initializer_list>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#define gamma 1.4
#define PI 3.141592653589793

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
	// template <typename T>
	// E_array& operator=(const T& other);
	// E_array operator*(double scalar) const;
	// friend E_array operator*(double scalar, const E_array& arr);
};

class F_array : public BaseArray {
public:
	F_array();
	F_array(const BaseArray& base);
	F_array(const std::initializer_list<double>& list);
	F_array& operator=(const BaseArray& other);
	// template <typename T>
	// F_array& operator=(const T& other);
	// F_array operator*(double scalar) const;
	// friend F_array operator*(double scalar, const F_array& arr);
};

class G_array : public BaseArray {
public:
	G_array();
	G_array(const BaseArray& base);
	G_array(const std::initializer_list<double>& list);
	G_array& operator=(const BaseArray& other);
	// template <typename T>
	// G_array& operator=(const T& other);
	// G_array operator*(double scalar) const;
	// friend G_array operator*(double scalar, const G_array& arr);
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
double r_b(double z);
double r_b_z(double z);
double xi_r(double r_s, double r_b);
double xi_theta(double xi, double r_s, double r_b, double r_s_theta);
double xi_z(double xi, double r_s, double r_b, double r_s_z, double r_b_z);
double q(double r, double theta, double z);
double lambda_r(double rho, double p, double u, double v, double w);
// double lambda_th(double r, double rho, double p, double u, double v, double w);
