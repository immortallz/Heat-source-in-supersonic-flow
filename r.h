#pragma once
#include <array>
#include <initializer_list>
#include <omp.h>
#include <vector>

constexpr double GAMMA = 1.4;
constexpr double PI = 3.141592653589793;

enum class BodyType {
    Cylindrical,
    Cone,
	Parabolic,
	DoubleCone
};

enum class FluxScheme {
	MacCormack,
	BeamWarming
};

struct FlowParams {
    double Mach_inf;
    double p_inf;
    double rho_inf;
    double a_inf;
    double V_inf;
    bool is_adiabatic;
};

struct BodyParams {
	double transitionPoint; // transition from cone to arbitrary body surface
	double bodyLength;
	BodyType bodyType;
};

struct NumericalParams {
	int N; // xi coordinate nodes count
	int M; // theta coordinate nodes count
	double CFL;
	int num_step_percent;
	int files_count;
	FluxScheme flux_scheme;
};

struct HeatSource {
    double x;
    double y;
    double z;
    double L;   // characteristic sourse length
    double Q;   // dimensionless intensity
};

extern FlowParams flowParams;
extern BodyParams bodyParams;
extern NumericalParams numericalParams;

std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b);

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b);

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, T scalar);

template <typename T>
std::vector<T> operator*(T scalar, const std::vector<T>& vec);

class BaseArray {
	protected:
		static constexpr int SIZE = 5;
	public:
		std::array<double, SIZE> data{};
		
		BaseArray();
		BaseArray(const std::initializer_list<double>& list);

		virtual BaseArray& operator=(const BaseArray& other);
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
	E_array& operator=(const BaseArray& other) override;
};

struct FluxPair {
	E_array E_left, E_right;
};

class F_array : public BaseArray {
public:
	F_array();
	F_array(const BaseArray& base);
	F_array(const std::initializer_list<double>& list);
	F_array& operator=(const BaseArray& other) override;
};

struct PhysicalParameters {
	double rho, p, u, v, w;
};

class G_array : public BaseArray {
public:
	G_array();
	G_array(const BaseArray& base);
	G_array(const std::initializer_list<double>& list);
	G_array& operator=(const BaseArray& other) override;
	PhysicalParameters get_physical_parameters(double r) const;
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
	R_array& operator=(const BaseArray& other) override;
};

G_array predictor(
	const E_array& Em, const E_array& Ep,
	const F_array& Fm, const F_array& Fp,
	const G_array &G,
	const R_array& R,
	double dr, double dth, double dz);

G_array corrector(
	const E_array& Em, const E_array& Ep,
	const F_array& Fm, const F_array& Fp,
	const G_array& Gm, const G_array& Gp,
	const R_array& R,
	double dr, double dth, double dz);

E_array get_E(const G_array& G, double r);
F_array get_F(const G_array& G, double r);
R_array get_R(const G_array& G, double r, double q);

double r_from_xi(double xi, double r_s, double r_b);
double r_b(double z);
double r_b_z(double z);
double xi_r(double r_s, double r_b);
double xi_theta(double xi, double r_s, double r_b, double r_s_theta);
double xi_z(double xi, double r_s, double r_b, double r_s_z, double r_b_z);
double q(double r, double theta, double z, const HeatSource &heatSource, bool is_adiabatic);
double lambda_r(double rho, double p, double u, double v, double w);
FluxPair get_fluxes(
	const std::vector<std::vector<E_array>>& E,
	const std::vector<std::vector<E_array>>& E_prev,
	int i, int j, FluxScheme scheme, bool is_predictor);

std::vector<double> solver(HeatSource heatSource);
