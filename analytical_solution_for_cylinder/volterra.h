#pragma once

#include <vector>
#include <string>

constexpr double GAMMA = 1.4;
constexpr double PI = 3.141592653589793;

struct VolterraParams {
    double Mach;        // Mach number
    double Q;           // heat source intensity
    double r_source;    // radial position of heat source
    double x_source;    // axial position of heat source (negative = upstream)
    double L;           // computation length
    int N;              // number of grid points
    double h_integral;  // integration step for numerical integrals
    std::string output_dir;
};

struct VolterraResult {
    std::vector<double> x;          // x coordinates
    std::vector<double> q;          // heat release distribution
    std::vector<double> dphi1_dx;   // d(phi1)/dx on body surface
    std::vector<double> dphiu1_dx;  // d(phi_u1)/dx perturbation potential
    std::vector<double> dphi1_dr;   // d(phi1)/dr radial derivative
    std::vector<double> sum;        // sum of potentials
    std::vector<double> derivative; // derivative of sum
};

// Initialize default parameters
VolterraParams get_default_volterra_params();

// Main solver function
VolterraResult solve_volterra(const VolterraParams& params);

// Individual computation functions
double compute_g_kernel(double x, const VolterraParams& params);
std::vector<double> solve_volterra_equation(const VolterraParams& params);
std::vector<double> compute_dphi1_dx(const std::vector<double>& q, const VolterraParams& params);
std::vector<double> compute_dphiu1_dx(const VolterraParams& params);
std::vector<double> compute_dphi1_dr(const std::vector<double>& q, const VolterraParams& params);

// Output functions
void save_results(const VolterraResult& result, const VolterraParams& params);
