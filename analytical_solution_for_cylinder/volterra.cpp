#include "volterra.h"
#include <cmath>
#include <fstream>
#include <algorithm>

VolterraParams get_default_volterra_params() {
    VolterraParams params;
    params.Mach = sqrt(2); // unused
    params.Q = pow(PI, 1.5) / 20.0 * 1.2255 * pow(2*1.4*101330/1.2255, 1.5) * pow(0.01, 2);
    params.r_source = 1.2;
    params.x_source = -1.2;  // default: xu = -ru
    params.L = 10.0;
    params.N = 10000;
    params.h_integral = 1e-3;
    params.output_dir = "output_results";
    return params;
}

/**
 * Compute the upper integration limit for g(x) and dphiu1_dx(x) integrals.
 * This determines the angular range of integration based on geometry.
 */
static double compute_upper_limit(double x, double r_source, double x_source) {
    const double r_body = 1.0;  // body radius (normalized)

    // y = intersection point of source cone with x=const plane
    double y = (r_source * r_source + r_body * r_body - (x - x_source) * (x - x_source)) / (2 * r_source);
    double z = std::sqrt(std::abs(1 - y * y));

    // Check if x is at the leading edge
    if (std::abs(x + 1) < 1e-15) {
        return 0;  // no contribution at leading edge
    }

    double upper_limit = PI;

    if (x < 1) {
        if (std::abs(y) < 1e-15) {
            upper_limit = PI / 2;
        } else {
            if (y > 0) {
                upper_limit = std::atan(z / y);
            } else {
                upper_limit = PI + std::atan(z / y);
            }
        }
    }

    return upper_limit;
}

/**
 * Compute g(x) - the kernel function for the Volterra integral equation.
 * This represents the influence of the heat source at position x.
 *
 * Old name: g(double x)
 */
double compute_g_kernel(double x, const VolterraParams& params) {
    double theta = 0;
    double Int = 0;
    double h = params.h_integral;

    double upper_limit = compute_upper_limit(x, params.r_source, params.x_source);

    if (upper_limit < 1e-10) {
        return 0;
    }

    const double r_body = 1.0;
    double ru = params.r_source;
    double xu = params.x_source;

    // Numerical integration using trapezoidal rule
    while (theta < upper_limit) {
        double denom = r_body * r_body + ru * ru - 2 * r_body * ru * std::cos(theta);
        double sqrt_term = std::sqrt((x - xu) * (x - xu) - denom);
        Int += h * std::cos(theta) * (r_body - ru * std::cos(theta)) * sqrt_term / denom;
        theta += h;
    }

    // Last partial step
    theta -= h;
    double denom = r_body * r_body + ru * ru - 2 * r_body * ru * std::cos(theta);
    double sqrt_term = std::sqrt((x - xu) * (x - xu) - denom);
    Int += (upper_limit - theta - h) * std::cos(theta) * (r_body - ru * std::cos(theta)) * sqrt_term / denom;

    Int *= params.Q / (PI * PI);

    return Int;
}

/**
 * Solve the Volterra integral equation to find q(x).
 * This is the heat release distribution along the body.
 *
 * Old name: qq(double *q)
 */
std::vector<double> solve_volterra_equation(const VolterraParams& params) {
    int N = params.N;
    double h = params.L / N;

    std::vector<double> q(N);

    // First point
    q[0] = compute_g_kernel(h - 1, params) / ((h + 1) * std::sqrt((h + 1) * (h + 1) - 1));

    // Subsequent points using Volterra equation
    for (int i = 2; i <= N; i++) {
        q[i - 1] = compute_g_kernel(i * h - 1, params);

        for (int j = 2; j <= i; j++) {
            double term1 = std::sqrt((j * h + 1) * (j * h + 1) - 1);
            double term2 = std::sqrt(((j - 1) * h + 1) * ((j - 1) * h + 1) - 1);
            q[i - 1] -= q[i - j] * (j * h + 1) * (term1 - term2);
        }

        q[i - 1] /= ((h + 1) * std::sqrt((h + 1) * (h + 1) - 1));
    }

    return q;
}

/**
 * Compute d(phi1)/dx - derivative of potential on the body surface.
 *
 * Old name: dphi1_dx_(double *q, double *dphi1_dx)
 */
std::vector<double> compute_dphi1_dx(const std::vector<double>& q, const VolterraParams& params) {
    int N = params.N;
    double h = params.L / N;

    std::vector<double> dphi1_dx(N, 0.0);

    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= i; j++) {
            double term1 = std::sqrt((j * h + 1) * (j * h + 1) - 1);
            double term2 = std::sqrt(((j - 1) * h + 1) * ((j - 1) * h + 1) - 1);
            dphi1_dx[i] += q[i - j] * (term1 - term2);
        }
        dphi1_dx[i] *= -1;
    }

    return dphi1_dx;
}

/**
 * Compute d(phi_u1)/dx at a single point - perturbation potential derivative.
 *
 * Old name: dphiu1_dx_(double x)
 */
static double compute_dphiu1_dx_at_point(double x, const VolterraParams& params) {
    double theta = 0;
    double Int = 0;
    double h = params.h_integral;

    double upper_limit = compute_upper_limit(x, params.r_source, params.x_source);

    if (upper_limit < 1e-10) {
        return 0;
    }

    const double r_body = 1.0;
    double ru = params.r_source;
    double xu = params.x_source;

    while (theta < upper_limit) {
        double denom = r_body * r_body + ru * ru - 2 * r_body * ru * std::cos(theta);
        Int += h * std::cos(theta) * std::acosh((x - xu) / std::sqrt(denom));
        theta += h;
    }

    theta -= h;
    double denom = r_body * r_body + ru * ru - 2 * r_body * ru * std::cos(theta);
    Int += (upper_limit - theta - h) * std::cos(theta) * std::acosh((x - xu) / std::sqrt(denom));

    Int *= params.Q / (PI * PI);

    return Int;
}

/**
 * Compute d(phi_u1)/dx for all points.
 *
 * Old name: dphiu1_dx_file(double *dphiu1_dx)
 */
std::vector<double> compute_dphiu1_dx(const VolterraParams& params) {
    int N = params.N;
    double h = params.L / N;

    std::vector<double> dphiu1_dx(N);

    for (int i = 0; i < N; i++) {
        dphiu1_dx[i] = compute_dphiu1_dx_at_point(i * h - 1, params);
    }

    return dphiu1_dx;
}

/**
 * Compute d(phi1)/dr - radial derivative of potential.
 *
 * Old name: dphi1_dr_(double *q, double *dphi1_dr)
 */
std::vector<double> compute_dphi1_dr(const std::vector<double>& q, const VolterraParams& params) {
    int N = params.N;
    double h = params.L / N;

    std::vector<double> dphi1_dr(N, 0.0);

    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= i; j++) {
            double term1 = std::sqrt((j * h + 1) * (j * h + 1) - 1);
            double term2 = std::sqrt(((j - 1) * h + 1) * ((j - 1) * h + 1) - 1);
            dphi1_dr[i] += q[i - j] * (j * h + 1) * (term1 - term2);
        }
    }

    return dphi1_dr;
}

/**
 * Main solver function - computes all results.
 */
VolterraResult solve_volterra(const VolterraParams& params) {
    VolterraResult result;

    int N = params.N;
    double h = params.L / N;

    // Generate x coordinates
    result.x.resize(N);
    for (int i = 0; i < N; i++) {
        result.x[i] = i * h - 1;
    }

    // Solve Volterra equation
    result.q = solve_volterra_equation(params);

    // Compute derivatives
    result.dphi1_dx = compute_dphi1_dx(result.q, params);
    result.dphiu1_dx = compute_dphiu1_dx(params);
    result.dphi1_dr = compute_dphi1_dr(result.q, params);

    // Compute sum and derivative
    result.sum.resize(N);
    result.derivative.resize(N);

    for (int i = 0; i < N; i++) {
        result.sum[i] = PI * (result.dphi1_dx[i] + result.dphiu1_dx[i]);
    }

    for (int i = 0; i < N - 1; i++) {
        result.derivative[i] = PI * (result.dphi1_dx[i + 1] + result.dphiu1_dx[i + 1]
                                    - result.dphi1_dx[i] - result.dphiu1_dx[i]) / h;
    }
    // Extrapolate last point
    result.derivative[N - 1] = result.derivative[N - 2] + (result.derivative[N - 2] - result.derivative[N - 3]);

    return result;
}

/**
 * Save all results to files.
 */
void save_results(const VolterraResult& result, const VolterraParams& params) {
    std::string dir = params.output_dir;
    int N = params.N;

    // q.txt
    {
        std::ofstream f(dir + "/q.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.q[i] << "\n";
        }
    }

    // dphi1_dx.txt
    {
        std::ofstream f(dir + "/dphi1_dx.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.dphi1_dx[i] << "\n";
        }
    }

    // dphiu1_dx.txt
    {
        std::ofstream f(dir + "/dphiu1_dx.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.dphiu1_dx[i] << "\n";
        }
    }

    // dphi1_dr.txt
    {
        std::ofstream f(dir + "/dphi1_dr.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.dphi1_dr[i] << "\n";
        }
    }

    // sum.txt
    {
        std::ofstream f(dir + "/sum.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.sum[i] << "\n";
        }
    }

    // derivative.txt
    {
        std::ofstream f(dir + "/derivative.txt");
        for (int i = 0; i < N; i++) {
            f << result.x[i] << " " << result.derivative[i] << "\n";
        }
    }
}
