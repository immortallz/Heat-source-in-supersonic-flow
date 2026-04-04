#include "r.h"

double r_from_xi(const double xi, const double r_s, const double r_b){
    return (r_s - r_b)*xi + r_b;
}

double xi_r(const double r_s, const double r_b) {
    return 1 / (r_s - r_b);
}

double xi_theta(const double xi, const double r_s, const double r_b, const double r_s_theta) {
    return -xi * r_s_theta / (r_s - r_b);
}

double xi_z(const double xi, const double r_s, const double r_b, const double r_s_z, const double r_b_z) {
    return -(r_b_z + xi * (r_s_z - r_b_z)) / (r_s - r_b);
}

double lambda_r(const double rho, const double p, const double u, const double v, const double w)
{
    const double lambda =
        (abs(rho*u*w) + sqrt(p*rho*(u*u + w*w) - p*p))
        / (rho*w*w - p);
    return lambda;
}
