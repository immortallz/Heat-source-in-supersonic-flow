#include "r.h"

double r_from_xi(double xi, double r_s, double r_b){
    return (r_s - r_b)*xi + r_b;
}

double xi_r(double r_s, double r_b) {
    return 1 / (r_s - r_b);
}

double xi_theta(double xi, double r_s, double r_b, double r_s_theta) {
    return -xi * r_s_theta / (r_s - r_b);
}

double xi_z(double xi, double r_s, double r_b, double r_s_z, double r_b_z) {
    return -(r_b_z + xi * (r_s_z - r_b_z)) / (r_s - r_b);
}

double lambda_r(double rho, double p, double u, double v, double w)
{
    double lambda = 
        (abs(rho*u*w) + sqrt(p*rho*(u*u + w*w) - p*p))
        / (rho*w*w - p);
    return lambda;
}
