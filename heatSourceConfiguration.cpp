#include "r.h"

double q(double r, double theta, double z, double x_q, double z_q, bool is_adiabatic) {
    if (is_adiabatic){
        return 0.0;
    }

    double
        y_q = 0,
        L_q = 0.02,
        x = r * cos(theta),
        y = r * sin(theta),
        Q = 1.0 / 10.0,
        q_0 = Q * params.V_inf * params.V_inf * params.V_inf / L_q;

    return q_0 * exp(-((x - x_q)*(x - x_q) + (y - y_q)*(y - y_q) + (z - z_q)*(z - z_q)) / L_q / L_q);
}
