#include "r.h"
#include <cmath>

double q(const double r, const double theta, const double z, const HeatSource &heatSource, const bool is_adiabatic) {
    if (is_adiabatic){
        return 0.0;
    }

    const double
        x = r * cos(theta),
        y = r * sin(theta),
        x_q = heatSource.x,
        y_q = heatSource.y,
        z_q = heatSource.z,
        L_q = heatSource.L,
        Q = heatSource.Q,
        q_0 = Q * flowParams.V_inf * flowParams.V_inf * flowParams.V_inf / L_q;

    return q_0 * exp(-((x - x_q)*(x - x_q) + (y - y_q)*(y - y_q) + (z - z_q)*(z - z_q)) / L_q / L_q);
}
