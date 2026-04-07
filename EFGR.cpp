#include "r.h"

// Реализация классов E_array и F_array
E_array::E_array() : BaseArray() {}
E_array::E_array(const BaseArray& base) : BaseArray(base) { }
E_array::E_array(const std::initializer_list<double>& list) : BaseArray(list) {}
E_array& E_array::operator=(const BaseArray& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

F_array::F_array() : BaseArray() {}
F_array::F_array(const BaseArray& base) : BaseArray(base) { }
F_array::F_array(const std::initializer_list<double>& list) : BaseArray(list) {}
F_array& F_array::operator=(const BaseArray& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

// Реализация класса G_array
G_array::G_array() : BaseArray() {}
G_array::G_array(const BaseArray& base) : BaseArray(base) { }
G_array::G_array(const std::initializer_list<double>& list) : BaseArray(list) {}
G_array& G_array::operator=(const BaseArray& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

double G_array::get_alpha() const {
    return -sqrt(
        (Gamma*Gamma - 1)*data[1]*data[1]
        + (Gamma*Gamma - 1)*data[2]*data[2]
        + Gamma*Gamma*data[3]*data[3]
        - 2*(Gamma*Gamma - 1)*data[0]*data[4]
    );
}

double G_array::get_rho(const double r) const {
    const double alpha = get_alpha();
    const double result =
        (data[0]*data[0]*(Gamma*data[3] + alpha))
        / ((Gamma - 1)*(2*data[0]*data[4] - data[1]*data[1] - data[2]*data[2]));
    return result / r;
}

double G_array::get_p(const double r) const {
    const double alpha = get_alpha();
    const double result = (data[3] + alpha) / (Gamma + 1);
    return result / r;
}

double G_array::get_u() const {
    return data[1] / data[0];
}

double G_array::get_v() const {
    return data[2] / data[0];
}

double G_array::get_w() const {
    const double alpha = get_alpha();
    const double result = (Gamma*data[3] - alpha) / ((Gamma + 1)*data[0]);
    return result;
}

// Реализация класса R_array
R_array::R_array() : BaseArray() {}
R_array::R_array(const BaseArray& base) : BaseArray(base) { }
R_array::R_array(const std::initializer_list<double>& list) : BaseArray(list) {}
R_array& R_array::operator=(const BaseArray& other) {
    if (this != &other) {
        data = other.data;
    }
    return *this;
}

// Функции вычислений
G_array predictor(
    const E_array& Em, const E_array& Ep,
    const F_array& Fm, const F_array& Fp,
    const G_array &G,
    const R_array& R,
    const double dr, const double dth, const double dz
) {
    G_array result = G;
    result = result - dz/dr * (Ep - Em) - dz/dth * (Fp - Fm) + dz * R;
    return result;
}

G_array corrector(
    const E_array& Em, const E_array& Ep,
    const F_array& Fm, const F_array& Fp,
    const G_array& Gm, const G_array& Gp,
    const R_array& R,
    const double dr, const double dth, const double dz
) {
    G_array result = (Gp + Gm) * 0.5;
    result = result - 0.5 * dz/dr * (Ep - Em) - 0.5 * dz/dth * (Fp - Fm) + 0.5 * dz * R;
    return result;
}

E_array get_E(const G_array& G, const double r) {
    const double rho = G.get_rho(r);
    const double p   = G.get_p(r);
    const double u   = G.get_u();
    const double v   = G.get_v();
    const double w   = G.get_w();
    const E_array result = {
        rho * u,
        rho * u * u + p,
        rho * u * v,
        rho * u * w,
        u * (Gamma/(Gamma - 1)*p + 0.5 * rho * (u*u + v*v + w*w))
    };
    return r * result;
}

F_array get_F(const G_array& G, const double r) {
    // Для функций get_F и get_R используем r = 1 при вычислении rho и p
    const double rho = G.get_rho(r);
    const double p   = G.get_p(r);
    const double u   = G.get_u();
    const double v   = G.get_v();
    const double w   = G.get_w();
    const F_array result = {
        rho * v,
        rho * u * v,
        rho * v * v + p,
        rho * v * w,
        v * (Gamma/(Gamma - 1)*p + 0.5 * rho * (u*u + v*v + w*w))
    };
    return result;
}

R_array get_R(const G_array& G, const double r, const double q) {
    const double rho = G.get_rho(r);
    const double p   = G.get_p(r);
    const double u   = G.get_u();
    const double v   = G.get_v();
    const R_array result = {
        0,
        rho * v * v + p,
        -rho * u * v,
        0,
        q
    };
    return result;
}

FluxPair get_fluxes(const std::vector<std::vector<E_array>>& E, int i, int j, const FluxScheme scheme, bool is_predictor) {
    FluxPair flux;

    if (scheme == FluxScheme::BeamWarming) {
        if (is_predictor && i == numericalParams.N - 1) {
            // Right boundary predictor - upwind from left
            flux.E_left = E[i][j];
            flux.E_right = 3*E[i][j] - 3*E[i-1][j] + E[i-2][j];
        } else if (!is_predictor && i == 0) {
            // Left boundary corrector - upwind from right
            flux.E_left = 3*E[i][j] - 3*E[i+1][j] + E[i+2][j];
            flux.E_right = E[i][j];
        } else {
            // Standard case
            if (is_predictor) {
                flux.E_left = E[i][j];
                flux.E_right = E[i + 1][j];
            } else {
                flux.E_left = E[i - 1][j];
                flux.E_right = E[i][j];
            }
        }
    } else {
        // Central difference
        if (is_predictor) {
            int idx = int(i == numericalParams.N - 1);
            flux.E_left = E[i - idx][j];
            flux.E_right = E[i + 1 - idx][j];
        } else {
            int idx = int(i == 0);
            flux.E_left = E[i - 1 + idx][j];
            flux.E_right = E[i + idx][j];
        }
    }

    return flux;
}
