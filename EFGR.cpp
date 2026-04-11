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
        (GAMMA*GAMMA - 1)*data[1]*data[1]
        + (GAMMA*GAMMA - 1)*data[2]*data[2]
        + GAMMA*GAMMA*data[3]*data[3]
        - 2*(GAMMA*GAMMA - 1)*data[0]*data[4]
    );
}

PhysicalParameters G_array::get_physical_parameters(const double r) const {
    const double alpha = get_alpha();
    PhysicalParameters prim{};

    prim.rho = (data[0]*data[0]*(GAMMA*data[3] + alpha))
             / ((GAMMA - 1)*(2*data[0]*data[4] - data[1]*data[1] - data[2]*data[2]))
             / r;
    prim.p = (data[3] + alpha) / (GAMMA + 1) / r;
    prim.u = data[1] / data[0];
    prim.v = data[2] / data[0];
    prim.w = (GAMMA*data[3] - alpha) / ((GAMMA + 1)*data[0]);

    return prim;
}

double G_array::get_rho(const double r) const {
    const double alpha = get_alpha();
    const double result =
        (data[0]*data[0]*(GAMMA*data[3] + alpha))
        / ((GAMMA - 1)*(2*data[0]*data[4] - data[1]*data[1] - data[2]*data[2]));
    return result / r;
}

double G_array::get_p(const double r) const {
    const double alpha = get_alpha();
    const double result = (data[3] + alpha) / (GAMMA + 1);
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
    const double result = (GAMMA*data[3] - alpha) / ((GAMMA + 1)*data[0]);
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
    const auto phys = G.get_physical_parameters(r);
    const E_array result = {
        phys.rho * phys.u,
        phys.rho * phys.u * phys.u + phys.p,
        phys.rho * phys.u * phys.v,
        phys.rho * phys.u * phys.w,
        phys.u * (GAMMA/(GAMMA - 1)*phys.p + 0.5 * phys.rho * (phys.u*phys.u + phys.v*phys.v + phys.w*phys.w))
    };
    return r * result;
}

F_array get_F(const G_array& G, const double r) {
    // Для функций get_F и get_R используем r = 1 при вычислении rho и p
    const auto phys = G.get_physical_parameters(r);
    const F_array result = {
        phys.rho * phys.v,
        phys.rho * phys.u * phys.v,
        phys.rho * phys.v * phys.v + phys.p,
        phys.rho * phys.v * phys.w,
        phys.v * (GAMMA/(GAMMA - 1)*phys.p + 0.5 * phys.rho * (phys.u*phys.u + phys.v*phys.v + phys.w*phys.w))
    };
    return result;
}

R_array get_R(const G_array& G, const double r, const double q) {
    const auto phys = G.get_physical_parameters(r);
    const R_array result = {
        0,
        phys.rho * phys.v * phys.v + phys.p,
        -phys.rho * phys.u * phys.v,
        0,
        q
    };
    return result;
}

FluxPair get_fluxes(
    const std::vector<std::vector<E_array>>& E,
    const std::vector<std::vector<E_array>>& E_prev,
    const int i, const int j, const FluxScheme scheme, const bool is_predictor)
{
    FluxPair flux;

    if (scheme == FluxScheme::BeamWarming) {
        if (is_predictor && i == numericalParams.N - 1) {
            // Right boundary predictor - upwind from left (Moretti scheme)
            flux.E_left = E[i][j];
            flux.E_right = 3*E[i][j] - 3*E[i-1][j] + E[i-2][j];
        } else if (!is_predictor && i == 0) {
            // Left boundary corrector - upwind from right (actual Beam-Warming)
            flux.E_left = E[i][j];
            flux.E_right = E[i+1][j] + (E_prev[i][j] - 2*E_prev[i+1][j] + E_prev[i+2][j]);
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
            const int idx = int(i == numericalParams.N - 1);
            flux.E_left = E[i - idx][j];
            flux.E_right = E[i + 1 - idx][j];
        } else {
            const int idx = int(i == 0);
            flux.E_left = E[i - 1 + idx][j];
            flux.E_right = E[i + idx][j];
        }
    }

    return flux;
}
