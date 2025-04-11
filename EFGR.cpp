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

double G_array::get_rho(double r) const {
    double alpha = -sqrt(
        (gamma*gamma - 1)*data[1]*data[1]
        + (gamma*gamma - 1)*data[2]*data[2]
        + gamma*gamma*data[3]*data[3]
        - 2*(gamma*gamma - 1)*data[0]*data[4]
    );
    double result =
        (data[0]*data[0]*(gamma*data[3] + alpha))
        / ((gamma - 1)*(2*data[0]*data[4] - data[1]*data[1] - data[2]*data[2]));
    return result / r;
}

double G_array::get_p(double r) const {
    double alpha = -sqrt(
        (gamma*gamma - 1)*data[1]*data[1]
        + (gamma*gamma - 1)*data[2]*data[2]
        + gamma*gamma*data[3]*data[3]
        - 2*(gamma*gamma - 1)*data[0]*data[4]
    );
    double result = (data[3] + alpha) / (gamma + 1);
    return result / r;
}

double G_array::get_u() const {
    return data[1] / data[0];
}

double G_array::get_v() const {
    return data[2] / data[0];
}

double G_array::get_w() const {
    double alpha = -sqrt(
        (gamma*gamma - 1)*data[1]*data[1]
        + (gamma*gamma - 1)*data[2]*data[2]
        + gamma*gamma*data[3]*data[3]
        - 2*(gamma*gamma - 1)*data[0]*data[4]
    );
    double result = (gamma*data[3] - alpha) / ((gamma + 1)*data[0]);
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
    E_array Em, 
    E_array Ep, 
    F_array Fm, 
    F_array Fp, 
    G_array G, 
    R_array R, 
    double dr, double dth, double dz
) {
    G_array result = G;
    result = result - dz/dr * (Ep - Em) - dz/dth * (Fp - Fm) + dz * R;
    return result;
}

G_array corrector(
    E_array Em, 
    E_array Ep, 
    F_array Fm, 
    F_array Fp, 
    G_array Gm, 
    G_array Gp, 
    R_array R, 
    double dr, double dth, double dz
) {
    G_array result;
    result = (Gp + Gm) * 0.5;
    result = result - 0.5 * dz/dr * (Ep - Em) - 0.5 * dz/dth * (Fp - Fm) + 0.5 * dz * R;
    return result;
}

E_array get_E(G_array G, double r) {
    double rho = G.get_rho(r);
    double p   = G.get_p(r);
    double u   = G.get_u();
    double v   = G.get_v();
    double w   = G.get_w();
    E_array result = {
        rho * u,
        rho * u * u + p,
        rho * u * v,
        rho * u * w,
        u * (gamma/(gamma - 1)*p + 0.5 * rho * (u*u + v*v + w*w))
    };
    return r * result;
}

F_array get_F(G_array G, double r) {
    // Для функций get_F и get_R используем r = 1 при вычислении rho и p
    double rho = G.get_rho(r);
    double p   = G.get_p(r);
    double u   = G.get_u();
    double v   = G.get_v();
    double w   = G.get_w();
    F_array result = {
        rho * v,
        rho * u * v,
        rho * v * v + p,
        rho * v * w,
        v * (gamma/(gamma - 1)*p + 0.5 * rho * (u*u + v*v + w*w))
    };
    return result;
}

R_array get_R(G_array G, double r, double q_val) {
    double rho = G.get_rho(r);
    double p   = G.get_p(r);
    double u   = G.get_u();
    double v   = G.get_v();
    R_array result = {
        0,
        rho * v * v + p,
        -rho * u * v,
        0,
        q_val
    };
    return result;
}

double r_from_xi(double xi, double r_s, double r_b)
{
    return (r_s - r_b)*xi + r_b;
}

double r_b(double z) {
    return tan(PI / 6) * z;
}

double r_b_z(double z, double dz){
    return (r_b(z + dz) - r_b(z - dz)) / dz * 0.5;
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

double q(double r, double theta, double z) {
    return 0;
}

double lambda_r(double rho, double p, double u, double v, double w)
{
    double lambda = 
        (abs(rho*u*w) + sqrt(p*rho*(u*u + w*w) - p*p))
        / (rho*w*w - p);
    return lambda;
}

double lambda_th(double r, double rho, double p, double u, double v, double w)
{
    double lambda =
        (abs(rho*v*w) + sqrt(p*rho*(v*v + w*w) - p*p))
        / (rho*w*w - p) / r;
    return lambda;
}
