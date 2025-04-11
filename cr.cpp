#include <iostream>
#include <array>
#include <initializer_list>
#include <cmath>
#include <vector>

#define gamma 1.4
#define PI 3.141592653589793

using namespace std;

class BaseArray {
protected:
    static const int SIZE = 5;

public:
    std::array<double, SIZE> data;
    BaseArray(){
        for (int i = 0; i < SIZE; ++i)
            data[i] = 0;
    }

    BaseArray(const std::initializer_list<double>& list) {
        int i = 0;
        for (auto it = list.begin(); it != list.end() && i < SIZE; ++it, ++i)
            data[i] = *it;
        for (; i < SIZE; ++i)
            data[i] = 0;
    }

    BaseArray& operator=(const BaseArray& other) {
        if (this != &other) {
            data = other.data; // Копируем данные
        }
        return *this;
    }

    BaseArray operator+(const BaseArray& other) const {
        BaseArray result;
        for (int i = 0; i < SIZE; i++) {
            result.data[i] = this->data[i] + other.data[i];
        }
        return result;
    }

    BaseArray operator-(const BaseArray& other) const {
        BaseArray result;
        for (int i = 0; i < SIZE; i++) {
            result.data[i] = this->data[i] - other.data[i];
        }
        return result;
    }

    BaseArray operator*(double scalar) const {
        BaseArray result;
        for (int i = 0; i < SIZE; i++) {
            result.data[i] = this->data[i] * scalar;
        }
        return result;
    }

    friend BaseArray operator*(double scalar, const BaseArray& arr) {
        return arr * scalar;
    }

    double& operator[](std::size_t index) {
        return data[index];
    }

    const double& operator[](std::size_t index) const {
        return data[index];
    }

    virtual BaseArray process() const {
        return *this;
    }

    virtual double get_rho(double r) const {
        return 0;
    }
    virtual double get_p(double r) const {
        return 0;
    }
    virtual double get_u() const {
        return 0;
    }
    virtual double get_v() const {
        return 0;
    }
    virtual double get_w() const {
        return 0;
    }

    virtual ~BaseArray() {}

    void print() const {
        std::cout << "[ ";
        for (int i = 0; i < SIZE; i++) {
            std::cout << data[i] << " ";
        }
        std::cout << "]" << std::endl;
    }
};

class E_array : public BaseArray {
public:
    E_array() : BaseArray() {}
    E_array(const std::initializer_list<double>& list) : BaseArray(list) {}
};

class F_array : public BaseArray {
public:
    F_array() : BaseArray() {}
    F_array(const std::initializer_list<double>& list) : BaseArray(list) {}

};

class G_array : public BaseArray {
public:
    G_array() : BaseArray() {}
    G_array(const std::initializer_list<double>& list) : BaseArray(list) {}

    double get_rho(double r) const override {
        double alpha = sqrt(
            (gamma*gamma - 1)*data[1]*data[1]
            + (gamma*gamma - 1)*data[2]*data[2]
            + gamma*gamma*data[3]*data[3]
            - 2*(gamma*gamma - 1)*data[0]*data[4]
        );
        double result =
            (data[0]*data[0]*(gamma*data[3] + alpha))
            / ((gamma - 1)*(2*data[0]*data[4] - data[1]*data[1] - data[2]*data[2]));
        return result * r;
    }

    double get_p(double r) const override {
        double alpha = sqrt(
            (gamma*gamma - 1)*data[1]*data[1]
            + (gamma*gamma - 1)*data[2]*data[2]
            + gamma*gamma*data[3]*data[3]
            - 2*(gamma*gamma - 1)*data[0]*data[4]
        );
        double result = (data[3] + alpha) / (gamma + 1);
        return result * r;
    }

    double get_u() const override {
        return data[1] / data[0];
    }

    double get_v() const override {
        return data[2] / data[0];
    }

    double get_w() const override {
        double alpha = sqrt(
            (gamma*gamma - 1)*data[1]*data[1]
            + (gamma*gamma - 1)*data[2]*data[2]
            + gamma*gamma*data[3]*data[3]
            - 2*(gamma*gamma - 1)*data[0]*data[4]
        );
        double result = (gamma*data[3] - alpha) / ((gamma + 1)*data[0]);
        return result;
    }
};

class R_array : public BaseArray {
public:
    R_array() : BaseArray() {}
    R_array(const std::initializer_list<double>& list) : BaseArray(list) {}
};

G_array predictor(E_array Ep, E_array Em, F_array Fp, F_array Fm, G_array G, R_array R, double dr, double dth, double dz){
    G_array result = G;
    result = result - dz/dr * (Ep - Em) - dz/dth * (Fp - Fm) + dz * R;
    return result;
}

G_array corrector(E_array Ep, E_array Em, F_array Fp, F_array Fm, G_array Gp, G_array Gm, R_array R, double dr, double dth, double dz){
    G_array result = (Gp + Gm)*0.5;
    result = result - 0.5*dz/dr * (Ep - Em) - 0.5*dz/dth * (Fp - Fm) + 0.5*dz * R;
    return result;
}

E_array get_E(G_array G, double r){
    double
        rho = G.get_rho(),
        p = G.get_p(),
        u = G.get_u(),
        v = G.get_v(),
        w = G.get_w();
    E_array result = {
        rho*u,
        rho*u*u + p,
        rho*u*v,
        rho*u*w,
        u*(gamma/(gamma - 1)*p + 0.5*rho*(u*u + v*v + w*w))
    };
    return r*result;
}

F_array get_F(G_array G){
    double
        rho = G.get_rho(),
        p = G.get_p(),
        u = G.get_u(),
        v = G.get_v(),
        w = G.get_w();
    F_array result = {
        rho*v,
        rho*u*v,
        rho*v*v + p,
        rho*v*w,
        v*(gamma/(gamma - 1)*p + 0.5*rho*(u*u + v*v + w*w))
    };
    return result;
}

R_array get_R(G_array G, double q){
    double
        rho = G.get_rho(),
        p = G.get_p(),
        u = G.get_u(),
        v = G.get_v();
    R_array result = {
        0,
        rho*v*v + p,
        -rho*u*v,
        0,
        q
    };
    return result;
}

double r_b(double z)
{
    return tan(PI / 6) * z;
}

double xi_r(double r_s, double r_b)
{
    return 1 / (r_s - r_b);
}

double xi_theta(double xi, double r_s, double r_b, double r_s_theta)
{
    return -xi * r_s_theta / (r_s - r_b);
}

double xi_z(double xi, double r_s, double r_b, double r_s_z, double r_b_z)
{
    return -(r_b_z + xi*(r_s_z - r_b_z)) / (r_s - r_b);
}

double q()
{
    return 0;
}

int main() {
    // Пример использования базового класса и производных
    G_array G = {1, 2, 3, 4, 5};

    G_array F;
    F = G;
    F.print();

    std::cout << G[0] << " " << G[4] << std::endl;
    std::cout << G.get_rho(3) << std::endl;
    std::cout << G.get_p(4) << std::endl;
    std::cout << G.get_u() << std::endl;


    int
        N = 10, // xi
        M = 10; // theta

    vector<vector<double>> r_s(M);

    double
        z0 = 1,
        r_b0, r_b_z0,
        r_s0, r_s_z0,
        h_xi = 1 / double(N - 1),
        h_th = PI / M,
        rho, p, u, v, w
        ;

    vector<vector<E_array>> E_prev(N, vector<E_array>(M)), E_next(N, vector<E_array>(M));
    vector<vector<F_array>> F_prev(N, vector<F_array>(M)), F_next(N, vector<F_array>(M));
    vector<vector<G_array>> G_prev(N, vector<G_array>(M)), G_next(N, vector<G_array>(M));
    vector<vector<R_array>> R_prev(N, vector<R_array>(M)), R_next(N, vector<R_array>(M));

    vector<double> phi_cone, rho_cone, p_cone, VR_cone, Vphi_cone;
    FILE *f_cone = fopen("rho.txt", "r");
    double a, b, c;
    while(!feof(f_cone))
    {
        fscanf(f_cone, "%lf %lf", &a, &b);
        phi_cone.push_back(a);
        rho_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("p.txt", "r");
    while(!feof(f_cone))
    {
        fscanf(f_cone, "%lf %lf", &a, &b);
        p_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("VR_Vtheta.txt", "r");
    while(!feof(f_cone))
    {
        fscanf(f_cone, "%lf %lf %lf", &a, &b, &c);
        VR_cone.push_back(b);
        Vphi_cone.push_back(c);
    }
    fclose(f_cone);

    double
        phi0 = phi_cone[phi_cone.size() - 1],
        phi1 = phi_cone[0],
        h_phi = (phi1 - phi0) / (phi_cone.size() - 1);

    r_s0 = z0 * tan(phi1);
    r_b0 = z0 * tan(phi0);
    r_s_z0 = tan(phi1);
    r_b_z0 = tan(phi0);

    int idx_phi;

    for(int i = 0, double r, xi = 0.0, phi; i < N; i++, xi += h_xi)
        for(int j = 0, double theta = 0.0; j < M; j++, theta += h_th)
        {
            r = xi*(r_s0 - r_b0) + r_b0;
            phi = atan(r / z0);
            idx_phi = phi_cone.size() - 1 - floor(phi / h_phi);

            rho = rho_cone[idx_phi];
            p = p_cone[idx_phi];
            u = VR_cone[idx_phi]*sin(phi_cone[idx_phi]) + Vphi_cone[idx_phi]*cos(phi_cone[idx_phi]);
            v = 0;
            w = VR_cone[idx_phi]*cos(phi_cone[idx_phi]) - Vphi_cone[idx_phi]*sin(phi_cone[idx_phi]);

            G_prev[i][j].data[0] = rho*w;
            G_prev[i][j].data[1] = rho*u*w;
            G_prev[i][j].data[2] = rho*v*w;
            G_prev[i][j].data[3] = rho*w*w + p;
            G_prev[i][j].data[4] = p / (gamma - 1) / rho + (u*u + v*v + w*w)*0.5;

            G_prev[i][j] = G_prev[i][j] * r;

            E_prev[i][j] = get_E(G_prev[i][j], r);
            F_prev[i][j] = get_F(G_prev[i][j]);
            R_prev[i][j] = get_R(G_prev[i][j], q());

            E_prev[i][j] =
                xi_r(r_s0, r_b0)*E_prev[i][j]
                + xi_theta(xi, r_s0, r_b0, 0)*F_prev[i][j]
                + xi_z(xi, r_s0, r_b0, r_s_z0, r_b_z0)*G_prev[i][j];
            
            R_prev[i][j] =
                R_prev[i][j]
                - 0/(r_s0 - r_b0) * F_prev[i][j]
                - (r_s_z0 - r_b_z0)/(r_s0 - r_b0) * G_prev[i][j];
        }

    return 0;
}
