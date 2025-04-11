#include "r.h"

int solver() {
    // Пример использования базового класса и производных
    // G_array G = {1, 2, 3, 4, 5};

    // G_array F;
    // F = G;
    // F.print();

    // std::cout << G[0] << " " << G[4] << std::endl;
    // std::cout << G.get_rho(3) << std::endl;
    // std::cout << G.get_p(4) << std::endl;
    // std::cout << G.get_u() << std::endl;
    ofstream rho_out("rho.txt");
    ofstream p_out("p.txt");
    ofstream u_out("u.txt");
    ofstream v_out("v.txt");
    ofstream w_out("w.txt");

    int
        N = 10, // xi
        M = 8; // theta

    double
        z0 = 1, z = z0, L = 3*z,
        r_b0, r_b_z0,
        r_s0, r_s_z0,
        dxi = 1 / double(N - 1),
        dth = PI / double(M - 1),
        rho, p, u, v, w,
        Mach = 3;

    double p0 = 101330, rho0 = 1.2255;
    double a0 = sqrt(gamma * p0 / rho0);
    double V_inf = Mach*a0;

    vector<vector<E_array>> E_prev(N, vector<E_array>(M)), E_next(N, vector<E_array>(M));
    vector<vector<F_array>> F_prev(N, vector<F_array>(M)), F_next(N, vector<F_array>(M));
    vector<vector<G_array>> G_prev(N, vector<G_array>(M)), G_next(N, vector<G_array>(M));
    vector<vector<R_array>> R_prev(N, vector<R_array>(M)), R_next(N, vector<R_array>(M));
    
    vector<vector<double>>
        rho_array(N, vector<double>(M)),
        p_array(N, vector<double>(M)),
        u_array(N, vector<double>(M)),
        v_array(N, vector<double>(M)),
        w_array(N, vector<double>(M));

    vector<vector<double>> r_s(M), r_s_theta(M), r_s_z(M);

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
        dphi = (phi1 - phi0) / (phi_cone.size() - 1);

    cout << phi0 << " " << phi1 << " " << dphi << endl;
    r_s0 = z0 * tan(phi1);
    r_b0 = z0 * tan(phi0);
    r_s_z0 = tan(phi1);
    r_b_z0 = tan(phi0);

    for(int j = 0; j < M; j++)
    {
        r_s[j].push_back(r_s0);
        r_s_theta[j].push_back(0);
        r_s_z[j].push_back(r_s_z0);
    }

    int idx_phi;

    double r = 0, xi = 0, phi, theta, delta_th;
    for(int i = 0; i < N; i++)
    {
        r = r_from_xi(xi, r_s0, r_b0);
        phi = atan(r / z0);
        idx_phi = max(int(phi_cone.size()) - 1 - int(round((phi - phi0) / dphi)), 0);
        cout << int(phi_cone.size()) << " " << int(floor((phi - phi0) / dphi)) << endl;
        cout << "real phi: " << phi << ", approximated phi:" << phi_cone[idx_phi] << endl;

        rho = rho_cone[idx_phi];
        p = p_cone[idx_phi];
        u = VR_cone[idx_phi]*sin(phi_cone[idx_phi]) + Vphi_cone[idx_phi]*cos(phi_cone[idx_phi]);
        v = 0;
        w = VR_cone[idx_phi]*cos(phi_cone[idx_phi]) - Vphi_cone[idx_phi]*sin(phi_cone[idx_phi]);
        printf(
            "rho = %lf, p = %lf, u = %lf, v = %lf, w = %lf\n",
            rho, p, u, v, w);
        for(int j = 0; j < M; j++)
        {
            theta = j * dth;
            G_prev[i][j].data[0] = rho*w;
            G_prev[i][j].data[1] = rho*u*w;
            G_prev[i][j].data[2] = rho*v*w;
            G_prev[i][j].data[3] = rho*w*w + p;
            G_prev[i][j].data[4] = w * (gamma * p / (gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
            
            G_prev[i][j] = G_prev[i][j] * r;

            E_prev[i][j] = get_E(G_prev[i][j], r);
            F_prev[i][j] = get_F(G_prev[i][j], r);
            R_prev[i][j] = get_R(G_prev[i][j], r, q(r, theta, z0));

            E_prev[i][j] =
                xi_r(r_s0, r_b0)*E_prev[i][j]
                + xi_theta(xi, r_s0, r_b0, 0)*F_prev[i][j]
                + xi_z(xi, r_s0, r_b0, r_s_z0, r_b_z0)*G_prev[i][j];
            
            R_prev[i][j] =
                R_prev[i][j]
                - 0/(r_s0 - r_b0) * F_prev[i][j]
                - (r_s_z0 - r_b_z0)/(r_s0 - r_b0) * G_prev[i][j];
        }
        xi += dxi;
    }

    xi = 0;
    for(int i = 0; i < N; i++)
    {
        theta = 0;
        r = r_from_xi(xi, r_s0, r_b0);
        for(int j = 0; j < M; j++)
        {
            theta = j * dth;
            cout << "\nxi:" << xi << ", r: " << r << ", theta: " << theta << endl;
            printf(
                "\nrho = %lf, p = %lf, u = %lf, v = %lf, w = %lf\n",
                G_prev[i][j].get_rho(r),
                G_prev[i][j].get_p(r),
                G_prev[i][j].get_u(),
                G_prev[i][j].get_v(),
                G_prev[i][j].get_w()
            );

            rho_array[i][j] = G_prev[i][j].get_rho(r);
            p_array[i][j] = G_prev[i][j].get_p(r);
            u_array[i][j] = G_prev[i][j].get_u();
            v_array[i][j] = G_prev[i][j].get_v();
            w_array[i][j] = G_prev[i][j].get_w();

            // Запись в файл!!!
            rho_out << rho_array[i][j] << " ";
            p_out << p_array[i][j] << " ";
            u_out << u_array[i][j] << " ";
            v_out << v_array[i][j] << " ";
            w_out << w_array[i][j] << " ";

            E_prev[i][j].print();
            F_prev[i][j].print();
            G_prev[i][j].print();
            R_prev[i][j].print();
            theta += dth;
        }
        xi += dxi;
        rho_out << "\n";
        p_out << "\n";
        u_out << "\n";
        v_out << "\n";
        w_out << "\n";
    }
    
    rho_out << "\n";
    p_out << "\n";
    u_out << "\n";
    v_out << "\n";
    w_out << "\n";

    double dr, dz;
    // int k = 0;
    vector<double> z_array;
    while(z < L){
        z_array.push_back(z);   
        // dz calculation from Spectral method
        for(int j = 0; j < M; j++)
            dr = min(dr, (r_s[j].back() - r_b(z))*dxi + r_b(z));
        for(int i = 0; i < N; i++)
            for(int j = 0; j < M; j++)
                dz = min(dz, 
                        lambda_r(
                            rho_array[i][j],
                            p_array[i][j],
                            u_array[i][j],
                            v_array[i][j],
                            w_array[i][j]
                            )/dr
                        + lambda_th(
                            r,
                            rho_array[i][j],
                            p_array[i][j],
                            u_array[i][j],
                            v_array[i][j],
                            w_array[i][j]
                            )/dth
                        );
        dz *= 0.9;
        
        xi = dxi;
        for(int i = 1; i < N - 1; i++){
            // Граница theta = 0
            theta = 0;
            r = r_from_xi(xi, r_s[0].back(), r_b(z));

            // predictor
            G_next[i][0] = predictor(
                E_prev[i][0],
                E_prev[i + 1][0],
                F_prev[i][0],
                F_prev[i][1],
                G_prev[i][0],
                R_prev[i][0],
                dr, dth, dz
            );
            E_next[i][0] = get_E(G_next[i][0], r);
            F_next[i][0] = get_F(G_next[i][0], r);
            R_next[i][0] = get_R(G_next[i][0], r, q(r, theta, z));

            //corrector
            G_next[i][0] = corrector(
                E_next[i - 1][0],
                E_next[i][0],
                F_next[i][1], // симметрия при theta -> -theta
                F_next[i][0],
                G_prev[i][0],
                G_next[i][0],
                R_next[i][0],
                dr, dth, dz
            );

            // Восстановление физических величин по вектору G
            rho_array[i][0] = G_next[i][0].get_rho(r);
            p_array[i][0] = G_next[i][0].get_p(r);
            u_array[i][0] = G_next[i][0].get_u();
            v_array[i][0] = G_next[i][0].get_v();
            w_array[i][0] = G_next[i][0].get_w();

            theta += dth;

            // Внутренние (по theta) узлы
            for(int j = 1; j < M - 1; j++){
                r = r_from_xi(xi, r_s[j].back(), r_b(z));

                // predictor
                G_next[i][j] = predictor(
                    E_prev[i][j],
                    E_prev[i + 1][j],
                    F_prev[i][j],
                    F_prev[i][j + 1],
                    G_prev[i][j],
                    R_prev[i][j],
                    dr, dth, dz
                );
                E_next[i][j] = get_E(G_next[i][j], r);
                F_next[i][j] = get_F(G_next[i][j], r);
                R_next[i][j] = get_R(G_next[i][j], r, q(r, theta, z));

                //corrector
                G_next[i][j] = corrector(
                    E_next[i - 1][j],
                    E_next[i][j],
                    F_next[i][j - 1],
                    F_next[i][j],
                    G_prev[i][j],
                    G_next[i][j],
                    R_next[i][j],
                    dr, dth, dz
                );

                // Восстановление физических величин по вектору G
                rho_array[i][j] = G_next[i][j].get_rho(r);
                p_array[i][j] = G_next[i][j].get_p(r);
                u_array[i][j] = G_next[i][j].get_u();
                v_array[i][j] = G_next[i][j].get_v();
                w_array[i][j] = G_next[i][j].get_w();

                theta += dth;
            }

            // Граница theta = PI
            r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));

            // predictor
            G_next[i][M - 1] = predictor(
                E_prev[i][M - 1],
                E_prev[i + 1][M - 1],
                F_prev[i][M - 1],
                F_prev[i][M - 2],
                G_prev[i][M - 1],
                R_prev[i][M - 1],
                dr, dth, dz
            );
            E_next[i][M - 1] = get_E(G_next[i][M - 1], r);
            F_next[i][M - 1] = get_F(G_next[i][M - 1], r);
            R_next[i][M - 1] = get_R(G_next[i][M - 1], r, q(r, theta, z));

            //corrector
            G_next[i][M - 1] = corrector(
                E_next[i - 1][M - 1],
                E_next[i][M - 1],
                F_next[i][M - 2], // симметрия при theta -> -theta
                F_next[i][M - 1],
                G_prev[i][M - 1],
                G_next[i][M - 1],
                R_next[i][M - 1],
                dr, dth, dz
            );

            // Восстановление физических величин по вектору G
            rho_array[i][M - 1] = G_next[i][M - 1].get_rho(r);
            p_array[i][M - 1] = G_next[i][M - 1].get_p(r);
            u_array[i][M - 1] = G_next[i][M - 1].get_u();
            v_array[i][M - 1] = G_next[i][M - 1].get_v();
            w_array[i][M - 1] = G_next[i][M - 1].get_w();

            xi += dxi;
        }
        
        // i = 0 and i = N - 1
        
        // Метод Аббета (поправка на поверхности тела)
        for(int j = 0; j < M; j++){
            double V_n = (u_array[0][j] - r_b_z(z, dz)*w_array[0][j])
                / sqrt(1 + r_b_z(z, dz)*r_b_z(z, dz)); // (V*, n)
            delta_th = asin(
                V_n / (
                    u_array[0][j]*u_array[0][j]
                    + v_array[0][j]*v_array[0][j]
                    + w_array[0][j]*w_array[0][j]
                )
            );
            theta += dth;

            // Поправка давления
            p_array[0][j] = p_array[0][j] * (
                1 + gamma*M*M*delta_th/sqrt(M*M - 1)
                + gamma*M*(((gamma + 1)*M*M*M*M - 4*(M*M - 1)) / (4*(M*M - 1)*(M*M - 1)))*delta_th*delta_th
            );

            //Поправка плотности
            rho_array[0][j] = rho0 * pow(p_array[0][j] / p0, 1/gamma);

            double H = gamma/(gamma - 1) * p0/rho0 + V_inf*V_inf;
            double V_abs = sqrt(2*(H - gamma/(gamma - 1) * p_array[0][j]/rho_array[0][j]));

            double Vr, Vth, Vz;
            // V_tau:
            Vr = u_array[0][j] - V_n / sqrt(1 + r_b_z(z, dz)*r_b_z(z, dz));
            Vth = v_array[0][j];
            Vz = w_array[0][j] + V_n * r_b_z(z, dz) / sqrt(1 + r_b_z(z, dz)*r_b_z(z, dz));

            double V_tau_abs = Vr*Vr + Vth*Vth + Vz*Vz;
            //Поправка скорости
            u_array[0][j] *= V_abs / V_tau_abs;
            v_array[0][j] *= V_abs / V_tau_abs;
            w_array[0][j] *= V_abs / V_tau_abs;
        }

        // Метод Томаса (поправка на поверхности ударной волны)
        for(int j = 0; j < M; j++){
            // Давление не меняем
            // Плотность из Р-Г
            rho_array[N - 1][j] = (
                rho0
                * ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p0)
                / ((gamma + 1)*p0 + (gamma - 1)*p_array[N - 1][j])
            );
            double V_inf_n = sqrt(
                ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p0)
                / (2 * rho0)
            );
            double V_n = rho0 / rho_array[N - 1][j] * V_inf_n;
            double beta = asin(V_inf_n / V_inf);
            r_s_z[j].push_back(tan(beta));
            r_s[j].push_back(r_s[j].back() + r_s_z[j].back()*dz); // ?????????
        }
        double r_s_theta_next;
        // Маккормак для r_s_theta
        r_s_theta_next = r_s_theta[0].back() + dz/dth * (r_s_z[0].back() - r_s_z[1].back()); // predictor
        // r_s_theta_next = 0.5*(r_s_theta_next + r_s_theta[0].back()) + 0.5*dz/dth * () ?????????????
        r_s_theta[0].push_back(r_s_theta_next);




        double nx, ny, nz, n_norm;
        for(int j = 0; j < M; j++){
            n_norm = (1 
                / sqrt(
                    1
                    + r_s_z[j].back()*r_s_z[j].back()
                    + (r_s_theta[j].back()/r_s[j].back())*(r_s_theta[j].back()/r_s[j].back())
                )
            );
            nx = 1 * n_norm;
            ny = -r_s_theta[j].back() / r_s[j].back() * n_norm;
            nz = -r_s_z[j].back() * n_norm;
            double V_inf_n = V_inf * nz;
            double V_inf_tau_x, V_inf_tau_y, V_inf_tau_z;
            V_inf_tau_x = -V_inf_n * nx;
            V_inf_tau_y = -V_inf_n * ny;
            V_inf_tau_z = V_inf - V_inf_n * nz;
            // double V_inf_n = sqrt(                                Проверить разницу V_inf_n
            //     ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p0)
            //     / (2 * rho0)
            // );
            double V_n = rho0 / rho_array[N - 1][j] * V_inf_n;
            u_array[N - 1][j] = V_inf_tau_x + V_n*nx;
            v_array[N - 1][j] = V_inf_tau_y + V_n*ny;
            w_array[N - 1][j] = V_inf_tau_z + V_n*nz;
        }
        //Поправка векторов E, F, G, R
        for(int i = 0; i < N; i += N - 1)
        {
            xi = i * dxi;
            for(int j = 0; j < M; j++){
                theta = j * dth;
                r = r_from_xi(xi, r_s[0].back(), r_b(z));

                G_prev[i][j].data[0] = rho_array[i][j]*w_array[i][j];
                G_prev[i][j].data[1] = rho_array[i][j]*u_array[i][j]*w_array[i][j];
                G_prev[i][j].data[2] = rho_array[i][j]*v_array[i][j]*w_array[i][j];
                G_prev[i][j].data[3] = rho_array[i][j]*w_array[i][j]*w_array[i][j] + p_array[i][j];
                G_prev[i][j].data[4] = (
                    w_array[i][j]
                    * (
                        gamma * p_array[i][j] / (gamma - 1) 
                        + rho_array[i][j]*(
                            u_array[i][j]*u_array[i][j]
                            + v_array[i][j]*v_array[i][j]
                            + w_array[i][j]*w_array[i][j]
                        )*0.5
                    )
                );
                
                E_prev[i][j] = get_E(G_prev[i][j], r);
                F_prev[i][j] = get_F(G_prev[i][j], r);
                R_prev[i][j] = get_R(G_prev[i][j], r, q(r, theta, z));

                // Нормализация
                E_prev[i][j] =
                    xi_r(r_s[j].back(), r_b(z))*E_prev[i][j]
                    + xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back())*F_prev[i][j]
                    + xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z, dz))*G_prev[i][j];
        
                R_prev[i][j] =
                    R_prev[i][j]
                    - r_s_theta[j].back()/(r_s[j].back() - r_b(z)) * F_prev[i][j]
                    - (r_s_z[j].back() - r_b_z(z, dz))/(r_s[j].back() - r_b(z)) * G_prev[i][j];
            }
        }
        // Запись в файл



        z += dz;
    }

    return 0;
}

void test()
{
    E_array E = {1,2,3,4,5};
    F_array F = {5,4,3,2,1};
    G_array G = {9,9,9,9,0};
    G_array R = {1,1,2,2,2};
    double rho = 1.2, p = 100000, u = 100, v = 200, w = 300, r = 3;
    rho = 1.615063, p = 149129.877770, u = 93.149968, v = 0.000000, w = 621.925338;

    G.data[0] = rho*w;
    G.data[1] = rho*u*w;
    G.data[2] = rho*v*w;
    G.data[3] = rho*w*w + p;
    G.data[4] = w * (gamma * p / (gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
    G = r * G;
    G.print();
    
    // G = 0.1*E + 1*G + 0.03*F;
    // G.print();
    // G = 1.0*G - 0.1*E - F*0.03;

    // (2*E).print();
    // (a*E + b*R).print();
    
    G.print();
    cout << G.get_rho(r) << endl;
    cout << G.get_p(r) << endl;
    cout << G.get_u() << endl;
    cout << G.get_v() << endl;
    cout << G.get_w() << endl;

    G.print();
    cout << G.get_rho(r) << endl;
}

int main()
{
    // solver();
    test();

    return 0;
}