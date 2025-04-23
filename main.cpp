#include "r.h"

int solver() {
    ofstream rho_out("rho_out.txt");
    ofstream p_out("p_out.txt");
    ofstream u_out("u_out.txt");
    ofstream v_out("v_out.txt");
    ofstream w_out("w_out.txt");

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
        Mach_inf = 3;

    double p_inf = 101330, rho_inf = 1.2255;
    double a0 = sqrt(gamma * p_inf / rho_inf);
    double V_inf = Mach_inf*a0;

    vector<vector<E_array>> E(N, vector<E_array>(M));
    vector<vector<F_array>> F(N, vector<F_array>(M));
    vector<vector<G_array>> G_prev(N, vector<G_array>(M)), G_next(N, vector<G_array>(M));
    vector<vector<R_array>> R(N, vector<R_array>(M));
    
    vector<vector<double>>
        rho_array(N, vector<double>(M)),
        p_array(N, vector<double>(M)),
        u_array(N, vector<double>(M)),
        v_array(N, vector<double>(M)),
        w_array(N, vector<double>(M));

    vector<vector<double>> r_s(M), r_s_theta(M), r_s_z(M);

    vector<double> phi_cone, rho_cone, p_cone, VR_cone, Vphi_cone;
    FILE *f_cone = fopen("rho.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        phi_cone.push_back(a);
        rho_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("p.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        p_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("VR_Vtheta.txt", "r");
    while(!feof(f_cone))
    {
        double a, b, c;
        fscanf(f_cone, "%lf %lf %lf", &a, &b, &c);
        VR_cone.push_back(b);
        Vphi_cone.push_back(c);
    }
    fclose(f_cone);

    double
        phi0 = phi_cone.back(),
        phi1 = phi_cone[0],
        dphi = (phi1 - phi0) / (phi_cone.size() - 1);

    std::cout << phi0 << " " << phi1 << " " << dphi << endl;
    r_b0 = z0 * tan(phi0);
    r_s0 = z0 * tan(phi1);
    r_b_z0 = tan(phi0);
    r_s_z0 = tan(phi1);

    for(int j = 0; j < M; j++)
    {
        r_s[j].push_back(r_s0);
        r_s_theta[j].push_back(0);
        r_s_z[j].push_back(r_s_z0);
    }

    double r, xi = 0, phi, theta, delta_th;
    for(int i = 0; i < N; i++)
    {
        int idx_phi;
        xi = i * dxi;
        r = r_from_xi(xi, r_s0, r_b0);
        phi = atan(r / z0);
        idx_phi = max(int(phi_cone.size()) - 1 - int(round((phi - phi0) / dphi)), 0);
        std::cout << int(phi_cone.size()) << " " << int(floor((phi - phi0) / dphi)) << endl;
        std::cout << "real phi: " << phi << ", approximated phi:" << phi_cone[idx_phi] << endl;

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

            E[i][j] = get_E(G_prev[i][j], r);
            F[i][j] = get_F(G_prev[i][j], r);
            R[i][j] = get_R(G_prev[i][j], r, q(r, theta, z0));

            // Normalization
            E[i][j] =
                xi_r(r_s0, r_b0)*E[i][j]
                + xi_theta(xi, r_s0, r_b0, 0)*F[i][j]
                + xi_z(xi, r_s0, r_b0, r_s_z0, r_b_z0)*G_prev[i][j];
            
            R[i][j] =
                R[i][j]
                - 0/(r_s0 - r_b0) * F[i][j]
                - (r_s_z0 - r_b_z0)/(r_s0 - r_b0) * G_prev[i][j];
        }
    }

    for(int i = 0; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s0, r_b0);
        for(int j = 0; j < M; j++){
            theta = j * dth;
            std::cout << "\nxi:" << xi << ", r: " << r << ", theta: " << theta << endl;
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

            E[i][j].print();
            F[i][j].print();
            G_prev[i][j].print();
            R[i][j].print();
        }
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

    double dz = 0.01, CFL = 0.9, lambda_xi, lambda_th;
    double MM, mm, xi_z_val, xi_r_val, xi_theta_val;
    vector<double> z_array;

    // main loop
    while(z < L){
        z_array.push_back(z);
        // dz calculation from Spectral method
        // for(int j = 0; j < M; j++)
        //     dr = min(dr, (r_s[j].back() - r_b(z))*dxi + r_b(z));
        for(int i = 0; i < N; i++){
            xi = i * dxi;
            for(int j = 0; j < M; j++){
                r = r_from_xi(xi, r_s[j].back(), r_b(z));
                xi_r_val = xi_r(r_s[j].back(), r_b(z));
                xi_theta_val = xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back());
                xi_z_val = xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z));
                double a = sqrt(gamma * p_array[i][j] / rho_array[i][j]);
                MM = w_array[i][j] / a;
                mm = (
                    u_array[i][j]*xi_r_val
                    + v_array[i][j]*xi_theta_val/r
                    + w_array[i][j]*xi_z_val
                ) / a;

                lambda_xi = (
                    abs(MM*mm - xi_z_val*xi_z_val)
                    + sqrt(
                        (mm - MM*xi_z_val)*(mm - MM*xi_z_val)
                        + (xi_r_val*xi_r_val + xi_theta_val*xi_theta_val/r/r)*(MM*MM - 1)
                    )
                ) / (MM*MM - 1);
                lambda_th = (
                    abs(
                        v_array[i][j]*w_array[i][j]/a/a
                    )
                    + sqrt(
                        (v_array[i][j]*v_array[i][j] + w_array[i][j]*w_array[i][j])/a/a - 1
                    )
                ) / (w_array[i][j]*w_array[i][j]/a/a - 1) / r;

                dz = min(dz, 1 / (lambda_xi/dxi + lambda_th/dth));
            }
        }
        // dz *= CFL;
        int idx;
        vector<double> r_s_theta_pred(M);
        for(int step = 0; step < 2; step++){ // step = 0: predictor, step = 1: corrector
            for(int i = 0; i < N; i++){
                xi = i*dxi;
                // Граница theta = 0
                theta = 0;
                r = r_from_xi(xi, r_s[0].back(), r_b(z));
                
                // predictor
                if(step == 0){
                    idx = int(i == N - 1);
                    G_next[i][0] = predictor(
                        E[i - idx][0],
                        E[i + 1 - idx][0],
                        F[i][0],
                        F[i][1],
                        G_prev[i][0],
                        R[i][0],
                        dxi, dth, dz
                    );
                }
                //corrector
                else{
                    F_array F_mirrored = (-1) * F[i][1];
                    F_mirrored[2] = (-1) * F_mirrored[2];
                    idx = int(i == 0);
                    G_next[i][0] = corrector(
                        E[i - 1 + idx][0],
                        E[i + idx][0],
                        F_mirrored, // отражение при theta -> -theta
                        F[i][0],
                        G_prev[i][0],
                        G_next[i][0],
                        R[i][0],
                        dxi, dth, dz
                    );
                }
                E[i][0] = get_E(G_next[i][0], r);
                F[i][0] = get_F(G_next[i][0], r);
                R[i][0] = get_R(G_next[i][0], r, q(r, theta, z));
    
    
                // Восстановление физических величин по вектору G
                rho_array[i][0] = G_next[i][0].get_rho(r);
                p_array[i][0] = G_next[i][0].get_p(r);
                u_array[i][0] = G_next[i][0].get_u();
                v_array[i][0] = G_next[i][0].get_v();
                w_array[i][0] = G_next[i][0].get_w();
    
                // Внутренние (по theta) узлы
                for(int j = 1; j < M - 1; j++){
                    theta = j*dth;
                    r = r_from_xi(xi, r_s[j].back(), r_b(z));
    
                    // predictor
                    if(step == 0){
                        idx = int(i == N - 1);
                        G_next[i][j] = predictor(
                            E[i - idx][j],
                            E[i + 1 - idx][j],
                            F[i][j],
                            F[i][j + 1],
                            G_prev[i][j],
                            R[i][j],
                            dxi, dth, dz
                        );
                    }
                    //corrector
                    else{
                        idx = int(i == 0);
                        G_next[i][j] = corrector(
                            E[i - 1 + idx][j],
                            E[i + idx][j],
                            F[i][j - 1],
                            F[i][j],
                            G_prev[i][j],
                            G_next[i][j],
                            R[i][j],
                            dxi, dth, dz
                        );
                    }
                    E[i][j] = get_E(G_next[i][j], r);
                    F[i][j] = get_F(G_next[i][j], r);
                    R[i][j] = get_R(G_next[i][j], r, q(r, theta, z));
    
                    // Восстановление физических величин по вектору G
                    rho_array[i][j] = G_next[i][j].get_rho(r);
                    p_array[i][j] = G_next[i][j].get_p(r);
                    u_array[i][j] = G_next[i][j].get_u();
                    v_array[i][j] = G_next[i][j].get_v();
                    w_array[i][j] = G_next[i][j].get_w();
                }
    
                // Граница theta = PI
                theta = PI;
                r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
    
                // predictor
                if(step == 0){
                    F_array F_mirrored = (-1) * F[i][M - 2];
                    F_mirrored[2] = (-1) * F_mirrored[2];
                    idx = int(i == N - 1);
                    G_next[i][M - 1] = predictor(
                        E[i - idx][M - 1],
                        E[i + 1 - idx][M - 1],
                        F[i][M - 1],
                        F_mirrored, // Отражение
                        G_prev[i][M - 1],
                        R[i][M - 1],
                        dxi, dth, dz
                    );
                }
                else{
                    idx = int(i == 0);
                    G_next[i][M - 1] = corrector(
                        E[i - 1 + idx][M - 1],
                        E[i + idx][M - 1],
                        F[i][M - 2],
                        F[i][M - 1],
                        G_prev[i][M - 1],
                        G_next[i][M - 1],
                        R[i][M - 1],
                        dxi, dth, dz
                    );
                }
                E[i][M - 1] = get_E(G_next[i][M - 1], r);
                F[i][M - 1] = get_F(G_next[i][M - 1], r);
                R[i][M - 1] = get_R(G_next[i][M - 1], r, q(r, theta, z));

                // Восстановление физических величин по вектору G
                rho_array[i][M - 1] = G_next[i][M - 1].get_rho(r);
                p_array[i][M - 1] = G_next[i][M - 1].get_p(r);
                u_array[i][M - 1] = G_next[i][M - 1].get_u();
                v_array[i][M - 1] = G_next[i][M - 1].get_v();
                w_array[i][M - 1] = G_next[i][M - 1].get_w();
            }
            
            // i = 0 and i = N - 1

            // Метод Аббета (поправка на поверхности тела)
            for(int j = 0; j < M; j++){
                theta = j * dth;
                double V_n = (u_array[0][j] - r_b_z(z)*w_array[0][j])
                    / sqrt(1 + r_b_z(z)*r_b_z(z)); // (V*, n)
                double Mach = sqrt(
                    (
                        u_array[0][j]*u_array[0][j]
                        + v_array[0][j]*v_array[0][j]
                        + w_array[0][j]*w_array[0][j]
                    )
                    /
                    (
                        gamma * p_array[0][j] / rho_array[0][j]
                    )
                );
                delta_th = asin(
                    V_n / sqrt(
                        u_array[0][j]*u_array[0][j]
                        + v_array[0][j]*v_array[0][j]
                        + w_array[0][j]*w_array[0][j]
                    )
                );

                // Поправка давления
                p_array[0][j] = p_array[0][j] * (
                    1 - gamma*Mach*Mach*delta_th/sqrt(Mach*Mach - 1)
                    + gamma
                        * Mach
                        * ((gamma + 1)*Mach*Mach*Mach*Mach - 4*(Mach*Mach - 1))
                        / (4*(Mach*Mach - 1)*(Mach*Mach - 1))
                        * delta_th * delta_th
                );
    
                //Поправка плотности
                rho_array[0][j] = rho_inf * pow(p_array[0][j] / p_inf, 1/gamma);
    
                double H = gamma/(gamma - 1) * p_inf/rho_inf + 0.5*V_inf*V_inf;
                double V_abs = sqrt(2*(H - gamma/(gamma - 1) * p_array[0][j]/rho_array[0][j]));
    
                double Vr_tau, Vth_tau, Vz_tau;
                // V_tau:
                Vr_tau = u_array[0][j] - V_n / sqrt(1 + r_b_z(z)*r_b_z(z));
                Vth_tau = v_array[0][j];
                Vz_tau = w_array[0][j] + V_n * r_b_z(z) / sqrt(1 + r_b_z(z)*r_b_z(z));
    
                double V_tau_abs = sqrt(Vr_tau*Vr_tau + Vth_tau*Vth_tau + Vz_tau*Vz_tau);
                //Поправка скорости
                u_array[0][j] = Vr_tau * V_abs / V_tau_abs;
                v_array[0][j] = Vth_tau * V_abs / V_tau_abs;
                w_array[0][j] = Vz_tau * V_abs / V_tau_abs;
            }
            
            // ПОПРАВКА ВЕКТОРОВ E F G R !!!!!!
            // ПОПРАВКА ВЕКТОРОВ E F G R !!!!!!
            // ПОПРАВКА ВЕКТОРОВ E F G R !!!!!!

            // Метод Томаса (поправка на поверхности ударной волны)
            for(int j = 0; j < M; j++){
                // Давление не меняем
                // Плотность из Р-Г
                rho_array[N - 1][j] = (
                    rho_inf
                    * ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)
                    / ((gamma + 1)*p_inf + (gamma - 1)*p_array[N - 1][j])
                );
                double V_inf_n = sqrt(
                    ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)
                    / (2 * rho_inf)
                );
                double V_n = rho_inf / rho_array[N - 1][j] * V_inf_n;
                // double beta = asin(V_inf_n / V_inf);
                // r_s_z[j].push_back(tan(beta));

                // r_s, r_s_theta, r_s_z calculation (pred-corr)
                if(step == 0){
                    r_s[j].push_back(r_s[j].back() + r_s_z[j].back()*dz);
                    r_s_theta_pred[j] = r_s_theta[j].back();
                    if(j == 0){
                        r_s_theta[j].push_back(
                            r_s_theta[j].back() + dz/dth*(r_s_z[j].back() - r_s_z[j + 1].back()) // симметрия j - 1 --> j + 1
                        );
                    }
                    else{
                        r_s_theta[j].push_back(
                            r_s_theta[j].back() + dz/dth*(r_s_z[j].back() - r_s_z[j - 1].back())
                        );
                    }
                    // new or old???????
                    r_s_z[j].push_back(
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (V_inf*V_inf - V_inf_n*V_inf_n)
                        )
                    );
                }
                else{
                    // Нужно ли корректировать r_s???
                    // r_s[j].back() = r_s[j].back() + r_s_z[j].back()*dz;
                    if(j < M - 1){
                        r_s_theta[j].back() = (
                            0.5 * (r_s_theta_pred[j] + r_s_theta[j].back())
                            + 0.5 * dz/dth * (r_s_z[j + 1].back() - r_s_z[j].back())
                        );
                    }
                    else{
                        r_s_theta[j].back() = (
                            0.5 * (r_s_theta_pred[j] + r_s_theta[j].back())
                            + 0.5 * dz/dth * (r_s_z[j - 1].back() - r_s_z[j].back()) // симметрия
                        );
                    }
                    r_s_z[j].back() = (
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (V_inf*V_inf - V_inf_n*V_inf_n)
                        )
                    );
                }

                // Нормаль к поверхности УВ
                double nx, ny, nz, n_norm; // n = (nx, ny, nz)
                n_norm = (1 / sqrt(
                        1
                        + r_s_z[j].back()*r_s_z[j].back()
                        + (r_s_theta[j].back()/r_s[j].back())*(r_s_theta[j].back()/r_s[j].back())
                    )
                );
                nx = 1 * n_norm;
                ny = -r_s_theta[j].back() / r_s[j].back() * n_norm;
                nz = -r_s_z[j].back() * n_norm;
                
                // Касательный вектор к поверхности УВ
                double V_inf_tau_x, V_inf_tau_y, V_inf_tau_z;
                V_inf_tau_x = 0 - V_inf_n * nx;
                V_inf_tau_y = 0 - V_inf_n * ny;
                V_inf_tau_z = V_inf - V_inf_n * nz;
                // double V_inf_n = sqrt(                                Проверить разницу V_inf_n
                //     ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)
                //     / (2 * rho_inf)
                // );
                u_array[N - 1][j] = V_inf_tau_x + V_n*nx;
                v_array[N - 1][j] = V_inf_tau_y + V_n*ny;
                w_array[N - 1][j] = V_inf_tau_z + V_n*nz;
            }
    
            //Поправка векторов E, F, G, R на границах
            for(int i = 0; i < N; i += N - 1){
                xi = i * dxi;
                for(int j = 0; j < M; j++){
                    theta = j * dth;
                    r = r_from_xi(xi, r_s[j].back(), r_b(z));

                    G_next[i][j].data[0] = rho_array[i][j]*w_array[i][j];
                    G_next[i][j].data[1] = rho_array[i][j]*u_array[i][j]*w_array[i][j];
                    G_next[i][j].data[2] = rho_array[i][j]*v_array[i][j]*w_array[i][j];
                    G_next[i][j].data[3] = rho_array[i][j]*w_array[i][j]*w_array[i][j] + p_array[i][j];
                    G_next[i][j].data[4] = (
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
                    G_next[i][j] = G_next[i][j] * r;
                    
                    E[i][j] = get_E(G_next[i][j], r);
                    F[i][j] = get_F(G_next[i][j], r);
                    R[i][j] = get_R(G_next[i][j], r, q(r, theta, z));
                }
            }
            // Нормализация
            for(int i = 0; i < N; i++){
                xi = i * dxi;
                for(int j = 0; j < M; j++){
                    theta = j * dth;
                    r = r_from_xi(xi, r_s[j].back(), r_b(z));
                    E[i][j] = (
                        xi_r(r_s[j].back(), r_b(z))*E[i][j]
                        + xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back())*F[i][j]
                        + xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z))*G_next[i][j]
                    );

                    R[i][j] = (
                        R[i][j]
                        - r_s_theta[j].back() / (r_s[j].back() - r_b(z)) * F[i][j]
                        - (r_s_z[j].back() - r_b_z(z)) / (r_s[j].back() - r_b(z)) * G_next[i][j]
                    );

                    // обновление G_prev <-- G_next после корректора (после окончания полного шага)
                    if(step == 1)
                        G_prev[i][j] = G_next[i][j];
                }
            }
        }
        // Запись в файл
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < M; j++)
            {
                rho_out << rho_array[i][j] << " ";
                p_out << p_array[i][j] << " ";
                u_out << u_array[i][j] << " ";
                v_out << v_array[i][j] << " ";
                w_out << w_array[i][j] << " ";
            }
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

        z += dz;

        cout << "dz = " << dz << endl;
    }
    return 0;
}

void test()
{
    // Пример использования базового класса и производных
    // G_array G = {1, 2, 3, 4, 5};

    // G_array F;
    // F = G;
    // F.print();

    // std::cout << G[0] << " " << G[4] << std::endl;
    // std::cout << G.get_rho(3) << std::endl;
    // std::cout << G.get_p(4) << std::endl;
    // std::cout << G.get_u() << std::endl;

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
    std::cout << G.get_rho(r) << endl;
    std::cout << G.get_p(r) << endl;
    std::cout << G.get_u() << endl;
    std::cout << G.get_v() << endl;
    std::cout << G.get_w() << endl;

    G.print();
    std::cout << G.get_rho(r) << endl;

    vector<double> vec(5, 1);
    std::cout << vec.back() << endl;
    vec.back() = 5;
    for(auto elem: vec)
        std::cout << elem << " ";
}

int main()
{
    solver();
    // test();

    return 0;
}