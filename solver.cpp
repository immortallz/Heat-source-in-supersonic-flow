#include "r.h"
#include <string>

std::vector<double> solver(HeatSource heatSource) {
    std::ofstream
        z_out("output_results/z_out.txt"),
        rho_out("output_results/rho_out.txt"),
        p_out("output_results/p_out.txt"),
        u_out("output_results/u_out.txt"),
        v_out("output_results/v_out.txt"),
        w_out("output_results/w_out.txt"),
        r_s_out("output_results/r_s_out.txt"),
        r_s_theta_out("output_results/r_s_theta_out.txt"),
        r_s_z_out("output_results/r_s_z_out.txt"),
        psi0_out("output_results/psi0_out.txt"),
        psi1_out("output_results/psi1_out.txt"),
        Fy_out("output_results/Fy_out.txt"),
        Mz_out("output_results/Mz_out.txt");
    
    int
        N = numericalParams.N, // xi
        M = numericalParams.M; // theta

    bool progress_bar = true;
    int num_step_percent = numericalParams.num_step_percent;
    std::vector<bool> progress_flag(num_step_percent, false);
    bool z_limit = true;
    int k = 1, files_count = numericalParams.files_count; // сколько записей в файл в единице по z

    double
        z0 = bodyParams.transitionPoint, z = z0, L = bodyParams.bodyLength,
        r_b0, r_b_z0,
        r_s0, r_s_z0,
        dxi = 1 / double(N - 1),
        dth = Pi / double(M - 1),
        rho, p, u, v, w;

    double Q = 0; // integral of q(r,theta,z)dr*dtheta*dz по всей расчетной области

    std::vector<std::vector<E_array>> E(N, std::vector<E_array>(M));
    std::vector<std::vector<F_array>> F(N, std::vector<F_array>(M));
    std::vector<std::vector<G_array>>
        G_prev(N, std::vector<G_array>(M)),
        G_next(N, std::vector<G_array>(M));
    std::vector<std::vector<R_array>> R(N, std::vector<R_array>(M));
    
    std::vector<std::vector<double>>
        rho_array(N, std::vector<double>(M)),
        p_array(N, std::vector<double>(M)),
        u_array(N, std::vector<double>(M)),
        v_array(N, std::vector<double>(M)),
        w_array(N, std::vector<double>(M));

    std::vector<std::vector<double>> r_s(M), r_s_theta(M), r_s_z(M);
    std::vector<double> psi0(N, 0), psi1(N, 0);

    std::vector<double> phi_cone, rho_cone, p_cone, VR_cone, Vphi_cone; // Векторы для начальных данных

    // Чтение (заранее заготовленных) начальных данных
    // из задачи об обтекании конуса
    FILE *f_cone = fopen("initial_cone_flow/output_results/rho_init.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        phi_cone.push_back(a);
        rho_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("initial_cone_flow/output_results/p_init.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        p_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("initial_cone_flow/output_results/VR_Vtheta_init.txt", "r");
    while(!feof(f_cone))
    {
        double a, b, c;
        fscanf(f_cone, "%lf %lf %lf", &a, &b, &c);
        VR_cone.push_back(b);
        Vphi_cone.push_back(c);
    }
    fclose(f_cone);

    double
        phi0 = phi_cone.back(), // Ввиду специфики решения уравнения обтекания конуса
        phi1 = phi_cone[0],     // угол phi отсчитывается наоборот: от УВ до поверхности тела
        dphi = (phi1 - phi0) / (phi_cone.size() - 1); // шаг постоянен
    
    // Начальные значения r_s, r_s_theta, r_s_z
    r_b0 = z0 * tan(phi0); // MUST BE EQUIVIVALENT TO r_b(z0)
    r_s0 = z0 * tan(phi1);
    r_b_z0 = tan(phi0); // MUST BE EQUIVIVALENT TO r_b_z(z0)
    r_s_z0 = tan(phi1);
    for(int j = 0; j < M; j++)
    {
        r_s[j].push_back(r_s0);
        r_s_theta[j].push_back(0);
        r_s_z[j].push_back(r_s_z0);
    }

    double r, xi, theta;
    double seconds;
    clock_t start_initial = clock();
    // Инициализация начального слоя
    for(int i = 0; i < N; i++)
    {
        int idx_phi;
        xi = i * dxi;
        r = r_from_xi(xi, r_s0, r_b0);

        // Синхронизация мелкой сетки (задача о конусе) с крупной сеткой
        double phi = atan(r / z0);
        // Ввиду обратного хода вычисления в задаче о конусе реверсивная индексация
        idx_phi = int(phi_cone.size()) - 1 - int((phi - phi0) / dphi);
        if(idx_phi < 0 || idx_phi >= int(phi_cone.size()))
            throw std::out_of_range("Индекс idx_phi = " + std::to_string(idx_phi) + " вне диапазона");

        // Инициализация/переход к переменным настоящей задачи
        rho = rho_cone[idx_phi];
        p = p_cone[idx_phi];
        u = VR_cone[idx_phi]*sin(phi_cone[idx_phi]) + Vphi_cone[idx_phi]*cos(phi_cone[idx_phi]);
        v = 0;
        w = VR_cone[idx_phi]*cos(phi_cone[idx_phi]) - Vphi_cone[idx_phi]*sin(phi_cone[idx_phi]);

        // Инициализация + Нормализация вектор-функций
        for(int j = 0; j < M; j++)
        {
            theta = j * dth;
            G_prev[i][j].data[0] = rho*w;
            G_prev[i][j].data[1] = rho*u*w;
            G_prev[i][j].data[2] = rho*v*w;
            G_prev[i][j].data[3] = rho*w*w + p;
            G_prev[i][j].data[4] = w * (Gamma * p / (Gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
            
            G_prev[i][j] = G_prev[i][j] * r;

            E[i][j] = get_E(G_prev[i][j], r);
            F[i][j] = get_F(G_prev[i][j], r);
            R[i][j] = get_R(G_prev[i][j], r, q(r, theta, z0, heatSource, flowParams.is_adiabatic));

            // Normalization
            E[i][j] =
                xi_r(r_s0, r_b0)*E[i][j]
                + xi_theta(xi, r_s0, r_b0, 0)*F[i][j]
                + xi_z(xi, r_s0, r_b0, r_s_z0, r_b_z0)*G_prev[i][j];
            
            R[i][j] =
                R[i][j]
                - 0/(r_s0 - r_b0) * F[i][j] // r_s_theta = 0 ввиду симметрии задачи о конусе
                - (r_s_z0 - r_b_z0)/(r_s0 - r_b0) * G_prev[i][j];
        }
    }
    clock_t finish_initial = clock();
    seconds = double(finish_initial - start_initial) / CLOCKS_PER_SEC;
    // std::cout << "Elapsed time (initial layer): " << seconds << "s" << std::endl;

    // Подъемная сила, поворачивающий момент
    double Fy = 0, Mz = 0;

    // Запись в файл
    z_out << z << " ";

    // #pragma omp parallel for private(xi, r, theta)
    for(int i = 0; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s0, r_b0);
        for(int j = 0; j < M; j++){
            theta = j * dth;

            rho_array[i][j] = G_prev[i][j].get_rho(r);
            p_array[i][j] = G_prev[i][j].get_p(r);
            u_array[i][j] = G_prev[i][j].get_u();
            v_array[i][j] = G_prev[i][j].get_v();
            w_array[i][j] = G_prev[i][j].get_w();
        }
    }
    // Запись в файлы
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
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
    Fy_out << Fy << " ";
    Mz_out << Mz << " ";

    for(int j = 0; j < M; j++){
        r_s_out << r_s0 << " ";
        r_s_theta_out << 0 << " ";
        r_s_z_out << r_s_z0 << " ";
    }
    r_s_out << "\n";
    r_s_theta_out << "\n";
    r_s_z_out << "\n";

    // Интегрирование уравнения d(psi)/dr = rho*w*r
    psi0[0] = 0;
    for(int i = 1; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s[0].back(), r_b(z));
        double dr = r_from_xi(xi, r_s[0].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[0].back(), r_b(z));
        psi0[i] = psi0[i - 1] + rho_array[i][0] * w_array[i][0] * r * dr;
    }
    psi1[0] = 0;
    for(int i = 1; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
        double dr = r_from_xi(xi, r_s[M - 1].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[M - 1].back(), r_b(z));
        psi1[i] = psi1[i - 1] + rho_array[i][M - 1] * w_array[i][M - 1] * r * dr;
    }

    for(int i = 0; i < N; i++){
        psi0_out << psi0[i] << " ";
        psi1_out << psi1[i] << " ";
    }
    psi0_out << "\n";
    psi1_out << "\n";

    double dz, CFL = 0.9, lambda_xi, lambda_th;
    double MM, mm, xi_z_val, xi_r_val, xi_theta_val;
    std::vector<double> z_array;

    // Главный цикл
    clock_t start = clock();
    while(z < L){
        z_array.push_back(z);
        if(progress_bar){
            for(int s = 0; s < num_step_percent; s++){
                if((z - z0) >= (s + 1)*(L - z0)/double(num_step_percent) && !progress_flag[s]){
                    std::cout << (s + 1)*100/double(num_step_percent) << "% completed" << std::endl;
                    progress_flag[s] = true;
                }
            }
        }
        dz = 0.01; // Начальное приближение шага интегрирования по z

        // Вычисление шага dz. Спектральный метод устойчивости
        // #pragma omp parallel for reduction(min: dz) private(xi, r, xi_r_val, xi_theta_val, xi_z_val, MM, mm, lambda_xi, lambda_th)
        for(int i = 0; i < N; i++){
            xi = i * dxi;
            for(int j = 0; j < M; j++){
                r = r_from_xi(xi, r_s[j].back(), r_b(z));
                xi_r_val = xi_r(r_s[j].back(), r_b(z));
                xi_theta_val = xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back());
                xi_z_val = xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z));
                
                double a = sqrt(Gamma * p_array[i][j] / rho_array[i][j]);
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

                // Обновление шага dz
                dz = std::min(dz, 1 / (lambda_xi/dxi + lambda_th/dth));
            }
        }
        // Запас устойчивости (коэффициент CFL)
        dz *= CFL;

        // r_s_pred, r_s_theta_pred - для сохранения предыдущего значения (фазы предиктор-корректор)
        // r_s_z_prev - копия r_s_z, защищенная от изменений по ходу алгоритма
        std::vector<double> r_s_pred(M), r_s_theta_pred(M), r_s_z_prev(M);
        int idx;
        for(int step = 0; step < 2; step++){ // step = 0: predictor, step = 1: corrector
            // #pragma omp parallel for private(xi, theta, idx)
            for(int i = 0; i < N; i++){
                xi = i * dxi;
                // Граница theta = 0
                theta = 0;
                if(step == 0) { // predictor
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
                else { //corrector
                    F_array F_mirrored = (-1) * F[i][1];
                    F_mirrored[2] = (-1) * F_mirrored[2];
                    idx = int(i == 0);
                    G_next[i][0] = corrector(
                        E[i - 1 + idx][0],
                        E[i + idx][0],
                        F_mirrored, // симметрия (у вектора F 1, 2, 4, 5 компоненты нечетны по theta, 3 - четна)
                        F[i][0],
                        G_prev[i][0],
                        G_next[i][0],
                        R[i][0],
                        dxi, dth, dz
                    );
                }

                // Внутренние (по theta) узлы
                for(int j = 1; j < M - 1; j++){
                    theta = j*dth;
                    if(step == 0){ // predictor
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
                    else{ //corrector
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
                }

                // Граница theta = Pi
                theta = Pi;
                if(step == 0){ // predictor
                    F_array F_mirrored = (-1) * F[i][M - 2];
                    F_mirrored[2] = (-1) * F_mirrored[2];
                    idx = int(i == N - 1);
                    G_next[i][M - 1] = predictor(
                        E[i - idx][M - 1],
                        E[i + 1 - idx][M - 1],
                        F[i][M - 1],
                        F_mirrored, // симметрия
                        G_prev[i][M - 1],
                        R[i][M - 1],
                        dxi, dth, dz
                    );
                }
                else{ // corrector
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
            }

            // r_s, r_s_theta, r_s_z calculation (pred-corr)
            for(int j = 0; j < M; j++)
                r_s_z_prev[j] = r_s_z[j].back(); // Создание защищенной копии r_s_z
            for(int j = 0; j < M; j++){
                if(step == 0){ // predictor
                    r_s_pred[j] = r_s[j].back();
                    r_s[j].push_back(r_s[j].back() + r_s_z_prev[j]*dz);
                    r_s_theta_pred[j] = r_s_theta[j].back();
                    if(j == 0){
                        r_s_theta[j].push_back(
                            r_s_theta[j].back() + dz/dth*(r_s_z_prev[j] - r_s_z_prev[j + 1])
                        );
                    }
                    else{
                        r_s_theta[j].push_back(
                            r_s_theta[j].back() + dz/dth*(r_s_z_prev[j] - r_s_z_prev[j - 1])
                        );
                    }
                    // r_s_theta[0].back() = 0; r_s_theta[M - 1].back() = 0; // ВЫБРАТЬ: либо зануление, либо симметрия
                    
                }
                else{ // corrector
                    r_s[j].back() = 0.5*(r_s_pred[j] + r_s[j].back()) + 0.5*dz*r_s_z_prev[j];
                    if(j < M - 1){
                        r_s_theta[j].back() = (
                            0.5 * (r_s_theta_pred[j] + r_s_theta[j].back())
                            + 0.5 * dz/dth * (r_s_z_prev[j + 1] - r_s_z_prev[j])
                        );
                    }
                    else{
                        r_s_theta[j].back() = (
                            0.5 * (r_s_theta_pred[j] + r_s_theta[j].back())
                            + 0.5 * dz/dth * (r_s_z_prev[j - 1] - r_s_z_prev[j]) // симметрия
                        );
                    }
                    // r_s_theta[0].back() = 0; r_s_theta[M - 1].back() = 0; // ВЫБРАТЬ: либо зануление, либо симметрия
                }
            }
            // Добавление шага dz только после предиктора,
            // так как на корректоре не происходит фактического шага
            if(step == 0)
                z += dz;

            // Восстановление векторов E, F, R и физических величин
            // #pragma omp parallel for private(xi, theta, r)
            for(int i = 0; i < N; i++){
                xi = i * dxi;
                // Граница theta = 0
                theta = 0;
                r = r_from_xi(xi, r_s[0].back(), r_b(z));
                E[i][0] = get_E(G_next[i][0], r);
                F[i][0] = get_F(G_next[i][0], r);
                R[i][0] = get_R(G_next[i][0], r, q(r, theta, z, heatSource, flowParams.is_adiabatic));
    
    
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
                    E[i][j] = get_E(G_next[i][j], r);
                    F[i][j] = get_F(G_next[i][j], r);
                    R[i][j] = get_R(G_next[i][j], r, q(r, theta, z, heatSource, flowParams.is_adiabatic));
    
                    // Восстановление физических величин по вектору G
                    rho_array[i][j] = G_next[i][j].get_rho(r);
                    p_array[i][j] = G_next[i][j].get_p(r);
                    u_array[i][j] = G_next[i][j].get_u();
                    v_array[i][j] = G_next[i][j].get_v();
                    w_array[i][j] = G_next[i][j].get_w();
                }
                
                // Граница theta = Pi
                theta = Pi;
                r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
                E[i][M - 1] = get_E(G_next[i][M - 1], r);
                F[i][M - 1] = get_F(G_next[i][M - 1], r);
                R[i][M - 1] = get_R(G_next[i][M - 1], r, q(r, theta, z, heatSource, flowParams.is_adiabatic));

                // Восстановление физических величин по вектору G
                rho_array[i][M - 1] = G_next[i][M - 1].get_rho(r);
                p_array[i][M - 1] = G_next[i][M - 1].get_p(r);
                u_array[i][M - 1] = G_next[i][M - 1].get_u();
                v_array[i][M - 1] = G_next[i][M - 1].get_v();
                w_array[i][M - 1] = G_next[i][M - 1].get_w();
            }
            
            // i = 0 and i = N - 1; поправки на границах

            // Метод Аббета (поправка на поверхности тела)
            for(int j = 0; j < M; j++){
                double delta_th;
                theta = j * dth;
                // (V*, n) - скалярное произведение
                double V_n = (u_array[0][j] - r_b_z(z)*w_array[0][j])
                    / sqrt(1 + r_b_z(z)*r_b_z(z));
                double Mach = sqrt(
                    (
                        u_array[0][j]*u_array[0][j]
                        + v_array[0][j]*v_array[0][j]
                        + w_array[0][j]*w_array[0][j]
                    )
                    /
                    (
                        Gamma * p_array[0][j] / rho_array[0][j]
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
                    1 - Gamma*Mach*Mach*delta_th/sqrt(Mach*Mach - 1)
                    + Gamma
                        * Mach
                        * ((Gamma + 1)*Mach*Mach*Mach*Mach - 4*(Mach*Mach - 1))
                        / (4*(Mach*Mach - 1)*(Mach*Mach - 1))
                        * delta_th * delta_th
                );

                //Поправка плотности
                rho_array[0][j] = flowParams.rho_inf * pow(p_array[0][j] / flowParams.p_inf, 1/Gamma);

                // Полная энтальпия и модуль скорости из интеграла Бернулли
                double H = Gamma/(Gamma - 1) * flowParams.p_inf/flowParams.rho_inf + 0.5*flowParams.V_inf*flowParams.V_inf;
                double V_abs = sqrt(2*(H - Gamma/(Gamma - 1) * p_array[0][j]/rho_array[0][j]));
    
                double Vr_tau, Vth_tau, Vz_tau;
                // Касательная компонента скорости V_tau:
                Vr_tau = u_array[0][j] - V_n / sqrt(1 + r_b_z(z)*r_b_z(z));
                Vth_tau = v_array[0][j];
                Vz_tau = w_array[0][j] + V_n * r_b_z(z) / sqrt(1 + r_b_z(z)*r_b_z(z));
    
                double V_tau_abs = sqrt(Vr_tau*Vr_tau + Vth_tau*Vth_tau + Vz_tau*Vz_tau);
                //Поправка скорости
                u_array[0][j] = Vr_tau * V_abs / V_tau_abs;
                v_array[0][j] = Vth_tau * V_abs / V_tau_abs;
                w_array[0][j] = Vz_tau * V_abs / V_tau_abs;

                // Подъемная сила, поворачивающий момент
                Fy += 2 * (-p_array[0][j] * cos(theta) * r_b(z) * dth * dz); // *2 - четность по углу
                 // TODO: describe mass center more precisely
                Mz += 2 * (-(z - 5.0/4.0 * bodyParams.bodyLength * 0.5) * p_array[0][j] * cos(theta) * r_b(z) * dth * dz);
            }

            // Метод Томаса (поправка на поверхности ударной волны)

            for(int j = 0; j < M; j++){
                // Давление не меняем
                // Плотность из условия Ранкина-Гюгонио
                rho_array[N - 1][j] = (
                    flowParams.rho_inf
                    * ((Gamma + 1)*p_array[N - 1][j] + (Gamma - 1)*flowParams.p_inf)
                    / ((Gamma + 1)*flowParams.p_inf + (Gamma - 1)*p_array[N - 1][j])
                );
                // Корень с отрицательным знаком, так как нормаль берется внешняя
                double V_inf_n = -sqrt(
                    ((Gamma + 1)*p_array[N - 1][j] + (Gamma - 1)*flowParams.p_inf)
                    / (2 * flowParams.rho_inf)
                );
                double V_n = flowParams.rho_inf / rho_array[N - 1][j] * V_inf_n;
                if(step == 0)
                    r_s_z[j].push_back(
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (flowParams.V_inf*flowParams.V_inf - V_inf_n*V_inf_n)
                        )
                    );
                else
                    r_s_z[j].back() = (
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (flowParams.V_inf*flowParams.V_inf - V_inf_n*V_inf_n)
                        )
                    );


                // Нормаль к поверхности УВ
                double nx, ny, nz, n_norm; // n = (nx, ny, nz)
                n_norm = sqrt(
                        1
                        + r_s_z_prev[j]*r_s_z_prev[j]
                        + (r_s_theta[j].back()/r_s[j].back()*r_s_theta[j].back()/r_s[j].back())
                    );
                nx = 1 / n_norm;
                ny = -r_s_theta[j].back() / r_s[j].back() / n_norm;
                nz = -r_s_z_prev[j] / n_norm;
                
                // Касательный вектор к поверхности УВ
                double V_inf_tau_x, V_inf_tau_y, V_inf_tau_z;
                V_inf_tau_x = 0 - V_inf_n * nx;
                V_inf_tau_y = 0 - V_inf_n * ny;
                V_inf_tau_z = flowParams.V_inf - V_inf_n * nz;

                u_array[N - 1][j] = V_inf_tau_x + V_n*nx;
                v_array[N - 1][j] = V_inf_tau_y + V_n*ny;
                w_array[N - 1][j] = V_inf_tau_z + V_n*nz;
            }

            // Функция тока в сечении theta = 0 (psi0) и theta = 1 (psi1)
            if(step == 1){
                // Интегрирование уравнения d(psi)/dr = rho*w*r
                psi0[0] = 0;
                for(int i = 1; i < N; i++){
                    xi = i * dxi;
                    r = r_from_xi(xi, r_s[0].back(), r_b(z));
                    double dr = r_from_xi(xi, r_s[0].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[0].back(), r_b(z));
                    psi0[i] = psi0[i - 1] + rho_array[i][0] * w_array[i][0] * r * dr;
                }
                psi1[0] = 0;
                for(int i = 1; i < N; i++){
                    xi = i * dxi;
                    r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
                    double dr = r_from_xi(xi, r_s[M - 1].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[M - 1].back(), r_b(z));
                    psi1[i] = psi1[i - 1] + rho_array[i][M - 1] * w_array[i][M - 1] * r * dr;
                }

                if(z_limit){
                    if(z > z0 + double(k)/double(files_count)){
                        for(int i = 0; i < N; i++){
                            psi0_out << psi0[i] << " ";
                            psi1_out << psi1[i] << " ";
                        }
                        psi0_out << "\n";
                        psi1_out << "\n";
                    }
                }
                else{
                    for(int i = 0; i < N; i++){
                        psi0_out << psi0[i] << " ";
                        psi1_out << psi1[i] << " ";
                    }
                    psi0_out << "\n";
                    psi1_out << "\n";
                }
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
                            Gamma * p_array[i][j] / (Gamma - 1) 
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
                    R[i][j] = get_R(G_next[i][j], r, q(r, theta, z, heatSource, flowParams.is_adiabatic));
                }
            }
            // Нормализация
            // #pragma omp parallel for private(xi, theta, r)
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
        // Q integration
        for(int i = 1; i < N; i++){
            xi = i * dxi;
            for(int j = 0; j < M; j++){
                r = r_from_xi(xi, r_s[j].back(), r_b(z));
                double dr = r_from_xi(xi, r_s[j].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[j].back(), r_b(z));
                theta = j * dth;
                Q += q(r, theta, z, heatSource, flowParams.is_adiabatic) * dr * r * dth * dz;
            }
        }

        // Запись в файл
        // rho_out << z << "\n";
        // p_out << z << "\n";
        // u_out << z << "\n";
        // v_out << z << "\n";
        // w_out << z << "\n";
        if(z_limit){
            if(z > z0 + double(k)/double(files_count)){
                z_out << z << " ";
                // std::cout << z << std::endl;
                for(int i = 0; i < N; i++){
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
                Fy_out << Fy << " ";
                Mz_out << Mz << " ";   
            }
        }
        else{
            z_out << z << " ";
            std::cout << z << std::endl;
            for(int i = 0; i < N; i++){
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
            Fy_out << Fy << " ";
            Mz_out << Mz << " ";
        }
        
        if(z_limit){
            if(z > z0 + double(k)/double(files_count)){
                for(int j = 0; j < M; j++){
                    r_s_out << r_s[j].back() << " ";
                    r_s_theta_out << r_s_theta[j].back() << " ";
                    r_s_z_out << r_s_z[j].back() << " ";
                }
                r_s_out << "\n";
                r_s_theta_out << "\n";
                r_s_z_out << "\n";
                k++;
            }
        }
        else{
            for(int j = 0; j < M; j++){
                r_s_out << r_s[j].back() << " ";
                r_s_theta_out << r_s_theta[j].back() << " ";
                r_s_z_out << r_s_z[j].back() << " ";
            }
            r_s_out << "\n";
            r_s_theta_out << "\n";
            r_s_z_out << "\n";
        }
    }
    if(progress_bar)
        std::cout << "100% completed" << std::endl;
    
    std::cout << "==================================\nLifting force: " << Fy << std::endl;
    std::cout << "Rotation momentum: " << Mz << std::endl;
    std::cout << "Q = " << Q << std::endl;

    clock_t finish = clock();
    seconds = double(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time (main loop): " << seconds << "s" << std::endl;

    std::vector<double> res = {Fy, Mz, Q};
    return res;
}
