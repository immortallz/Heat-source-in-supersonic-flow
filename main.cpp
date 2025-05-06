#include "r.h"
#include <string>
#include <ctime>

int solver() {
    ofstream z_out("z_out.txt");
    ofstream rho_out("rho_out.txt");
    ofstream p_out("p_out.txt");
    ofstream u_out("u_out.txt");
    ofstream v_out("v_out.txt");
    ofstream w_out("w_out.txt");
    ofstream r_s_out("r_s_out.txt");
    ofstream r_s_theta_out("r_s_theta_out.txt");
    ofstream r_s_z_out("r_s_z_out.txt");
    ofstream psi0_out("psi0_out.txt");
    ofstream psi1_out("psi1_out.txt");

    int
        N = 100, // xi
        M = 300; // theta

    bool debug = false, progress_bar = true;
    int num_step_percent = 100;
    vector<bool> progress_flag(num_step_percent, false);

    double
        z0 = 1.0, z = z0, L = z0 + 1,
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
    vector<double> psi0(N, 0), psi1(N, 0);

    vector<double> phi_cone, rho_cone, p_cone, VR_cone, Vphi_cone; // Векторы для начальных данных

    // Чтение (заранее заготовленных) начальных данных
    // из задачи об обтекании конуса
    FILE *f_cone = fopen("rho_init.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        phi_cone.push_back(a);
        rho_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("p_init.txt", "r");
    while(!feof(f_cone))
    {
        double a, b;
        fscanf(f_cone, "%lf %lf", &a, &b);
        p_cone.push_back(b);
    }
    fclose(f_cone);
    f_cone = fopen("VR_Vtheta_init.txt", "r");
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
    std::cout << phi0 << " " << phi1 << " " << dphi << endl;
    
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
                - 0/(r_s0 - r_b0) * F[i][j] // r_s_theta = 0 ввиду симметрии задачи о конусе
                - (r_s_z0 - r_b_z0)/(r_s0 - r_b0) * G_prev[i][j];
            
            if(debug && j == 0){ // debug output
                std::cout << "xi check initial" << endl;
                std::cout << xi_r(r_s0, r_b0) << " " << 1/z/(tan(phi1) - tan(phi0)) << endl;
                std::cout << xi_theta(xi, r_s0, r_b0, 0) << " " << 0.0 << endl;
                std::cout << xi_z(xi, r_s0, r_b0, r_s_z0, r_b_z0) << " " << -r/z/z/(tan(phi1) - tan(phi0)) << endl;
            }
        }
    }
    clock_t finish_initial = clock();
    seconds = double(finish_initial - start_initial) / CLOCKS_PER_SEC;
    // std::cout << "Elapsed time (initial layer): " << seconds << "s" << std::endl;

    // Подъемная сила, поворачивающий момент
    double Fy = 0, Mz = 0;

    // Запись в файл
    z_out << z << " ";
    // rho_out << z << "\n";
    // p_out << z << "\n";
    // u_out << z << "\n";
    // v_out << z << "\n";
    // w_out << z << "\n";
    #pragma omp parallel for private(xi, r, theta)
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
    // rho_out << "\n";
    // p_out << "\n";
    // u_out << "\n";
    // v_out << "\n";
    // w_out << "\n";

    // r_s_out << z << " ";
    // r_s_theta_out << z << " ";
    // r_s_z_out << z << " ";
    for(int j = 0; j < M; j++){
        r_s_out << r_s0 << " ";
        r_s_theta_out << 0 << " ";
        r_s_z_out << r_s_z0 << " ";
    }
    r_s_out << "\n";
    r_s_theta_out << "\n";
    r_s_z_out << "\n";

    // Интегрирование уравнения d(psi)/dr = rho*u*r
    psi0[0] = 0;
    for(int i = 1; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s[0].back(), r_b(z));
        double dr = r_from_xi(xi, r_s[0].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[0].back(), r_b(z));
        psi0[i] = psi0[i - 1] + rho_array[i][0] * u_array[i][0] * r * dr;
    }
    psi1[0] = 0;
    for(int i = 1; i < N; i++){
        xi = i * dxi;
        r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
        double dr = r_from_xi(xi, r_s[M - 1].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[M - 1].back(), r_b(z));
        psi1[i] = psi1[i - 1] + rho_array[i][M - 1] * u_array[i][M - 1] * r * dr;
    }

    // psi0_out << z << "\n";
    // psi1_out << z << "\n";
    for(int i = 0; i < N; i++){
        psi0_out << psi0[i] << " ";
        psi1_out << psi1[i] << " ";
    }
    psi0_out << "\n";
    psi1_out << "\n";

    double dz, CFL = 0.9, lambda_xi, lambda_th;
    double MM, mm, xi_z_val, xi_r_val, xi_theta_val;
    vector<double> z_array;

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
        #pragma omp parallel for reduction(min: dz) private(xi, r, xi_r_val, xi_theta_val, xi_z_val, MM, mm, lambda_xi, lambda_th)
        for(int i = 0; i < N; i++){
            xi = i * dxi;
            for(int j = 0; j < M; j++){
                r = r_from_xi(xi, r_s[j].back(), r_b(z));
                xi_r_val = xi_r(r_s[j].back(), r_b(z));
                xi_theta_val = xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back());
                xi_z_val = xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z));
                if(debug && j == 0){ // debug output
                    std::cout << "xi check before pred-corr" << endl;
                    std::cout << xi_r_val << " " << 1/z/(tan(phi1) - tan(phi0)) << endl;
                    std::cout << xi_theta_val << " " << 0.0 << endl;
                    std::cout << xi_z_val << " " << -r/z/z/(tan(phi1) - tan(phi0)) << endl;
                }
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

                // Обновление шага dz
                dz = min(dz, 1 / (lambda_xi/dxi + lambda_th/dth));
            }
        }
        // Запас устойчивости (коэффициент CFL)
        dz *= CFL;

        // r_s_pred, r_s_theta_pred - для сохранения предыдущего значения (фазы предиктор-корректор)
        // r_s_z_prev - копия r_s_z, защищенная от изменений по ходу алгоритма
        vector<double> r_s_pred(M), r_s_theta_pred(M), r_s_z_prev(M);
        int idx;
        for(int step = 0; step < 2; step++){ // step = 0: predictor, step = 1: corrector
            #pragma omp parallel for private(xi, theta, idx)
            for(int i = 0; i < N; i++){
                xi = i * dxi;
                // Граница theta = 0
                theta = 0;
                if(step == 0){ // predictor
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
                else{ //corrector
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

                // Граница theta = PI
                theta = PI;
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
                        if(debug) // debug output
                            std::cout << "r_s_theta[0] = " << r_s_theta[j].back() << endl;
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
                        if(debug)  // debug output
                            std::cout << "r_s_theta[M - 1] = " << r_s_theta[j].back() << endl;
                    }
                    // r_s_theta[0].back() = 0; r_s_theta[M - 1].back() = 0; // ВЫБРАТЬ: либо зануление, либо симметрия
                }
            }
            // Добавление шага dz только после предиктора,
            // так как на корректоре не происходит фактического шага
            if(step == 0)
                z += dz;
            if(debug) // debug output
                std::cout << "dz = " << dz << endl;

            // Восстановление векторов E, F, R и физических величин
            #pragma omp parallel for private(xi, theta, r)
            for(int i = 0; i < N; i++){
                xi = i * dxi;
                // Граница theta = 0
                theta = 0;
                r = r_from_xi(xi, r_s[0].back(), r_b(z));
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
                if(debug && j == 0) // debug output
                    std::cout << "\nABBET METHOD:\nV_n = " << V_n << "\nMach = " << Mach << "\ndelta_theta = " << delta_th << endl;

                // Поправка давления
                if(debug && j == 0) // debug output
                    std::cout << "old p = " << p_array[0][j] << " --> new p = " << p_array[0][j] * (1 - gamma*Mach*Mach*delta_th/sqrt(Mach*Mach - 1)+ gamma* Mach* ((gamma + 1)*Mach*Mach*Mach*Mach - 4*(Mach*Mach - 1))/ (4*(Mach*Mach - 1)*(Mach*Mach - 1))* delta_th * delta_th) << endl;
                p_array[0][j] = p_array[0][j] * (
                    1 - gamma*Mach*Mach*delta_th/sqrt(Mach*Mach - 1)
                    + gamma
                        * Mach
                        * ((gamma + 1)*Mach*Mach*Mach*Mach - 4*(Mach*Mach - 1))
                        / (4*(Mach*Mach - 1)*(Mach*Mach - 1))
                        * delta_th * delta_th
                );

                //Поправка плотности
                if(debug && j == 0) // debug output
                    std::cout << "old rho = " << rho_array[0][j] << " --> new rho = " << rho_inf * pow(p_array[0][j] / p_inf, 1/gamma) << endl;
                rho_array[0][j] = rho_inf * pow(p_array[0][j] / p_inf, 1/gamma);

                // Полная энтальпия и модуль скорости из интеграла Бернулли
                double H = gamma/(gamma - 1) * p_inf/rho_inf + 0.5*V_inf*V_inf;
                double V_abs = sqrt(2*(H - gamma/(gamma - 1) * p_array[0][j]/rho_array[0][j]));
    
                double Vr_tau, Vth_tau, Vz_tau;
                // Касательная компонента скорости V_tau:
                Vr_tau = u_array[0][j] - V_n / sqrt(1 + r_b_z(z)*r_b_z(z));
                Vth_tau = v_array[0][j];
                Vz_tau = w_array[0][j] + V_n * r_b_z(z) / sqrt(1 + r_b_z(z)*r_b_z(z));
    
                double V_tau_abs = sqrt(Vr_tau*Vr_tau + Vth_tau*Vth_tau + Vz_tau*Vz_tau);
                //Поправка скорости
                if(debug && j == 0){ // debug output
                    std::cout << "old u = " << u_array[0][j] << " --> new u = " << Vr_tau * V_abs / V_tau_abs << endl;
                    std::cout << "old v = " << v_array[0][j] << " --> new v = " << Vth_tau * V_abs / V_tau_abs << endl;
                    std::cout << "old w = " << w_array[0][j] << " --> new w = " << Vz_tau * V_abs / V_tau_abs << endl;
                }
                u_array[0][j] = Vr_tau * V_abs / V_tau_abs;
                v_array[0][j] = Vth_tau * V_abs / V_tau_abs;
                w_array[0][j] = Vz_tau * V_abs / V_tau_abs;

                // Подъемная сила, поворачивающий момент
                Fy += -p_array[0][j] * cos(theta) * r_b(z) * dth * dz;
                Mz += -z * p_array[0][j] * cos(theta) * r_b(z) * dth * dz;
            }

            // Метод Томаса (поправка на поверхности ударной волны)

            for(int j = 0; j < M; j++){
                // Давление не меняем
                // Плотность из условия Ранкина-Гюгонио
                if(debug && j == 0) // debug output
                    std::cout<<"\nTHOMAS METHOD:\nold rho = "<<rho_array[N - 1][j]<<" --> new rho = "<<(rho_inf* ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)/ ((gamma + 1)*p_inf + (gamma - 1)*p_array[N - 1][j])) << endl;
                rho_array[N - 1][j] = (
                    rho_inf
                    * ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)
                    / ((gamma + 1)*p_inf + (gamma - 1)*p_array[N - 1][j])
                );
                // Корень с отрицательным знаком, так как нормаль берется внешняя
                double V_inf_n = -sqrt(
                    ((gamma + 1)*p_array[N - 1][j] + (gamma - 1)*p_inf)
                    / (2 * rho_inf)
                );
                double V_n = rho_inf / rho_array[N - 1][j] * V_inf_n;
                if(step == 0)
                    r_s_z[j].push_back(
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (V_inf*V_inf - V_inf_n*V_inf_n)
                        )
                    );
                else
                    r_s_z[j].back() = (
                        sqrt(
                            V_inf_n * V_inf_n
                            * (1 + r_s_theta[j].back()*r_s_theta[j].back()/r_s[j].back()/r_s[j].back())
                            / (V_inf*V_inf - V_inf_n*V_inf_n)
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
                V_inf_tau_z = V_inf - V_inf_n * nz;

                if(debug && j == 0){ // debug output
                    std::cout << "old u = " << u_array[N - 1][j] << " --> new u = " << V_inf_tau_x + V_n*nx << endl;
                    std::cout << "old v = " << v_array[N - 1][j] << " --> new v = " << V_inf_tau_y + V_n*ny << endl;
                    std::cout << "old w = " << w_array[N - 1][j] << " --> new w = " << V_inf_tau_z + V_n*nz << endl;
                }
                u_array[N - 1][j] = V_inf_tau_x + V_n*nx;
                v_array[N - 1][j] = V_inf_tau_y + V_n*ny;
                w_array[N - 1][j] = V_inf_tau_z + V_n*nz;
            }

            // Функция тока в сечении theta = 0 (psi0) и theta = 1 (psi1)
            if(step == 1){
                // Интегрирование уравнения d(psi)/dr = rho*u*r
                psi0[0] = 0;
                for(int i = 1; i < N; i++){
                    xi = i * dxi;
                    r = r_from_xi(xi, r_s[0].back(), r_b(z));
                    double dr = r_from_xi(xi, r_s[0].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[0].back(), r_b(z));
                    psi0[i] = psi0[i - 1] + rho_array[i][0] * u_array[i][0] * r * dr;
                }
                psi1[0] = 0;
                for(int i = 1; i < N; i++){
                    xi = i * dxi;
                    r = r_from_xi(xi, r_s[M - 1].back(), r_b(z));
                    double dr = r_from_xi(xi, r_s[M - 1].back(), r_b(z)) - r_from_xi(xi - dxi, r_s[M - 1].back(), r_b(z));
                    psi1[i] = psi1[i - 1] + rho_array[i][M - 1] * u_array[i][M - 1] * r * dr;
                }

                // psi0_out << z << "\n";
                // psi1_out << z << "\n";
                for(int i = 0; i < N; i++){
                    psi0_out << psi0[i] << " ";
                    psi1_out << psi1[i] << " ";
                }
                psi0_out << "\n";
                psi1_out << "\n";
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
            #pragma omp parallel for private(xi, theta, r)
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
                    if(debug && j == 0){
                        std::cout << "xi check after A-T" << endl;
                        std::cout << xi_r(r_s[j].back(), r_b(z)) << " " << 1/z/(tan(phi1) - tan(phi0)) << endl;
                        std::cout << xi_theta(xi, r_s[j].back(), r_b(z), r_s_theta[j].back()) << " " << 0.0 << endl;
                        std::cout << xi_z(xi, r_s[j].back(), r_b(z), r_s_z[j].back(), r_b_z(z)) << " " << -r/z/z/(tan(phi1) - tan(phi0)) << endl;
                    }

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
        z_out << z << " ";
        // rho_out << z << "\n";
        // p_out << z << "\n";
        // u_out << z << "\n";
        // v_out << z << "\n";
        // w_out << z << "\n";
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
        // rho_out << "\n";
        // p_out << "\n";
        // u_out << "\n";
        // v_out << "\n";
        // w_out << "\n";

        // r_s_out << z << " ";
        // r_s_theta_out << z << " ";
        // r_s_z_out << z << " ";
        for(int j = 0; j < M; j++){
            r_s_out << r_s[j].back() << " ";
            r_s_theta_out << r_s_theta[j].back() << " ";
            r_s_z_out << r_s_z[j].back() << " ";
        }
        r_s_out << "\n";
        r_s_theta_out << "\n";
        r_s_z_out << "\n";
    }
    if(progress_bar)
        std::cout << "100% completed" << std::endl;
    
    std::cout << "\n=============================\nLifting force: " << Fy << std::endl;
    std::cout << "Rotation momentum: " << Mz << std::endl;

    clock_t finish = clock();
    seconds = double(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time (main loop): " << seconds << "s" << std::endl;
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
    R_array R = {1,1,2,2,2};
    double rho = 1.2, p = 100000, u = 100, v = 200, w = 300, r = 3;
    rho = 1.615063, p = 149129.877770, u = 93.149968, v = 1.000000, w = 621.925338;

    E.data[0] = rho*u;
    E.data[1] = rho*u*u + p;
    E.data[2] = rho*u*v;
    E.data[3] = rho*u*w;
    E.data[4] = u * (gamma * p / (gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
    E = r * E;

    F.data[0] = rho*v;
    F.data[1] = rho*u*v;
    F.data[2] = rho*v*v + p;
    F.data[3] = rho*v*w;
    F.data[4] = v * (gamma * p / (gamma - 1) + rho*(u*u + v*v + w*w)*0.5);

    G.data[0] = rho*w;
    G.data[1] = rho*u*w;
    G.data[2] = rho*v*w;
    G.data[3] = rho*w*w + p;
    G.data[4] = w * (gamma * p / (gamma - 1) + rho*(u*u + v*v + w*w)*0.5);
    G = r * G;

    R.data[0] = 0;
    R.data[1] = rho*v*v + p;
    R.data[2] = -rho*u*v;
    R.data[3] = 0;
    R.data[4] = 2.28;

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

    std::cout << "\n" << endl;

    E.print();
    get_E(G, r).print();
    F.print();
    get_F(G, r).print();
    R.print();
    get_R(G, r, 2.28).print();

    std::cout << "F_mirrored check:" << endl;
    F.print();
    F_array F_mirrored = (-1) * F;
    F_mirrored[2] = (-1) * F_mirrored[2];
    F_mirrored.print();
}

int main()
{
    clock_t start = clock();
    solver();
    clock_t finish = clock();
    double seconds = double(finish - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time: " << seconds << "s" << std::endl;
    // test();

    return 0;
}