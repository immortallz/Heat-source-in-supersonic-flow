#include "r.h"
#include <iostream>
#include <cmath>
#include <fstream>

int main()
{
	std::vector<double> v(2); // v[0] = V_R, v[1] = V_\theta
	constexpr double p0 = 101330, rho0 = 1.2255;
	double Mach0;
	const double a0 = sqrt(GAMMA * p0 / rho0);
	int theta0_int;
	std::cout << "enter theta0 (in degrees)" << std::endl;
	std::cin >> theta0_int;
	std::cout << "enter initial Mach" << std::endl;
	std::cin >> Mach0;
	const double theta0 = theta0_int * PI / 180.0;
	const double V1 = Mach0 * a0;
	const double C = 			//Bernoulli constant
		GAMMA/(GAMMA - 1) * p0 / rho0
		+ 0.5*(V1*V1),
		M = V1 / a0;

	const double beta = secant_method(p0, rho0, V1, theta0,
	                            PI / 6, PI / 5,
	                            0.00001, 1e-10
	);

	std::cout << "Calculated shockwave angle: " << beta << " rad = " << beta / PI * 180.0 << " degrees" << std::endl;

	const double rho1_rho2 = (GAMMA - 1)/(GAMMA + 1)
			+ 2 / ((GAMMA + 1) * M*M * sin(beta)*sin(beta));
	v[0] = V1 * cos(beta);
	v[1] = -rho1_rho2 * V1 * sin(beta);
	std::vector<double> sol(2);
	sol = RK(beta, theta0, -0.00001, v, "output_results/VR_Vtheta_init.txt", C);
	std::cout << "V_R, V_theta on body surface: " << sol[0] << ", " << sol[1] << std::endl;

	const double rho_s = rho0 / rho1_rho2;
	const double p_s = p0 * (2 * GAMMA * Mach0 * Mach0 * sin(beta) * sin(beta) / (GAMMA + 1) - (GAMMA - 1) / (GAMMA + 1));

	const double C_entr = p_s / pow(rho_s, GAMMA);

	FILE *f_rho = fopen("output_results/rho_init.txt", "w"), *f_p = fopen("output_results/p_init.txt", "w");

	std::ifstream input("output_results/VR_Vtheta_init.txt");
    if (!input) {
        std::cerr << "Не удалось открыть файл!" << std::endl;
        return 1;
    }

    // Три вектора для хранения чисел из каждой колонки
    double theta, VR, Vtheta;

    // Читаем по 3 числа в каждой строке
    double rho, p;
    while (input >> theta >> VR >> Vtheta) {
    	rho = pow((GAMMA-1)/GAMMA/C_entr*(C - 0.5*(VR*VR + Vtheta*Vtheta)), 1/(GAMMA-1));
    	p = C_entr * pow(rho, GAMMA);
    	fprintf(f_rho, "%.14lf %.14lf\n", theta, rho);
    	fprintf(f_p, "%.14lf %.14lf\n", theta, p);
    }
	std::cout << "Density and pressure on body surface: rho = " << rho << ", p = " << p << std::endl;

    fclose(f_rho);
    fclose(f_p);

	// RK_rho(rho0 / rho1_rho2, "VR_Vtheta.txt", "rho.txt", C);

	// FILE *beta_file = fopen(("beta_" + to_string(theta0_int) + ".txt").c_str(), "w");

	// for(double Mach = Mach0; Mach < 5.001; Mach += 0.01)
	// {
	// 	V1 = Mach*a0;
	// 	C = 			//Bernoulli constant
	// 		GAMMA/(GAMMA - 1) * p0 / rho0
	// 		+ 0.5*(V1*V1),
	// 	M = V1 / a0;
	// 	beta = newton_method(p0, rho0, V1, theta0,
	// 		PI/5,
	// 		0.001, 1e-7
	// 		);
	// 	std::cout << Mach << std::endl;
	// 	fprintf(beta_file, "%lf %lf\n", Mach, beta);
	// }
	// fclose(beta_file);
	return 0;
}