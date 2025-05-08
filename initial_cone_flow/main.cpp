#include "r.h"

int main()
{
	vector<double> v(2); //v[0] = V_R, v[1] = V_\theta
	double p0 = 101330, rho0 = 1.2255, V1, beta, theta0, Mach0;
	double a0 = sqrt(gamma * p0 / rho0);
	int theta0_int;
	cout << "enter theta0 (in degrees)" << endl;
	cin >> theta0_int;
	cout << "enter initial Mach" << endl;
	cin >> Mach0;
	theta0 = double(theta0_int * PI / 180.0);
	V1 = Mach0*a0;
	double C = 			//Bernoulli constant
		gamma/(gamma - 1) * p0 / rho0
		+ 0.5*(V1*V1),
		M = V1 / a0;
	
	beta = secant_method(p0, rho0, V1, theta0,
	PI/6, PI/5,
	0.00001, 1e-10
	);

	cout << "Calculated shockwave angle: " << beta << " rad = " << beta / PI * 180.0 << " degrees" << endl;

	double rho1_rho2 = (gamma - 1)/(gamma + 1)
			+ 2 / ((gamma + 1) * M*M * sin(beta)*sin(beta));
	v[0] = V1 * cos(beta);
	v[1] = -rho1_rho2 * V1 * sin(beta);
	vector<double> sol(2);
	sol = RK(beta, theta0, -0.00001, v, "VR_Vtheta_init.txt", C);
	cout << "V_R, V_theta on body surface: " << sol[0] << ", " << sol[1] << endl;

	double rho_s, p_s;
	rho_s = rho0 / rho1_rho2;
	p_s = p0 * (2*gamma*Mach0*Mach0*sin(beta)*sin(beta)/(gamma + 1) - (gamma-1)/(gamma+1));

	double C_entr = p0 / pow(rho0, gamma);

	FILE *f_rho = fopen("rho_init.txt", "w"), *f_p = fopen("p_init.txt", "w");

	ifstream input("VR_Vtheta_init.txt");
    if (!input) {
        cerr << "Не удалось открыть файл!" << endl;
        return 1;
    }

    // Три вектора для хранения чисел из каждой колонки
    double theta, VR, Vtheta;

    // Читаем по 3 числа в каждой строке
    double rho, p;
    while (input >> theta >> VR >> Vtheta) {
    	rho = pow((gamma-1)/gamma/C_entr*(C - 0.5*(VR*VR + Vtheta*Vtheta)), 1/(gamma-1));
    	p = C_entr * pow(rho, gamma);
    	fprintf(f_rho, "%.14lf %.14lf\n", theta, rho);
    	fprintf(f_p, "%.14lf %.14lf\n", theta, p);
    }
	cout << "Density and pressure on body surface: rho = " << rho << ", p = " << p << endl;

    fclose(f_rho);
    fclose(f_p);

	// RK_rho(rho0 / rho1_rho2, "VR_Vtheta.txt", "rho.txt", C);

	// FILE *beta_file = fopen(("beta_" + to_string(theta0_int) + ".txt").c_str(), "w");

	// for(double Mach = Mach0; Mach < 5.001; Mach += 0.01)
	// {
	// 	V1 = Mach*a0;
	// 	C = 			//Bernoulli constant
	// 		gamma/(gamma - 1) * p0 / rho0
	// 		+ 0.5*(V1*V1),
	// 	M = V1 / a0;
	// 	beta = newton_method(p0, rho0, V1, theta0,
	// 		PI/5,
	// 		0.001, 1e-7
	// 		);
	// 	cout << Mach << endl;
	// 	fprintf(beta_file, "%lf %lf\n", Mach, beta);
	// }
	// fclose(beta_file);
	return 0;
}