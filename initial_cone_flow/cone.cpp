#include "r.h"

vector<double> addVectors(const vector<double>& a, const vector<double>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Vectors must have the same size");
    }

    vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Vectors must have the same size");
    }

    vector<T> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template <typename T>
vector<T> operator*(const vector<T>& vec, T scalar) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

template <typename T>
vector<T> operator*(T scalar, const vector<T>& vec) {
    return vec * scalar; // Просто вызываем первую перегрузку
}

double a_sqaured(vector<double> v, double C)
{
	return (gamma - 1) * (C - 0.5*(v[0]*v[0] + v[1]*v[1]));
}

vector<double> func(double theta, vector<double> v, double C)
{
	// (v[0], v[1]) = (V_R, V_\theta)

	vector<double> result(v.size());
	double a2 = a_sqaured(v, C);
	result[0] = v[1];
	result[1] = -v[0] + a2*(v[0] + v[1]/tan(theta)) / (v[1]*v[1] - a2);
	return result;
}

vector<double> RK(double t0, double T, double h, vector<double> y0, string filename, double C)
{
	int
		s = 7, //num of stages
		dim = y0.size();
	vector<vector<double>>
		a(s + 1, vector<double>(s)),
		k(s + 1, vector<double>(dim));
	vector<double>
		b(s + 1),
		c(s + 1),
		nn(dim),
		y = y0,
		f;


	a[2][1]=0.2; a[3][1]=0.075; a[4][1]=44./45; a[5][1]=19372./6561; a[6][1]=9017./3168; a[7][1]=35./384;
	a[3][2]=0.225; a[4][2]=-56./15; a[5][2]=-25360./2187; a[6][2]=-355./33; a[7][2]=0;
	a[4][3]=32./9; a[5][3]=64448./6561; a[6][3]=46732./5247; a[7][3]=500./1113;
	a[5][4]=-212./729; a[6][4]=49./176; a[7][4]=125./192;
	a[6][5]=-5103./18656; a[7][5]=-2187./6784;
	a[7][6]=11./84;
	c[1]=0; c[2]=0.2; c[3]=0.3; c[4]=0.8; c[5]=8./9; c[6]=1; c[7]=1;
	b[1]=35./384; b[2]=0; b[3]=500./1113; b[4]=125./192; b[5]=-2187./6784; b[6]=11./84; b[7]=0;

	FILE *file = fopen(filename.c_str(),"w");

	int N = floor((T-t0)/h);

	for(int i = 0; i < N; i++)
	{
		fprintf(file, "%.14lf", t0 + i*h);
		for(double elem : y)
			fprintf(file, " %.14lf", elem);
		fprintf(file, "\n");

		for(int j = 1; j <= s; j++)
		{
			fill(nn.begin(), nn.end(), 0);
			for(int l = 1; l < j; l++)
			{
				for(int r = 0; r < dim; r++)
					nn[r] += a[j][l]*k[l][r];
			}

			//!!!!!
			//function is here
			f = func(t0 + i*h, y + (h*nn), C);
			for(int r = 0; r < f.size(); r++)
				k[j][r] = f[r];
		}
		for(int j = 1; j <= s; j++)
			for(int r = 0; r < dim; r++)
				y[r] += h*b[j]*k[j][r];
	}

	fprintf(file, "%.14lf", t0 + N*h);
	for(double elem : y)
		fprintf(file, " %.14lf", elem);
	fprintf(file, "\n");

	if(abs(T - t0 - N*h) > 1e-15){
		for(int j = 1; j <= s; j++)
		{
			fill(nn.begin(), nn.end(), 0);
			for(int l = 1; l < j; l++)
				for(int r = 0; r < dim; r++)
					nn[r] += a[j][l]*k[l][r];

			//!!!!!
			//function is here
			f = func(T, y + ((T - t0 - N*h)*nn), C);
			for(int r = 0; r < f.size(); r++)
				k[j][r] = f[r];
		}
		for(int j = 1; j <= s; j++)
			for(int r = 0; r < dim; r++)
				y[r] += (T - t0 - N*h)*b[j]*k[j][r];

		fprintf(file, "%.14lf", T);
		for(double elem : y)
			fprintf(file, " %.14lf", elem);
		fprintf(file, "\n");
	}
	fclose(file);

	return y; 
}

double secant_method(double p0, double rho0, double V1, double theta0,
	double beta0, double beta1,
	double h, double tol
	)
{
	vector<double> v(2), //v[0] = V_R, v[1] = V_\theta
		sol(2);
	double 
		C = 			//Bernoulli constant
			gamma/(gamma - 1) * p0 / rho0
			+ 0.5*(V1*V1),
		a0 = sqrt(gamma * p0 / rho0),
		M = V1 / a0,
		rho1_rho2;

	FILE *tmp;
	string tmp_filename = "tmp.txt";
	double tmp_beta;
	double err0, err1 = 1;
	while(abs(err1) > tol)
	{
		rho1_rho2 = (gamma - 1)/(gamma + 1)
			+ 2 / ((gamma + 1) * M*M * sin(beta0)*sin(beta0));
		v[0] = V1 * cos(beta0);
		v[1] = -rho1_rho2 * V1 * sin(beta0);
		sol = RK(beta0, theta0, -h, v, tmp_filename, C);
		err0 = sol[1];

		rho1_rho2 = (gamma - 1)/(gamma + 1)
			+ 2 / ((gamma + 1) * M*M * sin(beta1)*sin(beta1));
		v[0] = V1 * cos(beta1);
		v[1] = -rho1_rho2 * V1 * sin(beta1);
		sol = RK(beta1, theta0, -h, v, tmp_filename, C);
		err1 = sol[1];

		tmp_beta = beta1 - err1 * (beta1 - beta0) / (err1 - err0);
		beta0 = beta1;
		beta1 = tmp_beta;
	}
	return beta1;
}

double newton_method(double p0, double rho0, double V1, double theta0,
	double beta0,
	double h, double tol
	)
{
	vector<double> v(2), //v[0] = V_R, v[1] = V_\theta
		sol(2);
	double 
		C = 			//Bernoulli constant
			gamma/(gamma - 1) * p0 / rho0
			+ 0.5*(V1*V1),
		a0 = sqrt(gamma * p0 / rho0),
		M = V1 / a0,
		rho1_rho2;

	FILE *tmp;
	string tmp_filename = "tmp.txt";
	double tmp_beta;
	double err0 = 1, err1;
	while(abs(err0) > tol)
	{
		rho1_rho2 = (gamma - 1)/(gamma + 1)
			+ 2 / ((gamma + 1) * M*M * sin(beta0)*sin(beta0));
		v[0] = V1 * cos(beta0);
		v[1] = -rho1_rho2 * V1 * sin(beta0);
		sol = RK(beta0, theta0, -h, v, tmp_filename, C);
		err0 = sol[1];
		sol = RK(beta0 + h, theta0, -h, v, tmp_filename, C);
		err1 = sol[1];

		beta0 = beta0 - err0 * h / (err1 - err0);
	}
	return beta0;
}
