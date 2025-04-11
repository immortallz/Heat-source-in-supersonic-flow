#include "r.h"

void RK(double t0, double T, double h, vector<double> y0, string filename)
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
		y = y0;

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
		fprintf(file, "%lf", t0 + i*h);
		for(double elem : y)
			fprintf(file, " %lf", elem);
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
			k[j][0] = y[1] + h*nn[1];
			k[j][1] = -y[0] - h*nn[0];
			k[j][2] = -k[j][0]*k[j][1];
		}
		for(int j = 1; j <= s; j++)
			for(int r = 0; r < dim; r++)
				y[r] += h*b[j]*k[j][r];
	}

	fprintf(file, "%lf", t0 + N*h);
	for(double elem : y)
		fprintf(file, " %lf", elem);
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
			k[j][0] = y[1] + (T - t0 - N*h)*nn[1];
			k[j][1] = -y[0] - (T - t0 - N*h)*nn[0];
			k[j][2] = -k[j][0]*k[j][1];
		}
		for(int j = 1; j <= s; j++)
			for(int r = 0; r < dim; r++)
				y[r] += (T - t0 - N*h)*b[j]*k[j][r];

		fprintf(file, "%lf", T);
		for(double elem : y)
			fprintf(file, " %lf", elem);
		fprintf(file, "\n");
	}

	fclose(file);
}