#include "RK.h"

int main()
{
	vector<double> y0(3, 0);
	y0[0] = 0;
	y0[1] = 1;
	y0[2] = 0;

	RK(0, 2*PI, 0.001, y0, "output.txt");
	return 0;
}