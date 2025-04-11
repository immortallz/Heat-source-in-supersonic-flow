#include <iostream>

using namespace std;

int main()
{

	FILE *f = fopen("rho.txt", "r");
	double a, b;
	while(!feof(f))
	{
		fscanf(f, "%lf %lf", &a, &b);
		cout << a << endl;
	}
	return 0;
}