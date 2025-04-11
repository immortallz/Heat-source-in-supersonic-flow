#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#define PI 3.141592653589793

using namespace std;

void RK(double t0, double T, double h, vector<double> y0, string filename);