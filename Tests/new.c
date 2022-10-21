//"/Users/carlosalbertosalassalas/opt/anaconda3/include/python3.9/Python.h"
//#include "wctype.h"
#include <iostream>
#include <math.h>
#include <iomanip> 

using namespace std; 

double linear(double x, double a0, double a1){
	return a0 + a1*x;
}

// Potentials

double woodsSaxon(long double r, long double V0, long double a, long double r0){
	return (-1*V0)/(exp((r - r0)/a));
}



int main(int argc, char const *argv[])
{	
	cout << setprecision(16) << '\n';
	cout << woodsSaxon(1, 2, 3, 4) << '\n';
	cout << exp(1) << '\n';
	cout << M_PI << '\n';
	return 0;
}


