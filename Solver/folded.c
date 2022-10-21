#include <math.h>
#include <stdio.h>

double densityM(double r, double rho0M, double omega, double betaM){
	return rho0M*(1+omega*(pow(r, 2)))*exp(-betaM*(pow(r, 2)));
}

double densityAlpha(double r, double rho0Alpha, double betaAlpha){
	return rho0Alpha*exp(-betaAlpha*(pow(r, 2)));
}

double * sphericalVector(double r, double theta, double phi){
	
	double *vector; 

	vector[0] = r*sin(theta)*cos(phi);
	vector[1] = r*sin(theta)*sin(phi);
	vector[2] = r*cos(theta);

	return vector;
}						

double r12(double r, double theta, double phi, double r1, double theta1, double phi1, double r2, double theta2, double phi2){

	double *R;
	double *r1Vect;
	double *r2Vect;

	R 		= sphericalVector(r, theta, phi);
	r1Vect 	= sphericalVector(r1, theta1, phi1);
	r2Vect 	= sphericalVector(r2, theta2, phi2);

	double r12Vect[3]; 

	r12Vect[0] = R[0] - r1Vect[0] + r2Vect[0];
	r12Vect[1] = R[1] - r1Vect[1] + r2Vect[1];
	r12Vect[2] = R[2] - r1Vect[2] + r2Vect[2];

	return sqrt( pow(r12Vect[0], 2) + pow(r12Vect[1], 2) + pow(r12Vect[2], 2));
}

int main () {

	printf("Hello\n");
	//
}

/*

double doubleFoldingFunction(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma):

	def f(r1, theta1, phi1, r2, theta2, phi2): 

		r12Distance = r12(r, theta, phi, r1, theta1, phi1, r2, theta2, phi2)
		weight = (r1**2) * sin(theta1) * (r2**2) * sin(theta2)

		return densityM(r, rho0M, omega, betaM)*densityAlpha(r, rho0Alpha, betaAlpha)*gaussian(r12Distance, v0, gamma)*weight

	return vectorize(f)

def doubleFolding(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma):

	f = doubleFoldingFunction(r, theta, phi, rho0Alpha, rho0M, omega, betaAlpha, betaM, v0, gamma)

	return nquad(f, [	[0, 1], [0, pi], [0, 2*pi], 
						[0, 1], [0, pi], [0, 2*pi] ])


doubleFolding(r = 0, 
			 theta = 0, 
			 phi = 0, 
			 rho0Alpha = 0.4229, 
			 rho0M = 0.1644, 
			 omega = 0.4988, 
			 betaAlpha = 0.7024,
			 betaM = 0.3741,
			  v0 = 122.6225, 
			  gamma = 0.22)


def test(r1, theta1, phi1, r2, theta2, phi2):

	r1vec = array([	r1*sin(theta1)*cos(phi1),
						r1*sin(theta1)*sin(phi1), 
						r1*cos(theta1)])

	r2vec = array([	r2*sin(theta2)*cos(phi2),
						r2*sin(theta2)*sin(phi2), 
						r2*cos(theta2)])

	r12Distance = linalg.norm(r1vec - r2vec)

	return exp(-r12Distance**2)

nquad(test, [	[0, 1], [0, pi], [0, 2*pi], 
						[0, 1], [0, pi], [0, 2*pi] ])

*/