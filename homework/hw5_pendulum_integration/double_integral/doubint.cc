#include "nr3.h"
#include "quadrature.h"
#include "romberg.h"

using namespace std;

//Define function xmax of y
double xmax(double yy) {
	return 0.6 * yy + 0.2 * sqrt( 10.0 - 16.0 * yy * yy );
};

//define function xmin of y
double xmin(double yy) {
	return 0.6 * yy - 0.2 * sqrt( 10.0 - 16.0 * yy * yy );
};

struct Integrand {

	//Constructor
	Integrand(int a) {a=0;}
	//Destructor
	~Integrand() {cout << "Destructing Integrand . . .\n";}

	//Return value of integrand after analytically 
	//	evaluating inner integral
	double operator()(double y) {
		return (exp( -1.0 * xmin(y) * sin(y) ) - 
			exp( -1.0 * xmax(y) * sin(y) ) ) / sin(y);
	};
};

int main() {
	//Declare maximum y values
	const Doub ymax = sqrt(10.0) / 4.0;
	const Doub ymin = - sqrt(10.0) / 4.0;

	//Instantiate integrand
	Integrand integrand(0.0);

	//Instantiate Midpnt objects.
	Midpnt<Integrand> integratorU(integrand, 0.0, ymax);
	Midpnt<Integrand> integratorL(integrand, ymin, 0.0);

	//Romberg/midpoint integrate over y
	double integral = qromo(integratorU) + qromo(integratorL);

	cout << "Integral = " << setprecision(12) << integral << endl;

	return 0;
};
