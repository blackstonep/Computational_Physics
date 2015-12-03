#include "nr3.h"
#include "quadrature.h"
#include "romberg.h"

struct Integrand {
	double theta0;
	//constructor
	Integrand(double theta_in) {theta0 = theta_in; }
	//destructor
	~Integrand() {cout <<"destructing Integrand . . . "<<endl;};
	//Set Theta
	double set_theta(double theta) {return theta0 = theta; };
	//Actual Function:
	double operator()(double x) {return 1.0 / sqrt( cos(x) - cos(theta0) ); };
};


int main() {
	//Take arbitrary theta for targeted evaluation
	Doub theta0;
	cout << "Set theta0 (deg.) = \n";
	cin >> theta0;


	const Doub PI = 2.0 * acos(0.0);
	const Doub a = 0;
	const Doub b = theta0 * PI / 180.0;
	const Doub theta_max = PI / 2.0;

	//Testing Midsqu out
	Integrand integrand(theta0);
	Midsqu<Integrand> integrator(integrand, a, b);

	Doub integral = qromo(integrator);
	cout << "Integral for theta0 = " << theta0 << " is: " << integral << endl;

	//Creating data file
	ofstream data("data.dat");
	for (double i = 0.1 ; i < theta_max * 180.0 / PI ; i += 1.0) {
		
		//Define limits
		const Doub a = 0.0; 
		const Doub b = i * PI / 180.0;

		//Reset integral
		integrand.set_theta( b );
		Midsqu<Integrand> integrator_loop(integrand, a, b);

		Doub integral_loop = qromo(integrator_loop);

		//Output data: theta, integral, T / T0
		data << i << ", " << integral_loop << ", " 
			<< sqrt(2.0) * integral_loop / PI << endl;
	};

	return 0;
};