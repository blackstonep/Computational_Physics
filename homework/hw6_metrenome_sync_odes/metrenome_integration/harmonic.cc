#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;

const double pi = 2.0*acos(0.0);

//=========================================
//Standard struct for uncoupled metronome
//=========================================
struct Harmonic {
  double alpha;

  Harmonic(double alpha_in) {
  	alpha = alpha_in;
  }

  void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx) {
  	
  	const double epsilon=0.1;
  	const double theta0=0.1;

    const Doub theta = y[0];
    const Doub beta = y[1];   // = dtheta/dt
    dydx[0] = beta;
    dydx[1] = - alpha * sin(theta) - 
    	epsilon*( (theta / theta0)*(theta / theta0) - 1.0)*beta;
  }
};

//==================================
//Create coupling struct for Part E
//==================================

struct Coupling {
  double alpha;
  double eta;
  double otherTheta;

  Coupling(double alpha_in, double eta_in, double otherTheta_in) {
  	alpha = alpha_in;
  	eta = eta_in;
  	otherTheta = otherTheta_in;
  }

  void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx) {
  	
  	const double epsilon=0.1;
  	const double theta0=0.1;


    const Doub theta = y[0];
    const Doub beta = y[1];   // = dtheta/dt
    dydx[0] = beta;
    dydx[1] = - alpha * sin(theta) - 
    	epsilon*( (theta / theta0)*(theta / theta0) - 1)*beta 
    	- eta*(theta - otherTheta);
  }

  //Update theta2 within theta1 structure (visa versa)
  void changevar(double update) {
  	otherTheta = update;
  }
};


int main() {
	//================================================
	//
	//Create & initialize output file, alpha = 1.0
	//
	//================================================
	ofstream outfile;
	outfile.open("dataD_alpha1");
	outfile.setf(ios::left);
	outfile << setw(16) << "# times " << 
		setw(16) << "theta" << endl;
	outfile << "#====================================================" << endl;

	//
	//Constants and time values
	//
	const int nvar = 2;
	const double atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
	const double t_in = 0.0;
	const double t_final = 1000.0; 
	const int N_steps = 10000;
	const double delta_t = (t_final - t_in) / N_steps;
	VecDoub ystart(nvar);
	Output out;

	const Doub alpha1 = 1.0;

	//
	//Initial conditions
	//
	ystart[0] = 0.0;
	ystart[1] = 0.1;

	Harmonic oscillator1(alpha1);
	for (double t = t_in; t < t_final; t += delta_t) {
		Odeint<StepperDopr5<Harmonic> > ode(ystart, t, t+delta_t,
			atol, rtol, h1, hmin, out, oscillator1);
		ode.integrate();
		outfile << setw(16) << t << setw(16) << ystart[0] << endl;
	}

	outfile.close();

	//========================
	//
	//Now for alpha = 0.25
	//
	//========================

	outfile.open("dataD_alpha25");
	outfile.setf(ios::left);
	outfile << setw(16) << "# times " << 
		setw(16) << "theta" << endl;
	outfile << "#====================================================" << endl;

	const Doub alpha25 = 0.25;

	ystart[0] = 0.0;
	ystart[1] = 0.1;

	Harmonic oscillator25(alpha25);
	for (double t = t_in; t < t_final; t += delta_t) {
		Odeint<StepperDopr5<Harmonic> > ode(ystart, t, t+delta_t, 
			atol, rtol, h1, hmin, out, oscillator25);
		ode.integrate();
		outfile << setw(16) << t << setw(16) << ystart[0] << endl;
	}

	outfile.close();

	//================
	//
	//	  Part E
	// 
	//================

	outfile.open("dataE");
	outfile.setf(ios::left);
	outfile << setw(16) << "# times " << 
		setw(16) << "theta" << endl;
	outfile << "#====================================================" << endl;

	const double eta = 0.01;

	VecDoub ystart1(nvar);
	VecDoub ystart2(nvar);

	//
	//Adjust initial conditions
	//VERY SUPER AWESOME COOL
	//
	ystart1[0]=0.0;
	ystart1[1]=0.1;
	ystart2[0]=0.0;
	ystart2[1]=0.0;

	//
	//Instantiate coupled oscillators
	//
	Coupling coupled_1(alpha1, eta, ystart2[0]);
	Coupling coupled_2(alpha1, eta, ystart1[1]);

	for (double t; t < t_final; t += delta_t) {
		//
		//Make ODE's, update with new theta's.
		//
		Odeint<StepperDopr5<Coupling> > ode1(ystart1, t, t+delta_t, 
			atol, rtol, h1, hmin, out, coupled_1);
		Odeint<StepperDopr5<Coupling> > ode2(ystart2, t, t+delta_t,
			atol, rtol, h1, hmin, out, coupled_2); 

		//
		//Update local positions for theta1 and 2
		//
		ode1.integrate();
		ode2.integrate();

		//
		//Output data to file
		//
		outfile << t << ", " << ystart1[0] << 
			", " << ystart2[0] << endl;

		//
		//Update coupled position values
		//
		coupled_1.changevar(ystart2[0]);
		coupled_2.changevar(ystart1[0]);
	}

	outfile.close();




	return 0;
}