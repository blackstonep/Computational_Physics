#include "nr3.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "odeint.h"
#include "roots.h"

using namespace std; 

struct Legendre {
	double lambda;

	Legendre(double lambda_in) {lambda=lambda_in;}
	~Legendre() {;}

	void operator() (const double x, VecDoub_I &y, VecDoub_O &dydx) {
		
		//Set initial conditions
		const double yinit = y[0];
		const double beta = y[1];
		const double lambda_op = y[2];

		//
		//If statement to set initial conditions
		//
		if (x == 1) {
			dydx[0] = beta;
			dydx[1] = lambda_op / 4.0 * (lambda_op / 2.0 - 1.0);
			dydx[2] = 0.0;
		}
		else {
			dydx[0] = beta;						//dy/dx
			dydx[1] = 2.0*x*beta / ( 1.0 - x*x ) 
				- lambda_op*yinit/(1.0-x*x);	//d(beta)/dx
			dydx[2] = 0.0;  					//d(lambda)/dx
		}
	}
};

struct Shooting {
	//============================
	//Set constants and x ranges
	//============================
	const int nvar = 3;
	const double atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
	const double x_in=1.0, x_final= 0.0;
	const int N_steps = 1000;
	const double delta_x = (x_in - x_final) / double(N_steps);
	Output out;
	ofstream outfile;

	//====================================
	//Declare vector dimensions, vectors 
	//		for storing x, y and y', 
	//		and lambda. 
	//====================================
	const int vecdim = N_steps+1;

	double lambda;

	VecDoub *y_out;
	VecDoub *beta_out;
	VecDoub *x_out;

	//=======================================
	//Constructor: 
	//	1.) Initial conditions
	//	2.) Instantiate Legendre
	//	3.) Integrate Legendre, outputting
	//			x, y, and y' to data vectors.
	//========================================

	Shooting(double lambda_in) {
		lambda = lambda_in;
		VecDoub ystart(nvar);
		y_out = new VecDoub(vecdim);
		x_out = new VecDoub(vecdim);
		beta_out = new VecDoub(vecdim);
		
		//Initial Conditions
		
		ystart[0]=1.0;
		ystart[1]=lambda/2.0;
		ystart[2]=lambda;

		//Instantiate Legendre object 
		//		to be integrated

		Legendre legendre(lambda_in);

		//====================================
		//Integrate legendre from 1 to 0
		// 		
		//Output (x, y, y') to data vectors
		//====================================
		for (int i=0; i< vecdim; i++) {
			double x = x_in - double(i)*delta_x;
			(*x_out)[i] = x;
			(*y_out)[i] = ystart[0];
			(*beta_out)[i] = ystart[1];

			Odeint<StepperDopr5<Legendre> > ode(ystart, x, x-delta_x, 
				atol, rtol, h1, hmin, out, legendre);
			ode.integrate();
		}
	}
	~Shooting() {
		//cout << "Destructing Shooting for lambda = " << lambda << ". . .\n";
		delete y_out;
		delete x_out;
		delete beta_out;
	}

	//=====================================
	//Check if y or y' are zero (small), 
	//		i.e. is lambda an eigenvalue? 
	//			Even or odd? 
	//=====================================
	bool even() {
		return (abs((*beta_out)[vecdim-1]) < 1.e-10);
	}

	bool odd() {
		return (abs( (*y_out)[vecdim-1] ) < 1.e-10);
	}

	void print_data(Shooting &object) {
		//======================================
		//Print out a data file for when lambda 
		//		satisfies bdry cdtns and we 
		//		want to plot. 
		//======================================
		ostringstream filename;
		filename <<"data_lambda="<<floor(lambda)<<".dat"<<ends;
		outfile.open(filename.str().c_str());
		outfile.setf(ios::left);
		outfile << "# lambda = " << lambda << endl;
		outfile << "#===========================================\n";
		outfile << "#" << setw(16) << "x" << setw(16) << "y" 
					<< setw(16) << "y'" << endl;
		outfile << "#===========================================\n";

		for (int step = 0; step < vecdim; step++) {
			if (object.even()) {
				outfile << setw(16) << - (*x_out)[step] << setw(16) 
					<< (*y_out)[step] << setw(16) << - (*beta_out)[step]
					<< endl;
			}
			if (object.odd()) {
				outfile << setw(16) << - (*x_out)[step] << setw(16) 
					<< - (*y_out)[step] << setw(16) << (*beta_out)[step]
					<< endl;
			}
		}
		for (int step = vecdim - 1; step >= 0; step--) {
			outfile << setw(16) << (*x_out)[step] << setw(16) << (*y_out)[step] 
					<< setw(16) << (*beta_out)[step] << endl;
		}

		outfile.close();
	}



	double value() {
		//======================================
		//Want to return the value of y at 0
		//		to check bdry cdtns.
		//======================================
		return (*y_out)[vecdim - 1];
	};

	double deriv() {
		//======================================
		//Or, return the derivative to 
		// 		check bdry cdtns.
		//======================================
		return (*beta_out)[vecdim - 1];
	}
};

//=======================================
//Return y(0) as a function of lambda. 
//		Use to find lambda for 
//		ODD solutions. 
//=======================================
double integrate_value(double lambda) {
	Shooting shoot(lambda);

	return shoot.value();
}

//=======================================
//Return y'(0) as a function of lambda. 
//		Use to find lambda for 
//		EVEN solutions. 
//=======================================
double integrate_deriv(double lambda) {
	Shooting shoot(lambda);

	return shoot.deriv();
}

//===========================
//Create file containing
//	values for y(0) and y'(0)
//	WRT lambda. 
//===========================
void output_data(double lmax) {
	ofstream outfile;
	outfile.open("master_data.dat");
	outfile.setf(ios::left);
	outfile << "#" << setw(16) << "lambda" << setw(16) << "y" 
				<< setw(16) << "y'" << endl;
	outfile << "#===========================================\n";

	for (double lambda=0.0; lambda < lmax; lambda += 0.1) {
		Shooting shoot(lambda); 

		outfile << setw(16) << lambda << setw(16) 
				<< shoot.value() << setw(16) << shoot.deriv()
				<< endl;
	}
}



int main() {
	//=========================
	//Output y(0)(lambda) and 
	//	y'(0)(lambda) data for 
	//	visual inspection. 
	//=========================
	output_data(35.0);

	//
	//We want to find the roots of 
	//	y(0)(lambda), which will give
	//	us the lambdas for odd solutions.
	//
	//We use zbrak from -1 < lambda < 35. 
	//

	VecDoub odd_left(0);
	VecDoub odd_right(0);
	int nroot_odd;
	zbrak(integrate_value, -1.0, 35.0, 200, odd_left, 
		odd_right, nroot_odd);



	VecDoub odd_roots(nroot_odd);

	//Store odd-function lambdas
	for (int i=0; i < nroot_odd; i++) {
		odd_roots[i] = zbrent(integrate_value, odd_left[i], 
			odd_right[i], 1.e-10);
	}

	for (int i = 0; i<nroot_odd; i++) {
		cout<< "Odd Bracket " << i+1 << " = " << "[ " << odd_left[i] 
			<< " , " << odd_right[i] << " ]" << endl;
		cout << "Odd root " << i+1 << ": lambda = " << odd_roots[i] << endl << endl;
	}	

	//
	//Likewise, we find roots of y'(0)(lambda)
	//	to find lambda values for even functions
	//
	VecDoub even_left(0);
	VecDoub even_right(0);
	int nroot_even;
	zbrak(integrate_deriv, -1.0, 35.0, 100, even_left, 
		even_right, nroot_even);

	VecDoub even_roots(nroot_even);

	//Store even-function lambdas
	for (int i=0; i < nroot_even; i++) {
		even_roots[i] = zbrent(integrate_deriv, even_left[i], 
			even_right[i], 1.e-10);
	}

	for (int i = 0; i<nroot_even; i++) {
		cout<< "Even Bracket "<< i+1 << " = " << "[ " << even_left[i] 
			<< " , " << even_right[i] << " ]" << endl;
		cout << "Even root " << i+1 << ": lambda = " << even_roots[i] << endl << endl;
	}	

	for (int i =0;i<nroot_odd;i++) {
		for (int j = 0; j<nroot_even; j++) {
			Shooting odd_shoot(odd_roots[i]);
			Shooting even_shoot(even_roots[j]);
			odd_shoot.print_data(odd_shoot);
			even_shoot.print_data(even_shoot);
		}
	}


	return 0;
}