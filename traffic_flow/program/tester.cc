#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"

using namespace std;

//==================================//
//Density functor: 					//
//	Takes position x and returns	//
//	density rho for a set para-		//
//	meter lambda and center x.		//
//Constructor takes rho bar and 	//
//	del_rho. 						//
//==================================//
struct Rho_start {
	const double lambda = 0.1;
	const double center = 0.5;
	double rho_bar; 
	double del_rho;

	Rho_start(double rho_bar_in, double del_rho_in) {
		rho_bar = rho_bar_in;
		del_rho = del_rho_in;
	}
	~Rho_start() {}

	void set_coef(double rho_bar_in, double del_rho_in) {
		rho_bar = rho_bar_in;
		del_rho = del_rho_in;		
	}

	double operator()(double x) {
		return rho_bar + del_rho * exp(- (x-center)*(x-center) / (lambda*lambda));		
	}

};

//==================================================//
//Create structure to be integrated:				//
//													//
//	dt(rho_i) = (2 * rho_i - 1) * dx(rho_i)			//
//													//
//Where												//
//													//
//	dx(rho_i) = (rho_{i+1} - rho{i-1}) / 2*delta_x	//
//													//
//==================================================//
struct Traffic {
	int dim;
	double delta_x;
	double timer;

	Traffic(int gridpoints_in, double time_in) {
		dim = gridpoints_in;
		delta_x = 1.0 / double(dim - 1);
		timer = time_in;
	}
	~Traffic() {}

	void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
		//
		//Each current value of rho is exactly equal
		//	to the same entry of ystart. Odeint
		//	updates ystart, and *rho is never called
		//	or referenced so no need to update values. 
		//
		//Establish periodic boundary condition
		//	for derivative. 
		//
		for (int i=0; i<dim; i++) {
			if (i == 0 || i == dim-1) {
				dydx[i] = (2.0 * y[0] - 1.0) * (y[1] - y[dim-2]) / (2.0*delta_x);
			} else {
				dydx[i] = (2.0 * y[i] - 1.0) * (y[i+1] - y[i-1]) / (2.0*delta_x);
			}
		}


	}
};


//==========================================//
//Want to solve the differential equation:	//
//	ddt(rho) + ddx(f) = 0					//
//for rho. 									//
//											//
//Constructor: Rhobar, delrho, 				//
//	gridpoints, and end-time snapshot.		//
//============================================
struct TrafficIntegrator {
	VecDoub *ystart;

	int dim;
	const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;	
		const Doub t_in = 0.0;
		Doub t_final;
		const Doub x_in = 0.0;
		const Doub x_final = 1.0;
	Output out;
	int t_steps;
	int x_steps;
	double delta_t;
	double delta_x;

	TrafficIntegrator(double rho_bar_in, double del_rho_in, 
						double gridpoints_in, double t_final_in) {
		Rho_start rho_start(rho_bar_in, del_rho_in);
		ystart = new VecDoub(gridpoints_in);

		t_final = t_final_in;
		dim = gridpoints_in;
		t_steps = 1000;
		x_steps = gridpoints_in-1;
		delta_t = (t_final - t_in)/t_steps;
		delta_x = (x_final - x_in)/x_steps;

		//
		//Need to create vector with initial rho values. 
		// 
		for (double i=0.0; i<double(gridpoints_in); i++) {
			double x = i*delta_x;
			(*ystart)[i] = rho_start(x);
		}
	}
	~TrafficIntegrator() {
		delete ystart;
	}

	void print() {
		ofstream outfile;
		ostringstream filename;
		//
		//Create filename string. Gnuplot for loop only works
		//	on integers, so need the name to only contain 
		//	integer values. I.e., t=10*t_final. 
		//
		filename << "trash_data_t=" << int(10.0*t_final) << ".dat" << ends;

		outfile.open(filename.str().c_str());
		outfile.setf(ios::left);
		outfile << "#====================================================" << endl;		
		outfile << setw(8) << "#" << "t = " << t_final << endl;
		outfile << "#====================================================" << endl;
		outfile << setw(16) << "#x" << "rho" << endl;
		outfile << "#====================================================" << endl;

		for (int i=0; i<dim; i++) {
			double x = double(i) * delta_x;
			outfile << setw(16) << x << (*ystart)[i] << endl;
		}

		outfile.close();
	}

	void printMatrix() {
		ofstream outfile;
		ostringstream filename;
		filename << "matrix_data.dat" << ends;

		outfile.open(filename.str().c_str());
		outfile.setf(ios::left);
		outfile << "#====================================================" << endl;		
		outfile << setw(8) << "#" << "t = " << t_final << endl;
		outfile << "#====================================================" << endl;
		outfile << setw(16) << "#x" << "rho" << endl;
		outfile << "#====================================================" << endl;

		for (int i=0; i<dim; i++) {
			outfile << setw(6) << (*ystart)[i];
		}

		outfile << endl;

		outfile.close();

	}

	void operator()() {
		for (double t = t_in; t < t_final; t+=delta_t) {

			Traffic traffic(dim, t_final);

			Odeint<StepperDopr5<Traffic> > ode(*ystart, t, t+delta_t,
				atol, rtol, h1, hmin, out, traffic);
			ode.integrate();

			printMatrix();			
		}
	}	
};




int main() {
	double rho_bar = 0.6;
	double del_rho = 0.0005;

	int dim = 1001;

	double tmax = 5.0;
	double del_t = 0.2;

	for (double t_run=0.01; t_run<=tmax; t_run+=del_t) {

		TrafficIntegrator integrator(rho_bar, del_rho, dim, t_run);

		integrator();
		integrator.print();

		if (t_run == 0.01) {t_run =0;}
	}


/*	TrafficIntegrator tester(rho_bar, del_rho, dim, 0.01);
	tester();
	tester.print();
*/
	return 0;
}