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

	VecDoub *rho;

	Traffic(int gridpoints_in, double time_in) {
		dim = gridpoints_in;
		delta_x = 1.0 / double(dim - 1);
		timer = time_in;

		rho = new VecDoub(dim);
	}
	~Traffic() {
		delete rho;
	}

	void operator()(const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
		for (int i=0; i<dim; i++) {
			if (i==dim-1) {
				(*rho)[i] = (*rho)[0];
			} else{
				(*rho)[i] = y[i];
			}
		}	
			//
			//Establish periodic boundary condition.
			//
		for (int i=0; i<dim; i++) {
			if (i == 0 || i == dim-1) {
				dydx[i] = (2.0 * (*rho)[0] - 1.0) * ((*rho)[1] - (*rho)[dim-2]) / (2.0*delta_x);
			} else {
				dydx[i] = (2.0 * (*rho)[0] - 1.0) * ((*rho)[i+1] - (*rho)[i-1]) / (2.0*delta_x);
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
	Rho_start *rho_start;
	VecDoub *ystart;
	Traffic *traffic;

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
		rho_start = new Rho_start(rho_bar_in, del_rho_in);
		ystart = new VecDoub(gridpoints_in);
		traffic = new Traffic(gridpoints_in, t_final_in);

		t_final = t_final_in;
		dim = gridpoints_in;
		t_steps = gridpoints_in-1;
		x_steps = gridpoints_in-1;
		delta_t = (t_final - t_in)/t_steps;
		delta_x = (x_final - x_in)/x_steps;

		//
		//Need to create vector with initial rho values. 
		// 
		for (double i=0.0; i<double(gridpoints_in); i++) {
			double x = i*delta_x;
			(*ystart)[i] = (*rho_start)(x);
		}
	}
	~TrafficIntegrator() {
		delete rho_start;
		delete ystart;
		delete traffic;
	}

	void print() {
		//traffic -> printWave();
		ofstream outfile;
		ostringstream filename;
		filename << "a_data_t=" << int(10.0*t_final) << ".dat" << ends;

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

		cout << "Quarter Density = " << (*ystart)[dim/4] << endl;		
	}

	void operator()() {
		for (double t = t_in; t < t_final; t+=delta_t) {
			Odeint<StepperDopr5<Traffic> > ode(*ystart, t, t+delta_t,
				atol, rtol, h1, hmin, out, *traffic);
			ode.integrate();

			cout << "Time = " << t << endl;
			cout << "Quarter Density = " << (*ystart)[dim/4] << endl << endl;			
		}
	}	
};




int main() {
	double rho_bar = 0.5;
	double del_rho = 0.0005;

	int dim = 101;

	double tmax = 1.2;
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