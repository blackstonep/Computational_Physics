#include "nr3.h"
#include "roots.h"

using namespace std;

struct Proj {
	double v, theta, k, T0;
	Proj(double v_in, double theta_in, double k_in ) { //constructor?
		v = v_in;
		theta = theta_in;
		k = k_in;
		double g = 9.80665;
		T0 = 2.0*v_in*sin(theta_in)/g;
	}
	~Proj() {cout << "Destroying Proj..." << endl;} //desctructor
	double set_k(double coef) {return k = coef; }; //method set_k
	double operator()(double t) { 			//operator defn
		return (T0/2 + 1/k)*(1 - exp(-k*t))-t;
	}
};


template <class T> double functor_evaluator(T & func, double x) {
	return func(x);
};

int main() {
	double pi = 2*acos(0.0);
	double v0 = 40;
	double theta0 = pi/3;
	double k;
	double k0 = 0.0001;
	double kmax = 0.5;
	double x1 = 10.0;
	double x2 = 15.0;
	double xacc = .0000001;
	double Tapp;
	ofstream numerical("numerical.dat");
	ofstream approximate("approximate.dat");

	Proj proj(v0, theta0, k0);
	double T0;
	T0 = 2.0*v0*sin(theta0)/9.80665;
	for (k = k0; k <= kmax; k += .001) {
		proj.set_k(k);
		double root;
		if (zbrac(proj, x1, x2)) {
			root = rtbrent(proj, x1, x2, xacc);
			Tapp = T0-k*T0*T0/6;
			numerical << k << " " << root << endl;
			approximate << k << " " << Tapp << endl;
		}
		else {
			cout << "SOMETHING WENT WRONG" << endl;
		}
	}

	return 0;
}

