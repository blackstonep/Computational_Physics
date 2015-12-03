#include <iostream>
#include "nr3.h"
#include "interp_1d.h"

using namespace std;

struct InFall {

	int M;

	Poly_interp *data;
	VecDoub *rr, *tt;

	//constructor
	InFall(int M_in) { 
	
		M = M_in;
	
		const Doub PI = 2.0 * acos( 0.0 ) ;
		const Doub R0 = 149.6e+9 ;
		const Doub G = 6.67e-11;
		const Doub Msun = 1.989e+30 ;
		const Doub T0 = sqrt( ( R0 * R0 * R0) / ( 2.0 * G * Msun ) ) ;
		//VecDoub tt(M), rr(M); 
		tt = new VecDoub(M);
		rr = new VecDoub(M); 

		//construct data vectors
		for( int i=0 ; i < M ; i++) {
			double eta = (i + 1) * PI / M;
			(*rr)[i] = (R0 / 2.0) * (1.0 + cos(eta));
			(*tt)[i] = (T0 / 2.0) * (eta + sin(eta));
		};

		data = new Poly_interp(*tt, *rr, 4);	

	};
	//destructor
	~InFall() {
		cout<<"destructing InFall . . ."<<endl;
		delete data;
		delete tt;
		delete rr;
	};
	//Interpolator
	double operator()(double t) {
		//Convert t into seconds
		int tsec = t*24.0*3600.0;

		//Interpolate Radius value
		return data -> interp(tsec);
	};
};

int main() {
	int numpoints;
	ofstream rdata("rdata.dat");

	cout<< "Number of data points: ";
	cin>>numpoints;

	InFall infall(numpoints);

	for (int i = 10; i<64; i+=10) {
		cout<< "After " << i << " days, Radius = " << infall(double(i))<<endl;
		rdata<< i << ", "<<infall(double(i))<<endl;
	};

	return 0;
};