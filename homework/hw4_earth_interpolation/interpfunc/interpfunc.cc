#include "nr3.h"
#include "interp_1d.h"
#include "polcoef.h"

double func(double x) {
	return ( exp(- x * x )) / (1 + (25 * x * x ));
};



struct FuncInterp {

	int N;

	VecDoub *xx, *yy;
	Poly_interp *data;

	//constructor
	FuncInterp(int N_in) { 
		N = N_in; 

		//Allocating Data vectors
		xx = new VecDoub(N);
		yy = new VecDoub(N);
		

		double dx = 2.0 / (N-1);

		//Constructing Data Vectors
		for (int i = 0 ; i < N ; i++) {
			double ii = double(i);
			(*xx)[i] = -1.0 + ii*dx;
			(*yy)[i] = func( (*xx)[i] );
		};

		//Construct interpolator
		data = new Poly_interp(*xx, *yy, N);
	}; 
	//destructor
	~FuncInterp() {
		cout<<"Destructing FuncInterp . . ."<<endl;
		delete xx;
		delete yy;
		delete data;
	};

	//Return error for previous input
	double error() {return data -> dy;};	

	//Interpolate y value
	double operator()(double x) {
		return data -> interp(x);
	};
};

void output_data(int num) {
	int Nint = 10*num;
	double dxint = 2.0 / double(Nint - 1);

	FuncInterp interpolator(num);

	ofstream outfile;
	ostringstream filename;
	filename << "data_N=" << num << ".dat" << ends;
	outfile.open(filename.str().c_str());



	for (int i = 0 ; i < Nint ; i++) {
		double xint = -1.0 + i*dxint;
		outfile << xint << ", " << 
			interpolator(xint) << ", " << 
			func(xint) << ", " << 
			abs (interpolator(xint) - func(xint) ) << ", "<< 
			abs( interpolator.error() )<< endl;
	};

};


int main() {
	
	//Output data files for interpolation.
	output_data(5);
	output_data(20); 

	cout<<"Set Coef dimension?: ";
	int dim;
	cin >> dim;

	//Finding coefficients
	VecDoub xa(dim), ya(dim), coe(dim), cof(dim);
	double dx = 2.0 / double(dim-1);

	cout<<"For polcoe:" ;
	
	for (int i=0 ; i < dim ; i++) {
		xa[i] = -1.0 + i * dx;
		ya[i] = func( xa[i] );
	};

	polcoe(xa, ya, coe);
	for (int i=0 ; i < dim ; i++) {
		cout << "coe["<<i<<"] = " << coe[i] << endl;
	};

	cout<<"For polcof: " ;

	polcof(xa, ya, cof);
	for (int i=0 ; i < dim ; i++) {
		cout << "cof["<<i<<"] = " << cof[i] << endl;
	};




	return 0;
};