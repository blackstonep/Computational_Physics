#include "nr3.h"
#include "ludcmp.h"

using namespace std; 
const Doub pi = 2.0*acos(0.0);

//=======================
//Function for boundary 
//	cdtn at top of box.
//=======================
double top(double x) {
	if (x <= 0.5) {return 2.0*x;}
	else {return 2.0*(1.0-x);};
}

//=======================
//Function to return 
//	analytical value.
//	Order: 200
//=======================
double analytical(double x, double y) {
	const Doub nmax = 200.0;
	double ans = 0.0;
	for (double n=1.0; n<nmax; n++) {
		ans += sin(n*pi/2.0)/(n*n)*sin(n*pi*x)*sinh(n*pi*y)/sinh(n*pi);
	};
	return 8.0*ans/(pi*pi);
}

//========================
//Main structure to carry 
//	calculation. 
//========================
struct Temp {
	int L;
	int N;
	int runtime;
	double err; 
	VecDoub *b;
	VecDoub *T;
	MatDoub *a;
	ostringstream *filename;
	ofstream outfile;
	//
	//Constructor
	//
	Temp(int L_in) {
		L = L_in;
		N = (L+1)*(L+1);
		b = new VecDoub(N);
		T = new VecDoub(N);
		a = new MatDoub(N,N);
		filename = new ostringstream();
		(*filename) << "data_L=" << L << ".dat" << ends;
		//
		//Fill in matrix a and
		//	solution vector, b.
		//
		for (int j=0; j<=L; j++) {
			for (int ell=0; ell<=L; ell++) {
				int I = j*(L+1)+ell;
				//
				//Boundry Conditions
				//
				//Sides
				if (j==0||j==L||ell==0) {
					(*a)[I][I]=1.0; 
					(*b)[I]=0;
				}
				//Top
				else if (ell==L) {
					double x = double(j) / double(L);
					(*a)[I][I]=1.0;
					(*b)[I]=top(double(x));
				}
				//
				//Interior Points
				//
				else {
					(*a)[I][I]=-4.0;
					(*a)[I][I+1]=1.0;
					(*a)[I][I-1]=1.0;
					(*a)[I][I+L+1]=1.0;
					(*a)[I][I-(L+1)]=1.0;
					(*b)[I]=0.0;
				}
			}
		}

		//
		//Solve for the temperature
		//	supervector.
		//
		time_t tstart = time(0);
		LUdcmp alu((*a));
		alu.solve((*b),(*T));
		time_t tend = time(0);
		runtime = tend - tstart;
		err = (*T)[(N+1)/2] - analytical(0.5,0.5);
	}
	//
	//Destructor
	//
	~Temp() {
		delete b;
		delete T;
		delete a;
		delete filename;
		cout << "Destructing for L = " << L << endl;
	}
	//
	//Output the runtime of solution.
	//
	int running() {
		return runtime;
	}
	//
	//Calculate the error between numerical 
	//	and analytical at x = y = a/2.
	//	This has index 
	//		I = (L/2)*(L+2)
	//
	double error() {
		int I = (L/2)*(L+2);
		return (*T)[I] - analytical(0.5,0.5);
	}

	//
	//Create datafile and
	//	output result.
	//
	void output() {
		int width = 16;
		outfile.open((*filename).str().c_str());
		outfile.setf(ios::left);
		outfile << setw(2*width) << "#" << "L = " << L << endl;
		outfile << "#================================"
				<<  "================================\n";
		outfile << "#" << setw(width) << "x" << setw(width) << "y" 
					<< setw(width) << "Numerical" << setw(width) 
					<< "Analytical" << endl;
		outfile << "#================================"
				<<  "================================\n";

		for (int j=0; j<=L; j++) {
			for (int ell=0; ell<=L; ell++) {
				int I = j*(L+1) + ell;
				double del = 1.0 / double(L);
				double x = del*(double(j));
				double y = del*(double(ell));
				outfile << setw(width) << x << setw(width) << y 
							<< setw(width) << (*T)[I] << setw(width) 
							<< analytical(x,y) << endl << endl;
			}
		}
	}

	void altput() {
		int width = 16;
		outfile.open("altdata.dat");
		outfile.setf(ios::left);
		outfile << setw(2*width) << "#" << "L = " << L << ", y = a/2" << endl;
		outfile << "#================================"
				<<  "================================\n";
		outfile << "#" << setw(width) << "x" << setw(width) << "Numerical" << setw(width) 
					<< "Analytical" << endl;
		outfile << "#================================"
				<<  "================================\n";

		for (int j=0; j<=L; j++) {
				int ell = L / 2;
				int I = j*(L+1) + ell;
				double del = 1.0 / double(L);
				double x = del*(double(j));
				double y = del*(double(ell));
				outfile << setw(width) << x << setw(width) << (*T)[I] << setw(width) 
							<< analytical(x,y);
				if ()
		}
		outfile.close();		
	}
};


int main() {

	Temp temp(32);
	temp.output();

//	Temp alttemp(64);
//	alttemp.altput();

	//
	//Generate data for runtime vs. L
	//	(Takes a very long time to 
	//		compute! COMMENT OUT!!)
	//

/*	ofstream timefile;
	timefile.open("timedata.dat");
		timefile.setf(ios::left);
		timefile << "#" << "Run time vs. L" << endl;
		timefile << "#================================\n";
		timefile << "#" << setw(16) << "L" << "t\n" ;
		timefile << "#================================\n";

	for (int i=25; i<=80; i++) {
		Temp tempdata(i);
		timefile << setw(16) << i << tempdata.running() << endl;
	}
		timefile.close();
*/

	//
	//Calculating error at x = y = .5
	//

/*	ofstream errorfile;
	errorfile.open("errordata.dat");
		errorfile.setf(ios::left);
		errorfile << "#	Error vs. h = a/L" << endl;
		errorfile << "#==================================\n";
		errorfile << "#" << setw(16) << "h" << "delta T / T0\n";
		errorfile << "#==================================\n";

	for (int i = 16; i<=64; i+=2) {
		Temp errortemp(i);
		double h = 1.0 / double(i);
		errorfile << setw(16) << h << errortemp.error() << endl;

		cout << "Numerical = " << errortemp.error() << endl;

	}
		errorfile.close();
*/
	return 0;
}