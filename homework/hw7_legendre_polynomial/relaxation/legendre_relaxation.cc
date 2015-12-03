#include "nr3.h"
#include "stepper.h"
#include "stepperdopr5.h"
#include "odeint.h"
#include "roots.h"
#include "tridag.h"

using namespace std; 

struct TriMatrix {
	const int dim = 1000000;
	int ell;
	ofstream outfile;

	//============================//
	//	Instantiate vectors to 	  //
	//	solve tridag matrix eqn:  //
	//							  //
	//		[a:b:c].u = r 		  //
	//							  //
	//============================//

	VecDoub *a; //lower diag

	VecDoub *b;	//diag

	VecDoub *c;	//upper diag

	VecDoub *r; //source terms: all zero
					//		except for r[0]

	VecDoub *u; //y-values vector

	VecDoub *x;

	//===============================
	//Constructor: Fill in diagonal 
	//	etc. vectors. 
	//===============================

	TriMatrix(int ell_in) {
		ell = ell_in; 
		double lambda = double(ell*(ell+1));

		a = new VecDoub(dim);
		b = new VecDoub(dim);
		c = new VecDoub(dim); 
		r = new VecDoub(dim);
		u = new VecDoub(dim);
		x = new VecDoub(dim);

		//
		//Make the x-value vector
		//
		double x_in = 1.0;
		double x_final = 0.0;
		double dx = (x_in - x_final) / double(dim);
		for (int i=0;i<dim; i++) {
			(*x)[i] = x_in - double(i)*dx - (0.5*dx);
		}

		//
		//Define constants for convenience 
		//
		VecDoub alpha(dim);
			for (int i=0;i<dim;i++) {
				alpha[i] = (*x)[i] / (dx*(1.0-(*x)[i]*(*x)[i]));
			}
		VecDoub beta(dim); 
			for (int i=0;i<dim;i++) {
				beta[i] = lambda / (1.0 - (*x)[i]*(*x)[i]);
			}
		const double gamma = 1.0 / (dx*dx);

		//============================
		//Input boundary values into 
		//	vectors. 
		//============================

		(*a)[0] = 0.0;
		(*b)[0] = beta[0] - 2.0 * gamma;
		(*c)[0] = gamma - alpha[0];
		(*r)[0] = - (lambda*dx/4.0 +1.0) * (gamma +alpha[0]);

		//If we have an even solution...

		if (ell%2==0) {
			(*a)[dim-1] = gamma + alpha[dim-1];
			(*b)[dim-1] = beta[dim-1] - gamma - alpha[dim-1];
			(*c)[dim-1] = 0.0;
			(*r)[dim-1] = 0.0;
		}

		//If it's odd...
		else {
			(*a)[dim-1] = gamma + alpha[dim-1];
			(*b)[dim-1] = beta[dim-1] - 3.0 * gamma + alpha[dim-1];
			(*c)[dim-1] = 0.0;
			(*r)[dim-1] = 0.0;
		}

		//=====================
		//Now we fill in the 
		//	interior points!
		//=====================

		for (int i=1 ; i<(dim-1) ; i++) {
			(*a)[i] = gamma + alpha[i];
			(*b)[i] = beta[i] - 2.0 * gamma;
			(*c)[i] = gamma - alpha[i];
			(*r)[i] = 0.0;
		}

		//=====================
		//Use tridag to solve
		//	for the u vector!
		//=====================

		tridag(*a, *b, *c, *r, *u);
	}
	~TriMatrix() {
		delete a;
		delete b; 
		delete c; 
		delete r;
		delete u;
		delete x;
	}

	void printmatrix() {
		for (int i = 0; i < dim ; i++) {
			for (int j = 0; j < dim; j++) {
				if (i==j) {
					cout << (*b)[i] << "     ";
				}
				else if (j-i == 1) {
					cout << (*c)[i] << "     ";
				}
				else if (i-j == 1) {
					cout << (*a)[i] << "     ";
				}
				else {
					cout << 0 << "     ";
				}
				if (j==(dim-1)) {
					cout << endl;
				}
			}
		}
	}

	void output() {
		ostringstream filename; 
		filename << "relaxdata_lambda="<<ell*(ell+1)<<".dat"<<ends;
		outfile.open(filename.str().c_str());
		outfile.setf(ios::left);
		outfile << "# lambda = " << ell*(ell+1) << endl;
		outfile << "#======================\n";
		outfile << "#" << setw(16) << "x" << "y\n";
		outfile << "#======================\n";
		if (ell%2==1) {
			for (int i=0; i<dim; i++) {
				outfile << setw(16) << -(*x)[i] << setw(16) << -(*u)[i] << endl;
			}
			for (int i=dim-1;i>=0;i--) {
				outfile << setw(16) << (*x)[i] << setw(16) << (*u)[i] << endl;
			}
		}
		else {
			for (int i = 0; i<dim; i++) {
				outfile << setw(16) << -(*x)[i] << setw(16) << (*u)[i] << endl;
			}
			for (int i=dim-1; i>=0; i--) {
				outfile << setw(16) << (*x)[i] << setw(16) << (*u)[i] << endl;
			}
		}

	}

	void operator()(VecDoub &ans) {
		ans.resize(dim);
		for (int i = 0; i<dim ; i++) {
			ans[i] = (*u)[i];
		}
	}
};

int main() {
	for (int ell = 0; ell<5; ell++) {
		TriMatrix foo(ell);
		foo.output();
	}

/*	TriMatrix test(20);
	VecDoub ans(50);
	test(ans);
	for (int i =0; i<50; i++) {
		cout << "ans["<<i<<"] = " << ans[i] << endl;
	}
	test.printmatrix();
*/
	return 0;
}