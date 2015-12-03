#include "nr3.h"

using namespace std;

ofstream impdata("implicit_data.dat");
ofstream expdata("explicit_data.dat");
ofstream analyt("analytical_data.dat");

struct Euler {
	double lambda;
	int steps;
	double h;
	
	Euler(int lambda_in, int steps_in) {
		lambda=double(lambda_in);
		steps=steps_in;
		h=1.0/double(steps);
	};

	~Euler() {cout <<"Destructing Euler...\n";};

	//Collect explicit data
	void explic() {
		double x;
		double y = 1.0;

		for (double i=0;i<steps+1;i++) {
			x = i*h;
			expdata << x << ", " << y << endl;

			//Updata y = y + (y' * h)
			y += (-lambda * y) * (h);
		}
	}


	void implicit() {
		double x;
		double y = 1.0;
		for (double i=0;i<steps+1;i++) {
			x = i*h;
			impdata << x << ", " << y << endl;

			//Updata y = y / (1 + h * lambda)
			y = y / (1.0 + h * lambda);
		}
	}

	void analytical() {
		for (double i=0; i<steps+1; i++) {
			analyt << i*h << ", " << exp( - lambda * i*h) << endl;
		}
	}
};


int main() {
	cout<<"Enter number of steps: ";
	int n;
	cin >> n;
	Euler euler(100, n);

	euler.explic();
	euler.implicit();
	euler.analytical();

	return 0;
}