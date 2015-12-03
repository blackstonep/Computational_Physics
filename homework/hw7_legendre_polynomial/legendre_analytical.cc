#include "nr3.h"

using namespace std;

double legendre0(double x) {
	return 1.0;
}

double legendre1(double x) {
	return x;
}

double legendre2(double x) {
	return .5*(3.0*x*x - 1);
}

double legendre3(double x) {
	return .5* (5.0*x*x*x - 3.0*x);
}

double legendre4(double x) {
	return .125*(35.0*x*x*x*x - 30.0*x*x + 3.0);
}

int main() {
	ofstream data("analytical.dat");
	for (double i=-1.0;i<1.0;i+=0.0001) {
		data << i << " " << legendre0(i) << " " << legendre1(i) 
				<< " " << legendre2(i) << " " << legendre3(i) << " " 
				<< legendre4(i) << endl;
	}

	return 0;
}