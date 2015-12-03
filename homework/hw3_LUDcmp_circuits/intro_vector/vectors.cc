#include <iostream>
#include "nr3.h"

using namespace std;

int main() {
	VecDoub vec(3);
		double imax = 3.0;
		//Change Index: name[index] = newVal;
		for (double i=0.0; i < imax; i++) {
			vec[i]= i;
		};
		for (double i=0; i<imax; i++) {
			cout << "vec["<<i<<"] = "<< vec[i] << endl;
		};
		//Change Index: name[row][col] = Value
		//(Starting at index i=0)
	MatDoub mat(3,3);
		int kmax = 3;
		int jmax = 3;
		for (int k=0; k<kmax; k++) {
			for (int j=0; j<jmax; j++) {
				mat[k][j]=(k+1)*(j+1);
				cout << "mat["<<k<<"]["<<j<<"] = "<<mat[k][j]<< endl;
			};
		};


	return 0;
}
