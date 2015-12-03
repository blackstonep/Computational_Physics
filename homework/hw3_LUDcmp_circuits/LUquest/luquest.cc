#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

using namespace std;


int main() {
	VecDoub x(2), x1(2), x2(2), b(2);
	x[0]=1.0;
	x[1]=-1.0;
	x1[0]=0.999;
	x1[1]=-1.001;
	x2[0]=2.0;
	x2[1]=-2.40139;
	b[0]=0.22599;
	b[1]=0.338759;

	MatDoub a(2,2);
	a[0][0]=.789;
	a[0][1]=.56301;
	a[1][0]=1.182711;
	a[1][1]=.843952;

	LUdcmp alu(a);
	VecDoub x0(2);
	alu.solve(b,x0);
	cout<<"x0[0] = "<<setprecision(15) << x0[0]<<endl;
	cout<<"x0[1] = "<<setprecision(15)<< x0[1]<<endl;

	return 0;
};
