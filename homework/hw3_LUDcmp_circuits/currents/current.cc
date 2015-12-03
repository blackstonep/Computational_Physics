#include <iostream>
#include "nr3.h"
#include "ludcmp.h"

using namespace std;

int main() {
	int n;
	cout << "Choose number of resistors: ";
	cin >> n;
	MatDoub a(n,n);
	VecDoub resistor(n), voltage(n), rhs(n), current(n);
	a[0][0] = 1.0;
	rhs[0]=0.0;

//Setting resistor values and
//building matrix.
	for (int i=0;i<n;i++) {
		cout<<"Enter Resistor "<<i+1<<" = ";
		cin>>resistor[i];
		if (i>0) {
			a[i][0]=resistor[0];
			a[0][i]=-1.0;
			a[i][i]=resistor[i];
		};
	};
	cout<<"RESISTORS SET"<<endl;


//Setting voltage values and
//building RHS vector.
	for (int i=0;i<n;i++) {
		cout<<"Enter Voltage "<<i+1<<" = ";
		cin>>voltage[i];
		if (i>0) {
			rhs[i]=voltage[0]+voltage[i];
		};
	};

//Solving for currents.
	LUdcmp alu(a);
	alu.solve(rhs, current);
//Display currents.
	for (int i=0;i<n;i++) {
		cout<<"Current "<<i+1<<" = "<<current[i]<<endl;
	};
	cout<<"BEGIN PART D"<<endl;

//For part d:
	ofstream ddata("ddata.dat");
	for (int m=2;m<51;m++) {
		VecDoub curr(m), rs(m);
		MatDoub da(m,m);
		rs[0]=0.0;     //First
		da[0][0]=1.0;  //components
		for (int i=1;i<m;i++) {	//Set a and rs.
			for(int j=1;j<m;j++) { //Clean out da.
                                da[i][j]=0;
			};
			da[i][0]=1.0;
			da[0][i]=-1.0;
			da[i][i]=1.0;
			rs[i]=1.0;
		};

//Solve and output data.
		LUdcmp dalu(da);
		dalu.solve(rs,curr);
		ddata<<m<<", "<<curr[0]<<", "<< (m-1.0)/m <<endl;
	};

//For part e
	cout<<"STARTING PART E"<<endl;
	ofstream edata("edata.dat");
	for (int m=2; m<101;m++) {
		VecDoub curr(m), rs(m);
		MatDoub ea(m,m);
		rs[0]=0.0;
		ea[0][0]=1.0;
		for (int i=1;i<m;i++) {
			for (int j=1;j<m;j++) { //Clean out ea
				ea[i][j]=0.0;
			};
			rs[i]=2.0;
			ea[0][i]=-1.0;
			ea[i][0]=1.0;
			ea[i][i]=double(i);
		};
		LUdcmp ealu(ea);
		ealu.solve(rs,curr);
		if (m==100) {cout<<"For n=100, Current 1 = "<<curr[0]<<endl;};
		edata<<m<<", "<<curr[0]<<endl;
	};



	return 0;
};
