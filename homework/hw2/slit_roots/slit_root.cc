#include "nr3.h"
#include "roots.h"

double slit(double x) { 
  return (sin(x)*sin(x))/(x*x);
}

double dslit(double x) {
  return x-tan(x);
}

int main() {
	double  x1, x2, xacc, root;
	cout << "Set left limit: " << endl;
	cin >> x1;
	cout << "Set right limit: " << endl;
	cin >> x2;
	cout << "Set accuracy: " << endl;
	cin >> xacc;
	root = rtbis(dslit, x1, x2, xacc);
	cout << "Root at alpha = " << root << endl;
	cout << "Intensity Ratio = " << slit(root) << endl;


  return 0;
}
