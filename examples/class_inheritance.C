#include "nr3.h"
//
// base class
//
struct Parab {
  double A;
  // constructor
  Parab(double A_in) : A(A_in) { 
    cout << " constructing Parab..." << endl;
    }
  // destructor
  ~Parab() { cout << " destructing Parab..." << endl; }
  // set coefficient
  double set_coefficient(double coef) { return A = coef; };
  // here's where the parabola is:
  double operator()(double x) { return function(x); }
  double function(double x) { return A*x*x; }
  double virtual df(double x) = 0;    // virtual function - needs to be 
                                      // defined in derived class
};
//
// two derived classes: one computes derivative analytically...
//
struct Parab_analytical : Parab {
  Parab_analytical(double A_in) : Parab(A_in) {   // constructor
    cout << " constructing Parab_analytical..." << endl;
  };
  double df(double x) {
    return 2.0*A*x;
  }
};
//
// ... the other evaluates it numerically:
//
struct Parab_numerical : Parab {
  double dx;
  Parab_numerical(double A_in, double dx_in) : 
    Parab(A_in), dx(dx_in) {                     // constructor
    cout << " constructing Parab_numerical..." << endl;
  }; 
  double df(double x) {
    //    return (function(x+dx) - function(x-dx))/(2.0*dx);
    return (function(x+dx) - function(x))/dx;
  }
};
//
// Main code
// 
int main() {
  //
  // allocate derived classes:
  //
  double A = 0.5;
  double dx = 0.01;
  // can't compile base class directly because of virtual functions - try
  // uncommenting the next line:
  //  Parab functor(A);
  Parab_analytical analytical_deriv(A);
  Parab_numerical numerical_deriv(A,dx);
  double x = 1.0;
  //
  // change A...
  //
  A = 1.0;
  analytical_deriv.set_coefficient(A);
  numerical_deriv.set_coefficient(A);
  cout << " Analytical derivative at x = " << x << " : " 
       << analytical_deriv.df(x) << endl;
  cout << " Numerical  derivative at x = " << x << " : " 
       << numerical_deriv.df(x) << endl;
  //
  // this is how we can handle pointers to a class
  //
  // note: pointer to base class...
  Parab * parab_pointer;
  int answer = 0;
  cout << " Analytical [1] or numerical [2] derivative? ";
  cin >> answer;
  // ...but can allocate derived class:
  if (answer == 1)
    parab_pointer = new Parab_analytical(A); 
  else if (answer == 2)
    parab_pointer = new Parab_numerical(A,dx); 
  else {
    cout << " Don't understand answer: " << answer << endl;
    return 1;
  }
  A = 1.0;
  // use "->" when using pointer to class
  parab_pointer->set_coefficient(A); 
  cout << " Analytical derivative at x = " << x << " using pointer : " 
       << parab_pointer->df(x) << endl;
  // don't forget to delete allocated memory to free up memory 
  // (otherwise you may introduce "memory leaks")
  delete parab_pointer;    
};
