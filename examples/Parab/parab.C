#include "nr3.h"
//
// define function
//
double parabola(double x) {
  return 0.5*x*x;
};
//
// define functor
//

struct Parab {
  double A;
  // constructor
  Parab(double A_in) { A = A_in; }
  // destructor
  ~Parab() { cout << " destructing Parab..." << endl; }
  // set coefficient
  double set_coefficient(double coef) { return A = coef; };
  // here's where the parabola is:
  double operator()(double x) { return A*x*x; }
};
//
// define function evaluator
//
double function_evaluator(double func(double), double x) {
  return func(x);
};
//
// define functor evaluator
//
template <class T> double functor_evaluator(T & func, double x) {
  return func(x);
};
//
// main routine
//
int main() {
  double x = 2.0;
  // create instantiation of Parab, with coefficient A initialized to 0.5
  Parab parab(0.5);
  cout << " Function call:      " << parabola(x) << endl;
  cout << " Functor call:       " << parab(x) << endl;
  cout << " Function_evaluator: " << function_evaluator(parabola,x) << endl;
  cout << " Functor_evaluator of function:  " << functor_evaluator(parabola,x) << endl;
  cout << " Functor_evaluator of functor:   " << functor_evaluator(parab,x) << endl;
  // try changing coefficient
  parab.set_coefficient(1.0);
  cout << " Functor_evaluator of functor:   " << functor_evaluator(parab,x) << endl;



  return 0;
}
