#include "nr3.h"
#include "odeint.h"
#include "stepper.h"
#include "stepperdopr5.h"
//================================================================
//
// Example code for ODEINT
//
// Solves 
//    theta'' + alpha^2 * theta = 0
// as
//    theta' = kappa
//    kappa' = - alpha^2 * theta
//================================================================

struct Harmonic {
  Doub alpha;
  Doub alpha2;
  Harmonic(Doub alpha_in) : 
    alpha(alpha_in) { alpha2 = alpha*alpha; };
  void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
    const Doub theta = y[0];
    const Doub kappa = y[1];   // = dtheta/dt
    dydx[0] = kappa;
    dydx[1] = - alpha2 * theta;
  }
};

//================================================================
//
// Main code
//
//================================================================

int main() {
  //==============================================================
  // Create and initialize output file
  //==============================================================
  ofstream outfile;
  outfile.open("datafile");
  outfile.setf(ios::left);
  outfile << setw(16) << "# time " <<
    setw(16) << "theta" << endl;
  outfile << "#====================================================" << endl;
  const int nvar = 2;
  const Doub atol = 1.e-5, rtol = atol, h1=0.001, hmin=0.0;
  const Doub t_in = 0.0;
  const Doub t_final = 100.0;
  const int N_steps = 1000;
  const Doub delta_t = (t_final - t_in)/N_steps;
  VecDoub ystart(nvar);
  Output out;

  const Doub alpha = 1.0;
  //
  // initial conditions
  //
  ystart[0] = 1.0;
  ystart[1] = 0.0;
  //
  Harmonic oscillator(alpha);  
  for (Doub t = t_in; t < t_final; t += delta_t) {
    Odeint<StepperDopr5<Harmonic> > ode(ystart,t,t+delta_t,atol,
					rtol,h1,hmin,out,oscillator);
    ode.integrate();
    outfile << setw(16) << t << setw(16) << ystart[0] << endl;
    // cout << setw(16) << t << setw(16) << ystart[0] << endl;
  }
  //==============================================================
  // Close output file
  //==============================================================
  outfile.close();
}
