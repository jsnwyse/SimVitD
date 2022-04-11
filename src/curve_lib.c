#include "curve_lib.h"
//#include "bootstrap_required_libs.h"

// functions to simulate random cosines for 25OHD curves

void sim_curve_pars( double *ind_pars,  double *sim_pars, int type )
{
  // sim_pars: in order of entry: k.H, k.A, k.T, H.0, A.0, T.0, delta, delta.rate, tau, shape1, shape2
  // ind_pars: in order of entry: A, H, delta, omega, T, tau
  
  double a;
  a = sim_pars[4] / sim_pars[1] ; 
  a = a * a ;
  ind_pars[0] = rgamma( a, sim_pars[4] / a ) ; // amplitude
  
  ind_pars[1] = rnorm( sim_pars[3], sim_pars[0] ); // height
  
  if( type == 1 ) // for fixed-dose only
  {
    // simulation from truncated exponential
    a = runif(0.0,1.0);
    ind_pars[2] = sim_pars[6] + log( exp(-sim_pars[7] * sim_pars[6]) + ( 1.0 - exp( -sim_pars[7] * sim_pars[6] ) ) * a ) / sim_pars[7] ;
    ind_pars[3] = rbeta( sim_pars[9], sim_pars[10] ) ; // weighting (omega)
  }else{ 
    ind_pars[2] = 0.0; ind_pars[3] = 1.0;
  }
  
  if( type == 2 ) // for dynamic-dose only
  {
    a = sim_pars[5] / sim_pars[2] ;
    a = a * a ; 
    ind_pars[4] = rgamma( a, sim_pars[5] / a ) ; // threshold
  }else ind_pars[4] = 0.0;
  
  ind_pars[5] = sim_pars[8];
  
  return;
}

