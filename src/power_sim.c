#include "bootstrap_required_libs.h"
#include "bootstrap_lib.h"
#include "curve_lib.h"

// export 
SEXP SimVitD_powersim( SEXP call, SEXP op, SEXP args, SEXP rho);

double intens_function( double t) ;
double glf( double v );
double EXTRACT_REALSXP_ELEMENT( SEXP x, int i );
int EXTRACT_INTSXP_ELEMENT( SEXP x, int i );

double intens_function( double t ){ double a = 6.0 * ( t - floor( t / M_PI ) * M_PI ) ; if( a > M_PI && a < 4 * M_PI ) return( 0.1 ); else return(1.0);  }

double EXTRACT_REALSXP_ELEMENT( SEXP x, int i )
{
  PROTECT( x = coerceVector( x, REALSXP ) );
  double a = REAL(x)[i];
  UNPROTECT(1);
  return(a);
}

int EXTRACT_INTSXP_ELEMENT( SEXP x, int i )
{
  PROTECT( x = coerceVector( x, INTSXP ) );
  int a = INTEGER(x)[i];
  UNPROTECT(1);
  return(a);
}


SEXP SimVitD_powersim( SEXP call, SEXP op, SEXP args, SEXP rho)
{
  int k, a, i, j, n, c, exper,  N, trtype, tp, prot_count = 0, apply_hold = FALSE, 
    sig_count, B, len_n, hypothesis, verbose;
  int *sx, *ssx, *npl, *ntr, *x, *xpl, *xtr, *clt;
  double t, inten, *pars_ptr, *pars_ind, *pars_pl, *pars_tr, start0, end0, scale0, scalehold0, vlev, p, sig_level, z, v;
  double low_asy, high_asy, icept, slope, glf_at_vlev, grp0, grp1, grpdiff ; // glf
  double amplitude, height, delta, omega, threshdyn, mindet; 
  const char *type, *testty;
  SEXP t0, v0, intens, pr, fn_intens, R_fn_intens, fn_glog, pars, typ, endpts, 
  nvec, baserate, bootsim, numexper, numrep, group1, group2, effsize, cltapp;
  SEXP holdtime, app_holdtime, testtype, alpha, power, glfpar, res, vb;
  PROTECT_INDEX ipx;
  
  // pull arguments from .External2 call
  args = CDR(args); 
  fn_intens = CAR(args); args = CDR(args); // intensity function
  pars = CAR(args); args = CDR(args);      // parameters for curve generation
  typ = CAR(args); args = CDR(args);     // type of treatment -- fixed, dynamic etc
  endpts = CAR(args); args = CDR(args);   // start and end of study (in units of years)
  baserate = CAR(args); args = CDR(args); // base rate of exposures for Poisson process 
  app_holdtime = CAR(args); args = CDR(args); // apply a holding time?
  holdtime = CAR(args); args = CDR(args); // mean holding time after infection
  numexper = CAR(args); args = CDR(args); // number of experiments in n
  nvec = CAR(args); args = CDR(args);    // numbers in treatment and control grp
  numrep = CAR(args); args = CDR(args);  // number of repetitions (N)
  testtype = CAR(args); args = CDR(args); // test type- count/proportions
  alpha = CAR(args); args = CDR(args);    // significance level for rejection
  bootsim = CAR(args); args = CDR(args);    // number of bootstrap samples
  glfpar = CAR(args);  args = CDR(args);    // glf parameters
  cltapp = CAR(args); args = CDR(args);     //clt approximation
  vb = CAR(args);                         //verbose
  
  verbose = EXTRACT_INTSXP_ELEMENT(vb,0);
  
  //int clt_approx = TRUE; // need to create a vector of these
  sx = (int *) R_alloc( 2, sizeof(int) ); sx[0] = 0; sx[1] = 0;
  ssx = (int *) R_alloc( 2, sizeof(int) ); ssx[0] = 0; ssx[1] = 0;
  
  // PROTECT( t0 = allocVector(REALSXP, 1) ); prot_count++;
  //REAL(t0)[0] = 0.0;
  // set up R functions to be called for intensity and generalized logistic
  //R_fn_intens = PROTECT(lang2(fn_intens, R_NilValue)); prot_count++;
  // need to load args for initial call to both functions and then re-protect result
  //SETCADR( R_fn_intens, t0 );
  //PROTECT_WITH_INDEX(intens = eval(R_fn_intens, rho), &ipx); prot_count++;
  
  start0 = EXTRACT_REALSXP_ELEMENT(endpts, 0) * M_PI;
  end0 = EXTRACT_REALSXP_ELEMENT(endpts, 1) * M_PI;
  
  len_n = EXTRACT_INTSXP_ELEMENT( numexper, 0 );
  
  npl = (int *) R_alloc( len_n, sizeof(int) );
  ntr = (int *) R_alloc( len_n, sizeof(int) );
  
  PROTECT( nvec = coerceVector(nvec, INTSXP) );
  for( exper=0; exper<len_n; exper++ )
  { 
    npl[exper] = INTEGER(nvec)[exper]; 
    ntr[exper] = INTEGER(nvec)[len_n + exper]; 
  }
  UNPROTECT(1);
  
  // CLT approximation indicators
  clt = (int *) R_alloc( len_n, sizeof(int) );
  for( exper=0; exper<len_n; exper++ ) clt[exper] = EXTRACT_INTSXP_ELEMENT(cltapp, exper);
  
  N = EXTRACT_INTSXP_ELEMENT(numrep, 0);
  
  scale0 = 1.0 / ( EXTRACT_REALSXP_ELEMENT(baserate,0) * 52.0 / M_PI ) ;
  
  if( EXTRACT_INTSXP_ELEMENT(app_holdtime,0) == 1 ) apply_hold = TRUE;
  
  if( apply_hold ) scalehold0 = EXTRACT_REALSXP_ELEMENT(holdtime,0) / (52.0 * M_PI) ; else scalehold0 = INFINITY ;
  
  sig_level = EXTRACT_REALSXP_ELEMENT(alpha,0);
  
  B = EXTRACT_REALSXP_ELEMENT(bootsim,0); 
  
  low_asy = EXTRACT_REALSXP_ELEMENT(glfpar,0);
  high_asy = EXTRACT_REALSXP_ELEMENT(glfpar,1);
  icept = EXTRACT_REALSXP_ELEMENT(glfpar,2);
  slope = EXTRACT_REALSXP_ELEMENT(glfpar,3);
  
  type = CHAR(STRING_ELT(typ, 0));
  // reserve zero for placebo simulation
  if( strcmp(type, "fixed-dose") == 0 ) trtype = 1;
  if( strcmp(type, "dynamic-dose") == 0 ) trtype = 2;
  
  testty = CHAR(STRING_ELT(testtype, 0));
  // proportions or count
  if( strcmp(testty, "count") == 0 ) hypothesis = 1;
  if( strcmp(testty, "proportions") == 0 ) hypothesis = 2;
  
  // extract the parameters for 25OHD curve simulation
  pars_pl = (double *) R_alloc( 11, sizeof(double) ); pars_tr = (double *) R_alloc( 11, sizeof(double) );
  pars_ind = (double *) R_alloc( 6, sizeof(double) );
  PROTECT( pars = coerceVector(pars, REALSXP) );
  for( i=0; i<11; i++ )
  {
    pars_pl[i] = REAL(pars)[i]; 
    pars_tr[i] = REAL(pars)[11+i];  
  }
  UNPROTECT(1);
  
  // define holders once and then overwrite for each repetition
  i = 0;
  int max_n = 0, boot_flag = 0;
  while( i < len_n )
  {
    k = npl[i] < ntr[i] ? ntr[i] : npl[i] ;
    max_n = k < max_n ? max_n : k ;
    if( clt[i] == 0 ) boot_flag = 1;
    i++;
  }
  
  if( boot_flag ) // reassign x vectors if needed and reuse
  {
    xpl = (int *) R_alloc( max_n, sizeof(int) ); 
    xtr = (int *) R_alloc( max_n, sizeof(int) ) ; 
  }
  
  // storage of results
  PROTECT(res = allocVector(VECSXP, 4)); prot_count++;
  SEXP names;
  PROTECT(names = allocVector(STRSXP, 4)); prot_count++;
  SET_STRING_ELT(names, 0, mkChar("power"));
  SET_STRING_ELT(names, 1, mkChar("effsize"));
  SET_STRING_ELT(names, 2, mkChar("group1"));
  SET_STRING_ELT(names, 3, mkChar("group2"));
  setAttrib(res, R_NamesSymbol, names);
  
  PROTECT(power = allocVector(REALSXP, len_n)); prot_count++;
  PROTECT(effsize = allocVector(REALSXP, len_n)); prot_count++;
  PROTECT(group1 = allocVector(REALSXP, len_n)); prot_count++;
  PROTECT(group2 = allocVector(REALSXP, len_n)); prot_count++;
  
  GetRNGstate(); // RNG call
  
  for( exper=0; exper<len_n; exper++ )
  {
    sig_count = 0;
    grp0 = 0.0;
    grp1 = 0.0;
    grpdiff = 0.0;
    for( k=0; k<N; k++ )
    {
      R_CheckUserInterrupt(); // check to see if stop has been requested
      sx[0] = 0; sx[1] = 0; ssx[0] = 0; ssx[1] = 0;
      for( a=0; a<2; a++ ) // loop through each arm
      {
        // do each of the arms
        n = ( a == 0 )  ? npl[exper] : ntr[exper] ; // sample size
        pars_ptr = ( a == 0 ) ? pars_pl : pars_tr ; // parameters for arm
        x = (a == 0) ? xpl : xtr ; // x if needed
        tp = (a == 0) ? 0 : trtype ; // type placebo/fixed etc.
        for( i=0; i<n; i++ )
        {
          t = start0; // reset time      
          c = 0; // count infections
          sim_curve_pars( pars_ind, pars_ptr, tp ); // pars for individual i
          amplitude = pars_ind[0]; height = pars_ind[1];
          delta = pars_ind[2]; omega = pars_ind[3];
          threshdyn = pars_ind[4]; mindet = pars_ind[5];
          while( t < end0 )  // cycle through exposure times until end
          {
            t += rexp( scale0 ); // increment time
            if( t > end0 ) break;
            //REAL(t0)[0] = t;
            //REPROTECT( intens = eval( R_fn_intens, rho), ipx ) ;
            //REPROTECT( intens = coerceVector(intens, REALSXP), ipx);
            inten = intens_function( t );
            // accept-reject step
            if( runif(0.0,1.0) <  inten ) //REAL(intens)[0] )
            {
              vlev = height + amplitude * cos( 2.0 * t - M_PI ) ;
              if( tp == 1 ) vlev += delta * ( omega + 0.5 * (1.0 - omega) * ( 1.0 + sin( 2.0*t - 1.5 * M_PI ) ) ) ; // fixed dose
              if( tp == 2 ) vlev = ( vlev < threshdyn ) ? threshdyn : vlev ; // dynamic dose 
              // minimum detection threshold
              vlev = ( vlev < mindet ) ? mindet : vlev ;
              glf_at_vlev = low_asy + (high_asy - low_asy) / ( 1.0 + exp(icept + slope * vlev) );
              if( runif(0.0,1.0) < glf_at_vlev ) // generate 0-1 for infection
              {
                c += 1; 
                if( hypothesis == 2 ) break; // no need to keep iterating for proportions
                if( apply_hold ) t += rexp( scalehold0 ); // apply holding time if set
              }
            }
          }
          if( clt[exper] ){ sx[a] += c; ssx[a] += c * c; } else x[i] = c;
        }
      }
      
      /*monitor group means and effect size*/
      grp0 = ( grp0 * k + (double) sx[0] / npl[exper] ) / (k + 1.0) ; 
      grp1 = ( grp1 * k + (double) sx[1] / ntr[exper] ) / (k + 1.0) ;
      z = (double) sx[0] / npl[exper] - (double) sx[1] / ntr[exper] ; // estimated effect
      grpdiff = ( grpdiff * k + z ) / (k + 1.0);
      /*** CLT (app) ***/
      if( clt[exper] )
      {
        v = ( (double) ssx[0] - (double) sx[0] * sx[0] / npl[exper] ) / ( npl[exper] * (npl[exper] - 1) ); // variance estimate
        v += ( (double) ssx[1] - (double) sx[1] * sx[1] / ntr[exper] ) / ( ntr[exper] * (ntr[exper] - 1) );
        p = 1.0 - pnorm( z / sqrt(v), 0.0, 1.0, 1, 0 );
      }else{
        p = two_samp_boot( xpl, xtr, npl[exper], ntr[exper], B, 2 ) ; 
      }
      if( p < sig_level ) sig_count += 1;
    }
    REAL(power)[exper] = (double) sig_count / (double ) N ;
    REAL(effsize)[exper] = grpdiff ;
    REAL(group1)[exper] = grp0;
    REAL(group2)[exper] = grp1;
    // print an update message here
    if( verbose ) Rprintf("\n... placebo sample size %d complete...", npl[exper]);
  }
  
  PutRNGstate(); // RNG end
  
  SET_VECTOR_ELT(res, 0, power); SET_VECTOR_ELT(res, 1, effsize);
  SET_VECTOR_ELT(res, 2, group1); SET_VECTOR_ELT(res, 3, group2);
  
  UNPROTECT(prot_count);
  
  return res;
}

// method def and reg

static const R_ExternalMethodDef ExtEntries[] = {
  {"SimVitD_powersim", (DL_FUNC) &SimVitD_powersim, 16 },
  {NULL, NULL, 0}
};

void R_init_SimVitD( DllInfo *info )
{
  R_registerRoutines( info, NULL, NULL, NULL, ExtEntries );
  R_useDynamicSymbols( info, FALSE );
}
