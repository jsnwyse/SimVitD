#include "bootstrap_lib.h"

double two_samp_boot( int *x, int *y, int n1, int n2, int B, int type )
{
  int i, j, k, b, w, a, n, *st, *ptr ;
  double *del, *d, smd = 0.0, smd2 = 0.0, z, obs, r;
  double a0 = 0.0, b0, c0, s00 = 0.0, s0 = 0.0, s1 = 0.0, s2 = 0.0, p;
  
  del = (double *) R_alloc( 2, sizeof(double) );
  del[0] = 1.0 / n1 ; del[1] = 1.0 / n2 ;
  
  st = (int *) R_alloc( 2, sizeof(int) );
  st[0] = 0; st[1] = 0;
  for( a=0; a<2; a++ )
  {  
    n = (a == 0) ? n1 : n2 ;
    ptr = (a == 0) ? x : y ;
    for( k=0; k<n; k++ ) st[a] += ptr[k];
  }
  obs = (double) st[0] * del[0] - (double) st[1] * del[1]; // observed difference between groups
  
  // compute a0 using Jackknife
  d = (double *) R_alloc( n1 + n2, sizeof(double) );
  for( k=0; k<n1; k++ ){ d[k] = (double) ( st[0] - x[k] ) / (n1 - 1.0) - st[1] * del[1]; smd += d[k]; }
  for( k=0; k<n2; k++ ){ d[n1+k] = st[0] * del[0] - (double) ( st[1] - y[k] ) / (n2 - 1.0); smd += d[n1+k];  }
  smd /= (n1 + n2);
  for(k=0; k<n1+n2; k++ ){ d[k] -= smd; smd2 += d[k] * d[k]; }
  for( k=0; k<n1+n2; k++ ){ a0 += R_pow( d[k], 3.0 ); }
  a0 /= R_pow( smd2, 1.5);
  a0 /= -6.0 ;

  for( b=0; b<B; b++)
  {
    st[0] = 0; st[1] = 0;
    for( a=0; a<2; a++ )
    {
      n = (a == 0) ? n1 : n2 ;
      ptr = (a == 0) ? x : y ;
      for( k=0; k<n; k++ )
      { 
        r = floor( (double) n * runif(0.0,1.0) ) ; 
        j = (int) r; 
        st[a] += ptr[j]; 
      }
    }
    z = (double) st[0] * del[0] - (double) st[1] * del[1];
    if( z < obs ) s00 += 1.0;
    if( z < 0.0 ) s0 += 1.0; 
    if( z == 0.0 ) s1 += 1.0;
  }
  
  z = qnorm( s00 / B, 0.0, 1.0, 1, 0 ) ;
  r =  ( s0 + 0.5 * (s1 + 1.0) ) / (B + 1.0) ;
  b0 = qnorm( r, 0.0, 1.0, 1, 0 ) ;
  c0 = (2.0 - a0 * z + a0 * b0 ) * z - b0 ;
  c0 /= ( 1.0 - a0 * z + a0 * b0 ) ;

  p  =  pnorm( c0, 0.0, 1.0, 1, 0 ) ;

  return( 1.0 - p ); // should never get to here
}
