// coriolis.cpp

#include "Vector.h"
#include <iostream>
using namespace va;
using namespace std;

int main( void ) {

   const Vector ihat( 1., 0., 0. );   // points south
   const Vector jhat( 0., 1., 0. );   // points east
   const Vector khat( 0., 0., 1. );   // points up (radially from the earth)
   
   const double VX          = 0., VY = 0.5, VZ = 260., H = 0.;
   const double ALPHA       = atan( VZ / VY );
   const double V0          = sqrt( VX * VX + VY * VY + VZ * VZ );
   const double OMEGA_EARTH = 7.29e-5, G = 9.81, LAT = 30. * D2R;
   const double OMEGA       = OMEGA_EARTH;
   const double TWO_OMEGA   = 2. * OMEGA_EARTH;
   const double COS_LAT     = cos( LAT ), SIN_LAT = sin( LAT );
   const double COS_ALPHA   = cos( ALPHA ), SIN_ALPHA = sin( ALPHA );
   const Vector omega       = Vector( -COS_LAT, 0., SIN_LAT );   // unit vector
   const Vector g           = -G * khat;
   const Vector r0          = Vector( 0., 0., H );
   const Vector v0          = Vector( 0., V0 * COS_ALPHA, V0 * SIN_ALPHA );
   double T;
   if ( H == 0. )
      T  = 2. * V0 * SIN_ALPHA / ( G - ( 2. * OMEGA ) * COS_LAT * V0 * COS_ALPHA );
   else
      T = sqrt( 2. * H / G );
   
   double f1 = 0.5 * ( 1. - cos( 2. * OMEGA * T ) ) / OMEGA;
   double f2 = T - 0.5 * sin( 2. * OMEGA * T ) / OMEGA;
   
   Vector w1 = omega ^ g;
   Vector w2 = omega ^ w1;
   Vector g1 = ( omega ^ g ) / ( 2. * OMEGA );
   Vector g2 = omega ^ g1;
   Vector v1 = omega ^ v0;
   Vector v2 = omega ^ v1;
   
   Vector r, r1, r2, r_prime, g1_prime, g2_prime;
   double t2, f1_prime, f2_prime;
   cout << "t" << "\t" << "y0" << "\t" << "y1" << "\t" << "y2" << "\t" << "z" << endl;
   for ( double t = 0.; t <= T; t += 1. ) {
   
      t2 = t * t;
      f1 = 0.5 * ( 1. - cos( 2. * OMEGA * t ) ) / OMEGA;
      f2 = t - 0.5 * sin( 2. * OMEGA * t ) / OMEGA;
      f1_prime = 0.5 * ( 1. - cos( 2. * TWO_OMEGA * t ) ) / TWO_OMEGA;
      f2_prime = t - 0.5 * sin( 2. * TWO_OMEGA * t ) / TWO_OMEGA;
      g1_prime = ( omega ^ g ) / ( 2. * TWO_OMEGA );
      g2_prime = omega ^ g1_prime;
      r1 = r0 + v0 * t + 0.5 * g * t2;
      r2 = 0.5 * w2 * t2;
      r = r1 + r2 - f1 * ( v1 + g2 ) - f2 * ( g1 - v2 );
      r_prime = r1 + r2 - f1_prime * ( v1 + g2_prime ) - f2_prime * ( g1_prime - v2 );
      cout << t << "\t"                 // elapsed time
           << r1 * jhat << "\t"         // distance east with OMEGA = 0
           << r * jhat << "\t"          // distance east
           << r_prime * jhat << "\t"    // distance east with OMEGA = 2. * OMEGA_EARTH
           << r1 * khat << endl;        // height
   }
   return EXIT_SUCCESS;
}
