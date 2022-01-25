// centrifugal.cpp

#include "Vector.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace va;   // vector algebra namespace
using namespace std;

int main( void ) {

   Vector ihat( 1., 0., 0. ), jhat( 0., 1., 0. ), khat( 0., 0., 1. );   // unit vectors for cartesian coordinate system
   
   const double OMEGA  = 7.29e-5;     // Earth's rotational velocity
   const double G      = 9.81;        // acceleration due to gravity (m/s^2)
   const double H      = 100.;        // height (m)
   const double LAT    = 45. * D2R;   // northern latitude (deg converted to rad)
   const double OMEGA2 = OMEGA * OMEGA;
   const double T      = sqrt( 2. * H / G );   // time to fall
   const double T2      = T * T;
   const double T3      = T * T2;
   const double T4      = T2 * T2;
   const double THETA   = OMEGA * T;
   const double COS     = cos( THETA );
   const double SIN     = sin( THETA );
   const double COS_LAT = cos( LAT );
   const double SIN_LAT = sin( LAT );
   
   const Vector omega   = Vector( -COS_LAT, 0., SIN_LAT );   // unit vector for Earth's rotation velocity
   const Vector g       = -G * khat;
   const Vector r0      = Vector( 0., 0., H );
   const Vector r0_star = r0 + g / OMEGA2;
   const Vector v0      = Vector( 0., 0., 0. );
   
   const Vector R1 = r0 + v0 * T + 0.5 * g * T2;
   const Vector R2 = R1 + g / OMEGA2;
   const Vector R3 = omega ^ R2;
   const Vector R4 = omega ^ R3;
   const Vector R5 = omega ^ r0_star;
   const Vector R6 = omega ^ R5;
   
   const Vector W1 = omega ^ g;
   const Vector W2 = omega ^ W1;
   const Vector R0 = omega ^ ( omega ^ r0 );
   const Vector V1 = omega ^ v0;
   const Vector V2 = omega ^ V1;
   
   Vector r_exact  = R1 + R4 - SIN * R5 - COS * R6 - THETA * SIN * ( R6 + V1 / OMEGA ) + THETA * COS * ( R5 - V2 / OMEGA );
   Vector r_approx = R1 - OMEGA * ( V1 * T2 + W1 * T3 / 3. ) + OMEGA2 * ( -0.5 * R0 * T2 + 0.5 * V2 * T3 + 0.125 * W2 * T4 );
                        
   cout.precision(6);
   cout << "r_exact (cm)  = " << std::fixed << r_exact * 100. << endl;    // output in cm
   cout << "r_approx (cm) = " << std::fixed << r_approx * 100. << endl;   // output in cm

   return EXIT_SUCCESS;
}
