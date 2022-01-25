// lorentz.cpp: Solution to equation of motion for Lorentz force

#include "Vector.h"
#include <iostream>
#include <cmath>
using namespace va;   // vector algebra namespace
using namespace std;

int main( void ) {

   const Vector ihat( 1., 0., 0. );   // unit vector along x-axis
   const Vector jhat( 0., 1., 0. );   // unit vector along y-axis
   const Vector khat( 0., 0., 1. );   // unit vector along z-axis
   
   const double M     = 1.;
   const double Q     = 1.;
   const double EE    = 1.;
   const double BB    = 1.;
   const double OMEGA = Q * BB / M;
   const double T     = 8. * M_PI / ( OMEGA );
   double t, theta;
   
   const Vector r0( 0., 0., 0. );
   const Vector v0( 0., 0., 0. );
   const Vector E( 0., 0., EE );   // electric field along z-axis
   const Vector B( BB, 0., 0. );   // magnetic field along x-axis
   const Vector b = unit( B );     // unit vector along magnetic field
   
   Vector r;
   cout << "t" << "\t" << "y" << "\t" << "z" << endl;
   for ( t = 0; t <= T; t += 0.01 * T ) {
   
      theta = OMEGA * t;
      r = r0 + v0 * t + 0.5 * ( OMEGA * t * t / BB ) * b * ( b * E ) +
          ( ( 1. - cos( theta ) ) / OMEGA ) * ( -( b ^ v0 ) - ( b ^ ( b ^ E ) ) / BB ) +
          ( t - sin( theta ) / OMEGA ) * ( ( b ^ ( b ^ v0 ) ) - ( b ^ E ) / BB );
      cout << t << "\t" << r * jhat << "\t" << r * khat << endl;
   }
   return EXIT_SUCCESS;
}
