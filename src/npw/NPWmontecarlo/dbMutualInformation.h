/*
*  BASIX: Mutual Information
*   contains a function that calculates the positive Mutual
*   Information from a 2D histogram, in CImg format
*
*  author: Daan Broekhuizen
*  date: january 2010
*/

#include "../NPWopengl/NPWGLImage.h"

#include <math.h>

double dbGetMI( const NPWGLImage & hist, const bool verbose = false )
{
  int dimX, dimY, dimZ;
  hist.GetDimensions( dimX, dimY, dimZ );

  std::vector<double> xMarginal( dimX, 0.0 );
  std::vector<double> yMarginal( dimY, 0.0 );
  double Ha = 0.0;
  double Hb = 0.0;
  double Hab = 0.0;

  // get marginals and joint entropy: H(X,Y) = - SUM( p(x,y) * log(p(x,y)) )
  for ( int y = 0; y < dimY; ++y ) {
    for ( int x = 0; x < dimX; ++x ) {
      double tmp = static_cast<double>( hist(x,y) );
      yMarginal[x] += tmp;
      xMarginal[y] += tmp;
      if ( tmp > 0 )
      {
        Hab -= tmp * std::log(tmp);
      }
    }
  }

  // get the marginal entropies: H(X) = - SUM( p(x) * log(p(x)) )
  for ( int x = 0; x < dimX; ++x ) {
    double tmp = xMarginal[x];
    if ( tmp > 0 ) Ha -= tmp * std::log(tmp);
  }
  for ( int y = 0; y < dimY; ++y ) {
    double tmp = yMarginal[y];
    if ( tmp > 0 ) Hb -= tmp * std::log(tmp);
  }

  // Finally, get the MI
  double MI =  Ha + Hb - Hab ;

  if (verbose) {
    std::cout << "MI : " << Ha << " + " << Hb << " - " << Hab << " = " << MI << "\t"
              << "hist dimensions: " << dimX << " x " << dimY << std::endl;
  }
  return MI;

}
