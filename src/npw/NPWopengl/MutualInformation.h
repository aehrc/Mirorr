/*
 * MutualInformation.h
 *
 *  Created on: 16/06/2010
 *      Author: bro86j
 */

#ifndef MUTUAL_INFORMATION_H_
#define MUTUAL_INFORMATION_H_

//#include "NPWGLImage.h"
#include <cmath>
#include <vector>
#include <iostream>

namespace NPW
{

inline float getMutualInformationFromHistogram( const float * data, unsigned int bins )
{
    if ( bins < 1 )
    {
        std::cout << "getMutualInformationFromHistogramFaster() "
                  << "illegal number of bins : " << bins << std::endl;
        return -1.0f;
    }

    double * yMarginal = new double [bins];
    memset(yMarginal, 0, sizeof(double) * bins);

    double * xMarginal = new double [bins];
    memset(xMarginal, 0, sizeof(double) * bins);

    double Hx  = 0.0;
    double Hy  = 0.0;
    double Hxy = 0.0;

    // get marginals and combined entropy
    for ( unsigned int y = 0; y < bins; ++y )
    {
        int yOffset = y * bins;

        for (unsigned int x = 0; x < bins; ++x )
        {
            double val     = (double)data[ x + yOffset ];
            if ( val > 0.0 )
            {
                xMarginal[y] += val;
                yMarginal[x] += val;
                Hxy          -= val * std::log( val );
            }
        }
    }

    // get marginal entropies
    for ( unsigned int i = 0; i < bins; ++i )
    {
        double valX = xMarginal[i];
        if ( valX > 0.0f )
        {
            Hx    -= valX * std::log(valX);
        }
        double valY = yMarginal[i];
        if ( valY > 0.0f )
        {
            Hy    -= valY * std::log(valY);
        }
    }

    if (xMarginal) { delete [] xMarginal; xMarginal = 0; }
    if (yMarginal) { delete [] yMarginal; yMarginal = 0; }

    // return MI
    return ((float)(Hx + Hy - Hxy));
}

/*inline float getMutualInformationFromHistogram( const NPWGLImage & hist )
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
      if ( tmp > 0.0 )
      {
        Hab -= tmp * std::log(tmp);
      }
    }
  }

  // get the marginal entropies: H(X) = - SUM( p(x) * log(p(x)) )
  for ( int x = 0; x < dimX; ++x ) {
    double tmp = xMarginal[x];
    if ( tmp > 0.0 ) Ha -= tmp * std::log(tmp);
  }
  for ( int y = 0; y < dimY; ++y ) {
    double tmp = yMarginal[y];
    if ( tmp > 0.0 ) Hb -= tmp * std::log(tmp);
  }

  // Finally, get the MI
  return  ((float)(Ha + Hb - Hab));

}*/

} // end namespace NPW

#endif // MUTUAL_INFORMATION_H_
