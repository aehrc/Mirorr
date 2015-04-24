/*
 * MutualInformation.h
 *
 *  Created on: 16/06/2010
 *      Author: bro86j
 */

#ifndef MUTUAL_INFORMATION_H_
#define MUTUAL_INFORMATION_H_

#include <cmath>
#include <vector>
#include <iostream>

namespace NPW
{

inline float getMutualInformationFromHistogram( const float * data, int bins )
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
    for ( int y = 0; y < bins; ++y )
    {
        int yOffset = y * bins;

        for (int x = 0; x < bins; ++x )
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
    for ( int i = 0; i < bins; ++i )
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

} // end namespace NPW

#endif // MUTUAL_INFORMATION_H_
