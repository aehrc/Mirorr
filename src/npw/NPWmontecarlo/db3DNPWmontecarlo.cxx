/*=========================================================================
Program: MILX MixView
Module: db3DNPWmontecarlo.cxx
Author: Daan Broekhuizen
Modified by:
Language: C++
Created: Tue 25 May 2010 12:23:13 EST

Copyright: (c) 2009 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/
#include "db3DNPWmontecarlo.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "utils.h"
#include "sampling.h"

/*
*  3D NPW using monte carlo
*
*  AUTHOR: Daan Broekhuizen
*  DATE:   january 2010
*/


const unsigned int kMinimumNrOfSamples = 3;
const unsigned int kMaximumNrOfSamples = 1000;
const double kNrOfSamplesMultiplicationFactor = 5.0;

/* Adds a tetrahedron to the histogram
*  NOTE: maybe slow due to the necessary passing of all the variables?
*/
void dbAddTetrahedron( boost::shared_ptr< cimg::CImg<double> > hist,
                       const dbPoint3D<float> &p0,
                       const dbPoint3D<float> &p1,
                       const dbPoint3D<float> &p2,
                       const dbPoint3D<float> &p3,
                       double weight,
                       const float min1, const float min2,
                       const double binsize1, const double binsize2,
                       int NrOfSamples) {
  // Check if the tetrahedron is singlebinned
  int x = (int) floor( ( p0.getX() - min1 ) / binsize1 );
  int y = (int) floor( ( p0.getY() - min2 ) / binsize2 );
  if ( x == (int) floor( ( p1.getX() - min1 ) / binsize1 ) &&
       x == (int) floor( ( p2.getX() - min1 ) / binsize1 ) &&
       x == (int) floor( ( p3.getX() - min1 ) / binsize1 ) &&
       y == (int) floor( ( p1.getY() - min2 ) / binsize2 ) &&
       y == (int) floor( ( p2.getY() - min2 ) / binsize2 ) &&
       y == (int) floor( ( p3.getY() - min2 ) / binsize2 ) ) {
    // check boundaries to prevent access violations
    x = std::min( x, hist->width()-1 );
    y = std::min( y, hist->height()-1 );
    x = std::max( x, 0 );
    y = std::max( y, 0 );

    // add the weight to the correct bin
    (*hist)( x, y ) += weight * static_cast<double>(NrOfSamples);
    return;
  }

  // if NrOfSample == -1, then we adapt the Nr of Samples
  // on the volume of the tetrahedron
  if (NrOfSamples == -1) {
    // base Nr of samples on number of bins in bounding box
    // get min/max x and min/max y
    int minx = (int) (floor(p0.getX() - min1) / binsize1);
    int maxx = minx;
    int miny = (int) (floor(p0.getY() - min2) / binsize2);
    int maxy = miny;

    if( (int) (floor( p1.getX() - min1 ) / binsize1) < minx ) minx = (int) (floor( p1.getX() - min1 ) / binsize1);
    if( (int) (floor( p1.getX() - min1 ) / binsize1) > maxx ) maxx = (int) (floor( p1.getX() - min1 ) / binsize1);
    if( (int) (floor( p1.getY() - min2 ) / binsize2) < miny ) miny = (int) (floor( p1.getY() - min2 ) / binsize2);
    if( (int) (floor( p1.getY() - min2 ) / binsize2) > maxy ) maxy = (int) (floor( p1.getY() - min2 ) / binsize2);
    if( (int) (floor( p2.getX() - min1 ) / binsize1) < minx ) minx = (int) (floor( p2.getX() - min1 ) / binsize1);
    if( (int) (floor( p2.getX() - min1 ) / binsize1) > maxx ) maxx = (int) (floor( p2.getX() - min1 ) / binsize1);
    if( (int) (floor( p2.getY() - min2 ) / binsize2) < miny ) miny = (int) (floor( p2.getY() - min2 ) / binsize2);
    if( (int) (floor( p2.getY() - min2 ) / binsize2) > maxy ) maxy = (int) (floor( p2.getY() - min2 ) / binsize2);
    if( (int) (floor( p3.getX() - min1 ) / binsize1) < minx ) minx = (int) (floor( p3.getX() - min1 ) / binsize1);
    if( (int) (floor( p3.getX() - min1 ) / binsize1) > maxx ) maxx = (int) (floor( p3.getX() - min1 ) / binsize1);
    if( (int) (floor( p3.getY() - min2 ) / binsize2) < miny ) miny = (int) (floor( p3.getY() - min2 ) / binsize2);
    if( (int) (floor( p3.getY() - min2 ) / binsize2) > maxy ) maxy = (int) (floor( p3.getY() - min2 ) / binsize2);

    int NrOfBinsInBox = (maxx - minx) * (maxy - miny);
    NrOfSamples = NrOfBinsInBox * 5;


    // don't go over maximum or under minimum
    NrOfSamples = std::min( static_cast<unsigned int>(NrOfSamples), kMaximumNrOfSamples );
    NrOfSamples = std::max( static_cast<unsigned int>(NrOfSamples), kMinimumNrOfSamples );


    // weight is calculated in db3DNPWmontecarlo as:
    //    tetrahedronfraction / ( double(NrOfSamples) * NrOfNeighborhoods )
    //    and NrOfSamples was -1 in that case, therefore we divide by -NrOfSample here
    weight /= - static_cast<double>(NrOfSamples);
  }

  // now sample the tetrahedron
  for ( int i = 0; i < NrOfSamples; ++i ) {
    // get a point in the tetrahedron
    const dbPoint3D<float> samplePoint = dbSamplePointInTetrahedron(p0,p1,p2,p3);
    const float sampleX = samplePoint.getX();
    const float sampleY = samplePoint.getY();

    // convert it to histogram coordinates
    int x = (int) floor( (sampleX - min1) / binsize1 );
    int y = (int) floor( (sampleY - min2) / binsize2 );

    // check boundaries to prevent access violations
    x = std::min( x, hist->width()-1 );
    y = std::min( y, hist->height()-1 );
    x = std::max( x, 0 );
    y = std::max( y, 0 );

    // add the weight to the correct bin
    (*hist)( x, y ) += weight;
  }

}

/*-----------------------------------------------------------------------------------------------------------------------*/

void db3DNPWmontecarlo( boost::shared_ptr< cimg::CImg<double> > hist,
                        boost::shared_ptr< const cimg::CImg<float> > img1,
                        boost::shared_ptr< const cimg::CImg<float> > img2,
                        const int NrOfSamples,
                        bool verbosity ) {

  // check dimensions
  if( img1->width() != img2->width() || img1->height() != img2->height() || img1->depth() != img2->depth() ) {
    if (verbosity) {
      std::cerr << "  image 1: " << img1->width() << " x " << img1->height() << " x " << img1->depth() << std::endl;
      std::cerr << "  image 2: " << img2->width() << " x " << img2->height() << " x " << img2->depth() << std::endl;
    }
    throw std::invalid_argument("db3DNPWmontecarlo Error: images not of equal size!\n");
  }
  if( img1->width() < 2 || img1->height() < 2 || img1->depth() < 2 ) {
    throw std::invalid_argument("db3DNPWmontecarlo Error: images must be at least 2x2x2.\n");
  }
  if( 0 == hist->width() || 0 == hist->height() || 0 == hist->depth() ) {
    throw std::invalid_argument("db3DNPWmontecarlo Error: histogram has a dimension of zero.\n");
  }

  // get the min max of both images in one go, since they should be of equal size
  double binsize1,binsize2;
  float min1, max1, min2, max2;
  getBinsizes( img1, img2, hist->width(), hist->height(), binsize1, binsize2, min1, max1, min2, max2);

  // get the number of neighborhoods and the weight for each sample
  const unsigned int NrOfNeighborhoods = ( img1->width() - 1 ) * ( img1->height() - 1 ) * ( img1->depth() - 1 );
  const double innerSampleWeight = 1.0 / static_cast<double>( 3 * NrOfSamples * NrOfNeighborhoods );
  const double outerSampleWeight = 1.0 / static_cast<double>( 6 * NrOfSamples * NrOfNeighborhoods );
  dbPoint3D<float> alpha, beta, gamma, delta, epsilon, dzeta, eta, theta;

  // for each neighborhood...
  for( int iZ = 0; iZ < img1->depth() - 1; ++iZ ) {
    for( int iY = 0; iY < img1->height() - 1; ++iY ) {
      for( int iX = 0; iX < img1->width() - 1; ++iX ) {
        float alpha1   = (*img1)(iX  ,iY  ,iZ  );
        float beta1    = (*img1)(iX+1,iY  ,iZ  );
        float gamma1   = (*img1)(iX  ,iY+1,iZ  );
        float delta1   = (*img1)(iX+1,iY+1,iZ  );
        float epsilon1 = (*img1)(iX  ,iY  ,iZ+1);
        float dzeta1   = (*img1)(iX+1,iY  ,iZ+1);
        float eta1     = (*img1)(iX  ,iY+1,iZ+1);
        float theta1   = (*img1)(iX+1,iY+1,iZ+1);

        float alpha2   = (*img2)(iX  ,iY  ,iZ  );
        float beta2    = (*img2)(iX+1,iY  ,iZ  );
        float gamma2   = (*img2)(iX  ,iY+1,iZ  );
        float delta2   = (*img2)(iX+1,iY+1,iZ  );
        float epsilon2 = (*img2)(iX  ,iY  ,iZ+1);
        float dzeta2   = (*img2)(iX+1,iY  ,iZ+1);
        float eta2     = (*img2)(iX  ,iY+1,iZ+1);
        float theta2   = (*img2)(iX+1,iY+1,iZ+1);

        alpha.set(alpha1, alpha2, 0);
        beta.set(beta1, beta2, 0);
        gamma.set(gamma1, gamma2, 0);
        delta.set(delta1, delta2, 0);
        epsilon.set(epsilon1, epsilon2, 1);
        dzeta.set(dzeta1, dzeta2, 1);
        eta.set(eta1, eta2, 1);
        theta.set(theta1, theta2, 1);

        // tetrahedron 1 (containing 0,0,0) :
        dbAddTetrahedron(hist, alpha, beta, gamma, epsilon, outerSampleWeight,
                         min1, min2, binsize1, binsize2, NrOfSamples);

        // tetrahedron 2 (containing 1,1,0) :
        dbAddTetrahedron(hist, beta, gamma, delta, theta, outerSampleWeight,
                         min1, min2, binsize1, binsize2, NrOfSamples);

        // tetrahedron 3 (containing 1,0,1) :
        dbAddTetrahedron(hist, beta, epsilon, dzeta, theta, outerSampleWeight,
                         min1, min2, binsize1, binsize2, NrOfSamples);

        // tetrahedron 4 (containing 0,1,1) :
        dbAddTetrahedron(hist, gamma, epsilon, eta, theta, outerSampleWeight,
                         min1, min2, binsize1, binsize2, NrOfSamples);

        // tetrahedron 5 (the inner one):
        dbAddTetrahedron(hist, beta, gamma, epsilon, theta, innerSampleWeight,
                         min1, min2, binsize1, binsize2, NrOfSamples);
      }
    }
  }

  if ( verbosity ) {
    unsigned int NrOfTetrahedra = NrOfNeighborhoods * 5;
    std::cerr << "  NPW: succesfully filled a " << hist->width() << " x " << hist->width() << " histogram "
              << "with " << NrOfTetrahedra << " tetrahedra" << std::endl;
  }


}
