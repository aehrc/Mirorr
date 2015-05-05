/*=========================================================================
Program: MILX MixView
Module: utils.cxx
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
#include "utils.h"
#include "dbPoint3D.h"
#include <algorithm>

/*
*  BASIX: UTILS
*
*  handy functions
*
*  author: Daan Broekhuizen
*  date: january 2010
*/


void getBinsizes( boost::shared_ptr< const cimg::CImg<float> > img1,
                  boost::shared_ptr< const cimg::CImg<float> > img2,
                  const unsigned int NrOfBins1,
                  const unsigned int NrOfBins2,
                  double &binsize1,
                  double &binsize2,
                  float &min1, float &max1,
                  float &min2, float &max2 ) {

  min1 = (*img1)(0,0,0);
  max1 = min1;
  min2 = (*img2)(0,0,0);
  max2 = min2;

  for( int iZ = 0; iZ < img1->depth(); ++iZ ) {
    for( int iY = 0; iY < img1->height(); ++iY ) {
      for( int iX = 0; iX < img1->width(); ++iX ) {

        float val1 = (*img1)(iX,iY,iZ);
        float val2 = (*img2)(iX,iY,iZ);
        min1 = std::min( val1, min1 );
        max1 = std::max( val1, max1 );
        min2 = std::min( val2, min2 );
        max2 = std::max( val2, max2 );
      }
    }
  }

  binsize1 = static_cast<double>(max1 - min1) / static_cast<double>(NrOfBins1);
  binsize2 = static_cast<double>(max2 - min2) / static_cast<double>(NrOfBins2);

}

/*-----------------------------------------------------------------------------------------------------------------------*/

bool dbPointInTetrahedron( float xp, float yp, float zp,
                           float x1, float y1, float z1,
                           float x2, float y2, float z2,
                           float x3, float y3, float z3,
                           float x4, float y4, float z4 ) {
  // create a matrix
  cimg::CImg<float> m0( 4u,4u,1u,1u,   //size
                        x1,y1,z1,1.0,    //values
                        x2,y2,z2,1.0,
                        x3,y3,z3,1.0,
                        x4,y4,z4,1.0 );

  // get the sign of the determinant of the tetrahedron
  float D0 = (float) m0.det();
  int signD0 = dbSign( D0 );

  // if the sign is 0, the tetrahedron is degenerate
  if (signD0 == 0 ) return false;

  // now compare this to the sign of the other determinants,
  // if it is unequal, the point is not inside the tetrahedron
  cimg::CImg<float> m1( 4,4,1,1,
                        xp,yp,zp,1.0,
                        x2,y2,z2,1.0,
                        x3,y3,z3,1.0,
                        x4,y4,z4,1.0 );
  float D1 = (float) m1.det();
  if ( dbSign( D1 ) != signD0 ) return false;
  if ( dbSign( D1 ) == 0 ) return true;

  cimg::CImg<float> m2( 4,4,1,1,
                        x1,y1,z1,1.0,
                        xp,yp,zp,1.0,
                        x3,y3,z3,1.0,
                        x4,y4,z4,1.0 );
  float D2 = (float) m2.det();
  if ( dbSign( D2 ) != signD0 ) return false;
  if ( dbSign( D2 ) == 0 ) return true;

  cimg::CImg<float> m3( 4,4,1,1,
                        x1,y1,z1,1.0,
                        x2,y2,z2,1.0,
                        xp,yp,zp,1.0,
                        x4,y4,z4,1.0 );
  float D3 = (float) m3.det();
  if ( dbSign( D3 ) != signD0 ) return false;
  if ( dbSign( D3 ) == 0 ) return true;

  cimg::CImg<float> m4( 4,4,1,1,
                        x1,y1,z1,1.0,
                        x2,y2,z2,1.0,
                        x3,y3,z3,1.0,
                        xp,yp,zp,1.0 );
  float D4 = (float) m4.det();
  if ( dbSign( D4 ) != signD0 ) return false;
//  if ( dbSign( D4 ) == 0 ) return true;       // this check is not necessary, true is returned anyway

  // if all signs are equal, the point is within the tetrahedron
  return true;

}

/*-----------------------------------------------------------------------------------------------------------------------*/

float getTetrahedronVolume( const dbPoint3D<float> &p0,
                            const dbPoint3D<float> &p1,
                            const dbPoint3D<float> &p2,
                            const dbPoint3D<float> &p3 ) {
  // volume = ( dot(d-a, cross(c-a,b-a) ) ) / 6;
  dbPoint3D<float> tmp3 = p3 - p0;
  dbPoint3D<float> tmp2 = p2 - p0;
  dbPoint3D<float> tmp1 = p1 - p0;
  dbPoint3D<float> tmp4 = tmp2.cross(tmp1);

  return abs( tmp3.dot( tmp4 ) / 6.0f );
}

/*-----------------------------------------------------------------------------------------------------------------------*/
