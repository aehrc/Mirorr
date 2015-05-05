/*=========================================================================
Program: MILX MixView
Module: utils.h
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
#define cimg_display 0
#include "CImg.h"
#include <boost/shared_ptr.hpp>

/*
*  BASIX: UTILS
*
*  handy functions
*
*  author: Daan Broekhuizen
*  date: january 2010
*/

//#include "dbPoint3D.h"
template<typename T> class dbPoint3D;

namespace cimg = cimg_library;

// determines the sign of a float ( 0 if the value is 0.0 )
inline int dbSign( float x ) {
  return ( x > 0.0 ) - ( x < 0.0 );
}

// gets the binsizes of the joint histogram of two images
void getBinsizes( boost::shared_ptr< const cimg::CImg<float> > img1,
                  boost::shared_ptr< const cimg::CImg<float> > img2,
                  const unsigned int NrOfBins1,
                  const unsigned int NrOfBins2,
                  double &binsize1,
                  double &binsize2,
                  float &min1, float &max1,
                  float &min2, float &max2 );


/* determines whether a point is within a specified tetrahedron
*   OR on its edge!
*/
bool dbPointInTetrahedron( float xp, float yp, float zp,
                           float x1, float y1, float z1,
                           float x2, float y2, float z2,
                           float x3, float y3, float z3,
                           float x4, float y4, float z4 );

/* calculates the volume of a tetrahedron
*/
float getTetrahedronVolume( const dbPoint3D<float> &p0,
                            const dbPoint3D<float> &p1,
                            const dbPoint3D<float> &p2,
                            const dbPoint3D<float> &p3 );
