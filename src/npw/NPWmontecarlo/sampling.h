/*=========================================================================
Program: MILX MixView
Module: sampling.h
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
#include "dbPoint3D.h"
//template< typename T > class dbPoint3D;

/*
*  BASIX: SAMPLING
*   contains basic sampling functions
*
*  author: Daan Broekhuizen
*  date: january 2010
*/

/* returns a random float value between 0.0 and 1.0
*  NOTE:random 'enough'? Maybe some more sophisticated random number generator?
*/
inline float frand() {
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  //return ( float(rand()) / float(RAND_MAX));
}

/* Randomly samples a point in a tetrahedron
*
*  uses uniform sampling by folding the unit cube as described by Rocchini & Cignoni:
*     Generating Random Points in a Tetrahedron, 2001
*     http://jgt.akpeters.com/papers/RocchiniCignoni00/
*/
template<typename T>
dbPoint3D<T> dbSamplePointInTetrahedron( const dbPoint3D<T> &p0,
                                         const dbPoint3D<T> &p1,
                                         const dbPoint3D<T> &p2,
                                         const dbPoint3D<T> &p3 ) {
  T s = static_cast<T>(frand());
  T t = static_cast<T>(frand());
  T u = static_cast<T>(frand());
  if ( s + t > 1.0 ) { // cut'n fold the cube into a prism
    s = (T) (1.0) - s;
    t = (T) (1.0) - t;
  }
  if ( t + u > 1.0 ) { // cut'n fold the prism into a tetrahedron
    T tmp = u;
    u = (T) (1.0) - s - t;
    t = (T) (1.0) - tmp;
  } else if ( s + t + u > (T) (1.0) ) {
    T tmp = u;
    u = s + t + u - (T) (1.0);
    s = (T) (1.0) - t - tmp;
  }

  T a = 1 - s - t - u; // a,s,t,u are the barycentric coordinates of the random point.

  dbPoint3D<T> tmp = p0*a + p1*s + p2*t + p3*u;
  return tmp;
}
