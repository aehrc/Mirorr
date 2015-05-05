/*=========================================================================
Program: MILX MixView
Module: db3DNPWmontecarlo.h
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
*  3D NPW using monte carlo
*
*  AUTHOR: Daan Broekhuizen
*  DATE:   january 2010
*/

template< typename T > class dbPoint3D;

//namespace cimg = cimg_library;

/* adds a tetrahedron to the histogram
*  NOTE:
*  weight and NrOfSamples will be changed if we use adaptive monte carlo
*  therefor they cannot be const, but that doesn't matter, as they are
*  passed as a copy.
*/
void dbAddTetrahedron( boost::shared_ptr< cimg_library::CImg<double> > hist,
                       const dbPoint3D<float> &p0,
                       const dbPoint3D<float> &p1,
                       const dbPoint3D<float> &p2,
                       const dbPoint3D<float> &p3,
                       double weight,
                       const float min1, const float min2,
                       const float binsize1, const float binsize2,
                       int NrOfSamples );


/*-----------------------------------------------------------------------------------------------------------------------*/
/* fills a 2D joint histogram doing monte carlo 3D NPW on two images
*  set NrOfSamples to -1 for adaptive (on the volume o/t tetrahedron) samples
*/
void db3DNPWmontecarlo( boost::shared_ptr< cimg_library::CImg<double> > hist,
                        boost::shared_ptr< const cimg_library::CImg<float> > img1,
                        boost::shared_ptr< const cimg_library::CImg<float> > img2,
                        const int NrOfSamples,
                        bool verbosity);
