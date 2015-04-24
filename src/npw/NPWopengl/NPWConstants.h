/*
 * NPWConstants.h
 *
 *  Created on: 23/06/2010
 *      Author: bro86j
 */

#ifndef NPW_CONSTANTS_H_
#define NPW_CONSTANTS_H_

namespace npw
{

const unsigned int    nV               = 4; // number of vertices
const unsigned int    nMaxRenderPoints = 8;
        // even though this is actually 6, an array of 8 is faster than 6
        // http://www.eventhelix.com/realtimemantra/basics/optimizingcandcppcode.htm

const double EPSILON = 0.0001;

const double MINIMUM_AREA                             =   0.1;

const float  MINIMUM_LINE_LENGTH                      =   3.0f;
const float  MINIMUM_LINE_LENGTH_SQUARED              =   9.0f;

const float  MINIMUM_EDGE_DISTANCE                    =   2.0f;
const float  MINIMUM_EDGE_DISTANCE_SQUARED            =   4.0f;
const float  MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO  =   2.83f;
const float  MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO  =   8.0f;
const float  MINIMUM_EDGE_DISTANCE_TIMES_TWO          =   4.0f;
const float  MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_FOUR =  16.0f;

const float  VERTEX_TRIO_CORRECTION_DISTANCE          =   3.63f;
      // = MINIMUM_EDGE_DISTANCE / sin( pi/4 - arcsin( MINIMUM_EDGE_DISTANCE / MINIMUM_LINE_LENGTH ) )

const float  MINIMUM_POINT_DISTANCE_SQUARED           =   2.0f;
const float  MINIMUM_POINT_TO_LINE_DISTANCE_SQUARED   =   1.0f;
      // = 0.5 * MINIMUM_POINT_DISTANCE_SQUARED

// todo: a bit high?
const float  MAXIMUM_VERTEX_MOVEMENT_SQUARED          = 50.0f;

const float  MAXIMUM_VERTEX_CORRECTION                =   4.95f;
      // = 1.5 * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO + MINIMUM_POINT_DISTANCE / sqrt(2);

const float  SINUS_OF_120_DEGREES                     =  0.866025403784438646763723f;
const float  COSINUS_OF_120_DEGREES                   = -0.5f;
const float  SINUS_OF_240_DEGREES                     = -0.866025403784438646763723f;
const float  COSINUS_OF_240_DEGREES                   = -0.5f;

}

#endif // NPW_CONSTANTS_H_
