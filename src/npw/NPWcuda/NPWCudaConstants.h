/*
 * NPWConstants.h
 *
 *  Created on: 19/07/2010
 *      Author: bro86j
 */

#ifndef NPW_CONSTANTS_H_
#define NPW_CONSTANTS_H_

const unsigned int NPW_CUDA_N_CONVEX_HULL_VERTICES   = 19;

const unsigned int FLOATS_PER_GEO_UNIT               = 12;
const unsigned int FLOAT2S_PER_GEO_UNIT              =  6;
const unsigned int FLOAT4S_PER_GEO_UNIT              =  3;
const unsigned int NEIGHBORHOODS_PER_BLOCK_DIMENSION =  8;
const unsigned int VOXELS_PER_NEIGHBORHOOD           =  8;

//#define      FLOATS_PER_GEO_UNIT                12
//#define      FLOAT2S_PER_GEO_UNIT                6
//#define      FLOAT4S_PER_GEO_UNIT                3
//
//#define      NEIGHBORHOODS_PER_BLOCK_DIMENSION   8
//#define      VOXELS_PER_NEIGHBORHOOD             8


const float  NPW_CUDA_EPSILON                         =  1e-6f;
const float  NPW_CUDA_NEGATIVE_EPSILON                =  -NPW_CUDA_EPSILON;
const float  NPW_CUDA_ONE_PLUS_EPSILON                =  1.0f + NPW_CUDA_EPSILON;

const float  MINIMUM_LINE_LENGTH                      =   2.0f;
const float  MINIMUM_LINE_LENGTH_SQUARED              =   MINIMUM_LINE_LENGTH * MINIMUM_LINE_LENGTH;

const float  MINIMUM_AREA                             =   0.0f;
const float  MINIMUM_AREA_TIMES_TWO                   =   2.0f * MINIMUM_AREA;

const float  MINIMUM_EDGE_DISTANCE                    =   2.0f;
const float  MINIMUM_EDGE_DISTANCE_SQUARED            =   MINIMUM_EDGE_DISTANCE * MINIMUM_EDGE_DISTANCE;
const float  MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO  =   2.83f;
const float  MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO  =   2.0f * MINIMUM_EDGE_DISTANCE_SQUARED;
const float  MINIMUM_EDGE_DISTANCE_TIMES_TWO          =   2.0f * MINIMUM_EDGE_DISTANCE;
const float  MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_FOUR =   4.0f * MINIMUM_EDGE_DISTANCE_SQUARED;

const float  VERTEX_TRIO_CORRECTION_DISTANCE          =   3.63f;
      // = MINIMUM_EDGE_DISTANCE / sin( pi/4 - arcsin( MINIMUM_EDGE_DISTANCE / MINIMUM_LINE_LENGTH ) )

const float  MINIMUM_POINT_DISTANCE_SQUARED           =   1.0f;
const float  MINIMUM_POINT_TO_LINE_DISTANCE_SQUARED   =   0.5f * MINIMUM_POINT_DISTANCE_SQUARED;

// todo: a bit high?
const float  MAXIMUM_VERTEX_MOVEMENT_SQUARED          = 50.0f;

const float  MAXIMUM_VERTEX_CORRECTION                =   4.95f;
      // = 1.5 * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO + MINIMUM_POINT_DISTANCE / sqrt(2);

const float  SINUS_OF_120_DEGREES                     =  0.866025403784438646763723f;
const float  COSINUS_OF_120_DEGREES                   = -0.5f;
const float  SINUS_OF_240_DEGREES                     = -0.866025403784438646763723f;
const float  COSINUS_OF_240_DEGREES                   = -0.5f;

enum Shape
{
    QUAD = 0,
    TRIANGLE_WITH_THREE_SUBTRIANGLES = 1,
    TRIANGLE_WITH_ONE_VERTEX_ON_BORDER = 2,
    TRIANGLE_WITH_ONE_VERTEX_PAIR = 3,
    LINE_WITHOUT_VERTEX_PAIR = 4,
    LINE_WITH_VERTEX_PAIR_ON_END = 5,
    LINE_WITH_VERTEX_PAIR_ON_LINE = 6,
    LINE_WITH_TWO_VERTEX_PAIRS = 7,
    LINE_WITH_VERTEX_TRIO = 8,
    TRUE_POINT = 9,
    POINT_EQUIVALENT = 10,
    UNDEFINED = 11
};


#endif // NPW_CONSTANTS_H_
