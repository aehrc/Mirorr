/*
 * NPWCudaShapeExpansion.h
 *
 *  Created on: 19/07/2010
 *      Author: bro86j
 */

#ifndef NPW_CUDA_SHAPE_EXPANSION_H_
#define NPW_CUDA_SHAPE_EXPANSION_H_

#include "NPWCudaGeometryFunctions.h"
#include "NPWCudaConstants.h"

__device__ void correctLineWithVertexTrio(
        float2 * vertices, unsigned int n_Vertices,
        unsigned int * n_Duplicates_at_Vertices,
        float * geometryList )
{
    /*
     * determine what vertex is the trio
     */
    float2 move         = make_float2(0.0f,0.0f);
    float2 singleVertex = vertices[0];
    float2 vertexTrio   = vertices[1];
    if ( n_Duplicates_at_Vertices[0] == 2)
    {
        singleVertex  = vertices[1];
        vertexTrio    = vertices[0];
    }

    /*
     * move two vertices of the trio away from the line at 120 and 240 degrees
     */
    move        = Normalize( singleVertex - vertexTrio ) * MINIMUM_EDGE_DISTANCE_TIMES_TWO;
    move        = make_float2( move.x * COSINUS_OF_120_DEGREES - move.y * SINUS_OF_120_DEGREES,
                               move.x * SINUS_OF_120_DEGREES   + move.y * COSINUS_OF_120_DEGREES);
    vertices[2] = vertexTrio + move;

    move        = Normalize( singleVertex - vertexTrio) * MINIMUM_EDGE_DISTANCE_TIMES_TWO;
    move        = make_float2( move.x*COSINUS_OF_240_DEGREES - move.y*SINUS_OF_240_DEGREES,
                               move.x*SINUS_OF_240_DEGREES   + move.y*COSINUS_OF_240_DEGREES);
    vertices[3] = vertexTrio + move;

    /*
     * set geometry list
     */
    geometryList[2] = vertexTrio.x;
    geometryList[3] = vertexTrio.y;
    geometryList[4] = singleVertex.x;
    geometryList[5] = singleVertex.y;
    geometryList[6] = vertices[2].x;
    geometryList[7] = vertices[2].y;
    geometryList[8] = vertices[3].x;
    geometryList[9] = vertices[3].y;

    geometryList[0] = 3.0f;
}


__device__ void correctLineWithTwoVertexPairs(
        float2 * vertices, unsigned int n_Vertices,
        unsigned int * n_Duplicates_at_Vertices,
        float * geometryList )
{
    /*
     * calculate the intersect
     */
    float2 intersect = ( vertices[0] + vertices[1] ) * 0.5f;
    float2 move      = make_float2(0.0f,0.0f);

    /*
     * move pairs
     */
    move.x = vertices[0].y - vertices[1].y;
    move.y = vertices[1].x - vertices[0].x;
    move   = Normalize(move) * MINIMUM_EDGE_DISTANCE;

    vertices[2]  = vertices[1] + move;
    vertices[3]  = vertices[0] + move;
    vertices[0] -= move;
    vertices[1] -= move;

    /*
     * set render list
     */
    geometryList[ 2] = intersect.x;
    geometryList[ 3] = intersect.y;
    geometryList[ 4] = vertices[0].x;
    geometryList[ 5] = vertices[0].y;
    geometryList[ 6] = vertices[1].x;
    geometryList[ 7] = vertices[1].y;
    geometryList[ 8] = vertices[2].x;
    geometryList[ 9] = vertices[2].y;
    geometryList[10] = vertices[3].x;
    geometryList[11] = vertices[3].y;

    geometryList[0] = 4.0f;
}


__device__ void correctLineWithVertexPairOnEnd(
        float2 * vertices, unsigned int n_Vertices,
        unsigned int * n_Duplicates_at_Vertices,
        float * geometryList )
{
    /*
     * calculate the intersect and determine which is the pair end
     */
    float2 intersect = make_float2( 0.0f,0.0f );
    float2 move      = make_float2( 0.0f,0.0f );
    float2 pairEnd   = make_float2( vertices[0].x, vertices[0].y );
    float2 singleEnd = make_float2( vertices[2].x, vertices[2].y );

    // if pair is at vertices[0]
    if ( n_Duplicates_at_Vertices[0] == 1 )
    {
        if ( Abs( vertices[0].x - vertices[2].x ) >
             Abs( vertices[0].y - vertices[2].y ) )
        {
            intersect.x  = ( ( vertices[0].x * vertices[0].x ) -
                             ( vertices[1].x * vertices[2].x ) ) /
                           ( vertices[0].x * 2.0f - vertices[2].x  - vertices[1].x );
            float y_grad = ( vertices[2].y - vertices[0].y ) /
                           ( vertices[2].x - vertices[0].x );
            intersect.y  = ( intersect.x - vertices[0].x ) * y_grad + vertices[0].y;
        }
        else
        {
            intersect.y  = ( ( vertices[0].y * vertices[0].y ) -
                             ( vertices[1].y * vertices[2].y ) ) /
                           ( vertices[0].y * 2.0f - vertices[2].y - vertices[1].y );
            float x_grad = ( vertices[2].x - vertices[0].x ) /
                           ( vertices[2].y - vertices[0].y );
            intersect.x  = ( intersect.y - vertices[0].y ) * x_grad + vertices[0].x;
        }
    }
    // the pair is at vertices[2]
    else
    {
        pairEnd    = vertices[2];
        singleEnd  = vertices[0];

        if ( Abs( vertices[0].x - vertices[2].x ) >
             Abs( vertices[0].y - vertices[2].y ) )
        {
            intersect.x  = ( ( vertices[0].x * vertices[1].x ) -
                             ( vertices[2].x * vertices[2].x ) ) /
                           ( vertices[0].x  + vertices[1].x - vertices[2].x * 2.0f );
            float y_grad = ( vertices[2].y - vertices[0].y ) /
                           ( vertices[2].x - vertices[0].x );
            intersect.y  = ( intersect.x - vertices[0].x ) * y_grad + vertices[0].y;
        }
        else
        {
            intersect.y  = ( ( vertices[0].y * vertices[1].y ) -
                             ( vertices[2].y * vertices[2].y ) ) /
                           ( vertices[0].y + vertices[1].y - vertices[2].y * 2.0f );
            float x_grad = ( vertices[2].x - vertices[0].x ) /
                           ( vertices[2].y - vertices[0].y );
            intersect.x  = ( intersect.y - vertices[0].y ) * x_grad + vertices[0].x;
        }
    }

    /*
     * ensure the intersect is MINIMUM_EDGE_DISTANCE from the pair end
     */
    float distanceSquared = PointPointDistanceSquared( intersect, pairEnd );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED )
    {
        move  = Normalize( pairEnd - singleEnd );
        move *= ( MINIMUM_EDGE_DISTANCE - sqrtf(distanceSquared) );

        pairEnd += move;
    }

    /*
     * ensure the intersect is MINIMUM_EDGE_DISTANCE_TIMES_TWO from the single end
     */
    distanceSquared = PointPointDistanceSquared( intersect, singleEnd );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_FOUR )
    {
        move  = Normalize( singleEnd - pairEnd );
        move *= (MINIMUM_EDGE_DISTANCE_TIMES_TWO - sqrtf(distanceSquared) );

        singleEnd += move;
    }

    /*
     * move the pair perpendicular to the line with MINIMUM_EDGE_DISTANCE_TIMES_TWO
     */
    move.x = singleEnd.y - pairEnd.y;
    move.y = pairEnd.x   - singleEnd.x;
    move   = Normalize( move ) * MINIMUM_EDGE_DISTANCE_TIMES_TWO;

    vertices[3]  = pairEnd - move;
    pairEnd     += move;

    /*
     * set geometry list
     */
    geometryList[2] = intersect.x;
    geometryList[3] = intersect.y;
    geometryList[4] = pairEnd.x;
    geometryList[5] = pairEnd.y;
    geometryList[6] = vertices[3].x;
    geometryList[7] = vertices[3].y;
    geometryList[8] = singleEnd.x;
    geometryList[9] = singleEnd.y;

    geometryList[0] = 3.0f;
}


__device__ void correctLineWithVertexPairOnLine(
        float2 * vertices, unsigned int n_Vertices,
        unsigned int * n_Duplicates_at_Vertices,
        float * geometryList )
{
    /*
     * calculate the intersect
     */
    float2 move      = make_float2( 0.0f,0.0f );
    float2 intersect = make_float2( vertices[0].x, vertices[0].y );

    /*
      * check if the end points are a minimum distance of
      * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_2 from the pair
      */
     if ( PointPointDistanceSquared( vertices[0], vertices[1] )
             < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
     {
         move         = Normalize( vertices[0] - vertices[2] ) * MINIMUM_EDGE_DISTANCE;
         vertices[0] += move;
     }
     if ( PointPointDistanceSquared( vertices[1], vertices[2] )
             < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
     {
         move         = Normalize(vertices[2] - vertices[0]) * MINIMUM_EDGE_DISTANCE;
         vertices[2] += move;
     }

     /*
      * move vertex pair
      */
     move.x = vertices[0].y - vertices[2].y;
     move.y = vertices[2].x - vertices[0].x;
     move   = Normalize(move) * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO;

     vertices[3]  = vertices[1] + move;
     vertices[1] -= move;


     /*
      * set render list
      */
     geometryList[ 2] = intersect.x;
     geometryList[ 3] = intersect.y;
     geometryList[ 4] = vertices[1].x;
     geometryList[ 5] = vertices[1].y;
     geometryList[ 6] = vertices[0].x;
     geometryList[ 7] = vertices[0].y;
     geometryList[ 8] = vertices[3].x;
     geometryList[ 9] = vertices[3].y;
     geometryList[10] = vertices[2].x;
     geometryList[11] = vertices[2].y;

     geometryList[0] = 4.0f;
}

__device__ void correctLineWithoutVertexPair(
        float2 * vertices, unsigned int n_Vertices,
        unsigned int * n_Duplicates_at_Vertices,
        float * geometryList )
{
    float2 move      = make_float2(0.0f, 0.0f);
    float2 intersect = make_float2(0.0f,0.0f);

    /*
     * calculate the intersect
     */
    if ( Abs( vertices[0].x - vertices[3].x ) >
         Abs( vertices[0].y - vertices[3].y ) )
    {
        intersect.x  = ( ( vertices[0].x * vertices[1].x ) -
                         ( vertices[2].x * vertices[3].x ) ) /
                       ( vertices[0].x - vertices[3].x +
                         vertices[1].x - vertices[2].x );
        float y_grad = ( vertices[3].y - vertices[0].y ) /
                       ( vertices[3].x - vertices[0].x );
        intersect.y  = ( intersect.x - vertices[0].x ) * y_grad + vertices[0].y;
    }
    else
    {
        intersect.y  = ( ( vertices[0].y * vertices[1].y ) -
                         ( vertices[2].y * vertices[3].y ) ) /
                       ( vertices[0].y - vertices[3].y +
                         vertices[1].y - vertices[2].y );
        float x_grad = ( vertices[3].x - vertices[0].x ) /
                       ( vertices[3].y - vertices[0].y );
        intersect.x  = ( intersect.y - vertices[0].y ) * x_grad + vertices[0].x;
    }

    /*
     * ensure the two vertices on the line are at least MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_2
     * from their closest line end vertex
     */
    float distanceSquared = PointPointDistanceSquared( vertices[0], vertices[1] );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
    {
        move  = Normalize( vertices[0] - vertices[3] );
        move *= ( MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO - sqrtf(distanceSquared) );
        vertices[0] += move;
    }

    distanceSquared = PointPointDistanceSquared( vertices[3], vertices[2] );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
    {
        move  = Normalize( vertices[3] - vertices[0] );
        move *= ( MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO - sqrtf(distanceSquared) );
        vertices[3] += move;
    }

    move.x = vertices[0].y - vertices[3].y;
    move.y = vertices[3].x - vertices[0].x;
    move   = Normalize( move ) * MAXIMUM_VERTEX_CORRECTION;

    vertices[1] += move;
    vertices[2] -= move;

    /*
     * set render list and shape : 1->0->2->3
     */
    geometryList[ 2] = intersect.x;
    geometryList[ 3] = intersect.y;
    geometryList[ 4] = vertices[1].x;
    geometryList[ 5] = vertices[1].y;
    geometryList[ 6] = vertices[0].x;
    geometryList[ 7] = vertices[0].y;
    geometryList[ 8] = vertices[2].x;
    geometryList[ 9] = vertices[2].y;
    geometryList[10] = vertices[3].x;
    geometryList[11] = vertices[3].y;

    geometryList[0] = 4.0f;
}

/*
 * expands a triangle abc and stores the result in def
 */
__device__ bool expandTriangle( float2 a, float2 b, float2 c,
                                float & dx, float & dy,
                                float & ex, float & ey,
                                float & fx, float & fy,
                                unsigned int bins )
{
    float2 move = make_float2(0.0f,0.0f);

    /*
     * extend line ab
     */
    move.x      = a.y - b.y;
    move.y      = b.x - a.x;
    move        = Normalize( move ) * MINIMUM_EDGE_DISTANCE;
    float2 ab1  = a + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( PointsOnSameSideOfLine( a, b, ab1, c ) )
    {
        move    = -move;
        ab1     = a + move;
    }
    float2 ab2  = b + move;

    /*
     * extend line bc
     */
    move.x      = b.y - c.y;
    move.y      = c.x - b.x;
    move        = Normalize( move ) * MINIMUM_EDGE_DISTANCE;
    float2 bc1  = b + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( PointsOnSameSideOfLine( b, c, bc1, a) )
    {
        move    = -move;
        bc1     = b + move;
    }
    float2 bc2  = c + move;

    /*
     * extend line ca
     */
    move.x      = c.y - a.y;
    move.y      = a.x - c.x;
    move        = Normalize( move ) * MINIMUM_EDGE_DISTANCE;
    float2 ca1  = c + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( PointsOnSameSideOfLine( c, a, ca1, b ) )
    {
        move    = -move;
        ca1     = c + move;
    }
    float2 ca2  = a + move;

    /*
     * get the three intersections of the lines
     */
    float2 a_new = make_float2(0.0f,0.0f);
    float2 b_new = make_float2(0.0f,0.0f);
    float2 c_new = make_float2(0.0f,0.0f);

    if ( !LineSegmentLineSegmentIntersect( ab1,ab2, ca1,ca2, a_new ) )
    {
        return false;
    }
    if ( !LineSegmentLineSegmentIntersect( bc1,bc2, ab1,ab2, b_new ) )
    {
        return false;
    }
    if ( !LineSegmentLineSegmentIntersect( ca1,ca2, bc1,bc2, c_new ) )
    {
        return false;
    }

    /*
     * check if any point is moved over too large a distance
     */
    if ( PointPointDistanceSquared( a_new, a) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         PointPointDistanceSquared( b_new, b) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         PointPointDistanceSquared( c_new, c) > MAXIMUM_VERTEX_MOVEMENT_SQUARED )
    {
        return false;
    }


    /*
     * check if any of the new vertices is out of bounds
     */
    float maxVal = ((float)bins) + MAXIMUM_VERTEX_CORRECTION;
    float minVal = -MAXIMUM_VERTEX_CORRECTION;
    if ( ( a_new.x < minVal || a_new.x > maxVal || a_new.y < minVal || a_new.y > maxVal ) ||
         ( b_new.x < minVal || b_new.x > maxVal || b_new.y < minVal || b_new.y > maxVal ) ||
         ( c_new.x < minVal || c_new.x > maxVal || c_new.y < minVal || c_new.y > maxVal ) )
    {
        // at least one of the new vertices is out of bounds, cannot expand the triangle
        return false;
    }

    /*
     * set new points
     */
    dx = a_new.x;
    dy = a_new.y;
    ex = b_new.x;
    ey = b_new.y;
    fx = c_new.x;
    fy = c_new.y;

    return true; // success!

}



#endif // NPW_CUDA_SHAPE_EXPANSION_H_
