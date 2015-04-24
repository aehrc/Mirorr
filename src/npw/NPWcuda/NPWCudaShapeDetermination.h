/*
 * NPWCudaShapeDetermination.h
 *
 *  Created on: 19/07/2010
 *      Author: bro86j
 */

#ifndef NPW_CUDA_SHAPE_DETERMINATION_H_
#define NPW_CUDA_SHAPE_DETERMINATION_H_

//#define NPW_TRIANGLE_EXPANSION // NOTE: THIS IS NOT FULLY IMPLEMENTED, DO NOT USE!!!

#define GEOMETRY_LIST_SIZE_FOR_POINTS           4 // 1 vertex
#define GEOMETRY_LIST_SIZE_FOR_ONE_TRIANGLE     8 // 3 vertices
#define GEOMETRY_LIST_SIZE_FOR_TWO_TRIANGLES   10 // 4 vertices
#define GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES 10 // 4 vertices
#define GEOMETRY_LIST_SIZE_FOR_FOUR_TRIANGLES  12 // 5 vertices

#include "NPWCudaGeometryFunctions.h"
#include "NPWCudaShapeExpansion.h"
#include "NPWCudaConstants.h"


__device__ void sortVertices(
        float2 * unique_Vertices, unsigned int & n_Unique_Vertices,
        unsigned int * n_Duplicates_at_Unique_Vertices )
{
    switch ( n_Unique_Vertices )
    {

    /*
     * sorting is only necessary for 3 and 4 unique vertices
     */
    case 3:
    {
        float d01 = PointPointDistanceSquared(unique_Vertices[0], unique_Vertices[1]);
        float d02 = PointPointDistanceSquared(unique_Vertices[0], unique_Vertices[2]);
        float d12 = PointPointDistanceSquared(unique_Vertices[1], unique_Vertices[2]);

        if ( d02 > d01 && d02 > d12 )
        {
            // 012
            // already sorted
        }
        else if ( d01 > d12 && d01 > d02 )
        {
            // 120
            float2       tv = unique_Vertices[1];
            unsigned int td = n_Duplicates_at_Unique_Vertices[1];

            unique_Vertices[1]                 = unique_Vertices[2];
            n_Duplicates_at_Unique_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
            unique_Vertices[2]                 = unique_Vertices[0];
            n_Duplicates_at_Unique_Vertices[2] = n_Duplicates_at_Unique_Vertices[0];
            unique_Vertices[0]                 = tv;
            n_Duplicates_at_Unique_Vertices[0] = td;
        }
        else if ( d12 > d01 && d12 > d02 )
        {
            // 201
            float2       tv = unique_Vertices[0];
            unsigned int td = n_Duplicates_at_Unique_Vertices[0];

            unique_Vertices[0]                 = unique_Vertices[2];
            n_Duplicates_at_Unique_Vertices[0] = n_Duplicates_at_Unique_Vertices[2];
            unique_Vertices[2]                 = unique_Vertices[1];
            n_Duplicates_at_Unique_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
            unique_Vertices[1]                 = tv;
            n_Duplicates_at_Unique_Vertices[1] = td;
        }

        break;
    }

    case 4:
    {
        float d01 = PointPointDistanceSquared(unique_Vertices[0], unique_Vertices[1]);
        float d02 = PointPointDistanceSquared(unique_Vertices[0], unique_Vertices[2]);
        float d03 = PointPointDistanceSquared(unique_Vertices[0], unique_Vertices[3]);
        float d12 = PointPointDistanceSquared(unique_Vertices[1], unique_Vertices[2]);
        float d13 = PointPointDistanceSquared(unique_Vertices[1], unique_Vertices[3]);
        float d23 = PointPointDistanceSquared(unique_Vertices[2], unique_Vertices[3]);

        if ( d01 > d02 && d01 > d03 && d01 > d12 && d01 > d13 && d01 > d23 )
        {
            // 0xx1
            float2       tv = unique_Vertices[3];
            unsigned int td = n_Duplicates_at_Unique_Vertices[3];

            unique_Vertices[3]                 = unique_Vertices[1];
            n_Duplicates_at_Unique_Vertices[3] = n_Duplicates_at_Unique_Vertices[1];

            if ( d03 > d02 )
            {
                // 0231
                unique_Vertices[1]                 = unique_Vertices[2];
                n_Duplicates_at_Unique_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                unique_Vertices[2]                 = tv;
                n_Duplicates_at_Unique_Vertices[2] = td;
            }
            else
            {
                // 0321
                unique_Vertices[1]                 = tv;
                n_Duplicates_at_Unique_Vertices[1] = td;
            }
        }
        else if ( d02 > d01 && d02 > d03 && d02 > d12 && d02 > d13 && d02 > d23 )
        {
            // 0xx2
            float2       tv = unique_Vertices[3];
            unsigned int td = n_Duplicates_at_Unique_Vertices[3];

            unique_Vertices[3]                 = unique_Vertices[2];
            n_Duplicates_at_Unique_Vertices[3] = n_Duplicates_at_Unique_Vertices[2];

            if ( d03 > d01 )
            {
                // 0132
                unique_Vertices[2]                 = tv;
                n_Duplicates_at_Unique_Vertices[2] = td;
            }
            else
            {
                // 0312
                unique_Vertices[2]                 = unique_Vertices[1];
                n_Duplicates_at_Unique_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
                unique_Vertices[1]                 = tv;
                n_Duplicates_at_Unique_Vertices[1] = td;
            }
        }
        else if ( d03 > d01 && d03 > d02 && d03 > d12 && d03 > d13 && d03 > d23 )
        {
            // 0xx3

            if ( d02 > d01 )
            {
                // 0123
            }
            else
            {
                // 0213
                float2       tv = unique_Vertices[1];
                unsigned int td = n_Duplicates_at_Unique_Vertices[1];

                unique_Vertices[1]                 = unique_Vertices[2];
                n_Duplicates_at_Unique_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                unique_Vertices[2]                 = tv;
                n_Duplicates_at_Unique_Vertices[2] = td;
            }
        }
        else if ( d12 > d01 && d12 > d13 && d12 > d02 && d12 > d03 && d12 > d23 )
        {
            // 1xx2
            float2       tv0 = unique_Vertices[0];
            unsigned int td0 = n_Duplicates_at_Unique_Vertices[0];
            float2       tv3 = unique_Vertices[3];
            unsigned int td3 = n_Duplicates_at_Unique_Vertices[3];

            unique_Vertices[0]                 = unique_Vertices[1];
            n_Duplicates_at_Unique_Vertices[0] = n_Duplicates_at_Unique_Vertices[1];
            unique_Vertices[3]                 = unique_Vertices[2];
            n_Duplicates_at_Unique_Vertices[3] = n_Duplicates_at_Unique_Vertices[2];

            if ( d13 > d01 )
            {
                // 1032
                unique_Vertices[1]                 = tv0;
                n_Duplicates_at_Unique_Vertices[1] = td0;
                unique_Vertices[2]                 = tv3;
                n_Duplicates_at_Unique_Vertices[2] = td3;
            }
            else
            {
                // 1302
                unique_Vertices[1]                 = tv3;
                n_Duplicates_at_Unique_Vertices[1] = td3;
                unique_Vertices[2]                 = tv0;
                n_Duplicates_at_Unique_Vertices[2] = td0;
            }
        }
        else if ( d13 > d01 && d13 > d02 && d13 > d03 && d13 > d12 && d13 > d23 )
        {
            // 1xx3
            float2       tv = unique_Vertices[0];
            unsigned int td = n_Duplicates_at_Unique_Vertices[0];

            unique_Vertices[0]                 = unique_Vertices[1];
            n_Duplicates_at_Unique_Vertices[0] = n_Duplicates_at_Unique_Vertices[1];

            if ( d12 > d01 )
            {
                // 1023
                unique_Vertices[1]                 = tv;
                n_Duplicates_at_Unique_Vertices[1] = td;
            }
            else
            {
                // 1203
                unique_Vertices[1]                 = unique_Vertices[2];
                n_Duplicates_at_Unique_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                unique_Vertices[2]                 = tv;
                n_Duplicates_at_Unique_Vertices[2] = td;
            }
        }
        else if ( d23 > d01 && d23 > d02 && d23 > d03 && d23 > d12 && d23 > d13 )
        {
            // 2xx3
            float2       tv = unique_Vertices[0];
            unsigned int td = n_Duplicates_at_Unique_Vertices[0];

            unique_Vertices[0]                 = unique_Vertices[2];
            n_Duplicates_at_Unique_Vertices[0] = n_Duplicates_at_Unique_Vertices[2];

            if ( d12 > d02 )
            {
                // 2013
                unique_Vertices[2]                 = unique_Vertices[1];
                n_Duplicates_at_Unique_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
                unique_Vertices[1]                 = tv;
                n_Duplicates_at_Unique_Vertices[1] = td;
            }
            else
            {
                // 2103
                unique_Vertices[2]                 = tv;
                n_Duplicates_at_Unique_Vertices[2] = td;
            }
        }

        break;
    }

    default:
        break;

    }
}



__device__ void setVerticesAndDetectShape(
        float2 * unique_Vertices, unsigned int & n_Unique_Vertices,
        unsigned int * n_Duplicates_at_Unique_Vertices,
        unsigned int & index_of_maximum_intensity_vertex_in_triangle,
        unsigned int & index_of_opposite_vertex,
        Shape & shape,
        float2 v1, float2 v2, float2 v3, float2 v4 )
{
    bool all_points_are_on_a_line = true;
    shape                         = UNDEFINED;
    unsigned int index_of_middle_vertex_on_line;


    /*
     * add first vertex
     */
    unique_Vertices[0] = v1;
    n_Duplicates_at_Unique_Vertices[0] = 0;
    n_Unique_Vertices = 1;

    /*
     * add second vertex
     */
    if ( PointPointDistanceSquared(unique_Vertices[0], v2) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[0];
    }
    else
    {
        unique_Vertices[1] = v2;
        n_Duplicates_at_Unique_Vertices[1] = 0;
        ++n_Unique_Vertices;
    }

    /*
     * add third vertex
     */
    if ( PointPointDistanceSquared(unique_Vertices[0], v3) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[0];
    }
    else if ( n_Unique_Vertices > 1 &&
            PointPointDistanceSquared(unique_Vertices[1], v3) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[1];
    }
    else
    {
        /*
         * check if the third vertex is close to the existing line
         */
        if ( n_Unique_Vertices > 1 )
        {
            float2 closestPoint = make_float2(0.0f,0.0f);
            float t;
            if ( isVertexCloseToLine(
                    unique_Vertices[0], unique_Vertices[1], v3, closestPoint, t, false ) )
            {
                unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                ++n_Unique_Vertices;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 1 : 2 ) );
            }
            /*
             * also check if the existing vertices are close to the newly created lines
             * if they are, we have forced them to be on a line
             */
            // uV[0] to uV[1],v3
            else if ( isVertexCloseToLine(
                         unique_Vertices[1], v3, unique_Vertices[0], closestPoint, t, false ) )
            {
                unique_Vertices[0] = closestPoint;
                all_points_are_on_a_line = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 1 : ( t > 1.0f ? 2 : 0 ) );
                unique_Vertices[ n_Unique_Vertices ] = v3;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
            }
            // uV[1] to uV[0],v3
            else if ( isVertexCloseToLine(
                    unique_Vertices[0], v3, unique_Vertices[1], closestPoint, t, false ) )
            {
                unique_Vertices[1] = closestPoint;
                all_points_are_on_a_line = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 2 : 1 ) );
                unique_Vertices[ n_Unique_Vertices ] = v3;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
            }
            else
            {
                // it is a unique vertex and the three vertices are not on a line
                unique_Vertices[ n_Unique_Vertices ] = v3;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
                all_points_are_on_a_line = false;
            }
        }
        else
        {
            /*
             * there is only one unique vertex and v3 is the second unique vertex
             */
            unique_Vertices[ n_Unique_Vertices ] = v3;
            n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
            ++n_Unique_Vertices;
        }
    }


    /*
     * add fourth vertex
     */
    if ( PointPointDistanceSquared(unique_Vertices[0], v4) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[0];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 0;
        }
        else if ( n_Unique_Vertices == 1 )
        {
            shape = TRUE_POINT;
        }
    }
    else if ( n_Unique_Vertices > 1 &&
              PointPointDistanceSquared(unique_Vertices[1], v4) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[1];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 1;
        }
    }
    else if ( n_Unique_Vertices > 2 &&
              PointPointDistanceSquared(unique_Vertices[2], v4) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[2];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 2;
        }
    }
    else
    {
        if ( n_Unique_Vertices == 1 )
        {
            /*
             * v4 is the second unique point, so all others are identical
             */
            unique_Vertices[ n_Unique_Vertices ] = v4;

            n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
            ++n_Unique_Vertices;
            shape = LINE_WITH_VERTEX_TRIO;
        }
        /*
         * check if the fourth vertex is close to the existing line(s)
         * there are at least 2 other unique vertices
         */
        else
        {
            float2 closestPoint = make_float2(0.0f,0.0f);
            float t;

            bool putOnLine = false;
            // check if v4 is close to 0-1
            if ( isVertexCloseToLine(
                    unique_Vertices[0], unique_Vertices[1], v4, closestPoint, t, false ) )
            {
                unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
                putOnLine = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 1 : n_Unique_Vertices-1 ) );
            }
            /*
             * if there are already 3 points on a line, moving one of those points so v4 makes
             * a line with 2 other points breaks the original line.
             * in that case (i.e. a triangle), only try to move v4 onto the existing line
             */
            // check if 1 is close to v4-0
            else if ( ( n_Unique_Vertices < 3 || !all_points_are_on_a_line ) )
            {
                if ( isVertexCloseToLine(
                         unique_Vertices[0], v4, unique_Vertices[1], closestPoint, t, false ) )
                {
                    unique_Vertices[1] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? n_Unique_Vertices-1 : 1) );
                }
                // check if 0 is close to v4-1
                else if ( isVertexCloseToLine(
                        unique_Vertices[1], v4, unique_Vertices[0], closestPoint, t, false ) )
                {
                    unique_Vertices[0] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 1 : ( t > 1.0f ? n_Unique_Vertices-1 : 0 ) );
                }
            }

            if ( putOnLine )
            {
                if ( !all_points_are_on_a_line )
                {
                    shape = TRIANGLE_WITH_ONE_VERTEX_ON_BORDER;
                    index_of_maximum_intensity_vertex_in_triangle = index_of_middle_vertex_on_line;
                    index_of_opposite_vertex = 2;
                }
                else if ( n_Unique_Vertices == 4 )
                {
                    shape = LINE_WITHOUT_VERTEX_PAIR;
                }
            }

            /*
             * else, v4 is not close to an existing line or an existing vertex is not close to
             * a new line. I.e. v4 - uV[0] - uV[1] is not a line
             * but there may be another unique vertex!
             */
            else if ( n_Unique_Vertices > 2 )
            {
                /*
                 * check 6 lines
                 * we already checked:
                 *    v4 to 0-1
                 *    0  to 1-v4
                 *    1  to 0-v4
                 *
                 *    same as before, only move uV[0], uV[1] or uV[2]
                 *    if those points do not lie on a line
                 */
                // check if v4 is close to 0-2
                if ( isVertexCloseToLine(
                        unique_Vertices[0], unique_Vertices[2], v4,
                        closestPoint, t, false) )
                {
                    unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 0 : ( t > 1.0f ? 2 : n_Unique_Vertices-1 ) );
                    index_of_opposite_vertex = 1;
                }
                // check if v4 is close to 1-2
                else if ( isVertexCloseToLine(
                        unique_Vertices[1], unique_Vertices[2], v4,
                        closestPoint, t, false) )
                {
                    unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 1 : ( t > 1.0f ? 2 : n_Unique_Vertices-1 ) );
                    index_of_opposite_vertex = 0;
                }
                // check if 0 is close to 2-v4
                else if ( !all_points_are_on_a_line &&
                        isVertexCloseToLine(
                        unique_Vertices[2], v4, unique_Vertices[0],
                        closestPoint, t, false) )
                {
                    unique_Vertices[0] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 2 : ( t > 1.0f ? n_Unique_Vertices-1 : 0 ) );
                    index_of_opposite_vertex = 1;
                }
                // check if 1 is close to 2-v4
                else if ( !all_points_are_on_a_line &&
                        isVertexCloseToLine(
                        unique_Vertices[2], v4, unique_Vertices[1],
                        closestPoint, t, false) )
                {
                    unique_Vertices[1] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 2 : ( t > 1.0f ? n_Unique_Vertices-1 : 1 ) );
                    index_of_opposite_vertex = 0;
                }
                // check if 2 is close to 0-v4
                else if ( !all_points_are_on_a_line &&
                        isVertexCloseToLine(
                        unique_Vertices[0], v4, unique_Vertices[2],
                        closestPoint, t, false) )
                {
                    unique_Vertices[2] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 0 : ( t > 1.0f ? n_Unique_Vertices-1 : 2 ) );
                    index_of_opposite_vertex = 1;
                }
                // check if 2 is close to 1-v4
                else if ( !all_points_are_on_a_line &&
                        isVertexCloseToLine(
                        unique_Vertices[1], v4, unique_Vertices[2],
                        closestPoint, t, false) )
                {
                    unique_Vertices[2] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 1 : ( t > 1.0f ? n_Unique_Vertices-1 : 2 ) );
                    index_of_opposite_vertex = 0;
                }

                if ( putOnLine )
                {
                    if ( !all_points_are_on_a_line )
                    {
                        // v4 was put on a line, but the shape was already a triangle
                        shape = TRIANGLE_WITH_ONE_VERTEX_ON_BORDER;
                        index_of_maximum_intensity_vertex_in_triangle = index_of_middle_vertex_on_line;
                    }
                    else
                    {
                        // four unique points on a line
                        shape = LINE_WITHOUT_VERTEX_PAIR;
                    }
                }
                else if ( all_points_are_on_a_line )
                {
                    // v4 is not on a line, but uV[0], uv[1] and uv[2] are on a line
                    // and they are all unique vertices
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    all_points_are_on_a_line = false;

                    shape = TRIANGLE_WITH_ONE_VERTEX_ON_BORDER;

                    // maximum intensity is the middle point on the line uV[0] - uV[1] - uV[2]
                    index_of_maximum_intensity_vertex_in_triangle = index_of_middle_vertex_on_line;
                    index_of_opposite_vertex = 3;
                }
                else
                {
                    // there are 3 unique vertices, v4 is the fourth unique vertex
                    // and the shape is not a line
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    all_points_are_on_a_line = false;


                    // the shape is a TRIANGLE_WITH_THREE_SUBTRIANGLES or a QUAD
                    // determine which
                    if ( PointInTriangle(
                            unique_Vertices[0], unique_Vertices[1],
                            unique_Vertices[2], unique_Vertices[3]) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 3;
                    }
                    else if ( PointInTriangle(
                            unique_Vertices[1], unique_Vertices[2],
                            unique_Vertices[3], unique_Vertices[0]) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 0;
                    }
                    else if ( PointInTriangle(
                            unique_Vertices[2], unique_Vertices[3],
                            unique_Vertices[0], unique_Vertices[1]) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 1;
                    }
                    else if ( PointInTriangle(
                            unique_Vertices[3], unique_Vertices[0],
                            unique_Vertices[1], unique_Vertices[2]) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 2;
                    }
                    else
                    {
                        shape = QUAD;
                    }
                }
            }
            else
            {
                // there are 2 unique vertices and v4 is the third unique vertex
                // the shape is not a line
                unique_Vertices[ n_Unique_Vertices ] = v4;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
                all_points_are_on_a_line = false;

                shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
                // set index_of_maximum_intensity_vertex_in_triangle
                index_of_maximum_intensity_vertex_in_triangle =
//                    ( n_Duplicates_at_Unique_Vertices[1] == 1 ? 1 : 0 );
                        n_Duplicates_at_Unique_Vertices[1]; // look how smart!
            }
        }
    }

    /*
     * shapes that have already been determined:
     *
     * QUAD (intersect and order not yet determined)
     * TRIANGLE_WITH_ONE_VERTEX_PAIR ( pair index is known )
     * TRIANGLE_WITH_ONE_VERTEX_ON_BORDER ( border vertex index is known )
     * TRIANGLE_WITH_THREE_SUBTRIANGLES ( inner vertex index is known )
     * LINE_WITH_VERTEX_TRIO ( storing index is waste of space,
     *                         using n_Duplicates_at_Unique_Vertices works just as well )
     * LINE_WITHOUT_VERTEX_PAIR
     * TRUE_POINT
     *
     * So, all the renderer has to do is call setVerticesAndGetShape(), setRenderList() and
     * calculateHeight()
     * then the renderer can access the render list in this class
     * and render as many vertices as required
     */

    /*
     * if it is a line, we need to sort
     */
    if ( shape == UNDEFINED || shape == LINE_WITH_VERTEX_TRIO || shape == LINE_WITHOUT_VERTEX_PAIR )
    {
        sortVertices(unique_Vertices, n_Unique_Vertices, n_Duplicates_at_Unique_Vertices);


        if ( shape == UNDEFINED )
        {
            /*
             * detect the shapes we didn't detect yet
             */
            if ( n_Unique_Vertices == 2 )
            {
                shape = LINE_WITH_TWO_VERTEX_PAIRS;
            }
            else if ( n_Duplicates_at_Unique_Vertices[1] == 1 )
            {
                shape = LINE_WITH_VERTEX_PAIR_ON_LINE;
            }
            else
            {
                shape = LINE_WITH_VERTEX_PAIR_ON_END;
            }
        }
    }

    /*
     * if it's a line, it may be too short
     */
    if ( shape == LINE_WITHOUT_VERTEX_PAIR     || shape == LINE_WITH_TWO_VERTEX_PAIRS ||
         shape == LINE_WITH_VERTEX_PAIR_ON_END || shape == LINE_WITH_VERTEX_PAIR_ON_LINE ||
         shape == LINE_WITH_VERTEX_TRIO )
    {
        if ( PointPointDistanceSquared( unique_Vertices[0], unique_Vertices[n_Unique_Vertices-1] ) <
                MINIMUM_LINE_LENGTH_SQUARED )
        {
            unique_Vertices[0] = ( unique_Vertices[0] + unique_Vertices[n_Unique_Vertices-1] ) * 0.5f;
            shape = POINT_EQUIVALENT;
        }
    }

}

__device__ void getQuadIntersectAndSortVertices( float2 * vertices, float2 & intersect )
{
    // get the intersection between 02 and 13
    if ( LineSegmentLineSegmentIntersect(
            vertices[0], vertices[2], vertices[1], vertices[3], intersect ) ) // quad: 0123
    {
        return;
    }

    // get the intersection between 01 and 23
    if ( LineSegmentLineSegmentIntersect(
            vertices[0], vertices[1], vertices[2], vertices[3], intersect ) ) // quad: 0213
    {
        // swap 1 and 2
        float2 t    = vertices[1];
        vertices[1] = vertices[2];
        vertices[2] = t;
        return;
    }

    // get the intersection between 03 and 12
    if (LineSegmentLineSegmentIntersect(
            vertices[0], vertices[3], vertices[1], vertices[2], intersect ) ) // quad: 0132
    {
        // swap 2 and 3
        float2 t    = vertices[2];
        vertices[2] = vertices[3];
        vertices[3] = t;
        return;
    }
}

/*
 * sets the geometry list and returns how many floats contain relevant data
 */
__device__ unsigned int setGeometryList(
        float2 * unique_Vertices, unsigned int n_Unique_Vertices,
        unsigned int * n_Duplicates_at_Unique_Vertices,
        unsigned int index_of_maximum_intensity_vertex_in_triangle,
        unsigned int index_of_opposite_vertex,
        Shape shape,
        float * geometryList )
{
//    if ( //shape == TRIANGLE_WITH_ONE_VERTEX_ON_BORDER ||   // OK!
////         shape == TRIANGLE_WITH_ONE_VERTEX_PAIR ||        // OK!
////         shape == TRIANGLE_WITH_THREE_SUBTRIANGLES ||     // OK!
//         shape == QUAD )                                    // QUADS ARE WRONG!!!
//    {
//        return 0;
//    }

    switch ( shape )
    {

    case TRUE_POINT:
    case POINT_EQUIVALENT:
    {
        geometryList[2] = unique_Vertices[0].x;
        geometryList[3] = unique_Vertices[0].y;
        geometryList[0] = 0.0f;
        return GEOMETRY_LIST_SIZE_FOR_POINTS;
    }

    case LINE_WITH_VERTEX_TRIO:
    {
        correctLineWithVertexTrio(
                unique_Vertices,n_Unique_Vertices,n_Duplicates_at_Unique_Vertices,geometryList);
        return GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES;
    }

    case LINE_WITH_TWO_VERTEX_PAIRS:
    {
        correctLineWithTwoVertexPairs(
                unique_Vertices,n_Unique_Vertices,n_Duplicates_at_Unique_Vertices,geometryList);
        return GEOMETRY_LIST_SIZE_FOR_FOUR_TRIANGLES;
    }

    case LINE_WITH_VERTEX_PAIR_ON_END:
    {
        correctLineWithVertexPairOnEnd(
                unique_Vertices,n_Unique_Vertices,n_Duplicates_at_Unique_Vertices,geometryList);
        return GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES;
    }

    case LINE_WITH_VERTEX_PAIR_ON_LINE:
    {
        correctLineWithVertexPairOnLine(
                unique_Vertices,n_Unique_Vertices,n_Duplicates_at_Unique_Vertices,geometryList);
        return GEOMETRY_LIST_SIZE_FOR_FOUR_TRIANGLES;
    }

    case LINE_WITHOUT_VERTEX_PAIR:
    {
        correctLineWithoutVertexPair(
                unique_Vertices,n_Unique_Vertices,n_Duplicates_at_Unique_Vertices,geometryList);
        return GEOMETRY_LIST_SIZE_FOR_FOUR_TRIANGLES;
    }

    case TRIANGLE_WITH_ONE_VERTEX_PAIR:
    {
#ifdef NPW_TRIANGLE_EXPANSION
        geometryList[2] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle      ].x;
        geometryList[3] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle      ].y;

        if( !expandTriangle(
                unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle     ], // originals
                unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3],
                unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3],
                geometryList[4], geometryList[5], // output
                geometryList[6], geometryList[7],
                geometryList[8], geometryList[9],
                bins ) ) // number of bins, needed for out of bounds checking
        {
            geometryList[0] = 3.0f;

            return GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES;
        }
        else
        {
            geometryList[4] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3 ].x;
            geometryList[5] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3 ].y;
            geometryList[6] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3 ].x;
            geometryList[7] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3 ].y;

            return GEOMETRY_LIST_SIZE_FOR_ONE_TRIANGLE;
        }
#else
        geometryList[2] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle      ].x;
        geometryList[3] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle      ].y;
        geometryList[4] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3 ].x;
        geometryList[5] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3 ].y;
        geometryList[6] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3 ].x;
        geometryList[7] = unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3 ].y;

        geometryList[0] = 1.0f;

        return GEOMETRY_LIST_SIZE_FOR_ONE_TRIANGLE;
#endif

    }
    case TRIANGLE_WITH_ONE_VERTEX_ON_BORDER:
    {

#ifdef NPW_TRIANGLE_EXPANSION
        geometryList[2] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ].x;
        geometryList[3] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ].y;

        if( expandTriangle(
                unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%4], // originals
                unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%4],
                unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+3)%4],
                geometryList[4], geometryList[5], // output
                geometryList[6], geometryList[7],
                geometryList[8], geometryList[9],
                bins ) ) // number of bins, needed for out of bounds checking
        {
            geometryList[0] = 3.0f;

            return GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES;
        }
        else
        {
            unsigned int idx1 = (index_of_maximum_intensity_vertex_in_triangle+1)%4;
            unsigned int idx2 = (idx1+1)%4;

            while ( idx1 == index_of_opposite_vertex )
            {
                idx1 = (idx1+1)%4;
            }
            while ( idx2 == index_of_opposite_vertex || idx2 == idx1)
            {
                idx2 = (idx2+1)%4;
            }

            geometryList[4] = unique_Vertices[ index_of_opposite_vertex ].x;
            geometryList[5] = unique_Vertices[ index_of_opposite_vertex ].y;
            geometryList[6] = unique_Vertices[ idx1 ].x;
            geometryList[7] = unique_Vertices[ idx1 ].y;
            geometryList[8] = unique_Vertices[ idx2 ].x;
            geometryList[9] = unique_Vertices[ idx2 ].y;

            geometryList[0] = 2.0f;

            return GEOMETRY_LIST_SIZE_FOR_TWO_TRIANGLES;
        }
#else
        unsigned int idx1 = (index_of_maximum_intensity_vertex_in_triangle+1)%4;
        unsigned int idx2 = (idx1+1)%4;

        while ( idx1 == index_of_opposite_vertex )
        {
            idx1 = (idx1+1)%4;
        }
        while ( idx2 == index_of_opposite_vertex || idx2 == idx1)
        {
            idx2 = (idx2+1)%4;
        }

        geometryList[2] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ].x;
        geometryList[3] = unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ].y;
        geometryList[4] = unique_Vertices[ index_of_opposite_vertex ].x;
        geometryList[5] = unique_Vertices[ index_of_opposite_vertex ].y;
        geometryList[6] = unique_Vertices[ idx1 ].x;
        geometryList[7] = unique_Vertices[ idx1 ].y;
        geometryList[8] = unique_Vertices[ idx2 ].x;
        geometryList[9] = unique_Vertices[ idx2 ].y;

        geometryList[0] = 2.0f;

        return GEOMETRY_LIST_SIZE_FOR_TWO_TRIANGLES;
#endif
    }

    case TRIANGLE_WITH_THREE_SUBTRIANGLES:
    {
        geometryList[ 2] = unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ].x;
        geometryList[ 3] = unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ].y;
        geometryList[ 4] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ].x;
        geometryList[ 5] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ].y;
        geometryList[ 6] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ].x;
        geometryList[ 7] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ].y;
        geometryList[ 8] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ].x;
        geometryList[ 9] = unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ].y;

        geometryList[ 0] = 3.0f;

        return GEOMETRY_LIST_SIZE_FOR_THREE_TRIANGLES;
    }

    case QUAD:
    {
        float2 intersect = make_float2(0.0f,0.0f);
        getQuadIntersectAndSortVertices( unique_Vertices, intersect);

        geometryList[ 2] = intersect.x;
        geometryList[ 3] = intersect.y;
        geometryList[ 4] = unique_Vertices[0].x;
        geometryList[ 5] = unique_Vertices[0].y;
        geometryList[ 6] = unique_Vertices[1].x;
        geometryList[ 7] = unique_Vertices[1].y;
        geometryList[ 8] = unique_Vertices[2].x;
        geometryList[ 9] = unique_Vertices[2].y;
        geometryList[10] = unique_Vertices[3].x;
        geometryList[11] = unique_Vertices[3].y;

        geometryList[ 0] = 4.0f;

        return GEOMETRY_LIST_SIZE_FOR_FOUR_TRIANGLES;
    }

    default:
        geometryList[0] = 0.0f;
        return 1;
    }
}


__device__ float getVertexHeight(
        float * geometryList, float weight )
{
    float area2;

    switch ( (int) geometryList[0] )
    {

    case 0:
    {
        return weight;
    }

    case 1:
    {
        area2  = SignedTriangleAreaTimesTwo( geometryList[2], geometryList[3],
                                             geometryList[4], geometryList[5],
                                             geometryList[6], geometryList[7] );
        break;
    }

    case 2:
    case 3:
    {
        area2  = SignedTriangleAreaTimesTwo( geometryList[4], geometryList[5],
                                             geometryList[6], geometryList[7],
                                             geometryList[8], geometryList[9] );
        break;
    }

    case 4:
    {
        area2  = SignedQuadAreaTimesTwo( geometryList[ 4], geometryList[ 5],
                                         geometryList[ 6], geometryList[ 7],
                                         geometryList[ 8], geometryList[ 9],
                                         geometryList[10], geometryList[11] );
        break;
    }

    default:
    {
        return 0.0f;
    }

    }

    // geometryList[1] = 6.0f * weight / Max( Abs(area2), MINIMUM_AREA_TIMES_TWO );
    return 6.0f * weight / Abs(area2);

}


__device__ void getGeometryFromVertices(
                        float2 v1, float2 v2, float2 v3, float2 v4,
                        float * d_geometryData, unsigned int geo_unit_offset,
                        float weight )
{
    /*
     * initialize & declare variables
     */
    float2       unique_Vertices[4];
    unsigned int n_Unique_Vertices = 0;
    unsigned int n_Duplicates_at_Unique_Vertices[4];
    unsigned int index_of_maximum_intensity_vertex_in_triangle;
    unsigned int index_of_opposite_vertex;
    Shape        shape = UNDEFINED;

    for (int i = 0; i < 4; ++i)
    {
        unique_Vertices[i]                 = make_float2(0.0f,0.0f);
        n_Duplicates_at_Unique_Vertices[i] = 0;
    }

    /*
     * get the shape
     */
    setVerticesAndDetectShape(
            unique_Vertices, n_Unique_Vertices,
            n_Duplicates_at_Unique_Vertices,
            index_of_maximum_intensity_vertex_in_triangle,
            index_of_opposite_vertex,
            shape,
            v1,v2,v3,v4 );

    /*
     * set the geometry render list
     */
    float geometryList[FLOATS_PER_GEO_UNIT];
    for ( int i = 0; i < FLOATS_PER_GEO_UNIT; ++i )
    {
        geometryList[i] = 0.0f;
    }

    // number of floats defining this geometry unit
    unsigned int relevantSize = setGeometryList(
            unique_Vertices, n_Unique_Vertices,
            n_Duplicates_at_Unique_Vertices,
            index_of_maximum_intensity_vertex_in_triangle,
            index_of_opposite_vertex,
            shape,
            geometryList
            );

    /*
     * set the height
     */
    geometryList[1] = getVertexHeight(geometryList, weight);

    /*
     * save the geometry in the global structure
     */
    for ( unsigned int i = 0; i < relevantSize; ++i )
    {
        d_geometryData[geo_unit_offset+i] = geometryList[i];
    }

}

#endif // NPW_CUDA_SHAPE_DETERMINATION_H_

