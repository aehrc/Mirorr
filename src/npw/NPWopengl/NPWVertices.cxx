/*
 * NPWVertices.cxx
 *
 *  Created on: 23/06/2010
 *      Author: bro86j
 */

#include "NPWVertices.h"
#include "NPWGeometry.h"
#include <cmath>

NPWVertices::NPWVertices() :
    setRenderListFunction ( 0 )
{
#ifdef DB_COUNT
    for (int i = 0; i < (int)UNDEFINED+1; ++i)
    {
        shapeCounters[i] = 0;
    }
    triangleExpansionFails = 0;
    quadExpansionFails = 0;
#endif

    /*
     * shape expansion is off by default
     */
    turnShapeExpansionOff();
}

void NPWVertices::setParentWindow( NPWWindow * newParent )
{
    parentWindow = newParent;
}

void NPWVertices::setVerticesAndDetectShape(
        const Vector2 & v1, const Vector2 & v2, const Vector2 & v3, const Vector2 & v4 )
{
    all_points_are_on_a_line = true;
    shape                    = UNDEFINED;

#ifdef DB_DEBUG
    logStream.str("");
    isErrored = false;
    logStream << "setting vertices: \n"
              << "\tv1 : " << v1 << "\n"
              << "\tv2 : " << v2 << "\n"
              << "\tv3 : " << v3 << "\n"
              << "\tv4 : " << v4 << "\n"
              << std::endl;
#endif

    /*
     * add first vertex
     */
    unique_Vertices[0] = v1;
    n_Duplicates_at_Unique_Vertices[0] = 0;
    n_Unique_Vertices = 1;
#ifdef DB_DEBUG
    logStream << "v1 set" << std::endl;
#endif

    /*
     * add second vertex
     */
    if ( unique_Vertices[0].DistanceSquared( v2 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[0];
#ifdef DB_DEBUG
        logStream << "v2 is equal to v1" << std::endl;
#endif
    }
    else
    {
        unique_Vertices[1] = v2;
        n_Duplicates_at_Unique_Vertices[1] = 0;
        ++n_Unique_Vertices;
#ifdef DB_DEBUG
        logStream << "v2 is unique" << std::endl;
#endif
    }

    /*
     * add third vertex
     */
    if ( unique_Vertices[0].DistanceSquared( v3 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[0];
#ifdef DB_DEBUG
        logStream << "v3 is equal to v1" << std::endl;
#endif
    }
    else if ( n_Unique_Vertices > 1 &&
              unique_Vertices[1].DistanceSquared( v3 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Unique_Vertices[1];
#ifdef DB_DEBUG
        logStream << "v3 is equal to v2" << std::endl;
#endif
    }
    else
    {
        /*
         * check if the third vertex is close to the existing line
         */
        if ( n_Unique_Vertices > 1 )
        {
            Vector2 closestPoint;
            float t;
            if ( NPWGeometryMath::isVertexCloseToLine(
                    unique_Vertices[0], unique_Vertices[1], v3, closestPoint, t, false ) )
            {
                unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                ++n_Unique_Vertices;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 1 : 2 ) );
#ifdef DB_DEBUG
                logStream << "v3 is close to the line v1 - v2, setting it to be at "
                          << closestPoint << std::endl;
#endif
            }
            /*
             * also check if the existing vertices are close to the newly created lines
             * if they are, we have forced them to be on a line
             */
            // uV[0] to uV[1],v3
            else if ( NPWGeometryMath::isVertexCloseToLine(
                         unique_Vertices[1], v3, unique_Vertices[0], closestPoint, t, false ) )
            {
                unique_Vertices[0] = closestPoint;
                all_points_are_on_a_line = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 1 : ( t > 1.0f ? 2 : 0 ) );
#ifdef DB_DEBUG
                logStream << "uV[0] is close to the line uV[1] - v3, setting uV[0] to be "
                          << closestPoint << std::endl;
#endif
                unique_Vertices[ n_Unique_Vertices ] = v3;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
            }
            // uV[1] to uV[0],v3
            else if ( NPWGeometryMath::isVertexCloseToLine(
                    unique_Vertices[0], v3, unique_Vertices[1], closestPoint, t, false ) )
            {
                unique_Vertices[1] = closestPoint;
                all_points_are_on_a_line = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 2 : 1 ) );
#ifdef DB_DEBUG
                logStream << "uV[1] is close to the line uV[0] - v3, setting uV[1] to be "
                          << closestPoint << std::endl;
#endif
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
#ifdef DB_DEBUG
                logStream << "v3 is a unique point; there are now 3 unique points not a line" << std::endl;
#endif
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
#ifdef DB_DEBUG
            logStream << "v3 is the second unique vertex" << std::endl;
#endif
        }
    }

#ifdef DB_DEBUG
    logStream << "we have now placed 3 vertices:\n"
              << "\t" << n_Unique_Vertices << " unique vertices" << std::endl;
    for ( int i = 0; i < n_Unique_Vertices; ++i )
    {
        logStream << "\tuV[" << i << "] : " << unique_Vertices[i]
                  << "\tDuplicates[" << i << "] : " << n_Duplicates_at_Unique_Vertices[i] << std::endl;
    }
#endif

    /*
     * add fourth vertex
     */
    if ( unique_Vertices[0].DistanceSquared( v4 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
#ifdef DB_DEBUG
        logStream << "v4 is equal to uV[0]" << std::endl;
#endif
        ++n_Duplicates_at_Unique_Vertices[0];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 0;
#ifdef DB_DEBUG
            logStream << "shape is a triangle with one pair, pair is at "
                      << unique_Vertices[index_of_maximum_intensity_vertex_in_triangle] << std::endl;
#endif
        }
        else if ( n_Unique_Vertices == 1 )
        {
            shape = TRUE_POINT;
#ifdef DB_DEBUG
            logStream << "shape is a true point" << std::endl;
#endif
        }
    }
    else if ( n_Unique_Vertices > 1 &&
              unique_Vertices[1].DistanceSquared( v4 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
#ifdef DB_DEBUG
            logStream << "v4 is equal to uV[1]" << std::endl;
#endif
        ++n_Duplicates_at_Unique_Vertices[1];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 1;
#ifdef DB_DEBUG
            logStream << "shape is a triangle with one pair, pair is at "
                      << unique_Vertices[index_of_maximum_intensity_vertex_in_triangle] << std::endl;
#endif
        }
    }
    else if ( n_Unique_Vertices > 2 &&
              unique_Vertices[2].DistanceSquared( v4 ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
#ifdef DB_DEBUG
            logStream << "v4 is equal to uV[2]" << std::endl;
#endif
        ++n_Duplicates_at_Unique_Vertices[2];
        if ( !all_points_are_on_a_line )
        {
            shape = TRIANGLE_WITH_ONE_VERTEX_PAIR;
            index_of_maximum_intensity_vertex_in_triangle = 2;
#ifdef DB_DEBUG
            logStream << "shape is a triangle with one pair, pair is at "
                      << unique_Vertices[index_of_maximum_intensity_vertex_in_triangle] << std::endl;
#endif
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
#ifdef DB_DEBUG
            logStream << "v4 is the second unique point, "
                      << "the shape is a line with vertex trio" << std::endl;
#endif
        }
        /*
         * check if the fourth vertex is close to the existing line(s)
         * there are at least 2 other unique vertices
         */
        else
        {
            Vector2 closestPoint;
            float t;

            bool putOnLine = false;
            // check if v4 is close to 0-1
            if ( NPWGeometryMath::isVertexCloseToLine(
                    unique_Vertices[0], unique_Vertices[1], v4, closestPoint, t, false ) )
            {
                unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                ++n_Unique_Vertices;
                putOnLine = true;
                index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? 1 : n_Unique_Vertices-1 ) );
#ifdef DB_DEBUG
                logStream << "v4 is close to uV[0]-uV[1], the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
            }
            /*
             * if there are already 3 points on a line, moving one of those points so v4 makes
             * a line with 2 other points breaks the original line.
             * in that case (i.e. a triangle), only try to move v4 onto the existing line
             */
            // check if 1 is close to v4-0
            else if ( ( n_Unique_Vertices < 3 || !all_points_are_on_a_line ) )
            {
                if ( NPWGeometryMath::isVertexCloseToLine(
                         unique_Vertices[0], v4, unique_Vertices[1], closestPoint, t, false ) )
                {
                    unique_Vertices[1] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 0 : ( t > 1.0f ? n_Unique_Vertices-1 : 1) );
#ifdef DB_DEBUG
                    logStream << "uV[1] is close to uV[0]-v4, the vertex in the middle of the line is "
                              << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if 0 is close to v4-1
                else if ( NPWGeometryMath::isVertexCloseToLine(
                        unique_Vertices[1], v4, unique_Vertices[0], closestPoint, t, false ) )
                {
                    unique_Vertices[0] = closestPoint;
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                        ( t < 0.0f ? 1 : ( t > 1.0f ? n_Unique_Vertices-1 : 0 ) );
#ifdef DB_DEBUG
                    logStream << "uV[0] is close to uV[1]-v4, the vertex in the middle of the line is "
                              << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
            }

            if ( putOnLine )
            {
                if ( !all_points_are_on_a_line )
                {
                    shape = TRIANGLE_WITH_ONE_VERTEX_ON_BORDER;
                    index_of_maximum_intensity_vertex_in_triangle = index_of_middle_vertex_on_line;
#ifdef DB_DEBUG
                logStream << "v4 was put on a line, but the other points were not a line, "
                          << "it is a triangle with a vertex on the border, the border vertex is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                }
                else if ( n_Unique_Vertices == 4 )
                {
                    shape = LINE_WITHOUT_VERTEX_PAIR;
#ifdef DB_DEBUG
                    logStream << "it is a line without a vertex pair" << std::endl;
#endif
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
                if ( NPWGeometryMath::isVertexCloseToLine(
                        unique_Vertices[0], unique_Vertices[2], v4,
                        closestPoint, t, false) )
                {
                    unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 0 : ( t > 1.0f ? 2 : n_Unique_Vertices-1 ) );
#ifdef DB_DEBUG
                    logStream << "v4 is close to uV[0]-uV[2], the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if v4 is close to 1-2
                else if ( NPWGeometryMath::isVertexCloseToLine(
                        unique_Vertices[1], unique_Vertices[2], v4,
                        closestPoint, t, false) )
                {
                    unique_Vertices[ n_Unique_Vertices ] = closestPoint;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    putOnLine = true;
                    index_of_middle_vertex_on_line =
                            ( t < 0.0f ? 1 : ( t > 1.0f ? 2 : n_Unique_Vertices-1 ) );
#ifdef DB_DEBUG
                    logStream << "v4 is close to uV[1]-uV[2], the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if 0 is close to 2-v4
                else if ( !all_points_are_on_a_line &&
                        NPWGeometryMath::isVertexCloseToLine(
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
#ifdef DB_DEBUG
                    logStream << "uV[0] is close to uV[2]-v4, the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if 1 is close to 2-v4
                else if ( !all_points_are_on_a_line &&
                        NPWGeometryMath::isVertexCloseToLine(
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
#ifdef DB_DEBUG
                    logStream << "uV[1] is close to uV[2]-v4, the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if 2 is close to 0-v4
                else if ( !all_points_are_on_a_line &&
                        NPWGeometryMath::isVertexCloseToLine(
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
#ifdef DB_DEBUG
                    logStream << "uV[2] is close to uV[0]-v4, the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }
                // check if 2 is close to 1-v4
                else if ( !all_points_are_on_a_line &&
                        NPWGeometryMath::isVertexCloseToLine(
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
#ifdef DB_DEBUG
                    logStream << "uV[2] is close to uV[1]-v4, the vertex in the middle of the line is "
                          << unique_Vertices[index_of_middle_vertex_on_line] << std::endl;
#endif
                }

                if ( putOnLine )
                {
                    if ( !all_points_are_on_a_line )
                    {
                        // v4 was put on a line, but the shape was already a triangle
                        shape = TRIANGLE_WITH_ONE_VERTEX_ON_BORDER;
                        index_of_maximum_intensity_vertex_in_triangle = index_of_middle_vertex_on_line;
#ifdef DB_DEBUG
                        logStream << "v4 was put on a line, but the other points were not a line, "
                          << "it is a triangle with a vertex on the border, the border vertex is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                    }
                    else
                    {
                        // four unique points on a line
                        shape = LINE_WITHOUT_VERTEX_PAIR;
#ifdef DB_DEBUG
                        logStream << "it is a line without a vertex pair" << std::endl;
#endif
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
#ifdef DB_DEBUG
                    logStream << "v4 was not put on a line, but the other points were on a line, "
                          << "it is a triangle with a vertex on the border, the border vertex is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                }
                else
                {
                    // there are 3 unique vertices, v4 is the fourth unique vertex
                    // and the shape is not a line
                    unique_Vertices[ n_Unique_Vertices ] = v4;
                    n_Duplicates_at_Unique_Vertices[ n_Unique_Vertices ] = 0;
                    ++n_Unique_Vertices;
                    all_points_are_on_a_line = false;

#ifdef DB_DEBUG
                    logStream << "v4 is the fourth unique vertex, but no three points are on a line "
                              << "it is a quad or a triangle with three subtriangles" << std::endl;
#endif

                    // the shape is a TRIANGLE_WITH_THREE_SUBTRIANGLES or a QUAD
                    // determine which
                    if ( NPWGeometryMath::isVertexInTriangle(
                            &unique_Vertices[0], &unique_Vertices[1],
                            &unique_Vertices[2], &unique_Vertices[3] ) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 3;
#ifdef DB_DEBUG
                        logStream << "it is a triangle with three subtriangles, the fan center is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                    }
                    else if ( NPWGeometryMath::isVertexInTriangle(
                            &unique_Vertices[1], &unique_Vertices[2],
                            &unique_Vertices[3], &unique_Vertices[0] ) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 0;
#ifdef DB_DEBUG
                        logStream << "it is a triangle with three subtriangles, the fan center is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                    }
                    else if ( NPWGeometryMath::isVertexInTriangle(
                            &unique_Vertices[2], &unique_Vertices[3],
                            &unique_Vertices[0], &unique_Vertices[1] ) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 1;
#ifdef DB_DEBUG
                        logStream << "it is a triangle with three subtriangles, the fan center is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                    }
                    else if ( NPWGeometryMath::isVertexInTriangle(
                            &unique_Vertices[3], &unique_Vertices[0],
                            &unique_Vertices[1], &unique_Vertices[2] ) )
                    {
                        shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
                        index_of_maximum_intensity_vertex_in_triangle = 2;
#ifdef DB_DEBUG
                        logStream << "it is a triangle with three subtriangles, the fan center is "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
                    }
                    else
                    {
                        shape = QUAD;
#ifdef DB_DEBUG
                        logStream << "it is a quad" << std::endl;
#endif
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
#ifdef DB_DEBUG
                logStream << "it is a triangle with a vertex pair, the pair is not at v4 "
                          << "(because it is unique) but at "
                          << unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ]
                          << std::endl;
#endif
            }
        }
    }
#ifdef DB_DEBUG
    logStream << "we have now placed 4 vertices:\n"
              << "\t" << n_Unique_Vertices << " unique vertices" << std::endl;
    for ( int i = 0; i < n_Unique_Vertices; ++i )
    {
        logStream << "\tuV[" << i << "] : " << unique_Vertices[i]
                  << "\tDuplicates["<<i<<"] : " << n_Duplicates_at_Unique_Vertices[i] << std::endl;
    }
    if ( shape != UNDEFINED )
    {
        logStream << "we have also determined the shape is a " << getShapeString() << std::endl;
    }
#endif

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
        sortVertices();

#ifdef DB_DEBUG
        logStream << "before sorting :" << std::endl;
        for ( int i = 0; i < n_Unique_Vertices; ++i )
        {
            logStream << "\tuV[" << i << "] : "<< unique_Vertices[i]
                      << "\t duplicates[" << i << "] : " << n_Duplicates_at_Unique_Vertices[i]
                      << std::endl;
        }
        logStream << "after sorting :" << std::endl;
        for ( int i = 0; i < n_Unique_Vertices; ++i )
        {
            logStream << "\tsV[" << i << "] : "<< sorted_Vertices[i]
                      << "\t duplicates[" << i << "] : " << n_Duplicates_at_Sorted_Vertices[i]
                      << std::endl;
        }
#endif

        if ( shape == UNDEFINED )
        {
            /*
             * detect the shapes we didn't detect yet
             */
            if ( n_Unique_Vertices == 2 )
            {
                shape = LINE_WITH_TWO_VERTEX_PAIRS;
#ifdef DB_DEBUG
                logStream << "it is a line with two vertex pairs" << std::endl;
#endif
            }
            else if ( n_Duplicates_at_Sorted_Vertices[1] == 1 )
            {
                shape = LINE_WITH_VERTEX_PAIR_ON_LINE;
#ifdef DB_DEBUG
                logStream << "it is a line with a vertex pair on the line" << std::endl;
#endif
            }
            else
            {
                shape = LINE_WITH_VERTEX_PAIR_ON_END;
#ifdef DB_DEBUG
                logStream << "it is a line with a vertex pair on the end" << std::endl;
#endif
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
        if ( sorted_Vertices[0].DistanceSquared( sorted_Vertices[ n_Unique_Vertices-1 ] ) <
                MINIMUM_LINE_LENGTH_SQUARED )
        {
            shape = POINT_EQUIVALENT;
#ifdef DB_DEBUG
            logStream << "the line is too short, calculated the intersect to be at "
                      << intersect << " and rendering as point" << std::endl;
#endif
        }
    }

#ifdef DB_DEBUG
    logStream << "at the end of shape determination, it is a " << getShapeString() << std::endl;
#endif
#ifdef DB_COUNT
    ++shapeCounters[(int)shape];
#endif
}


void NPWVertices::sortVertices()
{
    switch ( n_Unique_Vertices )
    {

    /*
     * also sort for 1 and 2 unique vertices, because after calling this function,
     * other functions assume sorted_Vertices is set correctly
     */
    case 1:
    {
        sorted_Vertices[0] = unique_Vertices[0];
        n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
        break;
    }
    case 2:
    {
        sorted_Vertices[0] = unique_Vertices[0];
        sorted_Vertices[1] = unique_Vertices[1];
        n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
        n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[1];
        break;
    }

    case 3:
    {
        float d01 = unique_Vertices[0].DistanceSquared( unique_Vertices[1] );
        float d02 = unique_Vertices[0].DistanceSquared( unique_Vertices[2] );
        float d12 = unique_Vertices[1].DistanceSquared( unique_Vertices[2] );

        if ( d02 > d01 && d02 > d12 )
        {
            // 012
            sorted_Vertices[0] = unique_Vertices[0];
            sorted_Vertices[1] = unique_Vertices[1];
            sorted_Vertices[2] = unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[2];
        }
        else if ( d01 > d12 && d01 > d02 )
        {
            // 120
            sorted_Vertices[0] = unique_Vertices[1];
            sorted_Vertices[1] = unique_Vertices[2];
            sorted_Vertices[2] = unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[0];
        }
        else if ( d12 > d01 && d12 > d02 )
        {
            // 201
            sorted_Vertices[0] = unique_Vertices[2];
            sorted_Vertices[1] = unique_Vertices[0];
            sorted_Vertices[2] = unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
        }

        break;
    }

    case 4:
    {
        float d01 = unique_Vertices[0].DistanceSquared( unique_Vertices[1] );
        float d02 = unique_Vertices[0].DistanceSquared( unique_Vertices[2] );
        float d03 = unique_Vertices[0].DistanceSquared( unique_Vertices[3] );
        float d12 = unique_Vertices[1].DistanceSquared( unique_Vertices[2] );
        float d13 = unique_Vertices[1].DistanceSquared( unique_Vertices[3] );
        float d23 = unique_Vertices[2].DistanceSquared( unique_Vertices[3] );

        if ( d01 > d02 && d01 > d03 && d01 > d12 && d01 > d13 && d01 > d23 )
        {
            // 0xx1
            sorted_Vertices[0] = unique_Vertices[0];
            sorted_Vertices[3] = unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[1];

            if ( d03 > d02 )
            {
                // 0231
                sorted_Vertices[1] = unique_Vertices[2];
                sorted_Vertices[2] = unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[3];

            }
            else
            {
                // 0321
                sorted_Vertices[1] = unique_Vertices[3];
                sorted_Vertices[2] = unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[2];
            }
        }
        else if ( d02 > d01 && d02 > d03 && d02 > d12 && d02 > d13 && d02 > d23 )
        {
            // 0xx2
            sorted_Vertices[0] = unique_Vertices[0];
            sorted_Vertices[3] = unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[2];

            if ( d03 > d01 )
            {
                // 0132
                sorted_Vertices[1] = unique_Vertices[1];
                sorted_Vertices[2] = unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[3];
            }
            else
            {
                // 0312
                sorted_Vertices[1] = unique_Vertices[3];
                sorted_Vertices[2] = unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
            }
        }
        else if ( d03 > d01 && d03 > d02 && d03 > d12 && d03 > d13 && d03 > d23 )
        {
            // 0xx3
            sorted_Vertices[0] = unique_Vertices[0];
            sorted_Vertices[3] = unique_Vertices[3];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[0];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[3];


            if ( d02 > d01 )
            {
                // 0123
                sorted_Vertices[1] = unique_Vertices[1];
                sorted_Vertices[2] = unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[2];
            }
            else
            {
                // 0213
                sorted_Vertices[1] = unique_Vertices[2];
                sorted_Vertices[2] = unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
            }
        }
        else if ( d12 > d01 && d12 > d13 && d12 > d02 && d12 > d03 && d12 > d23 )
        {
            // 1xx2
            sorted_Vertices[0] = unique_Vertices[1];
            sorted_Vertices[3] = unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[2];


            if ( d13 > d01 )
            {
                // 1032
                sorted_Vertices[1] = unique_Vertices[0];
                sorted_Vertices[2] = unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[3];
            }
            else
            {
                // 1302
                sorted_Vertices[1] = unique_Vertices[3];
                sorted_Vertices[2] = unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[3];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[0];
            }
        }
        else if ( d13 > d01 && d13 > d02 && d13 > d03 && d13 > d12 && d13 > d23 )
        {
            // 1xx3
            sorted_Vertices[0] = unique_Vertices[1];
            sorted_Vertices[3] = unique_Vertices[3];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[1];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[3];


            if ( d12 > d01 )
            {
                // 1023
                sorted_Vertices[1] = unique_Vertices[0];
                sorted_Vertices[2] = unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[2];
            }
            else
            {
                // 1203
                sorted_Vertices[1] = unique_Vertices[2];
                sorted_Vertices[2] = unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[2];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[0];
            }
        }
        else if ( d23 > d01 && d23 > d02 && d23 > d03 && d23 > d12 && d23 > d13 )
        {
            // 2xx3
            sorted_Vertices[0] = unique_Vertices[2];
            sorted_Vertices[3] = unique_Vertices[3];
            n_Duplicates_at_Sorted_Vertices[0] = n_Duplicates_at_Unique_Vertices[2];
            n_Duplicates_at_Sorted_Vertices[3] = n_Duplicates_at_Unique_Vertices[3];


            if ( d12 > d02 )
            {
                // 2013
                sorted_Vertices[1] = unique_Vertices[0];
                sorted_Vertices[2] = unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[1];
            }
            else
            {
                // 2103
                sorted_Vertices[1] = unique_Vertices[1];
                sorted_Vertices[2] = unique_Vertices[0];
                n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Unique_Vertices[1];
                n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Unique_Vertices[0];
            }
        }

        break;
    }

    default:
        break;

    }
}

void NPWVertices::calculateIntersect()
{
    /*
     * very easy for line with pair on line, line with two pairs or line with trio.
     * for others, use formula described in report
     *
     * note: this is essentially the same formula:
     *   - get the x or y coordinate of the intersect depending on the gradient of the line
     *     formula Ix = ( ax bx - cx dx ) / ( ax - dx + bx - cx )      (same for y)
     *   - get the gradient of the line
     *   - get the other intersect coordinate by placing the obtained coordinate on the line
     */
    switch ( shape )
    {

    case LINE_WITH_TWO_VERTEX_PAIRS:
    {
        intersect = ( sorted_Vertices[0] + sorted_Vertices[1] ) / 2.0f;
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_LINE:
    {
        intersect = sorted_Vertices[1];
        break;
    }

    case LINE_WITH_VERTEX_TRIO:
    {
        intersect = n_Duplicates_at_Sorted_Vertices[0] == 2 ? sorted_Vertices[0] : sorted_Vertices[1];
        break;
    }

    case LINE_WITHOUT_VERTEX_PAIR:
    {
        if ( NPWGeometryMath::Abs( sorted_Vertices[0].x - sorted_Vertices[3].x ) >
             NPWGeometryMath::Abs( sorted_Vertices[0].y - sorted_Vertices[3].y ) )
        {
            intersect.x = (float)( (double)( sorted_Vertices[0].x * sorted_Vertices[1].x ) -
                                   (double)( sorted_Vertices[2].x * sorted_Vertices[3].x ) ) /
                                   (double)( sorted_Vertices[0].x - sorted_Vertices[3].x +
                                             sorted_Vertices[1].x - sorted_Vertices[2].x );
            double y_grad = (double)( sorted_Vertices[3].y - sorted_Vertices[0].y ) /
                            (double)( sorted_Vertices[3].x - sorted_Vertices[0].x );
            intersect.y = (float)( (double)( intersect.x - sorted_Vertices[0].x ) * y_grad +
                                   (double)sorted_Vertices[0].y );
        }
        else
        {
            intersect.y = (float)( (double)( sorted_Vertices[0].y * sorted_Vertices[1].y ) -
                                   (double)( sorted_Vertices[2].y * sorted_Vertices[3].y ) ) /
                                   (double)( sorted_Vertices[0].y - sorted_Vertices[3].y +
                                             sorted_Vertices[1].y - sorted_Vertices[2].y );
            double x_grad = (double)( sorted_Vertices[3].x - sorted_Vertices[0].x ) /
                            (double)( sorted_Vertices[3].y - sorted_Vertices[0].y );
            intersect.x = (float)( (double)( intersect.y - sorted_Vertices[0].y ) * x_grad +
                                   (double)sorted_Vertices[0].x );
        }
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_END:
    {
        /*
         * if pair is at sorted_Vertices[0]
         */
        if ( n_Duplicates_at_Sorted_Vertices[0] == 1 )
        {
            if ( NPWGeometryMath::Abs( sorted_Vertices[0].x - sorted_Vertices[2].x ) >
                 NPWGeometryMath::Abs( sorted_Vertices[0].y - sorted_Vertices[2].y ) )
            {
                intersect.x = (float)( (double)( sorted_Vertices[0].x * sorted_Vertices[0].x ) -
                                       (double)( sorted_Vertices[1].x * sorted_Vertices[2].x ) ) /
                                       (double)( sorted_Vertices[0].x * 2.0f -
                                                 sorted_Vertices[2].x  - sorted_Vertices[1].x );
                double y_grad = (double)( sorted_Vertices[2].y - sorted_Vertices[0].y ) /
                                (double)( sorted_Vertices[2].x - sorted_Vertices[0].x );
                intersect.y = (float)( (double)( intersect.x - sorted_Vertices[0].x ) * y_grad +
                                       (double)sorted_Vertices[0].y );
            }
            else
            {
                intersect.y = (float)( (double)( sorted_Vertices[0].y * sorted_Vertices[0].y ) -
                                       (double)( sorted_Vertices[1].y * sorted_Vertices[2].y ) ) /
                                       (double)( sorted_Vertices[0].y * 2.0f -
                                                 sorted_Vertices[2].y - sorted_Vertices[1].y );
                double x_grad = (double)( sorted_Vertices[2].x - sorted_Vertices[0].x ) /
                                (double)( sorted_Vertices[2].y - sorted_Vertices[0].y );
                intersect.x = (float)( (double)( intersect.y - sorted_Vertices[0].y ) * x_grad +
                                       (double)sorted_Vertices[0].x );
            }
        }
        /*
         * if pair is at sorted_Vertices[2]
         */
        else
        {
            if ( NPWGeometryMath::Abs( sorted_Vertices[0].x - sorted_Vertices[2].x ) >
                 NPWGeometryMath::Abs( sorted_Vertices[0].y - sorted_Vertices[2].y ) )
            {
                intersect.x = (float)( (double)( sorted_Vertices[0].x * sorted_Vertices[1].x ) -
                                       (double)( sorted_Vertices[2].x * sorted_Vertices[2].x ) ) /
                                       (double)( sorted_Vertices[0].x + sorted_Vertices[1].x -
                                                 sorted_Vertices[2].x * 2.0f );
                double y_grad = (double)( sorted_Vertices[2].y - sorted_Vertices[0].y ) /
                                (double)( sorted_Vertices[2].x - sorted_Vertices[0].x );
                intersect.y = (float)( (double)( intersect.x - sorted_Vertices[0].x ) * y_grad +
                                       (double)sorted_Vertices[0].y );
            }
            else
            {
                intersect.y = (float)( (double)( sorted_Vertices[0].y * sorted_Vertices[1].y ) -
                                       (double)( sorted_Vertices[2].y * sorted_Vertices[2].y ) ) /
                                       (double)( sorted_Vertices[0].y + sorted_Vertices[1].y -
                                                 sorted_Vertices[2].y * 2.0f );
                double x_grad = (double)( sorted_Vertices[2].x - sorted_Vertices[0].x ) /
                                (double)( sorted_Vertices[2].y - sorted_Vertices[0].y );
                intersect.x = (float)( (double)( intersect.y - sorted_Vertices[0].y ) * x_grad +
                                       (double)sorted_Vertices[0].x );
            }
        }

        break;
    }

    default:
        break;

    }

}


void NPWVertices::setRenderListWithoutShapeExpansion()
{
    switch ( shape )
    {

    case TRUE_POINT:
    {
        renderList[0] = &unique_Vertices[0];
        n_renderList_size = 1;
        break;
    }

    case POINT_EQUIVALENT:
    {
        // NOTE: we assume the intersect has been set
        renderList[0] = &intersect;
        n_renderList_size = 1;
        break;
    }

    case LINE_WITH_VERTEX_TRIO:
    {
        correctLineWithVertexTrio();
        break;
    }

    case LINE_WITH_TWO_VERTEX_PAIRS:
    {
        correctLineWithTwoVertexPairs();
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_END:
    {
        correctLineWithVertexPairOnEnd();
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_LINE:
    {
        correctLineWithVertexPairOnLine();
        break;
    }

    case LINE_WITHOUT_VERTEX_PAIR:
    {
        correctLineWithoutVertexPair();
        break;
    }

    case TRIANGLE_WITH_ONE_VERTEX_PAIR:
    {
        renderList[0] = &unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle      ];
        renderList[1] = &unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+1)%3 ];
        renderList[2] = &unique_Vertices[(index_of_maximum_intensity_vertex_in_triangle+2)%3 ];

        n_renderList_size = 3;

        break;
    }
    case TRIANGLE_WITH_ONE_VERTEX_ON_BORDER:
    {
        renderList[0] = &unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ];
        renderList[1] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ];
        renderList[2] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ];
        renderList[3] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ];

        n_renderList_size = 4;

        break;
    }

    case TRIANGLE_WITH_THREE_SUBTRIANGLES:
    {
        renderList[0] = &unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ];
        renderList[1] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ];
        renderList[2] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ];
        renderList[3] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ];
        renderList[4] = renderList[1];
        n_renderList_size = 5;

        break;
    }

    case QUAD:
    {
        getQuadIntersectAndSortVertices();

        renderList[0] = &intersect;
        renderList[1] = &sorted_Vertices[0];
        renderList[2] = &sorted_Vertices[1];
        renderList[3] = &sorted_Vertices[2];
        renderList[4] = &sorted_Vertices[3];
        renderList[5] = renderList[1];
        n_renderList_size = 6;

    }

    default:
        break;
    }
}


void NPWVertices::setRenderListWithShapeExpansion()
{
    switch ( shape )
    {

    case TRUE_POINT:
    {
        renderList[0] = &unique_Vertices[0];
        n_renderList_size = 1;
        break;
    }

    case POINT_EQUIVALENT:
    {
        // NOTE: we assume the intersect has been set
        renderList[0] = &intersect;
        n_renderList_size = 1;
        break;
    }

    case LINE_WITH_VERTEX_TRIO:
    {
        correctLineWithVertexTrio();
        break;
    }

    case LINE_WITH_TWO_VERTEX_PAIRS:
    {
        correctLineWithTwoVertexPairs();
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_END:
    {
        correctLineWithVertexPairOnEnd();
        break;
    }

    case LINE_WITH_VERTEX_PAIR_ON_LINE:
    {
        correctLineWithVertexPairOnLine();
        break;
    }

    case LINE_WITHOUT_VERTEX_PAIR:
    {
        correctLineWithoutVertexPair();
        break;
    }

    case TRIANGLE_WITH_ONE_VERTEX_PAIR:
    {
        intersect     = unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ];

        renderList[0] = &intersect;
        renderList[1] = &unique_Vertices[0];
        renderList[2] = &unique_Vertices[1];
        renderList[3] = &unique_Vertices[2];
        renderList[4] = renderList[1];

        n_renderList_size = 5;

        if ( !expandTriangle() )
        {
#ifdef DB_DEBUG
            logStream << "could not expand triangle, putting three points on a line..." << std::endl;
#endif
            putThreePointsOnALineAndSetRenderList();
        }

        break;
    }
    case TRIANGLE_WITH_ONE_VERTEX_ON_BORDER:
    {
        renderList[0] = &unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ];
        renderList[1] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ];
        renderList[2] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ];
        renderList[3] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ];
        renderList[4] = renderList[1];

        n_renderList_size = 5;

        if ( !expandTriangle() )
        {
#ifdef DB_DEBUG
            logStream << "could not expand triangle, putting four points on a line..." << std::endl;
#endif
            putFourPointsOnALineAndSetRenderList();
        }
        break;
    }

    case TRIANGLE_WITH_THREE_SUBTRIANGLES:
    {
        renderList[0] = &unique_Vertices[  index_of_maximum_intensity_vertex_in_triangle ];
        renderList[1] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+1)%4 ];
        renderList[2] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+2)%4 ];
        renderList[3] = &unique_Vertices[ (index_of_maximum_intensity_vertex_in_triangle+3)%4 ];
        renderList[4] = renderList[1];

        n_renderList_size = 5;

        if ( !checkIfRenderListIsRenderable() )
        {
#ifdef DB_DEBUG
            logStream << "triangle with three subtriangles is not renderable, expanding triangle"
                      << std::endl;
#endif
            if ( !expandTriangle() )
            {
#ifdef DB_DEBUG
            logStream << "could not expand triangle, putting four points on a line..." << std::endl;
#endif
                putFourPointsOnALineAndSetRenderList();
            }
        }

        break;
    }

    case QUAD:
    {
        getQuadIntersectAndSortVertices();

        renderList[0] = &intersect;
        renderList[1] = &sorted_Vertices[0];
        renderList[2] = &sorted_Vertices[1];
        renderList[3] = &sorted_Vertices[2];
        renderList[4] = &sorted_Vertices[3];
        renderList[5] = renderList[1];
        n_renderList_size = 6;

        if ( !checkIfRenderListIsRenderable() )
        {
#ifdef DB_DEBUG
            logStream << "quad is not renderable, expanding quad" << std::endl;
#endif
            if ( !expandQuad() )
            {
#ifdef DB_DEBUG
            logStream << "could not expand quad, putting four points on a line..." << std::endl;
#endif
                putFourPointsOnALineAndSetRenderList();
            }
        }
        break;
    }

    default:
        break;
    }

}

#ifdef DB_DEBUG
void NPWVertices::checkRenderList()
{
    logStream << "checking render list validity..." << std::endl;
    for ( unsigned int j = 0; j < n_renderList_size; ++j )
    {
        logStream << "\trenderList["<<j<<"] : " << *(renderList[j]) << std::endl;
    }
    if ( n_renderList_size == 5 )
    {
        /*
         * is centre in triangle?
         */
        if ( !NPWGeometryMath::isVertexInTriangle(
                renderList[1], renderList[2],
                renderList[3], renderList[0] ) )
        {
            logStream << "*** ERROR *** fan centre is not within triangle!" << std::endl;
            isErrored = true;
        }
    }
    else if ( n_renderList_size == 6 )
    {
        /*
         * is centre in quad?
         */
        for ( unsigned int i = 2; i < n_renderList_size; ++i )
        {
            unsigned int index = (i+1 > 5) ? 2 : i+1;
            if ( !NPWGeometryMath::arePointsOnSameSideOfLine(
                    *(renderList[i-1]), *(renderList[i]),
                    *(renderList[0])  , *(renderList[index]) ) )
            {
                logStream << "*** ERROR *** fan center is not within quad!" << std::endl;
                logStream << "   fan center " << *(renderList[0])
                          << " is not on the same side of "
                          << *(renderList[i-1]) << " <-> "
                          << *(renderList[i])
                          << " as " << *(renderList[index]) << std::endl;
                isErrored = true;
            }
        }
    }

    if ( isErrored )
    {
        std::cout << logStream.str() << std::endl;
    }
}
#endif


void NPWVertices::correctLineWithoutVertexPair()
{
    /*
     * calculate the intersect
     */
    calculateIntersect();

    /*
     * ensure the two vertices on the line are at least MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_2
     * from their closest line end vertex
     */
    float distanceSquared = sorted_Vertices[0].DistanceSquared( sorted_Vertices[1] );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
    {
        move = sorted_Vertices[0] - sorted_Vertices[3];
        move = move.Normalise() * ( MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO -
                                    std::sqrt(distanceSquared) );
#ifdef DB_DEBUG
         logStream << "moving " << sorted_Vertices[0] << " to " << sorted_Vertices[0]+move << std::endl;
#endif
        sorted_Vertices[0] += move;
    }

    distanceSquared = sorted_Vertices[3].DistanceSquared( sorted_Vertices[2] );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
    {
        move = sorted_Vertices[3] - sorted_Vertices[0];
        move = move.Normalise() * ( MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO -
                                    std::sqrt(distanceSquared) );
#ifdef DB_DEBUG
         logStream << "moving " << sorted_Vertices[3] << " to " << sorted_Vertices[3]+move << std::endl;
#endif
        sorted_Vertices[3] += move;
    }

    /*
     * let the points on the line be ordered A B C D
     * then move B perpendicular to the line over a distance of
     *     BD * tan ( arcsin( MINIMUM_EDGE_DISTANCE / CD ) )
     *
     * and move C perpendicular to the line over a distance of
     *     AC * tan ( arcsin( MINIMUM_EDGE_DISTANCE / AB ) )
     *
     * but don't move a point more than the maximum distance
     * even though the distance from the intersect to the line will no be
     * MINIMUM_EDGE_DISTANCE, it will be close
     *
     * NOTE: this is all very expensive and not worth it,
     *       just move over MAXIMUM_VERTEX_CORRECTION
     */
    move.x = sorted_Vertices[0].y - sorted_Vertices[3].y;
    move.y = sorted_Vertices[3].x - sorted_Vertices[0].x;
    move   = move.Normalise();

    /*
     * we'd like to keep the code though:
     */
//    float AB = sorted_Vertices[0].Distance( sorted_Vertices[1] );
//    float BC = sorted_Vertices[1].Distance( sorted_Vertices[2] );
//    float CD = sorted_Vertices[2].Distance( sorted_Vertices[3] );
//    float AC = AB + BC;
//    float BD = BC + CD;
//
//    float distance1 = NPWGeometryMath::Min(
//            BD * std::tan( std::asin( MINIMUM_EDGE_DISTANCE / CD ) ), MAXIMUM_VERTEX_CORRECTION );
//    float distance2 = NPWGeometryMath::Min(
//            AC * std::tan( std::asin( MINIMUM_EDGE_DISTANCE / AB ) ), MAXIMUM_VERTEX_CORRECTION );
//
//    sorted_Vertices[1] += (move * distance1);
//    sorted_Vertices[2] -= (move * distance2);
    // END OLD CODE

    sorted_Vertices[1] += (move * MAXIMUM_VERTEX_CORRECTION);
    sorted_Vertices[2] -= (move * MAXIMUM_VERTEX_CORRECTION);

    /*
     * set render list and shape
     */
    renderList[0] = &intersect;
    renderList[1] = &sorted_Vertices[1];
    renderList[2] = &sorted_Vertices[0];
    renderList[3] = &sorted_Vertices[2];
    renderList[4] = &sorted_Vertices[3];
    renderList[5] = renderList[1];

    n_renderList_size = 6;

    shape = QUAD;
}


void NPWVertices::correctLineWithVertexPairOnEnd()
{
    /*
     * calculate the intersect
     */
    calculateIntersect();

    /*
     * determine which vertex is the pair
     */
    Vector2 * pairEnd   = &sorted_Vertices[0];
    Vector2 * singleEnd = &sorted_Vertices[2];
    if ( n_Duplicates_at_Sorted_Vertices[0] != 1 )
    {
        pairEnd       = &sorted_Vertices[2];
        singleEnd     = &sorted_Vertices[0];
    }

    /*
     * ensure the intersect is MINIMUM_EDGE_DISTANCE from the pair end
     */
    float distanceSquared = intersect.DistanceSquared( *pairEnd );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED )
    {
        move = *pairEnd - *singleEnd;
        move = move.Normalise() * ( MINIMUM_EDGE_DISTANCE - std::sqrt(distanceSquared) );

        *pairEnd += move;
    }
    /*
     * ensure the intersect is MINIMUM_EDGE_DISTANCE_TIMES_TWO from the single end
     */
    distanceSquared = intersect.DistanceSquared( *singleEnd );
    if ( distanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_FOUR )
    {
        move = *singleEnd - *pairEnd;
        move = move.Normalise() * (MINIMUM_EDGE_DISTANCE_TIMES_TWO - std::sqrt(distanceSquared) );

        *singleEnd += move;
    }

    /*
     * move the pair perpendicular to the line with MINIMUM_EDGE_DISTANCE_TIMES_TWO
     */
    move.x = singleEnd->y - pairEnd->y;
    move.y = pairEnd->x   - singleEnd->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE_TIMES_TWO;

    sorted_Vertices[3] = *pairEnd - move;
    *pairEnd += move;

    /*
     * set render list and shape
     */
    renderList[0] = &intersect;
    renderList[1] = pairEnd;
    renderList[2] = &sorted_Vertices[3];
    renderList[3] = singleEnd;
    renderList[4] = renderList[1];
    n_renderList_size = 5;

    shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
}


void NPWVertices::correctLineWithVertexPairOnLine()
{
    /*
     * calculate the intersect
     */
    calculateIntersect();

#ifdef DB_DEBUG
    logStream << "intersect : " << intersect << std::endl;
#endif

    /*
      * check if the end points are a minimum distance of
      * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_2 from the pair
      */
     if ( sorted_Vertices[0].DistanceSquared(sorted_Vertices[1])
             < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
     {
         move = (sorted_Vertices[0] - sorted_Vertices[2]).Normalise() * MINIMUM_EDGE_DISTANCE;
#ifdef DB_DEBUG
         logStream << "moving " << sorted_Vertices[0] << " to " << sorted_Vertices[0]+move << std::endl;
#endif
         sorted_Vertices[0] += move;
     }
     if ( sorted_Vertices[2].DistanceSquared(sorted_Vertices[1])
             < MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO )
     {
         move = (sorted_Vertices[2] - sorted_Vertices[0]).Normalise() * MINIMUM_EDGE_DISTANCE;
#ifdef DB_DEBUG
         logStream << "moving " << sorted_Vertices[2] << " to " << sorted_Vertices[2]+move << std::endl;
#endif
         sorted_Vertices[2] += move;
     }

     /*
      * move vertex pair
      */
     move.x = sorted_Vertices[0].y - sorted_Vertices[2].y;
     move.y = sorted_Vertices[2].x - sorted_Vertices[0].x;
     move = move.Normalise() * MINIMUM_EDGE_DISTANCE_TIMES_SQRT_OF_TWO;

     sorted_Vertices[3]  = sorted_Vertices[1] + move;
     sorted_Vertices[1] -= move;

#ifdef DB_DEBUG
         logStream << "two points moved from the pair on the line:\n"
                   << "\t" << sorted_Vertices[3]
                   << "\t and " << sorted_Vertices[1] << std::endl;
#endif

     /*
      * set renderlist & shape
      */
     renderList[0] = &intersect;
     renderList[1] = &sorted_Vertices[1];
     renderList[2] = &sorted_Vertices[0];
     renderList[3] = &sorted_Vertices[3];
     renderList[4] = &sorted_Vertices[2];
     renderList[5] = renderList[1];

     n_renderList_size = 6;

     shape = QUAD;
}


void NPWVertices::correctLineWithTwoVertexPairs()
{
    /*
     * calculate the intersect
     */
    calculateIntersect();

    /*
     * move pairs
     */
    move.x = sorted_Vertices[0].y - sorted_Vertices[1].y;
    move.y = sorted_Vertices[1].x - sorted_Vertices[0].x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;

    sorted_Vertices[2]  = sorted_Vertices[1] + move;
    sorted_Vertices[3]  = sorted_Vertices[0] + move;
    sorted_Vertices[0] -= move;
    sorted_Vertices[1] -= move;

    /*
     * set render list and shape
     */
    renderList[0] = &intersect;
    renderList[1] = &sorted_Vertices[0];
    renderList[2] = &sorted_Vertices[1];
    renderList[3] = &sorted_Vertices[2];
    renderList[4] = &sorted_Vertices[3];
    renderList[5] = renderList[1];

    n_renderList_size = 6;

    shape = QUAD;
}


void NPWVertices::correctLineWithVertexTrio()
{
    /*
     * determine what vertex is the trio
     */
    Vector2 * singleVertex = &sorted_Vertices[0];
    Vector2 * vertexTrio   = &sorted_Vertices[1];
    if ( n_Duplicates_at_Unique_Vertices[0] == 2)
    {
        singleVertex  = &sorted_Vertices[1];
        vertexTrio    = &sorted_Vertices[0];
    }

    /*
     * move two vertices of the trio away from the line at 120 and 240 degrees
     */
    move  = (*singleVertex - *vertexTrio).Normalise();
    move  = Vector2(move.x * COSINUS_OF_120_DEGREES - move.y * SINUS_OF_120_DEGREES,
                    move.x * SINUS_OF_120_DEGREES   + move.y * COSINUS_OF_120_DEGREES)
            * MINIMUM_EDGE_DISTANCE_TIMES_TWO;
    sorted_Vertices[2] = *vertexTrio + move;

    move  = (*singleVertex - *vertexTrio).Normalise();
    move = Vector2(move.x*COSINUS_OF_240_DEGREES - move.y*SINUS_OF_240_DEGREES,
                   move.x*SINUS_OF_240_DEGREES   + move.y*COSINUS_OF_240_DEGREES);
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE_TIMES_TWO;
    sorted_Vertices[3] = *vertexTrio + move;

    /*
     * set render list and shape
     */
    renderList[0] = vertexTrio;
    renderList[1] = singleVertex;
    renderList[2] = &sorted_Vertices[2];
    renderList[3] = &sorted_Vertices[3];
    renderList[4] = renderList[1];

    n_renderList_size = 5;

    shape = TRIANGLE_WITH_THREE_SUBTRIANGLES;
}


bool NPWVertices::checkIfRenderListIsRenderable()
{
    /*
     * Check the distances from the center vertex to the edges
     */
    float minimumDistanceSquared = MINIMUM_EDGE_DISTANCE_SQUARED_TIMES_TWO;
    for ( unsigned int i = 2; i < n_renderList_size; ++i )
    {
        // check line renderList[i-1] -> renderList[i]
        Vector2 closestPoint;
        float unused;
        NPWGeometryMath::ClosestPointOnLine(*renderList[i-1], *renderList[i], *renderList[0],
                                            closestPoint, unused, true);
        float distanceSquared = renderList[0]->DistanceSquared( closestPoint );

        minimumDistanceSquared = NPWGeometryMath::Min( minimumDistanceSquared, distanceSquared );

        if ( minimumDistanceSquared < MINIMUM_EDGE_DISTANCE_SQUARED )
        {
            return false;
        }
    }

    return true;
}


bool NPWVertices::expandTriangle()
{
    /*
     * the triangle is renderList[1], renderList[2], renderList[3]
     * let's call these a, b and c for short
     */

    /*
     * extend line ab
     */
    move.x = renderList[1]->y - renderList[2]->y;
    move.y = renderList[2]->x - renderList[1]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 ab1 = *renderList[1] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(
            *renderList[1], *renderList[2], ab1, *renderList[3]) )
    {
        move *= -1.0f;
        ab1 = *renderList[1] + move;
    }
    Vector2 ab2 = *renderList[2] + move;

    /*
     * extend line bc
     */
    move.x = renderList[2]->y - renderList[3]->y;
    move.y = renderList[3]->x - renderList[2]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 bc1 = *renderList[2] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(
            *renderList[2], *renderList[3], bc1, *renderList[1]) )
    {
        move *= -1.0f;
        bc1 = *renderList[2] + move;
    }
    Vector2 bc2 = *renderList[3] + move;

    /*
     * extend line ca
     */
    move.x = renderList[3]->y - renderList[1]->y;
    move.y = renderList[1]->x - renderList[3]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 ca1 = *renderList[3] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(
            *renderList[3], *renderList[1], ca1, *renderList[2]) )
    {
        move *= -1.0f;
        ca1 = *renderList[3] + move;
    }
    Vector2 ca2 = *renderList[1] + move;

    /*
     * get the three intersections of the lines
     */
    bool unused;
    Vector2 a_new, b_new, c_new;

    NPWGeometryMath::LineLineIntersect2D( ab1,ab2, ca1,ca2, a_new, unused );
    NPWGeometryMath::LineLineIntersect2D( bc1,bc2, ab1,ab2, b_new, unused );
    NPWGeometryMath::LineLineIntersect2D( ca1,ca2, bc1,bc2, c_new, unused );

    /*
     * check if any point is moved over too large a distance
     */
    if ( a_new.DistanceSquared(*renderList[1]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         b_new.DistanceSquared(*renderList[2]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         c_new.DistanceSquared(*renderList[3]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED )
    {
#ifdef DB_COUNT
        ++triangleExpansionFails;
#endif
        return false;
    }

    /*
     * check if any of the new vertices is out of bounds
     */
    float maxX = (float)( parentWindow->getFrameBufferWidth()  -
                          parentWindow->getFrameBufferBorderSize() );
    float maxY = (float)( parentWindow->getFrameBufferHeight() -
                          parentWindow->getFrameBufferBorderSize() );
    float minX = -(float) parentWindow->getFrameBufferBorderSize();
    float minY = -(float) parentWindow->getFrameBufferBorderSize();
    if ( ( a_new.x < minX || a_new.x > maxX || a_new.y < minY || a_new.y > maxY ) ||
         ( b_new.x < minX || b_new.x > maxX || b_new.y < minY || b_new.y > maxY ) ||
         ( c_new.x < minX || c_new.x > maxX || c_new.y < minY || c_new.y > maxY ) )
    {
        // at least one of the new vertices is out of bounds, cannot expand the triangle
#ifdef DB_COUNT
        ++triangleExpansionFails;
#endif
        return false;
    }

    /*
     * set render list
     */
    *renderList[1] = a_new;
    *renderList[2] = b_new;
    *renderList[3] = c_new;
    // renderList[4] == renderList[1]

    return true; // success!

}


bool NPWVertices::expandQuad()
{
    /*
     * the quad is renderList[1], renderList[2], renderList[3], renderList[4]
     * let's call these a, b, c and d for short
     */

    /*
     * extend line ab
     */
    move.x = renderList[1]->y - renderList[2]->y;
    move.y = renderList[2]->x - renderList[1]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 ab1 = *renderList[1] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(*renderList[1], *renderList[2], ab1, intersect) )
    {
        move *= -1.0f;
        ab1 = *renderList[1] + move;
    }
    Vector2 ab2 = *renderList[2] + move;

    /*
     * extend line bc
     */
    move.x = renderList[2]->y - renderList[3]->y;
    move.y = renderList[3]->x - renderList[2]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 bc1 = *renderList[2] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(*renderList[2], *renderList[3], bc1, intersect) )
    {
        move *= -1.0f;
        bc1 = *renderList[2] + move;
    }
    Vector2 bc2 = *renderList[3] + move;

    /*
     * extend line cd
     */
    move.x = renderList[3]->y - renderList[4]->y;
    move.y = renderList[4]->x - renderList[3]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 cd1 = *renderList[3] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(*renderList[3], *renderList[4], cd1, intersect) )
    {
        move *= -1.0f;
        cd1 = *renderList[3] + move;
    }
    Vector2 cd2 = *renderList[4] + move;

    /*
     * extend line da
     */
    move.x = renderList[4]->y - renderList[1]->y;
    move.y = renderList[1]->x - renderList[4]->x;
    move = move.Normalise() * MINIMUM_EDGE_DISTANCE;
    Vector2 da1 = *renderList[4] + move;

    // if ab1 and c are not on the same side of ab, move is pointing in the wrong direction
    if ( NPWGeometryMath::arePointsOnSameSideOfLine(*renderList[4], *renderList[1], da1, intersect) )
    {
        move *= -1.0f;
        da1 = *renderList[4] + move;
    }
    Vector2 da2 = *renderList[1] + move;


    /*
     * get the three intersections of the line
     * we cannot pass the pointers directly, since they are const Vector2 pointers
     */
    bool unused;
    Vector2 a_new, b_new, c_new, d_new;

    NPWGeometryMath::LineLineIntersect2D( ab1,ab2, da1,da2, a_new, unused );
    NPWGeometryMath::LineLineIntersect2D( bc1,bc2, ab1,ab2, b_new, unused );
    NPWGeometryMath::LineLineIntersect2D( cd1,cd2, bc1,bc2, c_new, unused );
    NPWGeometryMath::LineLineIntersect2D( da1,da2, cd1,cd2, d_new, unused );

    /*
     * check if any of the new vertices is out of bounds
     */
    float maxX = (float)( parentWindow->getFrameBufferWidth()  -
                          parentWindow->getFrameBufferBorderSize() );
    float maxY = (float)( parentWindow->getFrameBufferHeight() -
                          parentWindow->getFrameBufferBorderSize() );
    float minX = -(float) parentWindow->getFrameBufferBorderSize();
    float minY = -(float) parentWindow->getFrameBufferBorderSize();
    if ( ( a_new.x < minX || a_new.x > maxX || a_new.y < minY || a_new.y > maxY ) ||
         ( b_new.x < minX || b_new.x > maxX || b_new.y < minY || b_new.y > maxY ) ||
         ( c_new.x < minX || c_new.x > maxX || c_new.y < minY || c_new.y > maxY ) ||
         ( d_new.x < minX || d_new.x > maxX || d_new.y < minY || d_new.y > maxY ) )
    {
        // at least one of the new vertices is out of bounds, cannot expand the triangle
#ifdef DB_COUNT
        ++quadExpansionFails;
#endif
        return false;
    }

    /*
     * check if any point is moved over too large a distance
     */
    if ( a_new.DistanceSquared(*renderList[1]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         b_new.DistanceSquared(*renderList[2]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         c_new.DistanceSquared(*renderList[3]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED ||
         d_new.DistanceSquared(*renderList[4]) > MAXIMUM_VERTEX_MOVEMENT_SQUARED )
    {
#ifdef DB_COUNT
        ++quadExpansionFails;
#endif
        return false;
    }

    /*
     * set renderlist
     */
    *renderList[1] = a_new;
    *renderList[2] = b_new;
    *renderList[3] = c_new;
    *renderList[4] = d_new;
    // renderList[5] == renderList[1]

    return true;
}


void NPWVertices::putThreePointsOnALineAndSetRenderList()
{
    /*
     * get the longest line
     */
    float distanceSquared;
    NPWGeometryMath::getLongestLineOfThreePoints(
            unique_Vertices[0], unique_Vertices[1], unique_Vertices[2],
            sorted_Vertices[0], sorted_Vertices[2], sorted_Vertices[1],
            distanceSquared);
    /*
     * figure out where the pair went
     */
    if ( unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ] == sorted_Vertices[0] )
    {
        n_Duplicates_at_Sorted_Vertices[0] = 1;
        n_Duplicates_at_Sorted_Vertices[1] = 0;
        n_Duplicates_at_Sorted_Vertices[2] = 0;
    }
    else if ( unique_Vertices[ index_of_maximum_intensity_vertex_in_triangle ] == sorted_Vertices[1] )
    {
        n_Duplicates_at_Sorted_Vertices[0] = 0;
        n_Duplicates_at_Sorted_Vertices[1] = 1;
        n_Duplicates_at_Sorted_Vertices[2] = 0;
    }
    else
    {
        n_Duplicates_at_Sorted_Vertices[0] = 0;
        n_Duplicates_at_Sorted_Vertices[1] = 0;
        n_Duplicates_at_Sorted_Vertices[2] = 1;
    }

    n_Unique_Vertices = 2;

    /*
     * project other vertex onto the line
     */
    float unused;
    Vector2 projected_new;
    NPWGeometryMath::ClosestPointOnLine( sorted_Vertices[0], sorted_Vertices[2],
            sorted_Vertices[1], projected_new, unused, false );

    /*
     * are the line segment end vertices far enough apart?
     */
    if ( sorted_Vertices[0].DistanceSquared( sorted_Vertices[3] ) < MINIMUM_LINE_LENGTH_SQUARED )
    {
        intersect = ( sorted_Vertices[0] + sorted_Vertices[3] ) / 2.0f;

        shape = POINT_EQUIVALENT;

        renderList[0] = &intersect;
        n_renderList_size = 1;

        return;
    }

    /*
     * is the projected vertex close to an end vertex?
     */
    if ( projected_new.DistanceSquared( sorted_Vertices[0] ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Sorted_Vertices[0];
    }
    else if ( projected_new.DistanceSquared( sorted_Vertices[2] ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Sorted_Vertices[2];
    }
    else
    {
        sorted_Vertices[1] = projected_new;
        ++n_Unique_Vertices;
    }

     /*
      * the end vertex may be at the wrong place in sorted_vertices
      */
     if ( n_Unique_Vertices == 2 )
     {
         sorted_Vertices[2] = sorted_Vertices[3];
         n_Duplicates_at_Sorted_Vertices[2] = n_Duplicates_at_Sorted_Vertices[3];
     }


     /*
      * we have all vertices sorted and guaranteed they're all far enough apart
      * now we just need to determine the line type and call the correct function
      */
     if ( n_Unique_Vertices == 2 )
     {
         if ( n_Duplicates_at_Sorted_Vertices[0] == 1 && n_Duplicates_at_Sorted_Vertices[1] == 1 )
         {
             shape = LINE_WITH_TWO_VERTEX_PAIRS;
             correctLineWithTwoVertexPairs();
         }
         else
         {
             shape = LINE_WITH_VERTEX_TRIO;
             correctLineWithVertexTrio();
         }
     }
     else
     {
         if ( n_Duplicates_at_Sorted_Vertices[1] == 1 )
         {
             shape = LINE_WITH_VERTEX_PAIR_ON_LINE;
             correctLineWithVertexPairOnLine();
         }
         else
         {
             shape = LINE_WITH_VERTEX_PAIR_ON_END;
             correctLineWithVertexPairOnEnd();
         }
     }

}

void NPWVertices::putFourPointsOnALineAndSetRenderList()
{
    /*
     * get the longest line
     */
    float distanceSquared;
    NPWGeometryMath::getLongestLineOfFourPoints(
            unique_Vertices[0], unique_Vertices[1], unique_Vertices[2], unique_Vertices[3],
            sorted_Vertices[0], sorted_Vertices[3], sorted_Vertices[1], sorted_Vertices[2],
            distanceSquared);
    n_Duplicates_at_Sorted_Vertices[0] = 0;
    n_Duplicates_at_Sorted_Vertices[3] = 0;
    n_Unique_Vertices = 2;

#ifdef DB_DEBUG
    logStream << "longest line is " << sorted_Vertices[0] << " <-> " << sorted_Vertices[3] << "\n"
              << "\tother 1 : " << sorted_Vertices[1] << "\n"
              << "\tother 2 : " << sorted_Vertices[2] << std::endl;
#endif

    /*
     * project other vertices onto the line
     */
    float unused;
    Vector2 projected1_new, projected2_new;
    NPWGeometryMath::ClosestPointOnLine( sorted_Vertices[0], sorted_Vertices[3],
            sorted_Vertices[1], projected1_new, unused, true );
    NPWGeometryMath::ClosestPointOnLine( sorted_Vertices[0], sorted_Vertices[3],
            sorted_Vertices[2], projected2_new, unused, true );

#ifdef DB_DEBUG
    logStream << "projected 1 : " << projected1_new << "\n"
              << "projected 2 : " << projected2_new << std::endl;;
#endif

    /*
     * are the line segment end vertices far enough apart?
     */
    if ( sorted_Vertices[0].DistanceSquared( sorted_Vertices[3] ) < MINIMUM_LINE_LENGTH_SQUARED )
    {
        intersect = ( sorted_Vertices[0] + sorted_Vertices[3] ) / 2.0f;

        shape = POINT_EQUIVALENT;

        renderList[0] = &intersect;
        n_renderList_size = 1;

#ifdef DB_DEBUG
        logStream << "longest line is not long enough, rendering center of line as point : "
                  << intersect << std::endl;
#endif

        return;
    }

    /*
     * is the first projected vertex close to an end vertex?
     */
    if ( projected1_new.DistanceSquared( sorted_Vertices[0] ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Sorted_Vertices[0];
#ifdef DB_DEBUG
        logStream << "projected 1 is close to line begin " << sorted_Vertices[0] << std::endl;
#endif
    }
    else if ( projected1_new.DistanceSquared( sorted_Vertices[3] ) < MINIMUM_POINT_DISTANCE_SQUARED )
    {
        ++n_Duplicates_at_Sorted_Vertices[3];
#ifdef DB_DEBUG
        logStream << "projected 1 is close to line end " << sorted_Vertices[3] << std::endl;
#endif
    }
    else
    {
        sorted_Vertices[1] = projected1_new;
        n_Duplicates_at_Sorted_Vertices[1] = 0;
        ++n_Unique_Vertices;
#ifdef DB_DEBUG
        logStream << "projected 1 is a unique vertex" << std::endl;
#endif
    }

    /*
     * is the second projected vertex close to an end vertex or the other projected vertex?
     */
     if ( projected2_new.DistanceSquared( sorted_Vertices[0] ) < MINIMUM_POINT_DISTANCE_SQUARED )
     {
         ++n_Duplicates_at_Sorted_Vertices[0];
#ifdef DB_DEBUG
        logStream << "projected 2 is close to line begin " << sorted_Vertices[0] << std::endl;
#endif
     }
     else if ( projected2_new.DistanceSquared( sorted_Vertices[3] ) < MINIMUM_POINT_DISTANCE_SQUARED )
     {
         ++n_Duplicates_at_Sorted_Vertices[3];
#ifdef DB_DEBUG
        logStream << "projected 2 is close to line end " << sorted_Vertices[1] << std::endl;
#endif
     }
     else if ( n_Unique_Vertices == 3 &&
               projected2_new.DistanceSquared( projected1_new ) < MINIMUM_POINT_DISTANCE_SQUARED)
     {
         ++n_Duplicates_at_Sorted_Vertices[1];
#ifdef DB_DEBUG
        logStream << "projected 2 is close to projected 1 " << projected1_new << std::endl;
#endif
     }
     else
     {
         sorted_Vertices[n_Unique_Vertices-1] = projected2_new;
         n_Duplicates_at_Sorted_Vertices[n_Unique_Vertices-1] = 0;
         ++n_Unique_Vertices;
#ifdef DB_DEBUG
        logStream << "projected 2 is a unique vertex" << std::endl;
#endif
     }

     /*
      * the end vertex may be at the wrong place in sorted_vertices
      */
     if ( n_Unique_Vertices != 4 )
     {
         sorted_Vertices[ n_Unique_Vertices-1 ] = sorted_Vertices[3];
         n_Duplicates_at_Sorted_Vertices[ n_Unique_Vertices-1 ] = n_Duplicates_at_Sorted_Vertices[3];
     }
     /*
      * projected 1 and 2 may be in the wrong order
      */
     else if ( sorted_Vertices[0].DistanceSquared( sorted_Vertices[1] ) >
               sorted_Vertices[0].DistanceSquared( sorted_Vertices[2] ) )
     {
         Vector2 tmpVertex                  = sorted_Vertices[1];
         int tmpDuplicates                  = n_Duplicates_at_Sorted_Vertices[1];
         sorted_Vertices[1]                 = sorted_Vertices[2];
         n_Duplicates_at_Sorted_Vertices[1] = n_Duplicates_at_Sorted_Vertices[2];
         sorted_Vertices[2]                 = tmpVertex;
         n_Duplicates_at_Sorted_Vertices[2] = tmpDuplicates;
     }

#ifdef DB_DEBUG
     logStream << "after putting four points on a line, the sorted vertices are : " << std::endl;
     for ( int i = 0; i < n_Unique_Vertices; ++i )
     {
         logStream << "\tsV[" << i << "] : " << sorted_Vertices[i]
                   << "\tduplicates[" << i << "] : " << n_Duplicates_at_Sorted_Vertices[i] << std::endl;
     }
#endif

     /*
      * we have sorted all vertices and guaranteed they're all far enough apart
      * now we just need to determine the line type and call the correct function
      */
     if ( n_Unique_Vertices == 2 )
     {
         if ( n_Duplicates_at_Sorted_Vertices[0] == 1 && n_Duplicates_at_Sorted_Vertices[1] == 1 )
         {
             shape = LINE_WITH_TWO_VERTEX_PAIRS;
#ifdef DB_DEBUG
             logStream << "the shape is now a line with two vertex pairs" << std::endl;
#endif
             correctLineWithTwoVertexPairs();
         }
         else
         {
             shape = LINE_WITH_VERTEX_TRIO;
#ifdef DB_DEBUG
             logStream << "the shape is now a line with a vertex trio" << std::endl;
#endif
             correctLineWithVertexTrio();
         }
     }
     else if ( n_Unique_Vertices == 3 )
     {
         if ( n_Duplicates_at_Sorted_Vertices[1] == 1 )
         {
             shape = LINE_WITH_VERTEX_PAIR_ON_LINE;
#ifdef DB_DEBUG
             logStream << "the shape is now a line with a vertex pair on the line" << std::endl;
#endif
             correctLineWithVertexPairOnLine();
         }
         else
         {
             shape = LINE_WITH_VERTEX_PAIR_ON_END;
#ifdef DB_DEBUG
             logStream << "the shape is now a line with a vertex pair on an end" << std::endl;
#endif
             correctLineWithVertexPairOnEnd();
         }
     }
     else
     {
         shape = LINE_WITHOUT_VERTEX_PAIR;
#ifdef DB_DEBUG
             logStream << "the shape is now a line without vertex pairs" << std::endl;
#endif
         correctLineWithoutVertexPair();
     }

}


void NPWVertices::getQuadIntersectAndSortVertices()
{
    bool intersectIsOnBothLines;
    bool solutionExists;

    // TODO: what's faster: if () else if () else if ()  OR   if () return if () return if () return

    // get the intersection between 01 and 23
    solutionExists = NPWGeometryMath::LineLineIntersect2D(
            unique_Vertices[0], unique_Vertices[1], unique_Vertices[2], unique_Vertices[3],
            intersect, intersectIsOnBothLines);
    if ( solutionExists && intersectIsOnBothLines ) // quad: 0213
    {
        sorted_Vertices[0] = unique_Vertices[0];
        sorted_Vertices[1] = unique_Vertices[2];
        sorted_Vertices[2] = unique_Vertices[1];
        sorted_Vertices[3] = unique_Vertices[3];
        return;
    }

    // get the intersection between 02 and 13
    solutionExists = NPWGeometryMath::LineLineIntersect2D(
            unique_Vertices[0], unique_Vertices[2], unique_Vertices[1], unique_Vertices[3],
            intersect, intersectIsOnBothLines);
    if ( solutionExists && intersectIsOnBothLines ) // quad: 0123
    {
        sorted_Vertices[0] = unique_Vertices[0];
        sorted_Vertices[1] = unique_Vertices[1];
        sorted_Vertices[2] = unique_Vertices[2];
        sorted_Vertices[3] = unique_Vertices[3];
        return;
    }

    // get the intersection between 03 and 12
    solutionExists = NPWGeometryMath::LineLineIntersect2D(
            unique_Vertices[0], unique_Vertices[3], unique_Vertices[1], unique_Vertices[2],
            intersect, intersectIsOnBothLines);
    if ( solutionExists && intersectIsOnBothLines ) // quad: 0132
    {
        sorted_Vertices[0] = unique_Vertices[0];
        sorted_Vertices[1] = unique_Vertices[1];
        sorted_Vertices[2] = unique_Vertices[3];
        sorted_Vertices[3] = unique_Vertices[2];
        return;
    }

}

#if defined( DB_DEBUG) || defined( DB_COUNT )
const char * NPWVertices::getShapeString( Shape s, bool padded )
{
    switch ( s )
        {
        case TRUE_POINT:
            return padded ? "TRUE_POINT                        " : "TRUE_POINT";
        case POINT_EQUIVALENT:
            return padded ? "POINT_EQUIVALENT                  " : "POINT_EQUIVALENT";
        case LINE_WITH_VERTEX_TRIO:
            return padded ? "LINE_WITH_VERTEX_TRIO             " : "LINE_WITH_VERTEX_TRIO";
        case LINE_WITH_TWO_VERTEX_PAIRS:
            return padded ? "LINE_WITH_TWO_VERTEX_PAIRS        " : "LINE_WITH_TWO_VERTEX_PAIRS";
        case LINE_WITHOUT_VERTEX_PAIR:
            return padded ? "LINE_WITHOUT_VERTEX_PAIR          " : "LINE_WITHOUT_VERTEX_PAIR";
        case LINE_WITH_VERTEX_PAIR_ON_END:
            return padded ? "LINE_WITH_VERTEX_PAIR_ON_END      " : "LINE_WITH_VERTEX_PAIR_ON_END";
        case LINE_WITH_VERTEX_PAIR_ON_LINE:
            return padded ? "LINE_WITH_VERTEX_PAIR_ON_LINE     " : "LINE_WITH_VERTEX_PAIR_ON_LINE";
        case TRIANGLE_WITH_ONE_VERTEX_PAIR:
            return padded ? "TRIANGLE_WITH_ONE_VERTEX_PAIR     " : "TRIANGLE_WITH_ONE_VERTEX_PAIR";
        case TRIANGLE_WITH_ONE_VERTEX_ON_BORDER:
            return padded ? "TRIANGLE_WITH_ONE_VERTEX_ON_BORDER" : "TRIANGLE_WITH_ONE_VERTEX_ON_BORDER";
        case TRIANGLE_WITH_THREE_SUBTRIANGLES:
            return padded ? "TRIANGLE_WITH_THREE_SUBTRIANGLES  " : "TRIANGLE_WITH_THREE_SUBTRIANGLES";
        case QUAD:
            return padded ? "QUAD                              " : "QUAD";
        default:
            return padded ? "UNDEFINED                         " : "UNDEFINED";
        }
}
#endif

#ifdef DB_COUNT
void NPWVertices::printCounters()
{
    for ( int i = 0; i < (int)UNDEFINED+1; ++i )
    {
        if ( shapeCounters[i] > 0 )
        {
            std::cout << "\t" << getShapeString( Shape(i), true ) << " : " << shapeCounters[i] << std::endl;
        }
    }

    std::cout << "\tfailed triangle expansions         : " << triangleExpansionFails << "\n"
              << "\tfailed quad expansions             : " << quadExpansionFails << std::endl;
}
#endif


void NPWVertices::turnShapeExpansionOn()
{
    setRenderListFunction = &NPWVertices::setRenderListWithShapeExpansion;
}


void NPWVertices::turnShapeExpansionOff()
{
    setRenderListFunction = &NPWVertices::setRenderListWithoutShapeExpansion;
}


bool NPWVertices::isShapeExpansionOn()
{
    return ( setRenderListFunction == &NPWVertices::setRenderListWithShapeExpansion );
}


void NPWVertices::setRenderList()
{
    (this->*setRenderListFunction)();
}
