/*
 * NPWVertices.h
 *
 *  Created on: 23/06/2010
 *      Author: bro86j
 */

#ifndef NPW_VERTICES_H_
#define NPW_VERTICES_H_

//#define DB_DEBUG
//#define DB_COUNT
//#define NPW_SAFE_BORDERS

#include "Vector2.h"
#include "NPWConstants.h"
#include "NPWWindow.h"

#ifdef DB_DEBUG
#include <iostream>
#endif


using namespace npw;

class NPWWindow;

class NPWVertices
{
public:

    /**
     * constructor
     */
    NPWVertices();

    /**
     * destructor
     */
    ~NPWVertices() {}

    /**
     * set the vertices.
     * this sets unique_Vertices, n_Unique_Vertices n_Duplicates_at_Unique_Vertices
     */
    void setVerticesAndDetectShape( const Vector2 & v1, const Vector2 & v2, const Vector2 & v3, const Vector2 & v4 );

    /**
     * sets the parent window
     */
    void setParentWindow( NPWWindow * newParent );

    /**
     * sets the render list
     */
    void setRenderList();

    /**
     * turn shape expansion on and off
     */
    void turnShapeExpansionOn();
    void turnShapeExpansionOff();
    bool isShapeExpansionOn();

    /**
     * the render list
     */
    Vector2 * renderList[nMaxRenderPoints];
    unsigned int n_renderList_size;


private:
    /**
      * sets the render list function
      */
     void setRenderListWithoutShapeExpansion();
     void setRenderListWithShapeExpansion();

     typedef void (NPWVertices::*setRenderListFnPt)();
     setRenderListFnPt setRenderListFunction;


    /**
     * calculates the intersect, only valid if points are on a line and already sorted
     */
    void calculateIntersect();

    /**
     * sort the vertices, only makes sense if they are on a line
     */
    void sortVertices();

    /**
     * for a quad, determine the ordering of the vertices, store those in sorted_Vertices
     * and also compute the intersect of the diagonals and store this in intersect
     */
    void getQuadIntersectAndSortVertices();


    /**
     * correct shapes
     */
    Vector2 move;

    void correctLineWithoutVertexPair();
    void correctLineWithVertexPairOnEnd();
    void correctLineWithVertexPairOnLine();
    void correctLineWithTwoVertexPairs();
    void correctLineWithVertexTrio();

    /**
     * checks whether the renderlist as it is now is renderable
     */
    bool checkIfRenderListIsRenderable();

    /**
     * expands the triangle defined by renderList[1] to renderList[4],
     * but doesn't move vertices too much
     */
    bool expandTriangle();

    /**
     * expands the quad defined by renderList[1] to renderList[5],
     * but doesn't move vertices too much
     */
    bool expandQuad();

    /**
     * If triangles or quads are really narrow, expansion will fail.
     * Treat them as a line instead by projecting other vertex/vertices
     */
    void putFourPointsOnALineAndSetRenderList();
    void putThreePointsOnALineAndSetRenderList();

    /**
     * shape enumerator
     */
    enum Shape {QUAD = 0,
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
                UNDEFINED = 11};
    Shape shape;
#if defined(DB_DEBUG ) || defined( DB_COUNT )
    const char * getShapeString( Shape s, bool padded = false);
    const char * getShapeString( bool padded = false ) { return getShapeString( shape, padded ); }
#endif

    /**
     * unique vertices, set by calling setVertices() with four vector2's
     * the maximum number of unique vertices is nV and the minimum is 1
     */
    Vector2 unique_Vertices[nV];

    /**
     * sorted vertices, calculated by calling sortVertices()
     */
    Vector2 sorted_Vertices[nV];

    /**
     * the number of unique vertices
     */
    int n_Unique_Vertices;

    /**
     * for every vertex in unique_Vertices, how many other vertices are at that position
     */
    int n_Duplicates_at_Unique_Vertices[nV];

    /**
     * for every vertex in sorted_Vertices, how many other vertices are at that position
     */
    int n_Duplicates_at_Sorted_Vertices[nV];

    /**
     * the intersect, i.e. point of maximum intensity
     */
    Vector2 intersect;

    /**
     * are all points on a single line?
     */
    bool all_points_are_on_a_line;

    /**
     * if it is a triangle,
     * this holds the index in unique_vertices of the vertex with maximum intensity
     */
    int index_of_maximum_intensity_vertex_in_triangle;

    /**
     * used internally during setting of the shape, if we put three points on a line,
     * this holds the index in unique_Vertices of the point in the middle
     */
    int index_of_middle_vertex_on_line;

    /**
     * we need to know about the window to get the framebuffer size and the border size
     */
    NPWWindow * parentWindow;

#ifdef DB_COUNT
    /**
     * counters
     */
    int shapeCounters[(int)UNDEFINED+1];
    int triangleExpansionFails;
    int quadExpansionFails;

public:
    void printCounters();

#endif

#ifdef DB_DEBUG
    std::stringstream logStream;
    bool isErrored;
#endif

};

#endif // NPW_VERTICES_H_
