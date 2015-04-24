/*
 * NPWGeometry.h
 *
 *  Created on: 17/05/2010
 *      Author: bro86j
 */

#include "Vector2.h"

#ifndef NPWGEOMETRY_H_
#define NPWGEOMETRY_H_

using namespace npw;

/*
 * some geometry functions, specialised for use in an NPW context
 */
namespace NPWGeometryMath
{

/**
 * constants
 */
static const float MINIMUM_POINT_DISTANCE_SQUARED = 0.25f; // how close does a point need to be from a line
                                                           // to be considered to lie on this line
static const double EPSILON = 0.0001;


/**
 * own versions of basic functions to prevent weird runtime errors
 */
inline float  Abs( float val )                  { return ( val < 0.0f ? -val : val ); }
inline double Abs( double val )                 { return ( val < 0.0 ? -val : val ); }
inline int    Floor( float val )                { return (int)val; }
inline int    Floor( double val )               { return (int)val; }
inline int    Ceil( float val )                 { return ((int)val) + 1; }
inline int    Ceil( double val )                { return ((int)val) + 1; }
inline float  Min( float val1, float val2 )     { return ( val1 < val2 ? val1 : val2 ); }
inline float  Max( float val1, float val2 )     { return ( val1 > val2 ? val1 : val2 ); }
inline double Min( double val1, double val2 )   { return ( val1 < val2 ? val1 : val2 ); }
inline double Max( double val1, double val2 )   { return ( val1 > val2 ? val1 : val2 ); }
inline int    Min( int val1, int val2 )         { return ( val1 < val2 ? val1 : val2 ); }
inline int    Max( int val1, int val2 )         { return ( val1 > val2 ? val1 : val2 ); }

/**
 * computes the area of a triangle
 */
inline float computeTriangleArea( const Vector2 * A, const Vector2 * B, const Vector2 * C )
{
	// area = 0.5 | (Ax - Cx) (By - Ay) - (Ax - Bx) (Cy - Ay) |
	// http://en.wikipedia.org/wiki/Triangle
	return 0.5f * Abs( (A->x - C->x) * (B->y - A->y) -
			           (A->x - B->x) * (C->y - A->y) );
}

/**
 * Determines whether p is on the line segment AB
 */
bool isVertexOnLine( const Vector2 * p, const Vector2 * A, const Vector2 * B );

/**
 *  Given 4 points, determines whether the first point p lies in the triangle ABC
 *    algorithm given at http://www.blackpawn.com/texts/pointinpoly/default.html
 */
bool isVertexInTriangle(
		const Vector2 * A, const Vector2 * B, const Vector2 * C, const Vector2 * p );

/**
 * Calculates the intersection p of the two rays p1p2 and p3p4 in 2D
 *  Returns false if no solution exists.
 *  Parameter t specifies whether p is on both lines
 *  i.e. whether the line segments intersect
 *
 *  calculates determinants as described on http://mathworld.wolfram.com/Line-LineIntersection.html
 */
bool LineLineIntersect2D(
		const Vector2 & p1,	const Vector2 & p2, const Vector2 & p3, const Vector2 & p4,
		Vector2 & p, bool & intersectIsOnBothLines );

/**
 * Finds the closest point on the line segment from a to b from c
 * which handles if a or b are actually the closest points on the line segment.
 * Updates the closestPoint and t parameters
 *
 * http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
 */
void ClosestPointOnLine(
        const Vector2 & a, const Vector2 & b, const Vector2 & c,
        Vector2 & closestPoint, float & t, bool restrictPointToLineSegment = true );

/**
 * determines if two points p1 and p2 are on the same side of the line a-b
 */
bool arePointsOnSameSideOfLine( const Vector2 & a, const Vector2 & b,
                                const Vector2 & p1, const Vector2 & p2 );

/**
 * returns the longest line possible for four given points.
 * updates the pointers to point to corresponding points (p1 and p2 are the two remaining points)
 * also updates the distance squared
 */
void getLongestLineOfFourPoints(
        const Vector2 & a, const Vector2 & b, const Vector2 & c, const Vector2 & d,
        Vector2 & lineBegin, Vector2 & lineEnd, Vector2 & p1, Vector2 & p2,
        float & distanceSquared);

/**
 * returns the longest line possible for four given points.
 * updates the pointers to point to corresponding points (p1 and p2 are the two remaining points)
 * also updates the distance squared
 */
void getLongestLineOfThreePoints(
        const Vector2 & a, const Vector2 & b, const Vector2 & c,
        Vector2 & lineBegin, Vector2 & lineEnd, Vector2 & p, float & distanceSquared );

/**
 * returns whether p is close to a b ( as defined by MINIMUM_POINT_TO_LINE_DISTANCE )
 * and updates closestPoint
 */
bool isVertexCloseToLine(
        const Vector2 & a, const Vector2 & b, const Vector2 & p, Vector2 & closestPoint,
        float & t, bool restrictPointToLineSegment = true );

} // end namespace NPWGeometryMath

#endif /* NPWGEOMETRY_H_ */
