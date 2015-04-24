/*
 * NPWGeometry.cxx
 *
 *  Created on: 18/05/2010
 *      Author: bro86j
 */

#include "NPWGeometry.h"
#include "NPWConstants.h"

bool NPWGeometryMath::isVertexOnLine(
        const Vector2 * p, const Vector2 * A, const Vector2 * B)
{
    /*
     * get the closest point on the line AB from p
     */
    Vector2 closestPoint;
    float unused;
    ClosestPointOnLine( *A, *B, *p, closestPoint, unused, true );

    return ( closestPoint.DistanceSquared( *p ) < MINIMUM_POINT_DISTANCE_SQUARED );
}

bool NPWGeometryMath::isVertexInTriangle(
		const Vector2 * A, const Vector2 * B, const Vector2 * C, const Vector2 * p )
{
    Vector2 v0 = *C - *A;
    Vector2 v1 = *B - *A;
    Vector2 v2 = *p - *A;

    double dot00 = (double)v0.Dot(v0);
    double dot01 = (double)v0.Dot(v1);
    double dot02 = (double)v0.Dot(v2);
    double dot11 = (double)v1.Dot(v1);
    double dot12 = (double)v1.Dot(v2);

    // barycentric coordinates
    double denom = (dot00 * dot11) - (dot01 * dot01);
    if ( denom < EPSILON ) {
        std::cout << "error determining if point " << *p
                  << " is in triangle " << *A << " - "
                  << *B << " - " << *C << std::endl;
        return false;
    }
    double u = ( (dot11 * dot02) - (dot01 * dot12) ) / denom;
    double v = ( (dot00 * dot12) - (dot01 * dot02) ) / denom;

    return (u > 0.0) && (v > 0.0) && (u + v < 1.0);
}

bool NPWGeometryMath::LineLineIntersect2D(
		const Vector2 & p1,	const Vector2 & p2, const Vector2 & p3, const Vector2 & p4,
		Vector2 & p, bool & intersectIsOnBothLines )
{
	// the denominator
	double d = ((double)(p1.x - p2.x) * (double)(p3.y - p4.y)) -
	           ((double)(p1.y - p2.y) * (double)(p3.x - p4.x));
	if ( d == 0.0 )
	{
		intersectIsOnBothLines = false; // there is no intersect, so this is neither true
										// nor false, but false makes more sense than true.
		return false; // lines are parallel
	}

	double pre  = (double)p1.x * (double)p2.y - (double)p1.y * (double)p2.x;
	double post = (double)p3.x * (double)p4.y - (double)p3.y * (double)p4.x;

	p.x = (float)( (pre * (p3.x - p4.x) - (p1.x - p2.x) * post) / d );
	p.y = (float)( (pre * (p3.y - p4.y) - (p1.y - p2.y) * post) / d );

	// check whether the intersect is on both lines by checking whether it is in the bounding box
	float minX = Min(p1.x, p2.x);
	float maxX = Max(p1.x, p2.x);
	float minY = Min(p1.y, p2.y);
	float maxY = Max(p1.y, p2.y);
	bool isOnLine1 = p.x <= maxX && p.x >= minX && p.y <= maxY && p.y >= minY;

	minX = Min(p3.x, p4.x);
	maxX = Max(p3.x, p4.x);
	minY = Min(p3.y, p4.y);
	maxY = Max(p3.y, p4.y);
	bool isOnLine2 = p.x <= maxX && p.x >= minX && p.y <= maxY && p.y >= minY;

	intersectIsOnBothLines = isOnLine1 && isOnLine2;
	return true;
}

void NPWGeometryMath::ClosestPointOnLine(
        const Vector2 & a, const Vector2 & b, const Vector2 & c,
        Vector2 & closestPoint, float & t, bool restrictPointToLineSegment )
{
    Vector2 ab = b - a;
    Vector2 ac = c - a;

    t = ac.Dot(ab) / ab.Dot();

    if ( restrictPointToLineSegment )
    {
        if (t < 0.0f)
        {
            t = 0.0f;
        }
        else if (t > 1.0f)
        {
            t = 1.0f;
        }
    }

    closestPoint = a + ab * t;
}

bool NPWGeometryMath::arePointsOnSameSideOfLine(
        const Vector2 & a, const Vector2 & b, const Vector2 & p1, const Vector2 & p2 )
{
    // get crossproduct of (b - a) and (c - a)
    // (technically, the z-component, but since they are Vector2's, x & y components are 0)
    Vector2 BminusA = b - a;
    Vector2 P1minusA = p1 - a;
    float cross1 = BminusA.x * P1minusA.y - BminusA.y * P1minusA.x;

    // get crossproduct of (b - a) and (d - a)
    Vector2 P2minusA = p2 - a;
    float cross2 = BminusA.x * P2minusA.y - BminusA.y * P2minusA.x;

    // if the sign of the cross products are equal, c & d are on the same side of ab
    return ( ( cross1 < 0.0f && cross2 < 0.0f ) || ( cross1 > 0.0f && cross2 > 0.0f ) );
}

void NPWGeometryMath::getLongestLineOfFourPoints(
        const Vector2 & a, const Vector2 & b, const Vector2 & c, const Vector2 & d,
        Vector2 & lineBegin, Vector2 & lineEnd, Vector2 & p1, Vector2 & p2,
        float & distanceSquared )
{
    // initialize with line ab
    lineBegin = a;
    lineEnd   = b;
    p1        = c;
    p2        = d;
    distanceSquared = lineBegin.DistanceSquared( lineEnd );

    // line ac
    float newDistanceSquared = lineBegin.DistanceSquared( c );
    if ( newDistanceSquared > distanceSquared )
    {
        lineEnd = c;
        p1      = b;
        distanceSquared = newDistanceSquared;
    }

    // line ad
    newDistanceSquared = lineBegin.DistanceSquared( d );
    if ( newDistanceSquared > distanceSquared )
    {
        lineEnd = d;
        p1      = b;
        p2      = c;
        distanceSquared = newDistanceSquared;
    }

    // line bc
    newDistanceSquared = b.DistanceSquared( c );
    if ( newDistanceSquared > distanceSquared )
    {
        lineBegin = b;
        lineEnd   = c;
        p1        = a;
        p2        = d;
        distanceSquared = newDistanceSquared;
    }

    // line bd
    newDistanceSquared = b.DistanceSquared( d );
    if ( newDistanceSquared > distanceSquared )
    {
        lineBegin = b;
        lineEnd   = d;
        p1        = a;
        p2        = c;
        distanceSquared = newDistanceSquared;
    }

    // line cd
    newDistanceSquared = c.DistanceSquared( d );
    if ( newDistanceSquared > distanceSquared )
    {
        lineBegin = c;
        lineEnd   = d;
        p1        = a;
        p2        = b;
        distanceSquared = newDistanceSquared;
    }
}

void NPWGeometryMath::getLongestLineOfThreePoints(
        const Vector2 & a, const Vector2 & b, const Vector2 & c,
        Vector2 & lineBegin, Vector2 & lineEnd, Vector2 & p, float & distanceSquared )
{
    // initialize with line ab
    lineBegin = a;
    lineEnd   = b;
    p         = c;
    distanceSquared = lineBegin.DistanceSquared( lineEnd );

    // line ac
    float newDistanceSquared = lineBegin.DistanceSquared( c );
    if ( newDistanceSquared > distanceSquared )
    {
        lineEnd = c;
        p       = b;
        distanceSquared = newDistanceSquared;
    }

    // line bc
    newDistanceSquared = b.DistanceSquared( c );
    if ( newDistanceSquared > distanceSquared )
    {
        lineBegin = b;
        lineEnd   = c;
        p         = a;
        distanceSquared = newDistanceSquared;
    }
}

bool NPWGeometryMath::isVertexCloseToLine(
        const Vector2 & a, const Vector2 & b, const Vector2 & p, Vector2 & closestPoint,
        float & t, bool restrictPointToLineSegment )
{
    ClosestPointOnLine( a, b, p, closestPoint, t, restrictPointToLineSegment );
    return ( closestPoint.DistanceSquared( p ) < MINIMUM_POINT_TO_LINE_DISTANCE_SQUARED );
}

