/*
 * NPWCudaGeometryFunctions.h
 *
 *  Created on: 20/07/2010
 *      Author: bro86j
 */

#ifndef NPW_CUDA_GEOMETRY_FUNCTIONS_H
#define NPW_CUDA_GEOMETRY_FUNCTIONS_H

#include "cuda_runtime.h"
#include "math_functions.h"
#include "NPWCudaConstants.h"


/*
 * basic float2 arithmetic
 */
inline __device__ float2 operator-(float2 &a)           { return make_float2(-a.x, -a.y); }
inline __device__ float2 floor(const float2 v)          { return make_float2(floor(v.x), floor(v.y)); }

inline __device__ float2 operator+(float2 a, float2 b)  { return make_float2(a.x + b.x, a.y + b.y); }
inline __device__ void operator+=(float2 &a, float2 b)  { a.x += b.x; a.y += b.y; }

inline __device__ float2 operator+(float2 a, float b)   { return make_float2(a.x + b, a.y + b); }
inline __device__ float2 operator+(float b, float2 a)   { return make_float2(a.x + b, a.y + b); }
inline __device__ void operator+=(float2 &a, float b)   { a.x += b; a.y += b; }

inline __device__ float2 operator-(float2 a, float2 b)  { return make_float2(a.x - b.x, a.y - b.y); }
inline __device__ void operator-=(float2 &a, float2 b)  { a.x -= b.x; a.y -= b.y; }

inline __device__ float2 operator-(float2 a, float b)   { return make_float2(a.x - b, a.y - b); }
inline __device__ float2 operator-(float b, float2 a)   { return make_float2(a.x - b, a.y - b); }
inline __device__ void operator-=(float2 &a, float b)   { a.x -= b; a.y -= b; }

inline __device__ float2 operator*(float2 a, float2 b)  { return make_float2(a.x * b.x, a.y * b.y); }
inline __device__ float2 operator*(float2 a, float s)   { return make_float2(a.x * s, a.y * s); }
inline __device__ float2 operator*(float s, float2 a)   { return make_float2(a.x * s, a.y * s); }

inline __device__ void operator*=(float2 &a, float s)   { a.x *= s; a.y *= s; }

inline __device__ float2 operator/(float2 a, float2 b)  { return make_float2(a.x / b.x, a.y / b.y); }
inline __device__ float2 operator/(float2 a, float s)   { float inv = 1.0f / s; return a * inv; }
inline __device__ float2 operator/(float s, float2 a)   { float inv = 1.0f / s; return a * inv; }

inline __device__ void operator/=(float2 &a, float s)   { float inv = 1.0f / s; a *= inv; }

inline __device__ bool operator==(float2 &a, float2 &b) { return a.x == b.x && a.y == b.y; }
inline __device__ bool operator!=(float2 &a, float2 &b) { return a.x != b.x || a.y != b.y; }

/*
 * my basic functions
 */
inline __device__ float Abs(float a)                     { return a < 0.0f ? -a : a; }
inline __device__ int   Abs(int a)                       { return a < 0    ? -a : a; }

inline __device__ float Min(float a, float b)            { return a < b ? a : b; }
inline __device__ float Min(float a, float b, float c)   { return Min(Min(a,b),c); }
inline __device__ float Max(float a, float b)            { return a > b ? a : b; }
inline __device__ float Max(float a, float b, float c)   { return Max(Max(a,b),c); }
inline __device__ int   Max(int a, int b)                { return a > b ? a : b; }
inline __device__ int   Min(int a, int b)                { return a < b ? a : b; }
inline __device__ float Min(int a, int b, int c)         { return Min(Min(a,b),c); }
inline __device__ float Max(int a, int b, int c)         { return Max(Max(a,b),c); }

inline __device__ int   Sign(float a) // returns 1 if positive, -1 if negative and 0 if equal to 0.0f
{
    return ( a > 0.0f ? 1 : (a < 0.0f ? -1 : 0) );
}

/*
 * magnitude of a float 2
 */
inline __device__ float MagnitudeSquared( float2 a )
{
    return a.x * a.x + a.y * a.y;
}
inline __device__ float Magnitude( float2 a )
{
    return sqrtf( MagnitudeSquared(a) );
}

/*
 * normalize a float2
 */
inline __device__ float2 Normalize( float2 a )
{
    float invDenom = 1.0f / Magnitude( a );
    return make_float2( a.x * invDenom, a.y * invDenom );
}

/*
 * convert multidimensional index to a linear index, in 2D and 3D
 */
inline __device__ int getArrayIndex2D(
        int x, int y, int xdim, int ydim )
{
    int clamped_x = ( x < 0 ? Abs(x) : (x >= xdim ? (2*xdim-x-1) : x) );
    int clamped_y = ( y < 0 ? Abs(y) : (y >= ydim ? (2*ydim-y-1) : y) );

    return clamped_x + __mul24( clamped_y, xdim );
}

inline __device__ unsigned int getArrayIndex3D(
        unsigned int x, unsigned int y, unsigned int z, unsigned int xDim, unsigned int yDim )
{
    // idx = x + y * xDim + z * xDim * yDim
    //     = x + xDim * ( z * yDim + y );
//    return x + __umul24( xDim, __umul24( z, yDim ) + y );
    return x + y * xDim + z * xDim * yDim;
}


/*
 * dot product of two float2's
 */
inline __device__ float  dotf(float ax, float ay, float bx, float by)
{
    return ax * bx + ay * by;
}
inline __device__ float  dotf(float2 a, float2 b)
{
    return a.x * b.x + a.y * b.y;
}


/*
 * pseudo cross product of two float2's
 */
inline __device__ float cross(float ax, float ay, float bx, float by)
{
    return ay * bx - ax * by;
}
inline __device__ float cross(float2 a, float2 b)
{
    return a.y * b.x - a.x * b.y;
}

/*
 * distance between two points
 */
inline __device__ float PointPointDistanceSquared( float ax, float ay, float bx, float by )
{
    return ( ax - bx ) * ( ax - bx ) + ( ay - by ) * ( ay - by );
}
inline __device__ float PointPointDistance( float ax, float ay, float bx, float by )
{
    return sqrtf( PointPointDistanceSquared(ax,ay,bx,by) );
}
inline __device__ float PointPointDistanceSquared( float2 a, float2 b )
{
    return ( a.x - b.x ) * ( a.x - b.x ) + ( a.y - b.y ) * ( a.y - b.y );
}
inline __device__ float PointPointDistance( float2 a, float2 b )
{
    return sqrtf( PointPointDistanceSquared(a,b) );
}

/*
 * returns the signed area of a triangle times two
 */
inline __device__  float SignedTriangleAreaTimesTwo( float2 a, float2 b, float2 c )
{
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}
inline __device__  float SignedTriangleAreaTimesTwo(
        float ax, float ay, float bx, float by, float cx, float cy )
{
    return (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
}

/*
 * returns the absolute area of a triangle
 */
inline __device__  float TriangleArea( float2 a, float2 b, float2 c )
{
    return 0.5f * Abs(SignedTriangleAreaTimesTwo(a,b,c));
}
inline __device__  float TriangleArea(
        float ax, float ay, float bx, float by, float cx, float cy )
{
    return 0.5f * Abs(SignedTriangleAreaTimesTwo(ax,ay,bx,by,cx,cy));
}

/*
 * returns the area of a quad, abcd should be ordered (either clockwise or counterclockwise)
 */
inline __device__ float SignedQuadAreaTimesTwo(
        float ax, float ay, float bx, float by, float cx, float cy, float dx, float dy )
{
    // area = 0.5 * | p x q | , where p and q are the diagonals
    //     note: z coordinate is 0.0
    //     ref: http://mathworld.wolfram.com/Quadrilateral.html
    // float2 p = c - a;
    // float2 q = d - b;
    // float abs_cross = abs( p.x * q.y - p.y * q.x );
    // area = 0.5f * abs_cross;
    return cross( cx-ax, cy-ay, dx-bx, dy-by );
}
inline __device__ float QuadArea(
        float ax, float ay, float bx, float by, float cx, float cy, float dx, float dy )
{
    // area = 0.5 * | p x q | , where p and q are the diagonals
    //     note: z coordinate is 0.0
    //     ref: http://mathworld.wolfram.com/Quadrilateral.html
    // float2 p = c - a;
    // float2 q = d - b;
    // float abs_cross = abs( p.x * q.y - p.y * q.x );
    // area = 0.5f * abs_cross;
    return 0.5f * Abs( SignedQuadAreaTimesTwo(ax,ay, bx,by, cx,cy, dx,dy) );
}

/*
 * test whether a triangle is clockwise oriented, based on
 */
inline __device__  bool IsTriangleClockwise( float2 a, float2 b, float2 c )
{
    return SignedTriangleAreaTimesTwo(a,b,c) < 0.0f;
}
inline __device__  bool IsTriangleClockwise(
        float ax, float ay, float bx, float by, float cx, float cy )
{
    return SignedTriangleAreaTimesTwo(ax,ay, bx,by, cx,cy) < 0.0f;
}

/*
 * cheap 2D point in triangle test for counter clockwise triangles
 */
inline __device__  bool IsPointInCounterClockwiseTriangle( float2 a, float2 b, float2 c, float2 p )
{
    if ( cross(p-a,b-a) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    if ( cross(p-b,c-b) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    if ( cross(p-c,a-c) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    return true;
}
inline __device__  bool IsPointInCounterClockwiseTriangle(
        float ax, float ay, float bx, float by, float cx, float cy, float px, float py )
{
    if ( cross(px-ax,py-ay, bx-ax,by-ay) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    if ( cross(px-bx,py-by, cx-bx,cy-by) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    if ( cross(px-cx,py-cy, ax-cx,ay-cy) < NPW_CUDA_NEGATIVE_EPSILON ) return false;
    return true;
}

/*
 * cheap 2D point in triangle test, where the triangle can be counterclockwise or clockwise
 */
inline __device__ bool PointInTriangle( float2 a, float2 b, float2 c, float2 p )
{
    return IsTriangleClockwise(a,b,c) ?
                IsPointInCounterClockwiseTriangle(a,c,b,p) :
                IsPointInCounterClockwiseTriangle(a,b,c,p);
}
inline __device__ bool PointInTriangle(
        float ax, float ay, float bx, float by, float cx, float cy, float px, float py )
{
    return IsTriangleClockwise(ax,ay, bx,by, cx,cy) ?
                IsPointInCounterClockwiseTriangle(ax,ay, cx,cy, bx,by, px,py) :
                IsPointInCounterClockwiseTriangle(ax,ay, bx,by, cx,cy, px,py);
}

/*
 * given a point p and a rectangle ab,
 * determines whether the point is inside the rectangle or on its edge
 */
inline __device__  bool PointInRectangle(
        float ax, float ay, float bx, float by, float px, float py )
{
    return ( px <= bx &&
             px >= ax &&
             py <= by &&
             py >= ay );
}
inline __device__  bool PointInRectangle( float2 a, float2 b, float2 p )
{
    return ( p.x <= b.x &&
             p.x >= a.x &&
             p.y <= b.y &&
             p.y >= a.y );
}

/*
 * determines whether p is on the same side of ab as s
 */
inline __device__ bool PointsOnSameSideOfLine( float2 a, float2 b, float2 p, float2 s )
{
    // get crossproduct of (b - a) and (p - a)
    float2 BminusA  = b - a;
    float2 PminusA  = p - a;
    float cross1    = cross( BminusA, PminusA);

    // get crossproduct of (b - a) and (s - a)
    float2 SminusA  = s - a;
    float cross2    = cross( BminusA, SminusA);

    // if the sign of the cross products are equal, c & d are on the same side of ab
    return ( Sign(cross1) == Sign(cross2) );
}

/*
 * find the intersection of two line segments
 */
inline __device__ bool LineSegmentLineSegmentIntersect(
        float2 a, float2 b, float2 c, float2 d,
        float2 & intersect )
{
    // the denominator
    float denom = ( (a.x - b.x) * (c.y - d.y) ) -
                  ( (a.y - b.y) * (c.x - d.x) );

    if ( denom == 0.0f )
    {
        return false; // lines are parallel
    }

    float invDenom = 1.0f / denom;

    float r = ( (a.y - c.y) * (d.x - c.x) - (a.x - c.x) * (d.y - c.y) ) * invDenom;
    float s = ( (a.y - c.y) * (b.x - a.x) - (a.x - c.x) * (b.y - a.y) ) * invDenom;

    if ( r < 0.0f || r > 1.0f || s < 0.0f || s > 1.0f )
    {
        return false;
    }

    intersect = a + r * ( b - a);

    return true;
}

inline __device__ bool LineSegmentLineSegmentIntersect(
        float ax, float ay, float bx, float by, float cx, float cy, float dx, float dy,
        float * ix, float * iy )
{
    // the denominator
    float denom = ( (ax - bx) * (cy - dy) ) -
                  ( (ay - by) * (cx - dx) );

    if ( denom == 0.0f )
    {
        return false; // lines are parallel
    }

    float invDenom = 1.0f / denom;

    float r = ( (ay - cy) * (dx - cx) - (ax - cx) * (dy - cy) ) * invDenom;
    float s = ( (ay - cy) * (bx - ax) - (ax - cx) * (by - ay) ) * invDenom;

    if ( r < 0.0f || r > 1.0f || s < 0.0f || s > 1.0f )
    {
        return false;
    }

    *ix = ax + r * ( bx - ax );
    *iy = ay + r * ( by - ay );

    return true;
}

/**
 * Finds the closest point on the line segment from a to b from c
 * which handles if a or b are actually the closest points on the line segment.
 * Updates the closestPoint and t parameters
 *
 * http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/
 */
inline __device__ void ClosestPointOnLine(
        float2 a, float2 b, float2 c, float2 & closestPoint,
        float & t, bool restrictPointToLineSegment = true )
{
    float2 ab = b - a;
    float2 ac = c - a;

    t = dotf(ac,ab) / dotf(ab,ab);

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

/**
 * returns whether p is close to a b ( as defined by MINIMUM_POINT_TO_LINE_DISTANCE )
 * and updates closestPoint
 * todo: change is into Is in function name
 */
inline __device__ bool isVertexCloseToLine(
        float2 a, float2 b, float2 p, float2 & closestPoint,
        float & t, bool restrictPointToLineSegment = true )
{
    ClosestPointOnLine( a, b, p, closestPoint, t, restrictPointToLineSegment );
    return ( PointPointDistanceSquared( closestPoint, p )
             < MINIMUM_POINT_TO_LINE_DISTANCE_SQUARED );
}

/*
 * takes a triangle, defined by a, b and c
 * where a is assumed to have a height of h and b and c a height of 0.0
 *
 * returns the interpolated value at p
 *
 * it is not checked whether p is in abc
 */
inline __device__ float InterpolateTriangle( float2 a, float2 b, float2 c, float2 p, float height )
{
    float2 v0 = c - a;
    float2 v1 = b - a;
    float2 v2 = p - a;

    float dot00 = dotf(v0,v0);
    float dot01 = dotf(v0,v1);
    float dot02 = dotf(v0,v2);
    float dot11 = dotf(v1,v1);
    float dot12 = dotf(v1,v2);

    float denom = (dot00 * dot11) - (dot01 * dot01);

    if ( denom == 0.0f )
    {
        // error determining if point is in triangle
        return 0.0f;
    }

    float invDenom = 1.0f / denom;

    float u = ( (dot11 * dot02) - (dot01 * dot12) ) * invDenom;
    float v = ( (dot00 * dot12) - (dot01 * dot02) ) * invDenom;

    float w_clamped = Max( 0.0f, Min( 1.0f, 1.0f - u - v) );

    return w_clamped * height;
}
inline __device__ float InterpolateTriangle(
        float ax, float ay, float bx, float by, float cx, float cy, float px, float py,
        float height )
{
    float v0x = cx - ax;
    float v0y = cy - ay;
    float v1x = bx - ax;
    float v1y = by - ay;
    float v2x = px - ax;
    float v2y = py - ay;

    float dot00 = dotf(v0x,v0y,v0x,v0y);
    float dot01 = dotf(v0x,v0y,v1x,v1y);
    float dot02 = dotf(v0x,v0y,v2x,v2y);
    float dot11 = dotf(v1x,v1y,v1x,v1y);
    float dot12 = dotf(v1x,v1y,v2x,v2y);

    float denom = (dot00 * dot11) - (dot01 * dot01);

    if ( denom == 0.0f )
    {
        // error determining if point is in triangle
        return 0.0f;
    }

    float invDenom = 1.0f / denom;

    float u = ( (dot11 * dot02) - (dot01 * dot12) ) * invDenom;
    float v = ( (dot00 * dot12) - (dot01 * dot02) ) * invDenom;

    float w_clamped = Max( 0.0f, Min( 1.0f, 1.0f - u - v) );

    return w_clamped * height;
}


#endif //NPW_CUDA_GEOMETRY_FUNCTIONS_H

