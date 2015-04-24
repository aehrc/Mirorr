/*
 * Vector2.cxx
 *
 *  Created on: 25/05/2010
 *      Author: bro86j
 */


#ifndef __VECTOR2_H
#define __VECTOR2_H

#include <memory.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace npw {

/**
 * \class Vector2
 *
 * \brief A 2D vector with standard set of operators.
 *
 * \ingroup MILXSimMath
 */
    class Vector2
    {
        public:
            float x, y;

            /**
            * Constructors
            */
            Vector2()                     : x(0.0f), y(0.0f)  { }
            Vector2(float x_n, float y_n) : x(x_n) , y(y_n)   { }
            Vector2(const Vector2 &v)     : x(v.x) , y(v.y)   { }

            /**
            * Operators
            */
            inline Vector2 operator +  (const Vector2 &v) const { return Vector2(x + v.x, y + v.y); }
            inline Vector2 operator -  (const Vector2 &v) const { return Vector2(x - v.x, y - v.y); }
            inline Vector2 operator *  (const Vector2 &v) const { return Vector2(x * v.x, y * v.y); }
            inline Vector2 operator /  (const Vector2 &v) const { return Vector2(x / v.x, y / v.y); }
            inline Vector2 operator += (const Vector2 &v)       { *this = *this + v; return *this; }
            inline Vector2 operator -= (const Vector2 &v)       { *this = *this - v; return *this; }
            inline Vector2 operator *= (const Vector2 &v)       { *this = *this * v; return *this; }
            inline Vector2 operator /= (const Vector2 &v)       { *this = *this / v; return *this; }

            /**
            * Float operators
            */
            inline Vector2 operator *  (float f) const          { return Vector2(x * f, y * f); }
            inline Vector2 operator /  (float f) const          { return Vector2(x / f, y / f); }
            inline Vector2 operator *= (float f)                { *this = *this * f; return *this; }
            inline Vector2 operator /= (float f)                { *this = *this / f; return *this; }

            /**
            * Array operators
            */
            inline float  &operator [] (unsigned int i)         { float *pData = ((float *) this + i); return *pData; }
            inline float   operator [] (unsigned int i) const   { return *((float *) this + i); }

            /**
            * Assignment operators
            */
            inline void operator = (const Vector2 &v)                  {  x = v.x; y = v.y; }

            /**
            * Comparison operators
            */
            bool operator == (const Vector2 &v) const           {
//            	return (memcmp(this, &v, sizeof(Vector2)) == 0);
            	return ( x == v.x && y == v.y );
            }
            bool operator != (const Vector2 &v) const           {
//            	return (memcmp(this, &v, sizeof(Vector2)) != 0);
            	return ( x != v.x || y != v.y );
            }

            /**
            * Complex operators
            */
            inline float Dot(const Vector2 &v) const            { return x * v.x + y * v.y; }
            inline float Dot() const                            { return x * x + y * y; }
            inline float Magnitude() const                      { return (float)sqrt(Dot()); }
            inline float MagnitudeSquared() const               { return Dot(); }

            friend std::ostream& operator << (std::ostream& os, const Vector2 & v)
            {
                os << "[ " << v.x << " " << v.y << " ]";
                return os;
            }

            /**
            * Normalise and return this vector
            */
            Vector2 Normalise()
            {
              float fMag = MagnitudeSquared();

              if (fMag > 0.0f)
                return (*this) * (1.0f / sqrtf(fMag));
              else
                return (*this);
            }

            /**
            * Return a rotated variant of this vector
            */
            Vector2 Rotate(float fAngle)
            {
              float fSin = (float)sin(fAngle);
              float fCos = (float)cos(fAngle);

              return Vector2(x*fCos - y*fSin, x*fSin + y*fCos);
            }

            /**
             * Return distance to another vector
             */
            float Distance(const Vector2 &v) const        { return (*this - v).Magnitude(); }
            float DistanceSquared(const Vector2 &v) const { return (*this - v).MagnitudeSquared(); }

    };

} // end namespace npw

#endif // __VECTOR2_H
