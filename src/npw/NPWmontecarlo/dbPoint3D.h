/*=========================================================================
Program: MILX MixView
Module: dbPoint3D.h
Author: Daan Broekhuizen
Modified by:
Language: C++
Created: Tue 25 May 2010 12:23:13 EST

Copyright: (c) 2009 CSIRO, Australia.

This software is protected by international copyright laws.
Any unauthorised copying, distribution or reverse engineering is prohibited.

Licence:
All rights in this Software are reserved to CSIRO. You are only permitted
to have this Software in your possession and to make use of it if you have
agreed to a Software License with CSIRO.

BioMedIA Lab: http://www.ict.csiro.au/BioMedIA/
=========================================================================*/

/*
*  BASIX: dbPoint3D
*
*  really basic point class
*
*  author: Daan Broekhuizen
*  date: january 2010
*/

template< typename T >
class dbPoint3D {
  private:
    T x;
    T y;
    T z;

  public:
    // CONSTRUCTORS
    dbPoint3D( const T x, const T y, const T z ) {
      this->x = x;
      this->y = y;
      this->z = z;
    }

    dbPoint3D( const dbPoint3D<T> &otherPoint ) {
      this->x = otherPoint.x;
      this->y = otherPoint.y;
      this->z = otherPoint.z;
    }

    dbPoint3D() {
      x = 0;
      y = 0;
      z = 0;
    }

    // BASIC ARITHMETIC

    // point addition
    dbPoint3D<T>& operator+=(const dbPoint3D<T> &otherPoint) {
      x += otherPoint.x;
      y += otherPoint.y;
      z += otherPoint.z;
      return *this;
    }
    dbPoint3D<T> operator+(const dbPoint3D<T> &rhs) const {
      return dbPoint3D<T>(*this) += rhs;
    }


    // scalar addition
    dbPoint3D<T>& operator+=(const T a) {
      x += a;
      y += a;
      z += a;
      return *this;
    }
    dbPoint3D<T> operator+(const T a) const {
      return dbPoint3D<T>(*this) += a;
    }

    // point subtraction
    dbPoint3D<T>& operator-=(const dbPoint3D<T> &otherPoint) {
      x -= otherPoint.x;
      y -= otherPoint.y;
      z -= otherPoint.z;
      return *this;
    }
    dbPoint3D<T> operator-(const dbPoint3D<T> &rhs) const {
      return dbPoint3D<T>(*this) -= rhs;
    }

    // scalar subtraction
    dbPoint3D<T>& operator-=(const T a) {
      x -= a;
      y -= a;
      z -= a;
      return *this;
    }
    dbPoint3D<T> operator-(const T a) const {
      return dbPoint3D<T>(*this) -= a;
    }


    // scalar multiplication
    dbPoint3D<T>& operator*=(const T a) {
      x *= a;
      y *= a;
      z *= a;
      return *this;
    }
    dbPoint3D<T> operator*(const T a) const {
      return dbPoint3D<T>(*this) *= a;
    }

    // scalar division
    dbPoint3D<T>& operator/=(const T a) {
      x /= a;
      y /= a;
      z /= a;
      return *this;
    }
    dbPoint3D<T> operator/(const T a) const {
      return dbPoint3D<T>(*this) /= a;
    }

    // elementwise point division
    dbPoint3D<T>& operator/=(const dbPoint3D<T> &p) {
      x /= p.x;
      y /= p.y;
      z /= p.z;
      return *this;
    }
    dbPoint3D<T>& operator/(const dbPoint3D<T> &rhs) const {
      return dbPoint3D<T>(*this) /= rhs;
    }

    // dot product
    T dot(dbPoint3D<T> &rhs) const {
      return (x*rhs.x) + (y*rhs.y) + (z*rhs.z);
    }

    // cross product
    dbPoint3D<T> cross(dbPoint3D<T> &rhs) const {
      return dbPoint3D<T>( y*rhs.z - z*rhs.y,
                           z*rhs.x - x*rhs.z,
                           x*rhs.y - y*rhs.x );
    }

    // ACCESS
    T getX() const {return x;}
    T getY() const {return y;}
    T getZ() const {return z;}

    void setX(const T x) { this->x = x; }
    void setY(const T y) { this->y = y; }
    void setZ(const T z) { this->z = z; }
    void set(const T x, const T y, const T z) {
      this->x = x;
      this->y = y;
      this->z = z;
    }

};