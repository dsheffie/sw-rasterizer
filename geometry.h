//[header]
// This program illustrates how the concept of vector and matrix can be implemented
// in C++. This is a light version of the implementation. It contains the most
// essential methods to manipulate vectors and matrices. It should be enough
// for most projects. Vectors and matrices are really the alphabet as we said
// in the lesson of any graphics application. It's really important you feel
// confortable with these techniques especially with the concepts of
// normalizing vectors, computing their length, computing the dot and cross products
// of two vectors, and the point- and vector-matrix multiplication (and knowing
// the difference between the two).
//[/header]
//[compile]
// c++ geometry.cpp  -o geometry -std=c++11
//[/compile]
//[ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//[/ignore]
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>

template<typename T>
class Vec2
{
public:
    Vec2() : x(0), y(0) {}
    Vec2(T xx) : x(xx), y(xx) {}
    Vec2(T xx, T yy) : x(xx), y(yy) {}
    Vec2 operator + (const Vec2 &v) const
    { return Vec2(x + v.x, y + v.y); }
    Vec2 operator / (const T &r) const
    { return Vec2(x / r, y / r); }
    Vec2 operator * (const T &r) const
    { return Vec2(x * r, y * r); }
    Vec2& operator /= (const T &r)
    { x /= r, y /= r; return *this; }
    Vec2& operator *= (const T &r)
    { x *= r, y *= r; return *this; }
    friend std::ostream& operator << (std::ostream &s, const Vec2<T> &v)
    {
        return s << '[' << v.x << ' ' << v.y << ']';
    }
    friend Vec2 operator * (const T &r, const Vec2<T> &v)
    { return Vec2(v.x * r, v.y * r); }
    T x, y;
};

typedef Vec2<float> Vec2f;
typedef Vec2<int> Vec2i;

//[comment]
// Implementation of a generic vector class - it will be used to deal with 3D points, vectors and normals.
// The class is implemented as a template. While it may complicate the code a bit, it gives us
// the flexibility later, to specialize the type of the coordinates into anything we want.
// For example: Vec3f if we want the coordinates to be floats or Vec3i if we want the coordinates to be integers.
//
// Vec3 is a standard/common way of naming vectors, points, etc. The OpenEXR and Autodesk libraries
// use this convention for instance.
//[/comment]
template<typename T>
class Vec3
{
public:
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3 operator + (const Vec3 &v) const
    { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator - (const Vec3 &v) const
    { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator - () const
    { return Vec3(-x, -y, -z); }
    Vec3 operator * (const T &r) const
    { return Vec3(x * r, y * r, z * r); }
    Vec3 operator * (const Vec3 &v) const
    { return Vec3(x * v.x, y * v.y, z * v.z); }
    T dotProduct(const Vec3<T> &v) const
    { return x * v.x + y * v.y + z * v.z; }
    Vec3& operator /= (const T &r)
    { x /= r, y /= r, z /= r; return *this; }
    Vec3& operator *= (const T &r)
    { x *= r, y *= r, z *= r; return *this; }
    Vec3 crossProduct(const Vec3<T> &v) const
    { return Vec3<T>(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
    T norm() const
    { return x * x + y * y + z * z; }
  T length() const {
    return std::sqrt(norm());
  }
    //[comment]
    // The next two operators are sometimes called access operators or
    // accessors. The Vec coordinates can be accessed that way v[0], v[1], v[2],
    // rather than using the more traditional form v.x, v.y, v.z. This useful
    // when vectors are used in loops: the coordinates can be accessed with the
    // loop index (e.g. v[i]).
    //[/comment]
    const T& operator [] (uint8_t i) const { return (&x)[i]; }
    T& operator [] (uint8_t i) { return (&x)[i]; }
    Vec3& normalize()
    {
        T n = norm();
        if (n > 0) {
            T factor = 1 / sqrt(n);
            x *= factor, y *= factor, z *= factor;
        }
        
        return *this;
    }

    friend Vec3 operator * (const T &r, const Vec3 &v)
    { return Vec3<T>(v.x * r, v.y * r, v.z * r); }
    friend Vec3 operator / (const T &r, const Vec3 &v)
    { return Vec3<T>(r / v.x, r / v.y, r / v.z); }

    friend std::ostream& operator << (std::ostream &s, const Vec3<T> &v)
    {
        return s << '[' << v.x << ' ' << v.y << ' ' << v.z << ']';
    }
    
    T x, y, z;
};

//[comment]
// Now you can specialize the class. We are just showing two examples here. In your code
// you can declare a vector either that way: Vec3<float> a, or that way: Vec3f a
//[/comment]
typedef Vec3<float> Vec3f;
typedef Vec3<int> Vec3i;

//[comment]
// Implementation of a generic 4x4 Matrix class - Same thing here than with the Vec3 class. It uses
// a template which is maybe less useful than with vectors but it can be used to
// define the coefficients of the matrix to be either floats (the most case) or doubles depending
// on our needs.
//
// To use you can either write: Matrix44<float> m; or: Matrix44f m;
//[/comment]
class Matrix44f
{
public:

    float x[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    Matrix44f() {}

    Matrix44f (float a, float b, float c, float d, float e, float f, float g, float h,
	       float i, float j, float k, float l, float m, float n, float o, float p)
    {
        x[0][0] = a;
        x[0][1] = b;
        x[0][2] = c;
        x[0][3] = d;
        x[1][0] = e;
        x[1][1] = f;
        x[1][2] = g;
        x[1][3] = h;
        x[2][0] = i;
        x[2][1] = j;
        x[2][2] = k;
        x[2][3] = l;
        x[3][0] = m;
        x[3][1] = n;
        x[3][2] = o;
        x[3][3] = p;
    }
    
    const float* operator [] (uint8_t i) const { return x[i]; }
    float* operator [] (uint8_t i) { return x[i]; }


    //[comment]
    // This method needs to be used for point-matrix multiplication. Keep in mind
    // we don't make the distinction between points and vectors at least from
    // a programming point of view, as both (as well as normals) are declared as Vec3.
    // However, mathematically they need to be treated differently. Points can be translated
    // when translation for vectors is meaningless. Furthermore, points are implicitly
    // be considered as having homogeneous coordinates. Thus the w coordinates needs
    // to be computed and to convert the coordinates from homogeneous back to Cartesian
    // coordinates, we need to divided x, y z by w.
    //
    // The coordinate w is more often than not equals to 1, but it can be different than
    // 1 especially when the matrix is projective matrix (perspective projection matrix).
    //[/comment]
  void multVecMatrix(const Vec3<float> &src, Vec3<float> &dst) const;
  

    //[comment]
    // Compute the inverse of the matrix using the Gauss-Jordan (or reduced row) elimination method.
    // We didn't explain in the lesson on Geometry how the inverse of matrix can be found. Don't
    // worry at this point if you don't understand how this works. But we will need to be able to
    // compute the inverse of matrices in the first lessons of the "Foundation of 3D Rendering" section,
    // which is why we've added this code. For now, you can just use it and rely on it
    // for doing what it's supposed to do. If you want to learn how this works though, check the lesson
    // on called Matrix Inverse in the "Mathematics and Physics of Computer Graphics" section.
    //[/comment]
  Matrix44f inverse() const;
  
  // \brief set current matrix to its inverse
  const Matrix44f & invert() {
    *this = inverse();
    return *this;
  }
  
    friend std::ostream& operator << (std::ostream &s, const Matrix44f &m)
    {
        std::ios_base::fmtflags oldFlags = s.flags();
        int width = 12; // total with of the displayed number
        s.precision(5); // control the number of displayed decimals
        s.setf (std::ios_base::fixed);
        
        s << "[" << std::setw (width) << m[0][0] <<
             " " << std::setw (width) << m[0][1] <<
             " " << std::setw (width) << m[0][2] <<
             " " << std::setw (width) << m[0][3] << "\n" <<
            
             " " << std::setw (width) << m[1][0] <<
             " " << std::setw (width) << m[1][1] <<
             " " << std::setw (width) << m[1][2] <<
             " " << std::setw (width) << m[1][3] << "\n" <<
            
             " " << std::setw (width) << m[2][0] <<
             " " << std::setw (width) << m[2][1] <<
             " " << std::setw (width) << m[2][2] <<
             " " << std::setw (width) << m[2][3] << "\n" <<
            
             " " << std::setw (width) << m[3][0] <<
             " " << std::setw (width) << m[3][1] <<
             " " << std::setw (width) << m[3][2] <<
             " " << std::setw (width) << m[3][3] << "]";
        
        s.flags (oldFlags);
        return s;
    }
};



