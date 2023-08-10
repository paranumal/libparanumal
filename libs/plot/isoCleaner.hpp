/*

The MIT License (MIT)

Copyright (c) 2017-2023 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


// NBN: concats 4 files + iso_1_base.hpp:
// ... /libs/plot/isoSurf/iso_2_vertex.hpp
// ... /libs/plot/isoSurf/iso_3_face.hpp
// ... /libs/plot/isoSurf/iso_4_mesh.hpp
// ... /libs/plot/isoSurf/iso_5_quadric.hpp

#pragma once

#include <algorithm>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <typeindex>
#include <functional>
#include <map>

#include <cassert>
#include <cstring>

// define arg types
#define ISO_COORD_TYPE_FP32   1
#define ISO_SCALAR_TYPE_FP32  1
#define ISO_QUALITY_IS_FP32   1

// simplify loop syntax
// #define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )
// #define loopj(start_l,end_l) for ( int j=start_l;j<end_l;++j )
// #define loopk(start_l,end_l) for ( int k=start_l;k<end_l;++k )

typedef unsigned int uint;
extern int g_procid;


#ifndef _MSC_VER

// TODO: expose simple logging scheme
// void nnLOG(int n, const char* format_str, ...);
// void nnMSG(int n, const char* format_str, ...);
// void nnTRC(int n, const char* format_str, ...);
// 
// void nnLOG(const stdS & msg, int n);
// void nnMSG(const stdS & msg, int n);
// void nnTRC(const stdS & msg, int n);

void nnMSG(int n, const char* format_str, ...) {

  int g_MSG_FLAG = 5;
  static char bufM[1024];
  if (n <= g_MSG_FLAG) {
    va_list arglist;
    va_start(arglist, format_str);
    int nUsed = -1;
    nUsed = vsnprintf  (bufL, 1023,       format_str, arglist);
  //nUsed = vsnprintf_s(bufL, 1023, 1000, format_str, arglist);
    assert(nUsed >= 0);
    va_end(arglist);

    fprintf(stdout, "%s", bufM);
    fflush(stdout);
  }
}

void nnMSG(const stdS& msg, int n) {
  nnMSG(n, "%s", msg.c_str());
}

#define nnLOG nnMSG
#define nnTRC nnMSG

#endif


namespace iso {

  class   MyVertex2;
  typedef MyVertex2  VertexType;
  typedef MyVertex2* VertexPointer;

  class   MyFace2;
  typedef MyFace2  FaceType;
  typedef MyFace2* FacePointer;

  typedef const MyVertex2* ConstVertexPointer;
  typedef const MyFace2* ConstFacePointer;

  typedef int FlagType;

  template <class P3ScalarType> class Point3;

  typedef Point3<int>	    Point3i;
  typedef Point3<float>   Point3f;
  typedef Point3<double>  Point3d;

  // scalar values used for penalizing the removal of edges
  typedef float ColorType;    // scalars from isosurf kernel
  typedef float QualityType;  // scalars for quadric area, error, ...

  //-------------------------------------------------------
  // various types used by MyVertex2 / MyFace2 / MyMesh2
  //-------------------------------------------------------
#if (ISO_COORD_TYPE_FP32)
  typedef iso::Point3f  CoordType;
#else
  typedef iso::Point3d  CoordType;
#endif

#if (ISO_SCALAR_TYPE_FP32)
  typedef float   ScalarType;
#else
  typedef double  ScalarType;
#endif


  // forward
  template <class T> class Box3;

  // Templated class for representing a point in 3D space.
  template <class P3ScalarType> class Point3
  {
  protected:

    P3ScalarType _v[3];   // The only data member

  public:

    typedef P3ScalarType ScalarType;

    Point3() { _v[0] = _v[1] = _v[2] = (P3ScalarType)0; }
    Point3(const P3ScalarType nx, const P3ScalarType ny, const P3ScalarType nz) { _v[0] = nx; _v[1] = ny; _v[2] = nz; }
    Point3(Point3 const& p) = default;

    // Copy from Point with different template
    template<class Q>
    Point3(Point3<Q> const& p) { _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2]; }
    Point3(const P3ScalarType nv[3]) { _v[0] = nv[0]; _v[1] = nv[1]; _v[2] = nv[2]; }

    // Default copy assignment
    Point3& operator = (Point3 const& p) {        // = default;
      _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2];
      return *this;
    }

    // Copy assignment from Point with different template
    template<class Q>
    Point3& operator= (Point3<Q> const& p) {
      _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2];
      return *this;
    }

    void SetZero() { _v[0] = 0; _v[1] = 0; _v[2] = 0; }
    void SetCoord(const P3ScalarType nx, const P3ScalarType ny, const P3ScalarType nz) { _v[0] = nx; _v[1] = ny; _v[2] = nz; }

    // Pad values outside [0..2] range with 0.
    P3ScalarType Ext(const int i) const {
      if (i >= 0 && i <= 2) return _v[i];
      else return 0;
    }

    template <class Q>
    void Import(const Point3<Q>& b) {
      _v[0] = P3ScalarType(b[0]); _v[1] = P3ScalarType(b[1]); _v[2] = P3ScalarType(b[2]);
    }

    // convert from Point3 of different scalar type {int, fp32, fp64}
    template <class Q>
    static inline Point3 Construct(const Point3<Q>& b) {
      return Point3(P3ScalarType(b[0]), P3ScalarType(b[1]), P3ScalarType(b[2]));
    }

    static inline Point3 Zero() { return Point3(0, 0, 0); }
    static inline Point3 One() { return Point3(1, 1, 1); }

    inline P3ScalarType& operator [] (const int i) {
      assert(i >= 0 && i < 3);
      return _v[i];
    }
    inline const P3ScalarType& operator [] (const int i) const {
      assert(i >= 0 && i < 3);
      return _v[i];
    }
    inline const P3ScalarType& X() const { return _v[0]; }
    inline const P3ScalarType& Y() const { return _v[1]; }
    inline const P3ScalarType& Z() const { return _v[2]; }
    inline P3ScalarType& X() { return _v[0]; }
    inline P3ScalarType& Y() { return _v[1]; }
    inline P3ScalarType& Z() { return _v[2]; }
    inline const P3ScalarType* V() const { return _v; }
    inline P3ScalarType* V() { return _v; }

    inline P3ScalarType& V(const int i) {
      assert(i >= 0 && i < 3);
      return _v[i];
    }
    inline const P3ScalarType& V(const int i) const {
      assert(i >= 0 && i < 3);
      return _v[i];
    }

    inline Point3 operator + (Point3 const& p) const { return Point3<P3ScalarType>(_v[0]+p._v[0], _v[1]+p._v[1], _v[2]+p._v[2]); }
    inline Point3 operator - (Point3 const& p) const { return Point3<P3ScalarType>(_v[0]-p._v[0], _v[1]-p._v[1], _v[2]-p._v[2]); }
    inline Point3 operator * (const P3ScalarType s) const { return Point3<P3ScalarType>(_v[0]*s, _v[1]*s, _v[2]*s); }
    inline Point3 operator / (const P3ScalarType s) const { return Point3<P3ScalarType>(_v[0]/s, _v[1]/s, _v[2]/s); }

    // Dot product
    inline P3ScalarType operator * (Point3 const& p) const { return (_v[0] * p._v[0] + _v[1] * p._v[1] + _v[2] * p._v[2]); }
    inline P3ScalarType dot(const Point3& p) const { return (*this) * p; }

    // Cross product
    inline Point3 operator ^ (Point3 const& p) const {
      return Point3 <P3ScalarType>(
        _v[1] * p._v[2] - _v[2] * p._v[1],
        _v[2] * p._v[0] - _v[0] * p._v[2],
        _v[0] * p._v[1] - _v[1] * p._v[0]);
    }

    inline Point3& operator += (Point3 const& p) { _v[0] += p._v[0]; _v[1] += p._v[1]; _v[2] += p._v[2]; return *this; }
    inline Point3& operator -= (Point3 const& p) { _v[0] -= p._v[0]; _v[1] -= p._v[1]; _v[2] -= p._v[2]; return *this; }
    inline Point3& operator *= (const P3ScalarType s) { _v[0] *= s; _v[1] *= s; _v[2] *= s; return *this; }
    inline Point3& operator /= (const P3ScalarType s) { _v[0] /= s; _v[1] /= s; _v[2] /= s; return *this; }

    inline P3ScalarType Norm()        const { return sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]); }
    inline P3ScalarType SquaredNorm() const { return     (_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]); }

    inline Point3& Scale(const P3ScalarType sx, const P3ScalarType sy, const P3ScalarType sz) {
      _v[0] *= sx; _v[1] *= sy; _v[2] *= sz; return *this;
    }
    inline Point3& Scale(const Point3& p) { 
      _v[0] *= p._v[0]; _v[1] *= p._v[1]; _v[2] *= p._v[2]; return *this;
    }

    inline Point3& Normalize() {
      P3ScalarType n = P3ScalarType(sqrt(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]));
      if (n > P3ScalarType(0)) { _v[0] /= n; _v[1] /= n; _v[2] /= n; }
      return *this;
    }
    inline void normalize() {
      this->Normalize();
    }
    inline Point3 normalized() const {
      Point3<P3ScalarType> p = *this;
      p.normalize();
      return p;
    }

    // defined below after Box3<T> declaration
    Box3<P3ScalarType> GetBBox(Box3<P3ScalarType>& bb) const;

    size_t MaxCoeffId() const {
      if (_v[0] > _v[1]) return _v[0] > _v[2] ? 0 : 2;
      else               return _v[1] > _v[2] ? 1 : 2;
    }

    // Comparison Operators 1:
    inline bool operator == (Point3 const& p) const {
      return _v[0] == p._v[0] && _v[1] == p._v[1] && _v[2] == p._v[2];
    }
    inline bool operator != (Point3 const& p) const {
      return _v[0] != p._v[0] || _v[1] != p._v[1] || _v[2] != p._v[2];
    }

    // Comparison Operators 2: note reversed z-priority ordering
    inline bool operator <  (Point3 const& p) const {
      return (_v[2] != p._v[2]) ? (_v[2] < p._v[2]) :
             (_v[1] != p._v[1]) ? (_v[1] < p._v[1]) :
             (_v[0] < p._v[0]);
    }
    inline bool operator >  (Point3 const& p) const {
      return (_v[2] != p._v[2]) ? (_v[2] > p._v[2]) :
             (_v[1] != p._v[1]) ? (_v[1] > p._v[1]) :
             (_v[0] > p._v[0]);
    }
    inline bool operator <= (Point3 const& p) const {
      return (_v[2] != p._v[2]) ? (_v[2] < p._v[2]) :
             (_v[1] != p._v[1]) ? (_v[1] < p._v[1]) :
             (_v[0] <= p._v[0]);
    }
    inline bool operator >= (Point3 const& p) const {
      return (_v[2] != p._v[2]) ? (_v[2] > p._v[2]) :
             (_v[1] != p._v[1]) ? (_v[1] > p._v[1]) :
             (_v[0] >= p._v[0]);
    }
    inline Point3 operator - () const {
      return Point3<P3ScalarType>(-_v[0], -_v[1], -_v[2]);
    }

  }; // end class Point3


  //-------------------------------------
  // miscellaneous Point3<T> routines:
  //-------------------------------------

  template <class P3ScalarType>
  inline P3ScalarType Norm(Point3<P3ScalarType> const& p) { return p.Norm(); }

  template <class P3ScalarType>
  inline P3ScalarType SquaredNorm(Point3<P3ScalarType> const& p) { return p.SquaredNorm(); }

  template <class P3ScalarType>
  inline P3ScalarType Distance(Point3<P3ScalarType> const& p1, Point3<P3ScalarType> const& p2) { return (p1 - p2).Norm(); }

  //-------------------------------------
  // miscellaneous 3D routines:
  //-------------------------------------

  // 2D plane in 3D space, stored as a (Point3) direction + scalar offset
  template <class T, bool NORM = true> class Plane3 {
  public:
    typedef T ScalarType;
    typedef Point3<T> PointType;

  private:
    ScalarType _offset;   // Distance from origin
    PointType _dir;       // Direction (not normalized unless NORM is true)

  public:
    Plane3() {}
    Plane3(const ScalarType& dist, const PointType& dir) { Set(dist, dir); }

    template <class Q>
    inline void Import(const Plane3<Q, false>& b) {
      _offset = ScalarType(b.Offset());
      _dir = Point3<T>::Construct(b.Direction());
    }

    const ScalarType& Offset() const { return _offset; }
    ScalarType& Offset() { return _offset; }
    void SetOffset(const ScalarType& o) { _offset = o; }

    const PointType& Direction() const { return _dir; }
    void SetDirection(const PointType& dir) {
      _dir = dir;
      if (NORM) _dir.Normalize();
    }

    void Set(const ScalarType& off, const PointType& dir) {
      if (NORM) {
        const ScalarType normFactor = dir.Norm();
        this->_dir = dir / normFactor;
        this->_offset = off / normFactor;
      }
      else {
        this->_offset = off;
        this->_dir = dir;
      }
    }
    void Set(const PointType& dir, const ScalarType& off) { Set(off, dir); }

    bool operator == (Plane3 const& p) const { return _offset == p._offset && _dir == p._dir; }
    bool operator != (Plane3 const& p) const { return _offset != p._offset || _dir != p._dir; }

    // Project a point on the plane
    PointType Projection(const PointType& p) const {
      ScalarType k = p.dot(_dir) - _offset;
      return p - _dir * k;
    }

    // Mirror the point wrt the plane
    PointType Mirror(const PointType& p) const {
      PointType mirr = Projection(p);
      mirr += mirr - p;
      return mirr;
    }

    void Normalize() { _dir.Normalize(); }

    // define the plane passing through three points
    void Init(const PointType& p0, const PointType& p1, const PointType& p2) {
      _dir = (p2 - p0) ^ (p1 - p0);
      if (NORM) Normalize();
      _offset = p0.dot(_dir);
    }

    // define the plane passing through a point with normal
    void Init(const PointType& p0, const PointType& norm) {
      _dir = norm;
      if (NORM) Normalize();
      _offset = p0.dot(_dir);
    }
  };	// end class Plane3

  typedef Plane3<float>  Plane3f;
  typedef Plane3<double> Plane3d;


  // axis-aligned 3D bounding box, stored as two Point3<T>. 
  template <class BoxScalarType> class Box3
  {
  public:

    typedef BoxScalarType ScalarType;

    Point3<BoxScalarType> min;
    Point3<BoxScalarType> max;

    inline Box3() { this->SetNull(); }
    inline Box3(const Point3<BoxScalarType>& mi, const Point3<BoxScalarType>& ma) { min = mi; max = ma; }
    inline ~Box3() { }

    inline bool operator == (const Box3<BoxScalarType>& p) const { return min == p.min && max == p.max; }
    inline bool operator != (const Box3<BoxScalarType>& p) const { return min != p.min || max != p.max; }

    // Expand the box by a vector delta.
    void Offset(const BoxScalarType s) { Offset(Point3<BoxScalarType>(s, s, s)); }
    void Offset(const Point3<BoxScalarType>& delta) { min -= delta; max += delta; }

    void Set(const Point3<BoxScalarType>& p) { min = max = p; }

    // Set the bounding box to a "null" value
    void SetNull() {
      min.X() = 1; max.X() = -1;
      min.Y() = 1; max.Y() = -1;
      min.Z() = 1; max.Z() = -1;
    }

    // Expand the current bbox to include new point
    void Add(const Point3<BoxScalarType>& p) {
      if (IsNull()) Set(p);
      else {
        if (min.X() > p.X()) min.X() = p.X();
        if (min.Y() > p.Y()) min.Y() = p.Y();
        if (min.Z() > p.Z()) min.Z() = p.Z();

        if (max.X() < p.X()) max.X() = p.X();
        if (max.Y() < p.Y()) max.Y() = p.Y();
        if (max.Z() < p.Z()) max.Z() = p.Z();
      }
    }

    void Intersect(const Box3<BoxScalarType>& b) {
      if (min.X() < b.min.X()) min.X() = b.min.X();
      if (min.Y() < b.min.Y()) min.Y() = b.min.Y();
      if (min.Z() < b.min.Z()) min.Z() = b.min.Z();

      if (max.X() > b.max.X()) max.X() = b.max.X();
      if (max.Y() > b.max.Y()) max.Y() = b.max.Y();
      if (max.Z() > b.max.Z()) max.Z() = b.max.Z();

      if (min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z()) {
        SetNull();    // no intersection with b
      }
    }

    BoxScalarType Diag() const { return Distance(min, max); }
    bool IsNull() const { return min.X() > max.X() || min.Y() > max.Y() || min.Z() > max.Z(); }

    // Check if vertex lies inside a bounding box.
    // Note: point is NOT inside if it's ON a bounding face.
    bool Contains(const CoordType& p) const {
      return (min.X() < p.X() && max.X() > p.X() &&
              min.Y() < p.Y() && max.Y() > p.Y() &&
              min.Z() < p.Z() && max.Z() > p.Z());
    }

    inline BoxScalarType DimX() const { return max.X() - min.X(); }
    inline BoxScalarType DimY() const { return max.Y() - min.Y(); }
    inline BoxScalarType DimZ() const { return max.Z() - min.Z(); }

  }; // end class Box3


  template <class T>
  Box3<T> Point3<T>::GetBBox(Box3<T>& bb) const {
    bb.Set(*this);
    return bb;
  }

  typedef Box3<short>   Box3s;
  typedef Box3<int>     Box3i;
  typedef Box3<float>   Box3f;
  typedef Box3<double>  Box3d;


  // Simple class abstracting a grid of voxels in a 3d space:
  // - bbox defines the real extent of the grid in 3d space;
  // - siz is the number of cells for each side
  //
  // OBJTYPE:      type of the indexed objects.
  // SCALARTYPE:   type of scalars for structure's internal data 
  //               (this may differ from the object's scalar type).

  template <class SCALARTYPE>
  class BasicGrid
  {
  public:

    typedef SCALARTYPE            ScalarType;
    typedef Box3<ScalarType>      Box3x;
    typedef Point3<ScalarType>    CoordType;
    typedef BasicGrid<SCALARTYPE> GridType;

    Box3x bbox;

    CoordType dim;    // Spatial extent of the bounding box
    Point3i   siz;    // Number of voxels in each dimension of the grid
    CoordType voxel;  // Dimensions of a single cell/voxel

    // Calculate 3d extent of both the entire grid and a single voxel
    void ComputeDimAndVoxel() {
      this->dim = this->bbox.max - this->bbox.min;
      this->voxel[0] = this->dim[0] / this->siz[0];
      this->voxel[1] = this->dim[1] / this->siz[1];
      this->voxel[2] = this->dim[2] / this->siz[2];
    }

    // Given a 3D point, return integer coords of cell containing the point
    void PToIP(const CoordType& p, Point3i& pi) const {
      CoordType t = p - bbox.min;
      pi[0] = int(t[0] / voxel[0]);
      pi[1] = int(t[1] / voxel[1]);
      pi[2] = int(t[2] / voxel[2]);
    }

    // Given a 3D point, return integer coords of cell containing the point
    Point3i GridP(const Point3<ScalarType>& p) const {
      Point3i pi;
      PToIP(p, pi);
      return pi;
    }

    // Given a cell index return the lower corner of the cell
    template <class OtherScalarType>
    void IPiToPf(const Point3i& pi, Point3<OtherScalarType>& p) const {
      p[0] = bbox.min[0] + ((OtherScalarType)pi[0]) * voxel[0];
      p[1] = bbox.min[1] + ((OtherScalarType)pi[1]) * voxel[1];
      p[2] = bbox.min[2] + ((OtherScalarType)pi[2]) * voxel[2];
    }

    // Given a cell index return the corresponding box
    void IPiToBox(const Point3i& pi, Box3x& b) const {
      CoordType p;
      p[0] = ((ScalarType)pi[0]) * voxel[0];
      p[1] = ((ScalarType)pi[1]) * voxel[1];
      p[2] = ((ScalarType)pi[2]) * voxel[2];
      p += bbox.min;
      b.min = p;
      b.max = (p + voxel);
    }

    // Given a cell index return the center of the cell itself
    void IPiToBoxCenter(const Point3i& pi, CoordType& c) const {
      CoordType p;
      IPiToPf(pi, p);
      c = p + voxel / ScalarType(2.0);
    }

    // Same as IPiToPf but handles different types of Point3<T>.
    template <class OtherScalarType>
    void IPfToPf(const Point3<OtherScalarType>& pi, Point3<OtherScalarType>& p) const {
      p[0] = ((OtherScalarType)pi[0]) * voxel[0] + bbox.min[0];
      p[1] = ((OtherScalarType)pi[1]) * voxel[1] + bbox.min[1];
      p[2] = ((OtherScalarType)pi[2]) * voxel[2] + bbox.min[2];
    }

    // Given a cell in <ScalarType> coordinates, compute 
    // the corresponding cell in integer coordinates
    void BoxToIBox(const Box3x& b, Box3i& ib) const {
      PToIP(b.min, ib.min);
      PToIP(b.max, ib.max);
    }

    // Given a cell in integer coordinates, compute the 
    // corresponding cell in <ScalarType> coordinates
    void IBoxToBox(const Box3i& ib, Box3x& b) const {
      IPiToPf(ib.min, b.min);
      IPiToPf(ib.max + Point3i(1, 1, 1), b.max);
    }
  };  // end class BasicGrid


  template<class scalar_type>
  void BestDim(const Box3<scalar_type> box, const scalar_type voxel_size, Point3i& dim)
  {
    Point3<scalar_type> box_size = box.max - box.min;
    int64_t elem_num = (int64_t)(box_size[0] / voxel_size + 0.5) *
                       (int64_t)(box_size[1] / voxel_size + 0.5) *
                       (int64_t)(box_size[2] / voxel_size + 0.5);
    BestDim(elem_num, box_size, dim);
  }


  template<class scalar_type>
  void BestDim(const int64_t elems, const Point3<scalar_type>& size, Point3i& dim)
  {
    // Calculates the size of the grid as a function of 
    // the bounding box ratio and the number of elements

    const int64_t mincells = 1;   // minimum number of cells
    const double GFactor = 1;     // #cells = #elems * GFactor
    double diag = size.Norm();    // box diagonal
    double eps = diag * 1e-4;     // tolerance factor

    assert(elems > 0);
    assert(size[0] >= 0.0);
    assert(size[1] >= 0.0);
    assert(size[2] >= 0.0);


    int64_t ncell = (int64_t)(elems * GFactor);	// Calc number of voxels
    if (ncell < mincells) ncell = mincells;

    dim[0] = dim[1] = dim[2] = 1;

    if (size[0] > eps) {
      if (size[1] > eps) {
        if (size[2] > eps) {
          double k = pow((double)(ncell / (size[0]*size[1]*size[2])), (1.0/3.0));
          dim[0] = int(size[0] * k);
          dim[1] = int(size[1] * k);
          dim[2] = int(size[2] * k);
        }
        else {
          dim[0] = int(sqrt(ncell * size[0] / size[1]));
          dim[1] = int(sqrt(ncell * size[1] / size[0]));
        }
      }
      else {
        if (size[2] > eps) {
          dim[0] = int(sqrt(ncell * size[0] / size[2]));
          dim[2] = int(sqrt(ncell * size[2] / size[0]));
        }
        else {
          dim[0] = int(ncell);
        }
      }
    }
    else {
      if (size[1] > eps) {
        if (size[2] > eps) {
          dim[1] = int(sqrt(ncell * size[1] / size[2]));
          dim[2] = int(sqrt(ncell * size[2] / size[1]));
        }
        else {
          dim[1] = int(ncell);
        }
      }
      else if (size[2] > eps) {
        dim[2] = int(ncell);
      }
    }
    dim[0] = std::max(dim[0], 1);
    dim[1] = std::max(dim[1], 1);
    dim[2] = std::max(dim[2], 1);
  }


  template <class VertHash, class VertCtnr>
  uint GridGetInBox(VertHash& _Si,
                    const iso::Box3<typename VertHash::ScalarType>& _bbox,
                    VertCtnr& _objectPtrs)
  {
    // return a list of pointers to all vertices 
    // that lie within the specified bbox

    typename VertHash::CellIterator first, last, l;
    _objectPtrs.clear();
    iso::Box3i ibbox;
    Box3i Si_ibox(Point3i(0, 0, 0), _Si.siz - Point3i(1, 1, 1));
    _Si.BoxToIBox(_bbox, ibbox);
    ibbox.Intersect(Si_ibox);
    if (ibbox.IsNull()) {
      return 0;
    }
    else {
      int ix, iy, iz;
      for (ix = ibbox.min[0]; ix <= ibbox.max[0]; ix++) {
        for (iy = ibbox.min[1]; iy <= ibbox.max[1]; iy++) {
          for (iz = ibbox.min[2]; iz <= ibbox.max[2]; iz++) {
            _Si.Grid(ix, iy, iz, first, last);
            for (l = first; l != last; ++l) {
              if (!(**l).IsD()) {

                typename VertHash::ObjPtr elem = &(**l);

                if (_bbox.Contains(elem->cP())) {
                  _objectPtrs.push_back(elem);
                }
              }
            }
          }
        }
      }
      return (static_cast<uint>(_objectPtrs.size()));
    }
  }

  // (1) HashFunctor 
  // (2) see also std::hash<CellType*> in class SimpleTri

  struct HashFunctor {

    size_t operator()(const Point3i& p) const {
      const size_t _HASH_P0 = 73856093u;
      const size_t _HASH_P1 = 19349663u;
      const size_t _HASH_P2 = 83492791u;

      return size_t(p.V(0)) * _HASH_P0 ^ size_t(p.V(1)) * _HASH_P1 ^ size_t(p.V(2)) * _HASH_P2;
    }

    bool operator()(const Point3i& s1, const Point3i& s2) const {
      return (s1 < s2);   // test if s1 ordered before s2
    }
  };

  // Spatial Hash Table for Hashing as described in:
  // "Optimized Spatial Hashing for Collision Detection of Deformable Objects", 
  // Matthias Teschner, Bruno Heidelberger, Matthias Muller, 
  // Danat Pomeranets, Markus Gross

  template <typename ObjType, class FLT = double>
  class SpatialHashTable : public BasicGrid<FLT>
  {
  public:
    typedef SpatialHashTable                SpatialHashType;
    typedef ObjType* ObjPtr;
    typedef Point3<ScalarType>              CoordType;
    typedef typename BasicGrid<FLT>::Box3x  Box3x;

    // Hash table definition.
    // The hash indexes into a virtual grid structure, using a
    // std::multimap to store multiple faces in each grid cell.

    typedef typename std::unordered_multimap<Point3i, ObjType*, HashFunctor> HashType;
    typedef typename HashType::iterator     HashIterator;
    typedef typename HashType::value_type   HashValue;

    HashType hash_table;

    bool    Empty() const { return hash_table.empty(); }
    size_t  CellSize(const Point3i& cell) { return hash_table.count(cell); }
    bool    EmptyCell(const Point3i& cell) const { return hash_table.find(cell) == hash_table.end(); }

    template <class OBJITER>
    void Set(const OBJITER& _oBegin, const OBJITER& _oEnd, const Box3x& _bbox = Box3x())
    {
      // Insert a 3d mesh into a grid of voxels.
      Box3x& bbox = this->bbox;
      CoordType& dim = this->dim;
      Point3i& siz = this->siz;
      CoordType& voxel = this->voxel;

      Box3x b; OBJITER i;
      int _size = (int)std::distance<OBJITER>(_oBegin, _oEnd);
      if (!_bbox.IsNull()) {
        this->bbox = _bbox;
      }
      else {
        for (i = _oBegin; i != _oEnd; ++i) {
          this->bbox.Add((*i).cP());        // add each point to bounding box
        }
        bbox.Offset(bbox.Diag() / 100.0);   // inflate the calculated bounding box
      }

      dim = bbox.max - bbox.min;
      BestDim(_size, dim, siz);

      voxel[0] = dim[0] / siz[0];           // find voxel size
      voxel[1] = dim[1] / siz[1];
      voxel[2] = dim[2] / siz[2];

      for (i = _oBegin; i != _oEnd; ++i) {
        // convert {x,y,z} coords to {i,j,k} voxel
        Point3i cell; BasicGrid<FLT>::PToIP(i->cP(), cell);
        hash_table.insert(HashValue(cell, &(*i)));  // add vertex to voxel
      }
    }

    // Class to abstract a HashIterator
    struct CellIterator {
      CellIterator() {}
      HashIterator t;
      ObjPtr& operator *() { return (t->second); }
      ObjPtr  operator *() const { return (t->second); }
      bool operator != (const CellIterator& p) const { return t != p.t; }
      void operator ++() { t++; }
    };

    // return the simplexes on a specified cell
    void Grid(int x, int y, int z, CellIterator& first, CellIterator& last) {
      this->Grid(iso::Point3i(x, y, z), first, last);
    }

    // return the simplexes on a specified cell
    void Grid(const Point3i& _c, CellIterator& first, CellIterator& end) {
      std::pair<HashIterator, HashIterator> CellRange = hash_table.equal_range(_c);
      first.t = CellRange.first;
      end.t = CellRange.second;
    }

  }; // end class SpatialHashTable


  /////////////////////////////////////////////////////////
  // 
  // start: isoSurf/iso_2_vertex.hpp
  // 
  /////////////////////////////////////////////////////////

  // forward
  bool CholSolve(double a[6], double b[3], Point3d& x);

  class Quadric
  {
    // Encode a quadric function: f(x) = xAx + bx + c
    // where A is a symmetric 3x3 matrix, 
    //       b a vector and 
    //       c a scalar constant.
    //
    // Note: quadric calculations done in double precision, 
    //       but args can be either float or double

  public:

    double a[6];  // Symmetric Matrix 3x3 : a11 a12 a13 a22 a23 a33
    double b[3];  // Vector r3
    double c;     // Scalar (-1 means null/un-initialized quadric)
    double m_val; // store value calculated by Apply()

    Quadric() : c(-1.0), m_val(0) {}

    bool IsValid() const { return c >= 0; }
    void SetInvalid() { c = -1.0; }

    // Initialize quadric with squared distance from Plane p
    template <class PlaneType>
    void ByPlane(const PlaneType& p) {
      const PlaneType::PointType& pD = p.Direction();

      a[0] = pD[0] * pD[0];	// a11
      a[1] = pD[1] * pD[0];	// a12 (=a21)
      a[2] = pD[2] * pD[0];	// a13 (=a31)
      a[3] = pD[1] * pD[1];	// a22
      a[4] = pD[2] * pD[1];	// a23 (=a32)
      a[5] = pD[2] * pD[2];	// a33

      double pOff = p.Offset();
      b[0] = (-2.0) * pOff * pD[0];
      b[1] = (-2.0) * pOff * pD[1];
      b[2] = (-2.0) * pOff * pD[2];
      c = pOff * pOff;
    }

    void SetZero() {
      a[0] = a[1] = a[2] = a[3] = a[4] = a[5] = 0.0;
      b[0] = b[1] = b[2] = 0.0;
      c = 0.0;
    }

    void operator = (const Quadric& q) {
      assert(q.IsValid());
      a[0] = q.a[0]; a[1] = q.a[1]; a[2] = q.a[2];
      a[3] = q.a[3]; a[4] = q.a[4]; a[5] = q.a[5];
      b[0] = q.b[0]; b[1] = q.b[1]; b[2] = q.b[2];
      c = q.c;
    }

    void operator += (const Quadric& q) {
      assert(IsValid() && q.IsValid());
      a[0] += q.a[0]; a[1] += q.a[1]; a[2] += q.a[2];
      a[3] += q.a[3]; a[4] += q.a[4]; a[5] += q.a[5];
      b[0] += q.b[0]; b[1] += q.b[1]; b[2] += q.b[2];
      c += q.c;
    }

    void operator *= (const double w) {
      assert(IsValid());
      a[0] *= w; a[1] *= w; a[2] *= w;
      a[3] *= w; a[4] *= w; a[5] *= w;
      b[0] *= w; b[1] *= w; b[2] *= w;
      c *= w;
    }

    // Evaluate quadric at point p.
    template <class T>
    T Apply(const Point3<T>& p) {
      assert(IsValid());

      m_val = (p[0]*p[0]*a[0]
             + 2*p[0]*p[1]*a[1] + 2*p[0]*p[2]*a[2] + p[0]*b[0] + p[1]*p[1]*a[3]
             + 2*p[1]*p[2]*a[4] + p[1]*b[1] + p[2]*p[2]*a[5] + p[2]*b[2]
             + c);

      return m_val;
    }

    static double& RelativeErrorThr() { static double _err = 1.0e-6; return _err; }

    // return determinant of symmetric A(3,3)
    double sym_det() {
      double dd = (a[0]*a[3]*a[5]) + (a[1]*a[4]*a[2]) + (a[2]*a[1]*a[4])
                - (a[2]*a[3]*a[2]) - (a[4]*a[4]*a[0]) - (a[5]*a[1]*a[1]);
      return dd;
    }

    bool Minimum(Point3d& x) {

      // Optimal vertex location is given by minimizing the quadric:
      // 
      //      xAx + bx + c 
      // 
      // A(3,3) is PSD/invertible, so we can minimize the 
      // quadric by solving the first derivative for x:
      // 
      //     2 Ax + b = 0 
      //       Ax     = -b/2
      // 
      // Note: solving with b2 = (-b/2.0);

      double b2[3] = { -b[0]/2.0, -b[1]/2.0, -b[2]/2.0 };
      bool ok = CholSolve(this->a, b2, x);
      return ok;
    }

  };


  class MyVertex2
  {
  public:

    iso::Quadric m_q;
    iso::Quadric& Qd() { return m_q; }  // return quadric for this vertex

    MyVertex2() {
      m_flags = 0;      // bitflags: track if (Deleted/Visited/Border) etc. 
      m_fp = nullptr;   // VFAdj: nullptr ==> face pointer not initialized
      m_zp = -1;        //         -1     ==> edge index not initialized
      m_color = 0;
      m_cost = 0;
      m_imark = 0;
      m_coord.SetZero();
    }

    void CopyVertex(const MyVertex2& R) {

      // used by PermutateVertexVector / CompactVertexVector
      // TODO: check which data needs to be copied

      IMark() = R.cIMark();     // copy Mark
      Flags() = R.cFlags();     // copy Flags
      C() = R.cC();             // copy scalar ("Color")
      Q() = R.cQ();             // copy variance/cost ("Quality")
      P() = R.cP();             // copy {x,y,z} coords

      if (IsVFInitialized()) {
        VFp() = R.cVFp();       // copy Vector-Face info
        VFi() = R.cVFi();
      } 
      else {
        VFClear();
      }
    }

  public:
    // BitFlags
          int& Flags()       { return m_flags; }
    const int& Flags() const { return m_flags; }
    const int cFlags() const { return m_flags; }
    static bool HasFlags()   { return true; }
  private:
    int  m_flags;

  public:
    // Per-vertex Vertex-Face adjacency ("VFAdj").
    // Stores a pointer to first face of a list of faces.
    // Assumes Face uses corresponding face::VFAdj component.
    FacePointer& VFp()       { return m_fp; }
    FacePointer cVFp() const { return m_fp; }

    char& VFi()       { return m_zp; }
    char  VFi() const { return m_zp; }
    char cVFi() const { return m_zp; }

    bool IsNull() const { return m_zp == -1; }
    bool IsVFInitialized() const { return cVFi() != -1; }
    void VFClear() { m_fp = nullptr;  m_zp = -1; }

  private:
    FacePointer m_fp;   // nullptr ==> face pointer not initialized
    char        m_zp;   //  -1     ==> edge index not initialized


  public:
    // Geometric Position of the vertex (Point3f or Point3d)
          CoordType& P()        { return m_coord; }
    const CoordType& P() const  { return m_coord; }
    const CoordType& cP() const { return m_coord; }

  private:
    CoordType m_coord;

  public:
    // "Color" : a scalar value for this vertex
    const ColorType& C() const { return m_color; }
          ColorType& C()       { return m_color; }
          ColorType cC() const { return m_color; }

    static bool HasColor() { return true; }
    bool IsColorEnabled() const { return true; }

  private:
    ColorType m_color;    // e.g.: rho, |v|, qw

  public:
    // Per vertex "quality": stores the "cost" of removing an edge.
    const QualityType& Q() const { return m_cost; }
          QualityType& Q()       { return m_cost; }
          QualityType cQ() const { return m_cost; }
  private:
    QualityType m_cost;   // typedef'd as <short/fp32/fp64>

  public:
    // Per vertex Incremental Mark
    const int& IMark() const { return m_imark; }
          int& IMark()       { return m_imark; }
    int       cIMark() const { return m_imark; }
    void   InitIMark()       { m_imark = 0; }

  private:
    int m_imark;

  public:

    // Manage bitflags (MyVertex2)
    enum {                  // flags used:
      DELETED  = 0x0001,    // Vertex is deleted from the mesh
      NOTREAD  = 0x0002,    // Vertex of the mesh is not readable
      NOTWRITE = 0x0004,    // Vertex is not writable
      MODIFIED = 0x0008,    // Vertex is modified
      VISITED  = 0x0010,    // Vertex has been visited
      BORDER   = 0x0100,    // Vertex is on Border
      USER0    = 0x0200     // First user bit
    };

    bool IsD() const { return (this->cFlags() & DELETED) != 0; }      //  is vertex deleted?
    bool IsR() const { return (this->cFlags() & NOTREAD) == 0; }      //  is vertex readable?
    bool IsW() const { return (this->cFlags() & NOTWRITE) == 0; }     //  is vertex modifiable?
    bool IsRW()const { return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0; } // is vertex read/write?
    bool IsB() const { return (this->cFlags() & BORDER) != 0; }       //  is vertex on a border?
    bool IsV() const { return (this->cFlags() & VISITED) != 0; }      //  has vertex been visited?

    void SetFlags(int flagp) { this->Flags() = flagp; }
    void ClearFlags() { this->Flags() = 0; }

    void SetD()   { this->Flags() |= DELETED; }     // delete vertex from mesh
    void ClearW() { this->Flags() |= NOTWRITE; }    // mark vertex as writable
    void SetW()   { this->Flags() &= (~NOTWRITE); } // mark vertex as not writable
    void SetB()   { this->Flags() |= BORDER; }      // mark vertex as border
    void ClearB() { this->Flags() &= ~BORDER; }     // mark vertex as not border
    void SetV()   { this->Flags() |= VISITED; }     // mark vertex as visited
    void ClearV() { this->Flags() &= ~VISITED; }    // mark vertex as note visited

    // Unused flags (MyVertex2)
    //void ClearD() { this->Flags() &= (~DELETED); }  // un-delete vertex
    //void SetR()   { this->Flags() &= (~NOTREAD); }  // mark vertex as readable
    //void ClearR() { this->Flags() |= NOTREAD; }     // mark vertex as not readable


    // Return first "unused" user-bit
    static int& FirstUnusedBitFlag() { static int b = USER0; return b; }

    // Allocate a bit that can be used by user. Updates FirstUnusedBitFlag.
    static int NewBitFlag() {
      int bitForTheUser = FirstUnusedBitFlag();
      FirstUnusedBitFlag() = FirstUnusedBitFlag() << 1;
      return bitForTheUser;
    }

    // De-allocate a pre allocated bit. Updates FirstUnusedBitFlag.
    // Note: deallocate bits in reverse order of allocation (as in a stack)
    static bool DeleteBitFlag(int bitval) {
      if (FirstUnusedBitFlag() >> 1 == bitval) {
        FirstUnusedBitFlag() = FirstUnusedBitFlag() >> 1;
        return true;
      }
      assert(0);        // should not get here
      return false;
    }

    bool IsUserBit(int userBit)    { return (this->Flags() & userBit) != 0; }
    void SetUserBit(int userBit)   {         this->Flags() |= userBit; }
    void ClearUserBit(int userBit) {         this->Flags() &= (~userBit); }

    void GetBBox(iso::Box3<ScalarType>& bb) const { bb.Set(this->cP()); }

  };  // end class MyVertex2


  bool CholSolve(double a[6], double b[3], Point3d& x)
  {
    // Cholesky decomposition: M[3][3] = LL'
    // Off-diag elements in L overwrite lower tri of M,
    // diagonal elements of L are stored in vector p

    double M[3][3] = {{ a[0], a[1], a[2] },
                      {  0.0, a[3], a[4] },
                      {  0.0,  0.0, a[5] }};

    double p[3] = { 0.0 }; // store diagonal of L
    double sum = 0.0;

    //-------------------------------------
    // Cholesky decomposition: A = LL'
    //-------------------------------------
    {
      if (M[0][0] <= 0.0) return false;     // i = 0; check SPD
      p[0] = sqrt(M[0][0]);
      M[1][0] = M[0][1] / p[0];
      M[2][0] = M[0][2] / p[0];

      sum = M[1][1] - (M[1][0]*M[1][0]);
      if (sum <= 0.0) return false;         // i = 1; check SPD
      p[1] = sqrt(sum);
      sum = M[1][2] - (M[1][0]*M[2][0]);
      M[2][1] = sum / p[1];

      sum = M[2][2] - (M[2][0]*M[2][0]) - (M[2][1] * M[2][1]);
      if (sum <= 0.0) return false;         // i = 2; check SPD
      p[2] = sqrt(sum);
    }

    //-------------------------------------
    // Solve Ax = b
    //-------------------------------------
    {
      // 1. solve Ly = b
      //------------------
      x[0] = b[0] / p[0];
      sum  = b[1] - (M[1][0]*x[0]);                   // i = 0;
      x[1] = sum / p[1];                              // i = 1;
      sum  = b[2] - (M[2][0]*x[0]) - (M[2][1]*x[1]);  // i = 2;
      x[2] = sum / p[2];

      // 2. solve L'x = y
      //------------------
      x[2] /= p[2];
      sum  = x[1] - (M[2][1]*x[2]);                   // i = 2;
      x[1] = sum / p[1];                              // i = 1;
      sum  = x[0] - (M[1][0]*x[1]) - (M[2][0]*x[2]);  // i = 0;
      x[0] = sum / p[0];
    }

    // Check:
    // double relative_error = (M*x - b).norm() / b.norm();
    // if (relative_error > RelativeErrorThr()) return false;

    return true;
  }


  /////////////////////////////////////////////////////////
  // 
  // start: isoSurf/iso_3_face.hpp
  // 
  /////////////////////////////////////////////////////////


  class MyFace2
  {
  public:

    MyFace2() {
      m_flags = 0;                              // bitflags: track if (Deleted/Visited/Border) etc. 
      m_vfp[0] = m_vfp[1] = m_vfp[2] = nullptr; // VFAdj: nullptr ==> face pointers not initialized
      m_vfi[0] = m_vfi[1] = m_vfi[2] = -1;      //         -1     ==> edge ids not initialized
      m_ffp[0] = m_ffp[1] = m_ffp[2] = nullptr; // FFAdj: nullptr == not initialized
      m_ffi[0] = m_ffi[1] = m_ffi[2] = -1;
      v[0] = v[1] = v[2] = nullptr;         // 3 vertex pointers: nullptr ==> not initialized
      m_imark = 0;
    }

    void CopyFace(const MyFace2& R) {

      // used by CompactFaceVector()
      // TODO: check which data needs to be copied

      IMark() = R.cIMark();         // copy Mark
      Flags() = R.cFlags();         // copy Flags
      for (int i=0; i<3; ++i) {
        V(i) = R.cV(i);             // copy Vertex pointers
      }
      for (int i=0; i<3; ++i) {
        if (IsVFInitialized(i)) {
          VFp(i) = R.cVFp(i);       // copy Vector-Face info
          VFi(i) = R.cVFi(i);
        }
        else {
          VFClear(i);
        }
      }
      for (int i=0; i<3; ++i) {
        if (IsFFInitialized(i)) {
          FFp(i) = R.cFFp(i);       // copy Face-Face info
          FFi(i) = R.cFFi(i);
        }
        else {
          FFClear(i);
        }
      }
    }

  public:
    // BitFlags
    int& Flags() { return m_flags; }
    const int& Flags() const { return m_flags; }
    const int cFlags() const { return m_flags; }
    static bool HasFlags() { return true; }
  private:
    int  m_flags;

  public:
    // Per Face Vertex-Face adjacency ("VFAdj")
    // Stores a pointer to next face of list of faces incident on a vertex.
    // Assumes Vertex uses corresponding vertex::VFAdj component.
    FacePointer& VFp(const int j)       { assert(j >= 0 && j < 3); return m_vfp[j]; }
    FacePointer  VFp(const int j) const { assert(j >= 0 && j < 3); return m_vfp[j]; }
    FacePointer cVFp(const int j) const { assert(j >= 0 && j < 3); return m_vfp[j]; }

    char& VFi(const int j)       { return m_vfi[j]; }
    char  VFi(const int j) const { return m_vfi[j]; }
    char cVFi(const int j) const { return m_vfi[j]; }

    bool IsVFInitialized(const int j) const { return cVFi(j) != -1; }
    void VFClear(int j) { VFp(j) = nullptr; VFi(j) = -1; }

  private:
    FacePointer m_vfp[3];
    char        m_vfi[3];

  public:
    // Per Face "Face-Face" adjacency ("FFAdj")
    // 
    // Encodes the adjacency of faces through edges; for 2-manifold 
    // edges it just points to the other face. For non manifold edges 
    // (where more than 2 faces share the same edge) stores a pointer 
    // to the next face of the ring of faces incident on a edge.
    // 
    // Note: border faces point to themselves.

    FacePointer& FFp(const int j)       { assert(j >= 0 && j < 3); return m_ffp[j]; }
    FacePointer  FFp(const int j) const { assert(j >= 0 && j < 3); return m_ffp[j]; }
    FacePointer cFFp(const int j) const { assert(j >= 0 && j < 3); return m_ffp[j]; }

    char& FFi(const int j) { return m_ffi[j]; }
    char  FFi(const int j) const { return m_ffi[j]; }
    char cFFi(const int j) const { return m_ffi[j]; }

    FacePointer& FFp1(const int j)       { return FFp((j + 1) % 3); }
    FacePointer& FFp2(const int j)       { return FFp((j + 2) % 3); }
    FacePointer  FFp1(const int j) const { return FFp((j + 1) % 3); }
    FacePointer  FFp2(const int j) const { return FFp((j + 2) % 3); }
    FacePointer cFFp1(const int j) const { return FFp((j + 1) % 3); }
    FacePointer cFFp2(const int j) const { return FFp((j + 2) % 3); }

    bool IsFFInitialized(const int j) const { return cFFp(j) != nullptr; }
    void FFClear(int j) { FFp(j) = nullptr; FFi(j) = -1; } // set to known "cleared" state

  private:
    FacePointer m_ffp[3];
    char        m_ffi[3];

  public:
    // Three vertices for a triangular face
    // Stored as 3 pointers to the VertexType

          VertexType*&   V(const int j)       { assert(j >= 0 && j < 3); return v[j]; }
    const VertexType*    V(const int j) const { assert(j >= 0 && j < 3); return v[j]; }
    const VertexPointer cV(const int j) const { assert(j >= 0 && j < 3); return v[j]; }

    CoordType& P(const int j) { assert(j >= 0 && j < 3); return v[j]->P(); }
    const CoordType& P(const int j) const { assert(j >= 0 && j < 3); return v[j]->P(); }
    CoordType cP(const int j) const { assert(j >= 0 && j < 3); return v[j]->cP(); }

    VertexType*& V0(const int j) { return V(j); }
    VertexType*& V1(const int j) { return V((j + 1) % 3); }
    VertexType*& V2(const int j) { return V((j + 2) % 3); }

    const VertexType* V0(const int j) const { return V(j); }
    const VertexType* V1(const int j) const { return V((j + 1) % 3); }
    const VertexType* V2(const int j) const { return V((j + 2) % 3); }
    const VertexType* cV0(const int j) const { return cV(j); }
    const VertexType* cV1(const int j) const { return cV((j + 1) % 3); }
    const VertexType* cV2(const int j) const { return cV((j + 2) % 3); }

    CoordType& P0(const int j) { return V(j)->P(); }
    CoordType& P1(const int j) { return V((j + 1) % 3)->P(); }
    CoordType& P2(const int j) { return V((j + 2) % 3)->P(); }

    const CoordType& P0(const int j) const { return V(j)->P(); }
    const CoordType& P1(const int j) const { return V((j + 1) % 3)->P(); }
    const CoordType& P2(const int j) const { return V((j + 2) % 3)->P(); }

    const CoordType& cP0(const int j) const { return cV(j)->P(); }
    const CoordType& cP1(const int j) const { return cV((j + 1) % 3)->P(); }
    const CoordType& cP2(const int j) const { return cV((j + 2) % 3)->P(); }

  private:
    VertexType* v[3];

  public:
    // Per vertex Incremental Mark
    int&     IMark()       { return m_imark; }
    int      IMark() const { return m_imark; }
    int     cIMark() const { return m_imark; }
    void InitIMark()       { m_imark = 0; }

  private:
    int m_imark;

  public:

    // pointers to the three vertices defining this face
    int VN()  const { return 3; }
    int Prev(const int& i) const { return (i + (3 - 1)) % 3; }
    int Next(const int& i) const { return (i + 1) % 3; }

    // Manage bitflags (MyFace2)
    enum {                    // flags used:
      DELETED  = 0x00000001,  // Face is deleted from the mesh
      NOTREAD  = 0x00000002,  // Face of the mesh is not readable
      NOTWRITE = 0x00000004,  // Face of the mesh is not writable
      VISITED  = 0x00000010,  // Face has been visited

      BORDER0  = 0x00000040,  // Border flags, 
      BORDER1  = 0x00000080,  // assumes BORDERi = BORDER0<<i
      BORDER2  = 0x00000100,
      BORDER012 = BORDER0 | BORDER1 | BORDER2,

      NORMX = 0x00000200,     // Face Orientation Flags, used for 
      NORMY = 0x00000400,     // computing point-face distance
      NORMZ = 0x00000800,
      USER0 = 0x00200000      // First user bit
    };

    bool IsD() const  { return (this->cFlags() & DELETED) != 0; }   //  if Face is deleted
    bool IsR() const  { return (this->cFlags() & NOTREAD) == 0; }   //  if Face is readable
    bool IsW() const  { return (this->cFlags() & NOTWRITE) == 0; }  //  if Face is modifiable
    bool IsRW() const { return (this->cFlags() & (NOTREAD | NOTWRITE)) == 0; }  // if face is readable + modifiable
    bool IsV() const  { return (this->cFlags() & VISITED) != 0; }   //  if Face is Modified

    void SetFlags(int flagp) { this->Flags() = flagp; }
    void ClearFlags() { this->Flags() = 0; }

    void SetD()   { this->Flags() |= DELETED; }       // set deleted
    void SetW()   { this->Flags() &= (~NOTWRITE); }   // marks the Face as writable
    void ClearW() { this->Flags() |= NOTWRITE; }      // marks the Face as notwritable
    void SetV()   { this->Flags() |= VISITED; }       // set as visited the Face
    void ClearV() { this->Flags() &= ~VISITED; }      // set as unvisited the Face

    // Unused flags (MyFace2)
    // void ClearD() { this->Flags() &= (~DELETED); }    // undelete the Face
    // void SetR()   { this->Flags() &= (~NOTREAD); }    // marks Face as readable
    // void ClearR() { this->Flags() |= NOTREAD; }       // marks the Face as not readable


    bool IsB(int i) const { return (this->cFlags() & (BORDER0 << i)) != 0; }
    void SetB(int i)      {         this->Flags() |= (BORDER0 << i); }
    void ClearB(int i)    {         this->Flags() &= (~(BORDER0 << i)); }

    // Return first "unused" user-bit
    static int& FirstUnusedBitFlag() { static int b = USER0; return b; }

    // Allocate a bit that can be used by user. Updates FirstUnusedBitFlag.
    static inline int NewBitFlag() {
      int bitForTheUser = FirstUnusedBitFlag();
      FirstUnusedBitFlag() = FirstUnusedBitFlag() << 1;
      return bitForTheUser;
    }

    // De-allocate a pre allocated bit. Updates FirstUnusedBitFlag.
    // Note: deallocate bits in reverse order of allocation (as in a stack)
    static inline bool DeleteBitFlag(int bitval) {
      if (FirstUnusedBitFlag() >> 1 == bitval) {
        FirstUnusedBitFlag() = FirstUnusedBitFlag() >> 1;
        return true;
      }
      assert(0);        // should not get here
      return false;
    }

    bool IsUserBit(int userBit)    { return (this->Flags() & userBit) != 0; }
    void SetUserBit(int userBit)   {         this->Flags() |= userBit; }
    void ClearUserBit(int userBit) {         this->Flags() &= (~userBit); }

    void GetBBox(iso::Box3<ScalarType>& bb) const {
      if (this->IsD()) { bb.SetNull(); return; }
      bb.Set(cP(0)); bb.Add(cP(1)); bb.Add(cP(2));
    }

  };  // end class MyFace2


  /////////////////////////////////////////////////////////
  // 
  // start: isoSurf/iso_4_mesh.hpp
  // 
  /////////////////////////////////////////////////////////


  class MyMesh2;
  typedef typename MyMesh2  MeshType;

  typedef typename std::vector<VertexType>  VertContainer;
  typedef typename std::vector<FaceType>    FaceContainer;

  typedef typename VertContainer::iterator        VertexIterator;
  typedef typename VertContainer::const_iterator  ConstVertexIterator;

  typedef typename FaceContainer::iterator        FaceIterator;
  typedef typename FaceContainer::const_iterator  ConstFaceIterator;


  class MyMesh2
  {
  public:
    MyMesh2() : m_vn(0), m_fn(0), m_imark(0) {}
    ~MyMesh2() { Clear(); }

    VertContainer   vert;     // Container of vertices
    FaceContainer   face;     // Container of faces

    int VN() const { return m_vn; }
    int FN() const { return m_fn; }

    int& FaceNumber() { return m_fn; }
    int& VertexNumber() { return m_vn; }

    bool IsEmpty() const {
      return vert.empty() && face.empty();
    }

    void Clear() {
      vert.clear(); face.clear(); m_vn = 0; m_fn = 0; m_imark = 0;
    }

    int m_vn;       // Current number of non-[D]eleted vertices
    int m_fn;       // Current number of non-[D]eleted faces
    int m_imark;    // Current incremental mark

    iso::Box3<ScalarType> m_bbox;   // Bounding box of the mesh

  private:
    // disable copying of a mesh. Use Append()
    MyMesh2 operator = (const MyMesh2&) = delete;
    MyMesh2(const MyMesh2&) = delete;
  };

  size_t Index(const MeshType& m, const VertexType& v) {
    return ((&v) - (&*m.vert.begin()));
  }

  size_t Index(const MeshType& m, const VertexType* vp) {
    return (vp - (&*m.vert.begin()));
  }

  // Initialize the imark-system of the faces
  void InitFaceIMark(MeshType& m) {
    FaceIterator f;
    for (f = m.face.begin(); f != m.face.end(); ++f)
      if (!(*f).IsD() && (*f).IsR() && (*f).IsW())
        (*f).InitIMark();
  }

  // Initialize the imark-system of the vertices
  void InitVertexIMark(MeshType& m) {
    VertexIterator vi;
    for (vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if (!(*vi).IsD() && (*vi).IsRW())
        (*vi).InitIMark();
  }


  // Marking and unmarking vertices and faces:
  // Elements of a mesh are considered "marked" if their 
  // internal imark matches the mesh's internal counter.

  // Access the incremental mark. Generally, use IsMarked() and Mark()
  int& IMark(MeshType& m) { return m.m_imark; }
  // Does vertex mark match mesh mark?
  bool IsMarked(const MeshType& m, ConstVertexPointer  v) { return v->cIMark() == m.m_imark; }
  // Does face mark match mesh mark?.
  bool IsMarked(const MeshType& m, ConstFacePointer f) { return f->cIMark() == m.m_imark; }
  // Make element mark match mesh mark.
  void Mark(MeshType& m, VertexPointer v) { v->IMark() = m.m_imark; }
  void Mark(MeshType& m, FacePointer   f) { f->IMark() = m.m_imark; }
  // Incrementing mesh mark "un-marks" all elements:
  void UnMarkAll(MeshType& m) { ++m.m_imark; }

  // edge[j] of face f is a border if FFp(j) points to face j
  bool IsBorder(FaceType const& f, const int j) { return f.cFFp(j) == &f; }

  void RequireVertexCompactness(const MeshType& m) {
    if (m.vert.size() != size_t(m.m_vn))
      throw std::string("Vertex Vector Contains deleted elements");
  }
  void RequireFaceCompactness(const MeshType& m) {
    if (m.face.size() != size_t(m.m_fn))
      throw std::string("Face Vector Contains deleted elements");
  }
  void RequireCompactness(const MeshType& m) {
    RequireVertexCompactness(m);
    RequireFaceCompactness(m);
  }

  // compute or update the bounding box of a mesh..
  void UpdateBounding(MeshType& m) {
    m.m_bbox.SetNull();
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
      if (!(*vi).IsD()) {
        m.m_bbox.Add((*vi).cP());
      }
    }
  };

  // simple data structure for computing face-face adjacency.
  class PEdge {
  public:
    VertexPointer  v[2];  // the two Vertex pointer are ordered!
    FacePointer    f;     // the face where this edge belong
    int            z;     // index in [0..2] of the edge of the face
    bool isBorder;

    PEdge() { v[0] = v[1] = nullptr; f = nullptr; z = -1; isBorder = false; }
    PEdge(FacePointer  pf, const int nz) { this->Set(pf, nz); }

    void Set(FacePointer  pf, const int nz) {
      assert((pf != nullptr) && (nz >= 0) && (nz < 3));
      v[0] = pf->V(nz);
      v[1] = pf->V(pf->Next(nz));
      assert(v[0] != v[1]); // coincident verts ==> face f is degenerate 

      if (v[0] > v[1]) std::swap(v[0], v[1]);
      f = pf;
      z = nz;
    }

    inline bool operator < (const PEdge& pe) const {
      if (v[0] < pe.v[0]) return true;
      else if (v[0] > pe.v[0]) return false;
      else return v[1] < pe.v[1];
    }
    inline bool operator == (const PEdge& pe) const {
      return v[0] == pe.v[0] && v[1] == pe.v[1];
    }
  };


  // Fill a vector with all edges in the mesh. Each edge is stored 
  // the number of times it appears in mesh, with the referring face. 
  void FillEdgeVector(MeshType& m, std::vector<PEdge>& edgeVec) {
    edgeVec.reserve(3 * m.m_fn);
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int j = 0; j < 3; ++j) {
          edgeVec.push_back(PEdge(&*fi, j));
        }
      }
    }
  }

  // Update Face-Face topological relations to allow each face 
  // to find other faces sharing their edges.
  void FaceFace(MeshType& m) {
    if (m.m_fn == 0) return;

    std::vector<PEdge> e;
    FillEdgeVector(m, e);
    sort(e.begin(), e.end());               // sort the vertices

    int ne = 0;                             // number of real edges
    typename std::vector<PEdge>::iterator pe, ps;
    ps = e.begin(); pe = e.begin();
    do
    {
      if (pe == e.end() || !(*pe == *ps)) {
        typename std::vector<PEdge>::iterator q, q_next;
        for (q = ps; q < pe - 1; ++q) {     // process a block of equal edges
          assert(((*q).z >= 0) && ((*q).z < 3));
          q_next = q;
          ++q_next;
          assert(((*q_next).z >= 0) && ((*q_next).z < 3));
          (*q).f->FFp(q->z) = (*q_next).f;  // connect list of faces
          (*q).f->FFi(q->z) = (*q_next).z;
        }
        assert(((*q).z >= 0) && ((*q).z < 3));
        (*q).f->FFp((*q).z) = ps->f;
        (*q).f->FFi((*q).z) = ps->z;
        ps = pe;
        ++ne;                               // increment number of edges
      }
      if (pe == e.end()) {
        break;
      }
      ++pe;
    } while (true);
  }


  // Update Vertex-Face connectivity.
  //
  // For each vertex, enables retrieval of all faces sharing this vertex.
  // Note: isolated vertices have a null list of faces.
  void VertexFace(MeshType& m) {
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
      (*vi).VFp() = nullptr;  // note that (null,-1) means uninitiazlied 
      (*vi).VFi() = 0;        // while (null,0) is valid for isolated verts.
    }
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int j = 0; j < 3; ++j) {
          (*fi).VFp(j) = (*fi).V(j)->VFp();
          (*fi).VFi(j) = (*fi).V(j)->VFi();
          (*fi).V(j)->VFp() = &(*fi);
          (*fi).V(j)->VFi() = j;
        }
      }
    }
  }

  // Manage per-vertex and per-face flags (e.g. "border" flags).
  void VertexClear(MeshType& m, unsigned int FlagMask = 0xffffffff) {
    int andMask = ~FlagMask;
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
      if (!(*vi).IsD()) (*vi).Flags() &= andMask;
  }

  void FaceClear(MeshType& m, unsigned int FlagMask = 0xffffffff) {
    int andMask = ~FlagMask;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
      if (!(*fi).IsD()) (*fi).Flags() &= andMask;
  }

  void VertexClearV(MeshType& m) { VertexClear(m, VertexType::VISITED); }
  void FaceClearV(MeshType& m) { FaceClear(m, FaceType::VISITED); }
  void FaceClearB(MeshType& m) { FaceClear(m, FaceType::BORDER012); }


  class EdgeSorter
  {
  public:
    VertexPointer v[2];
    FacePointer f;
    int z;

    EdgeSorter() { v[0] = v[1] = nullptr; f = nullptr; z = -1; }
    void Set(const FacePointer pf, const int nz) {
      assert((pf != nullptr) && (nz >= 0) && (nz < 3));
      v[0] = pf->V(nz);
      v[1] = pf->V((nz + 1) % 3);
      assert(v[0] != v[1]);
      if (v[0] > v[1]) std::swap(v[0], v[1]);
      f = pf; z = nz;
    }
    inline bool operator <  (const EdgeSorter& pe) const {
      if (v[0] < pe.v[0]) return true;
      else if (v[0] > pe.v[0]) return false;
      else return v[1] < pe.v[1];
    }
    inline bool operator == (const EdgeSorter& pe) const {
      return v[0] == pe.v[0] && v[1] == pe.v[1];
    }
    inline bool operator != (const EdgeSorter& pe) const {
      return v[0] != pe.v[0] || v[1] != pe.v[1];
    }
  };

  // Computes per-face border flags without requiring topology info
  void FaceBorderFromNone(MeshType& m)
  {
    std::vector<EdgeSorter> e;
    typename FaceIterator pf;
    typename std::vector<EdgeSorter>::iterator p;

    for (VertexIterator v = m.vert.begin(); v != m.vert.end(); ++v) (*v).ClearB();
    if (m.m_fn == 0) return;

    FaceIterator fi; int n_edges = 0;
    for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD()) n_edges += 3;
    e.resize(n_edges);

    p = e.begin();                        // load face data
    for (pf = m.face.begin(); pf != m.face.end(); ++pf) {
      if (!(*pf).IsD()) {
        for (int j = 0; j < 3; ++j) {
          (*p).Set(&(*pf), j);
          (*pf).ClearB(j);
          ++p;
        }
      }
    }
    assert(p == e.end());
    sort(e.begin(), e.end());							// sort by vertex

    typename std::vector<EdgeSorter>::iterator pe, ps;
    ps = e.begin(); pe = e.begin();
    do {
      if (pe == e.end() || *pe != *ps) {  // find block of equal edges
        if (pe - ps == 1) {     // if (blocksize==1), this edge has no 
          ps->f->SetB(ps->z);   // neighboring face, so mark as "border"
        }
        ps = pe;
      }
      if (pe == e.end()) break;
      ++pe;
    } while (true);
  }

  class VFIterator
  {
  public:
    FaceType* f;  // Pointer to the face of the half-edge
    int z;        // Index of the vertex

    VFIterator() : f(0), z(-1) {}
    VFIterator(FaceType* _f, const int& _z) { f = _f; z = _z;  assert(z >= 0); }
    VFIterator(VertexType* _v) { f = _v->VFp(); z = _v->VFi(); assert(z >= 0); }

    FaceType*& F() { return f; }
    int& I() { return z; }

    VertexType* V()  const { return f->V(z); }
    VertexType* const& V0() const { return f->V0(z); }
    VertexType* const& V1() const { return f->V1(z); }
    VertexType* const& V2() const { return f->V2(z); }

    bool End() const { return f == nullptr; }
    void operator++() {
      FaceType* t = f;
      f = t->VFp(z);
      z = t->VFi(z);
    }
    void operator++(int) {
      ++(*this);
    }
  };

  void FaceBorderFromVF(MeshType& m) {
    FaceClearB(m);
    int visitedBit = VertexType::NewBitFlag();

    // edge calculation:
    // for each vertex look for adjacent vertices touched 
    // by only one face (or by an odd number of faces)
    const int BORDERFLAG[3] = { FaceType::BORDER0, FaceType::BORDER1, FaceType::BORDER2 };

    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
      if (!(*vi).IsD()) {
        for (VFIterator vfi(&*vi); !vfi.End(); ++vfi) {
          vfi.f->V1(vfi.z)->ClearUserBit(visitedBit);
          vfi.f->V2(vfi.z)->ClearUserBit(visitedBit);
        }
        for (VFIterator vfi(&*vi); !vfi.End(); ++vfi) {
          if (vfi.f->V1(vfi.z)->IsUserBit(visitedBit))  vfi.f->V1(vfi.z)->ClearUserBit(visitedBit);
          else vfi.f->V1(vfi.z)->SetUserBit(visitedBit);
          if (vfi.f->V2(vfi.z)->IsUserBit(visitedBit))  vfi.f->V2(vfi.z)->ClearUserBit(visitedBit);
          else vfi.f->V2(vfi.z)->SetUserBit(visitedBit);
        }
        for (VFIterator vfi(&*vi); !vfi.End(); ++vfi) {
          if (vfi.f->V(vfi.z) < vfi.f->V1(vfi.z) && vfi.f->V1(vfi.z)->IsUserBit(visitedBit))
            vfi.f->Flags() |= BORDERFLAG[vfi.z];
          if (vfi.f->V(vfi.z) < vfi.f->V2(vfi.z) && vfi.f->V2(vfi.z)->IsUserBit(visitedBit))
            vfi.f->Flags() |= BORDERFLAG[(vfi.z + 2) % 3];
        }
      }
    }
    VertexType::DeleteBitFlag(visitedBit);
  }

} // end namespace iso


// compiler complains if defined inside iso namespace
namespace std {
  template<>
  struct hash<iso::Point3i> {
    std::size_t operator()(const iso::Point3i& s) const {
      return std::hash<int>()(s[0]) ^ std::hash<int>()(s[1]) ^ std::hash<int>()(s[2]);
    }
  };
}


namespace iso {

  class AverageColorCell
  {
  public:

    AverageColorCell()
      : cnt(0), id(0), p(0, 0, 0), c(0), q(0) {}  // n(0, 0, 0), 

    ColorType   Col() const  { return (c / (float)cnt); }
    QualityType Cost() const { return (q / (float)cnt); }
    CoordType   Pos() const  { return (p / (float)cnt); }

    // Un-normalized face normals: want to drop small faces, and
    // weigh average with the size of faces falling in this cell.
    void AddFaceVertex(FaceType& f, int i) {
      p += f.cV(i)->cP();
      c += f.cV(i)->cC();
      q += f.cV(i)->cQ();
      cnt++;
    }

    int cnt, id;    // count verts added to this cell
    CoordType p;    // average position of verts in cell
    ColorType c;    // average of colors
    QualityType q;  // average of (quadric) "costs"
  };


  //-----------------------------------
  // Clustering routine
  //-----------------------------------
  template <class CellType> class Clustering
  {
  public:

    class SimpleTri {
    public:
      CellType* v[3];  // three cells define vertices for a face.

      int ii(int i) const { return *((int*)(&(v[i]))); }
      bool operator < (const SimpleTri& p) const {
        return	(v[2] != p.v[2]) ? (v[2] < p.v[2]) :
                (v[1] != p.v[1]) ? (v[1] < p.v[1]) :
                (v[0] < p.v[0]);
      }

      void sort() {
        if (v[0] > v[1]) std::swap(v[0], v[1]); // now v0 < v1
        if (v[0] > v[2]) std::swap(v[0], v[2]); // now v0 is the minimum
        if (v[1] > v[2]) std::swap(v[1], v[2]); // sorted!
      }

      bool operator == (const SimpleTri& pt) const {
        return ((pt.v[0] == v[0]) && (pt.v[1] == v[1]) && (pt.v[2] == v[2]));
      }

      // Hashing Function;
      size_t operator()(const SimpleTri& pt) const {
        return std::hash<CellType*>()(pt.v[0]) ^ std::hash<CellType*>()(pt.v[1]) ^ std::hash<CellType*>()(pt.v[2]);
      }
    };


    // Init() parameters:
    // 
    // _numcells is an approximate total number of cells composing the 
    //        grid surrounding the objects (usually a large number)
    //        e.g. [_numcells==1,000,000] means a 100x100x100 grid
    // 
    // _cellsize is the length of each edge of each grid cell.
    //        e.g. [_cellsize==2.0] means that all vertexes in a 
    //             2x2x2 cell will be clustered togheter
    //
    // Note:
    // _numcells is used only if _cellsize is zero.
    // _cellsize implies a measure of the maximum error introduced
    //           during the simplification (half cell edge length)

    void Init(iso::Box3<ScalarType> _mbb, int _numcells, ScalarType _cellsize = 0) {
      GridCell.clear();
      TriSet.clear();
      Grid.bbox = _mbb;
      // inflate the bb calculated
      ScalarType infl = (_cellsize == (ScalarType)0) ? (Grid.bbox.Diag() / _numcells) : (_cellsize);
      Grid.bbox.min -= CoordType(infl, infl, infl);
      Grid.bbox.max += CoordType(infl, infl, infl);
      Grid.dim = Grid.bbox.max - Grid.bbox.min;
      if (0 == _cellsize)
        iso::BestDim(_numcells, Grid.dim, Grid.siz);
      else
        Grid.siz = iso::Point3i::Construct(Grid.dim / _cellsize);

      // find voxel size
      Grid.voxel[0] = Grid.dim[0] / Grid.siz[0];
      Grid.voxel[1] = Grid.dim[1] / Grid.siz[1];
      Grid.voxel[2] = Grid.dim[2] / Grid.siz[2];
    }

    iso::BasicGrid<ScalarType> Grid;

    std::unordered_set<SimpleTri, SimpleTri> TriSet;
    typedef typename std::unordered_set<SimpleTri, SimpleTri>::iterator TriHashSetIterator;
    std::unordered_map<iso::Point3i, CellType> GridCell;

    void AddMesh(MeshType& m) {
      FaceIterator fi;
      for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD()) {
        iso::Point3i pi; SimpleTri st;
        for (int i = 0; i < 3; ++i) {
          Grid.PToIP((*fi).cV(i)->cP(), pi); // map point coords to voxel index
          st.v[i] = &(GridCell[pi]);
          st.v[i]->AddFaceVertex((*fi), i);
        }
        if ((st.v[0] != st.v[1]) && (st.v[0] != st.v[2]) && (st.v[1] != st.v[2])) {
          // (three distinct voxels) ==> add non-degenerate SimpleTri
          st.sort();
          TriSet.insert(st);
        }
      }
    }


    void AddMesh_Clip(MeshType& m, const dfloat(&cb)[2][3]) {

      FaceIterator fi;
      for (fi = m.face.begin(); fi != m.face.end(); ++fi) if (!(*fi).IsD()) {

        bool doAdd = true;
        for (int i = 0; i < 3; ++i) {

          // if ANY vertex in this face is outside the clip-box,
          // mark the face as "[D]eleted" and skip to next face:

          const CoordType& p = (*fi).cV(i)->cP();
          if ((cb[0][0] > p.X()) || (cb[1][0] < p.X())) { fi->SetD(); doAdd = false; break; }
          if ((cb[0][1] > p.Y()) || (cb[1][1] < p.Y())) { fi->SetD(); doAdd = false; break; }
          if ((cb[0][2] > p.Z()) || (cb[1][2] < p.Z())) { fi->SetD(); doAdd = false; break; }
        }

        if (doAdd) {
          // all vertices in this face are inside clip region
          iso::Point3i pi; SimpleTri st;
          for (int i = 0; i < 3; ++i) {
            Grid.PToIP((*fi).cV(i)->cP(), pi);
            st.v[i] = &(GridCell[pi]);
            st.v[i]->AddFaceVertex((*fi), i);
          }
          if ((st.v[0] != st.v[1]) && (st.v[0] != st.v[2]) && (st.v[1] != st.v[2])) {
            st.sort();
            TriSet.insert(st);
          }
        }
      }
    }


    void ExtractMesh(MeshType& m) {
      m.Clear();
      if (GridCell.empty())  return;

      m.vert.resize(GridCell.size());   // new size of vertex array
      m.m_vn = GridCell.size();         // number of "non-Deleted" verts

      typename std::unordered_map<iso::Point3i, CellType>::iterator gi;
      int i = 0;
      for (gi = GridCell.begin(); gi != GridCell.end(); ++gi) {
        m.vert[i].P() = (*gi).second.Pos();
        m.vert[i].C() = (*gi).second.Col();   // "color"
        m.vert[i].Q() = (*gi).second.Cost();  // "cost" (for quadric edge collapse)
        (*gi).second.id = i;
        ++i;
      }

      m.face.resize(TriSet.size());   // new size of face array
      m.m_fn = TriSet.size();         // number of "non-[D]eleted" faces

      TriHashSetIterator ti;
      i = 0;
      for (ti = TriSet.begin(); ti != TriSet.end(); ++ti) {
        m.face[i].V(0) = &(m.vert[(*ti).v[0]->id]);
        m.face[i].V(1) = &(m.vert[(*ti).v[1]->id]);
        m.face[i].V(2) = &(m.vert[(*ti).v[2]->id]);
        i++;
      }

    }
  }; // end class Clustering


  //-------------------------------------------------------
  // Utility cleaning routines
  //-------------------------------------------------------

  // Helper class to update pointers, used when allocation 
  // operations invalidate the pointers to mesh elements.
  // Typically used after allocation of new vertexes, edges
  // and faces, or when compacting containers to get rid of 
  // [D]eleted elements.
  // 
  // Notes: 
  // THis object allows the updating of a pointer immediately 
  // after it becomes invalid. It can also be used to PREVENT 
  // update of the internal pointers (e.g. when building all
  // the internal connections when importing a mesh).

  template<class SimplexPointerType>
  class PointerUpdater {
  public:
    PointerUpdater() 
      : newBase(0), oldBase(0), newEnd(0), oldEnd(0), preventUpdateFlag(false) {}

    void Clear() { newBase = oldBase = newEnd = oldEnd = 0; remap.clear(); }

    // Update a pointer to an element of a mesh after a reallocation
    // The updating is done correctly ONLY if this PointerUpdater has 
    // been passed to the corresponding allocation call.
    void Update(SimplexPointerType& vp) {
      if (vp<oldBase || vp>oldEnd) return;
      assert(vp >= oldBase);
      assert(vp < oldEnd);
      vp = newBase + (vp - oldBase);
      if (!remap.empty()) {
        vp = newBase + remap[vp - newBase];
      }
    }

    // Return true if the allocation operation that initialized 
    // this PointerUpdater has caused a reallocation
    bool NeedUpdate() {
      if ((oldBase && (newBase != oldBase) && (!preventUpdateFlag)) || (!remap.empty())) {
        return true;
      }
      return false;
    }

    SimplexPointerType newBase, oldBase, newEnd, oldEnd;

    // This vector keeps the new position of an element. 
    // Uninitialized elements have max_int value to denote 
    // an element that has not to be remapped.
    std::vector<size_t> remap;

    // when true no update is considered necessary.
    bool preventUpdateFlag;
  };


  // Mark vertex as "[D]eleted", update mesh's vertex counter
  void DeleteVertex(MeshType& m, VertexType& v) {
    assert(&v >= &m.vert.front() && &v <= &m.vert.back());
    assert(!v.IsD());
    v.SetD();
    --m.m_vn;
  }

  // Remove vertices that are not referenced by any face or edge.
  // If (!DeleteVertexFlag), just count the unreferenced verts.
  int RemoveUnreferencedVertex(MeshType& m, bool DeleteVertexFlag = true) {
    std::vector<bool> referredVec(m.vert.size(), false);
    int deleted = 0;
    for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
      if (!(*fi).IsD())
        for (auto j = 0; j < 3; ++j)
          referredVec[Index(m, (*fi).V(j))] = true;

    if (!DeleteVertexFlag) {
      return std::count(referredVec.begin(), referredVec.end(), false);
    }
    for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
      if ((!(*vi).IsD()) && (!referredVec[Index(m, *vi)])) {
        DeleteVertex(m, *vi);
        ++deleted;
      }
    }
    return deleted;
  }

  class SortedTriple {
  public:
    SortedTriple() { v[0] = v[1] = v[2] = 0; fp = nullptr; }
    SortedTriple(uint v0, uint v1, uint v2, FacePointer _fp) {
      v[0] = v0; v[1] = v1; v[2] = v2; fp = _fp;
      std::sort(v, v + 3);
    }
    bool operator < (const SortedTriple& p) const {
      return (v[2] != p.v[2]) ? (v[2] < p.v[2]) :
             (v[1] != p.v[1]) ? (v[1] < p.v[1]) :
             (v[0] < p.v[0]);
    }
    bool operator == (const SortedTriple& s) const {
      if ((v[0] == s.v[0]) && (v[1] == s.v[1]) && (v[2] == s.v[2])) return true;
      return false;
    }
    uint v[3];
    FacePointer fp;
  };

  void DeleteFace(MeshType& m, FaceType& f) {
    assert(&f >= &m.face.front() && &f <= &m.face.back());
    assert(!f.IsD());
    f.SetD();
    --m.m_fn;
  }

  // Removes duplicate faces in the mesh (same vertex refs).
  int RemoveDuplicateFace(MeshType& m) {
    std::vector<SortedTriple> fvec;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        fvec.push_back(SortedTriple(Index(m, (*fi).V(0)), Index(m, (*fi).V(1)), Index(m, (*fi).V(2)), &*fi));
      }
    }
    std::sort(fvec.begin(), fvec.end());
    int total = 0;
    for (int i = 0; i<int(fvec.size()) - 1; ++i) {
      if (fvec[i] == fvec[i + 1]) {
        DeleteFace(m, *(fvec[i].fp));
        total++;
      }
    }
    return total;
  }


  void PermutateVertexVector(MeshType& m, PointerUpdater<VertexPointer>& pu)
  {
    // Permute the vertex vector. After calling this function,
    //  m.vert[ newVertIndex[i] ] = m.vert[i];
    // i.e. newVertIndex[i] is the new index of the vertex i

    if (m.vert.empty()) return;
    for (size_t i = 0; i < m.vert.size(); ++i) {
      if (pu.remap[i] < size_t(m.m_vn)) {
        assert(!m.vert[i].IsD());
        m.vert[pu.remap[i]].CopyVertex(m.vert[i]);
      }
    }

    // setup the pointer updater
    pu.oldBase = &m.vert[0];
    pu.oldEnd = &m.vert.back() + 1;

    // resize to fit the number of "non-Deleted" vertices
    m.vert.resize(m.m_vn);

    // setup the pointer updater
    pu.newBase = (m.vert.empty()) ? 0 : &m.vert[0];
    pu.newEnd = (m.vert.empty()) ? 0 : &m.vert.back() + 1;

    // For each face, update its Vertex pointers
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int i = 0; i < 3; ++i) {
          size_t oldIndex = (*fi).V(i) - pu.oldBase;
          assert(pu.oldBase <= (*fi).V(i) && oldIndex < pu.remap.size());
          (*fi).V(i) = pu.newBase + pu.remap[oldIndex];
        }
      }
    }
  }

  // Compact vector of vertices, removing deleted elements.
  // Deleted elements are put to the end of the vector and the 
  // vector is resized. Order between elements is preserved 
  // but not their position (hence the PointerUpdater)

  void CompactVertexVector(MeshType& m, PointerUpdater<VertexPointer>& pu) {
    // If already compacted, return
    if (m.m_vn == (int)m.vert.size()) return;

    // newVertIndex [<old_vert_position>] gives the 
    // new position of the vertex in the vector;
    pu.remap.resize(m.vert.size(), std::numeric_limits<size_t>::max());

    size_t pos = 0;
    size_t i = 0;
    for (i = 0; i < m.vert.size(); ++i) {
      if (!m.vert[i].IsD()) {
        pu.remap[i] = pos;
        ++pos;
      }
    }

    assert((int)pos == m.m_vn);
    PermutateVertexVector(m, pu);
  }

  void CompactVertexVector(MeshType& m) {
    PointerUpdater<VertexPointer>  pu;
    CompactVertexVector(m, pu);
  }

  int ClusterVertex(MeshType& m, const float rad) {
    if (m.m_vn == 0) return 0;
    int mergedCnt = 0;

    CompactVertexVector(m);

    iso::SpatialHashTable<VertexType, ScalarType> sht;
    sht.Set(m.vert.begin(), m.vert.end());

    VertexClearV(m);
    std::vector<VertexType*> closests;
    for (VertexIterator viv = m.vert.begin(); viv != m.vert.end(); ++viv) {
      if (!(*viv).IsD() && !(*viv).IsV()) {
        (*viv).SetV();

        CoordType p = viv->cP();
        iso::Box3<ScalarType> bb(p - CoordType(rad, rad, rad),
          p + CoordType(rad, rad, rad));

        // find all points within bbox of point p
        GridGetInBox(sht, bb, closests);

        for (size_t i = 0; i < closests.size(); ++i) {
          ScalarType dist = Distance(p, closests[i]->cP());
          if (dist < rad && !closests[i]->IsV()) {
            mergedCnt++;
            closests[i]->SetV();
            closests[i]->P() = p;
          }
        }
      }
    }
    return mergedCnt;
  }

  int RemoveDegenerateFace(MeshType& m) {
    // Delete faces with less than three distinct vertices
    int count_fd = 0;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        if ((*fi).V(0) == (*fi).V(1) || (*fi).V(0) == (*fi).V(2) || (*fi).V(1) == (*fi).V(2)) {
          DeleteFace(m, *fi);
          count_fd++;
        }
      }
    }
    return count_fd;
  }

  // comparison function for removing duplicate vertices
  struct RemoveDupVert_Compare {
    inline bool operator()(VertexPointer const& a, VertexPointer const& b) {
      return ((*a).cP() == (*b).cP()) ? (a < b) : ((*a).cP() < (*b).cP());
    }
  };

  // Remove all duplicate vertices in the mesh (same spatial positions).
  // Call this function BEFORE building any topology information.
  int RemoveDuplicateVertex(MeshType& m, bool RemoveDegenerateFlag = true) {
    if (m.vert.size() == 0 || m.m_vn == 0) return 0;

    std::map<VertexPointer, VertexPointer> mp;
    size_t i, j;
    VertexIterator vi;
    int deleted = 0;
    int k = 0;
    size_t num_vert = m.vert.size();
    std::vector<VertexPointer> perm(num_vert);
    for (vi = m.vert.begin(); vi != m.vert.end(); ++vi, ++k)
      perm[k] = &(*vi);

    RemoveDupVert_Compare c_obj;
    std::sort(perm.begin(), perm.end(), c_obj);

    j = 0;
    i = j;
    mp[perm[i]] = perm[j];
    ++i;
    for (; i != num_vert; ) {
      if ((!(*perm[i]).IsD()) && (!(*perm[j]).IsD()) && (*perm[i]).P() == (*perm[j]).cP()) {
        VertexPointer t = perm[i];
        mp[perm[i]] = perm[j];
        ++i;
        DeleteVertex(m, *t);
        deleted++;
      }
      else {
        j = i;
        ++i;
      }
    }

    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (k = 0; k < 3; ++k) {
          if (mp.find((VertexPointer)(*fi).V(k)) != mp.end()) {
            (*fi).V(k) = &*mp[(*fi).V(k)];
          }
        }
      }
    }

    if (RemoveDegenerateFlag) {
      RemoveDegenerateFace(m);
    }
    return deleted;
  }

  // forward
  void CompactFaceVector(MeshType& m);

  // Merge all vertices that are closer than the given radius.
  int MergeCloseVertex(MeshType& m, const float radius) {
    int mergedCnt = 0;
    mergedCnt = ClusterVertex(m, radius);   // merged verts will be duplicates
    RemoveDuplicateVertex(m, true);         // mark duplicates as Deleted

#if (1)
    // TODO: Do we NEED to clean up everything?
    RemoveUnreferencedVertex(m);
    RemoveDuplicateFace(m);
    CompactFaceVector(m);
    CompactVertexVector(m);
#endif

    return mergedCnt;
  }

  // Compact face vector by removing deleted elements.
  void CompactFaceVector(MeshType& m, PointerUpdater<FacePointer>& pu) {
    // If already compacted, return
    if (m.m_fn == (int)m.face.size()) return;

    // newFaceIndex [<old_face_pos>] gives the new position 
    // of the face in the vector;
    pu.remap.resize(m.face.size(), std::numeric_limits<size_t>::max());

    size_t pos = 0;
    for (size_t i = 0; i < m.face.size(); ++i) {
      if (!m.face[i].IsD()) {
        if (pos != i) {
          m.face[pos].CopyFace(m.face[i]);
        }
        pu.remap[i] = pos;
        ++pos;
      }
    }
    assert((int)pos == m.m_fn);

    FacePointer fbase = &m.face[0];
    // Loop on the vertices to correct VF relation
    for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
      if (!(*vi).IsD()) {
        if ((*vi).IsVFInitialized() && (*vi).VFp() != nullptr) {
          size_t oldIndex = (*vi).cVFp() - fbase;
          assert(fbase <= (*vi).cVFp() && oldIndex < pu.remap.size());
          (*vi).VFp() = fbase + pu.remap[oldIndex];
        }
      }
    }

    // Loop on the faces to correct VF and FF relations
    pu.oldBase = &m.face[0];
    pu.oldEnd = &m.face.back() + 1;

    m.face.resize(m.m_fn);
    pu.newBase = (m.face.empty()) ? 0 : &m.face[0];
    pu.newEnd = (m.face.empty()) ? 0 : &m.face.back() + 1;

    // Finally, update the various face pointers (if initialized)
    // which define the VF and FF connectivity relations
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!(*fi).IsD()) {
        for (int i=0; i<3; ++i) {
          if ((*fi).IsVFInitialized(i) && (*fi).VFp(i) != nullptr) {
            size_t oldIndex = (*fi).VFp(i) - fbase;
            assert(fbase <= (*fi).VFp(i) && oldIndex < pu.remap.size());
            (*fi).VFp(i) = fbase + pu.remap[oldIndex];
          }
        }

        for (int i=0; i<3; ++i) {
          if ((*fi).cFFp(i) != nullptr) {
            size_t oldIndex = (*fi).FFp(i) - fbase;
            assert(fbase <= (*fi).FFp(i) && oldIndex < pu.remap.size());
            (*fi).FFp(i) = fbase + pu.remap[oldIndex];
          }
        }
      }
    }
  }

  void CompactFaceVector(MeshType& m) {
    PointerUpdater<FacePointer>  pu;
    CompactFaceVector(m, pu);
  }

  //-------------------------------------------------------
  // manage "Connected Components"
  //-------------------------------------------------------
  class ConnectedComponentIterator
  {
  public:
    void operator ++() {
      FacePointer fpt = sf.top();
      sf.pop();
      for (int j = 0; j < fpt->VN(); ++j) {
        if (!IsBorder(*fpt, j)) {
          FacePointer l = fpt->FFp(j);
          if (!IsMarked(*mp, l)) {
            Mark(*mp, l);
            sf.push(l);
          }
        }
      }
    }

    void start(MeshType& m, FacePointer p) {
      mp = &m;
      while (!sf.empty()) sf.pop();
      UnMarkAll(m);
      Mark(m, p);
      sf.push(p);
    }

    bool completed() { return sf.empty(); }
    FacePointer operator *() { return sf.top(); }

  private:
    std::stack<FacePointer> sf;
    MeshType* mp;
  };

  // Compute the set of connected components in a mesh.
  // For each connected component, store a pair <int,faceptr> 
  // with its size and a pointer into the connected component.

  typedef typename std::vector<std::pair<int, FacePointer> >  VecPiF;

  int ConnectedComponents(MeshType& m, VecPiF& CCV) {
    CCV.clear();
    FaceClearV(m);
    std::stack<FacePointer> sf;
    FacePointer fpt = &*(m.face.begin());
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!((*fi).IsD()) && !(*fi).IsV())
      {
        (*fi).SetV();
        CCV.push_back(std::make_pair(0, &*fi));
        sf.push(&*fi);
        while (!sf.empty()) {
          fpt = sf.top();
          ++CCV.back().first;
          sf.pop();
          for (int j = 0; j < fpt->VN(); ++j) {
            if (!IsBorder(*fpt, j)) {
              FacePointer l = fpt->FFp(j);
              if (!(*l).IsV()) {
                (*l).SetV();
                sf.push(l);
              }
            }
          }
        }
      }
    }
    return int(CCV.size());
  }

  int CountConnectedComponents(MeshType& m) {
    VecPiF CCV;
    return ConnectedComponents(m, CCV);
  }


  // Delete small / isolated clusters of faces
  std::pair<int, int> RemoveSmallConnectedComponents(MeshType& m, int maxCCSize) {
    VecPiF CCV;
    int TotalCC = ConnectedComponents(m, CCV);
    int DeletedCC = 0;

    ConnectedComponentIterator ci;
    for (unsigned int i = 0; i < CCV.size(); ++i) {
      std::vector<FacePointer> FPV;
      if (CCV[i].first < maxCCSize) {
        DeletedCC++;
        for (ci.start(m, CCV[i].second); !ci.completed(); ++ci)
          FPV.push_back(*ci);

        typename std::vector<FacePointer>::iterator fpvi;
        for (fpvi = FPV.begin(); fpvi != FPV.end(); ++fpvi) {
          DeleteFace(m, (**fpvi));
        }
      }
    }
    return std::make_pair(TotalCC, DeletedCC);
  }


  //-------------------------------------------------------
  // manage non-manifold faces
  //-------------------------------------------------------
  bool IsManifold(FaceType const& f, const int j) {
    assert(f.cFFp(j) != nullptr);           // assumes Face-Face info is up-to-date
    return ((f.cFFp(j) == &f) ||            // edge is on border OR
      (&f == f.cFFp(j)->cFFp(f.cFFi(j))));  // neighbor of neighbor is self
  }

  class Pos
  {
    // Stores a "position" over a face in a mesh containing:
    //  - a pointer to the current face,
    //  - the index of one edge and a pointer to 
    //    one of the vertices of the edge.

  public:
    typedef Pos  PosType;

    FaceType* f;      // Pointer to the face of the half-edge
    int z;            // Index of the edge
    VertexType* v;    // Pointer to the vertex

    Pos() : f(0), z(-1), v(0) {}

    Pos(FaceType* fp, int zp, VertexType* vp) {
      f = fp; z = zp; v = vp;
      assert((vp == fp->V0(zp)) || (vp == fp->V1(zp)));
    }

    Pos(FaceType* fp, int zp) { f = fp; z = zp; v = f->V(zp); }

    Pos(FaceType* fp, VertexType* vp) {
      f = fp;
      v = vp;
      for (int i = 0; i < f->VN(); ++i)
        if (f->V(i) == v) { z = f->Prev(i); break; }
    }

    VertexType*& V() { return v; }
    VertexType* V() const { return v; }
    int& E() { return z; }
    int          E() const { return z; }
    FaceType*& F() { return f; }
    FaceType* F() const { return f; }

    // Returns the face index of the vertex inside the face.
    // Note that this is DIFFERENT from using the z member 
    // that denotes the edge index inside the face.
    // 
    // Expect: (Vind != (z+1)%3) && (Vind == z || Vind = (z+2)%3)

    int VInd() const {
      for (int i = 0; i < 3; ++i) if (v == f->V(i)) return i;
      assert(0);
      return -1;
    }
    inline bool operator == (PosType const& p) const {
      return (f == p.f && z == p.z && v == p.v);
    }
    inline bool operator != (PosType const& p) const {
      return (f != p.f || z != p.z || v != p.v);
    }

    // compare: (1) face pointers, (2) index of edge, (3) vertex pointers
    inline bool operator <= (PosType const& p) const {
      return	(f != p.f) ? (f < p.f) : (z != p.z) ? (z < p.z) : (v <= p.v);
    }
    inline bool operator < (PosType const& p) const {
      if ((*this) == p) return false;
      return ((*this) <= p);
    }

    void SetNull() { f = 0; v = 0; z = -1; }
    bool IsNull() const { return f == 0 || v == 0 || z < 0; }

    void NextF() {        // Change face along z 
      FaceType* t = f;    // Note: handles non-manifold edge/face.
      f = t->FFp(z);
      z = t->FFi(z);
    }

    bool pos_IsBorder() const { return IsBorder(*f, z); }
    bool pos_IsManifold() { return IsManifold(*f, z); }

    void Set(FaceType* const fp, int const zp, VertexType* const vp) {
      f = fp; z = zp; v = vp;
      assert(f->V(f->Prev(z)) != v && (f->V(f->Next(z)) == v || f->V(z) == v));
    }

    void Set(FaceType* const pFace, VertexType* const pVertex) {
      f = pFace; v = pVertex;
      for (int i = 0; i < f->VN(); ++i) if (f->V(i) == v) { z = f->Prev(i); break; }
    }
  };

  // Count number of faces incident upon a complex edge
  int ComplexSize(FaceType& f, const int e) {
    if (IsBorder(f, e))  return 1;
    if (IsManifold(f, e)) return 2;

    // Non 2-manifold case
    Pos fpos(&f, e);
    int cnt = 0;
    do {
      fpos.NextF();
      assert(!fpos.pos_IsBorder());
      assert(!fpos.pos_IsManifold());
      ++cnt;
    } while (fpos.f != &f);
    assert(cnt > 2);
    return cnt;
  }

  // Detach face f from adjacent face connected via edge e.
  // Designed to be used when faces are not two-manifold.
  // Note: assumes Face-Face adjacency info is up-to-date.

  void FFDetach(FaceType& f, const int e) {
    assert(!IsBorder(f, e));                // don't detach border edges
    int complexity = ComplexSize(f, e);
    assert(complexity > 0);
    Pos FirstFace(&f, e);
    Pos LastFace(&f, e);
    FirstFace.NextF();
    LastFace.NextF();
    int cnt = 0;

    // In the case of non 2-manifold faces, we continue to 
    // advance LastFace until it becomes the one preceeding 
    // face to be erased (TODO: consider deleting "fins")

    while (LastFace.f->FFp(LastFace.z) != &f) {
      assert(ComplexSize(*LastFace.f, LastFace.z) == complexity);
      assert(!LastFace.pos_IsManifold());   // enter here only if edge is non-manifold 
      assert(!LastFace.pos_IsBorder());
      LastFace.NextF();
      cnt++;
      assert(cnt < 100);
    }

    assert(LastFace.f->FFp(LastFace.z) == &f);
    assert(f.FFp(e) == FirstFace.f);

    // Link last face to first, skipping the face to be detached;
    LastFace.f->FFp(LastFace.z) = FirstFace.f;
    LastFace.f->FFi(LastFace.z) = FirstFace.z;
    assert(ComplexSize(*LastFace.f, LastFace.z) == complexity - 1);

    // Finally, self-connect target edge (make it a "border").
    f.FFp(e) = &f;
    f.FFi(e) = e;
    assert(ComplexSize(f, e) == 1);
  }

  ScalarType DoubleArea(const FaceType& t) {
    return Norm((t.cP(1) - t.cP(0)) ^ (t.cP(2) - t.cP(0)));
  }

  // Helper for sorting non-manifold faces by area. 
  // Used in RemoveNonManifoldFace()
  struct CompareAreaFP {
    bool operator ()(FacePointer const& f1, FacePointer const& f2) const {
      // TODO: sort "fins" first?
      return DoubleArea(*f1) < DoubleArea(*f2);
    }
  };

  // Removal of faces that were incident on a non manifold edge.
  int RemoveNonManifoldFace(MeshType& m) {
    // assumes caller has called FaceFace(m);
    FaceIterator fi;
    int count_fd = 0;
    std::vector<FacePointer> ToDelVec;

    for (fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!fi->IsD()) {
        if ((!IsManifold(*fi, 0)) ||
          (!IsManifold(*fi, 1)) ||
          (!IsManifold(*fi, 2)))
        {
          ToDelVec.push_back(&*fi);
        }
      }
    }

    std::sort(ToDelVec.begin(), ToDelVec.end(), CompareAreaFP());

    for (size_t i = 0; i < ToDelVec.size(); ++i) {
      if (!ToDelVec[i]->IsD()) {
        FaceType& ff = *ToDelVec[i];
        if ((!IsManifold(ff, 0)) ||
          (!IsManifold(ff, 1)) ||
          (!IsManifold(ff, 2)))
        {
          for (int j = 0; j < 3; ++j)
            if (!IsBorder(ff, j))
              FFDetach(ff, j);

          DeleteFace(m, ff);
          count_fd++;
        }
      }
    }

    return count_fd;
  }


  // Count number of non-manifold EDGES in a mesh 
  // (i.e. edges with > 2 incident faces)
  // 
  // Note: does not check for non-manifold VERTICES.

  int CountNonManifoldEdgeFF(MeshType& m) {
    int nmfBit[3];
    nmfBit[0] = FaceType::NewBitFlag();
    nmfBit[1] = FaceType::NewBitFlag();
    nmfBit[2] = FaceType::NewBitFlag();

    FaceClear(m, nmfBit[0] + nmfBit[1] + nmfBit[2]);

    int edgeCnt = 0;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!fi->IsD()) {
        for (int i = 0; i < 3; ++i) {
          if (!IsManifold(*fi, i)) {
            if (!(*fi).IsUserBit(nmfBit[i])) {
              ++edgeCnt;
              // follow the ring of faces incident on edge i;
              Pos nmf(&*fi, i);
              do {
                nmf.F()->SetUserBit(nmfBit[nmf.E()]);
                nmf.NextF();
              } while (nmf.f != &*fi);
            }
          }
        }
      }
    }
    return edgeCnt;
  }


  // Change orientation of a face by inverting index of two verts.
  void SwapEdge(FaceType& f, const int z) {

    std::swap(f.V0(z), f.V1(z));    // swap V0(z) with V1(z)

    // store info to preserve topology
    int z1 = (z + 1) % 3;
    int z2 = (z + 2) % 3;
    FaceType* g1p = f.FFp(z1);
    FaceType* g2p = f.FFp(z2);
    int g1i = f.FFi(z1);
    int g2i = f.FFi(z2);

    // g0 face topology is not affected by the swap
    if (g1p != &f) {
      g1p->FFi(g1i) = z2; f.FFi(z2) = g1i;
    }
    else {
      f.FFi(z2) = z2;
    }

    if (g2p != &f) {
      g2p->FFi(g2i) = z1;  f.FFi(z1) = g2i;
    }
    else {
      f.FFi(z1) = z1;
    }

    // finalize swap
    f.FFp(z1) = g2p;
    f.FFp(z2) = g1p;
  }


  // Check if given face is oriented consistently 
  // with face adjacent to the specified edge.
  bool CheckOrientation(FaceType& f, int z) {
    if (IsBorder(f, z)) {
      return true;
    }
    else {
      FaceType* g = f.FFp(z);
      int gi = f.FFi(z);
      if (f.V0(z) == g->V1(gi))
        return true;
      else
        return false;
    }
  }

  void OrientCoherentlyMesh(MeshType& m, bool& _IsOriented, bool& _IsOrientable)
  {
    //######################################
    // See OrientConnectedComponents() below
    //######################################

    // assumes caller has called FaceFace(m);
    bool IsOrientable = true, IsOriented = true;

    FaceClearV(m);
    std::stack<FacePointer> faces;
    for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi)
    {
      if (!fi->IsD() && !fi->IsV()) {
        // each face put in the stack is selected (and oriented)
        fi->SetV();
        faces.push(&(*fi));
        while (!faces.empty())
        {
          FacePointer fp = faces.top();
          faces.pop();

          // make consistently oriented the adjacent faces
          for (int j = 0; j < 3; j++) {
            if (!IsBorder(*fp, j) && IsManifold(*fp, j)) {
              FacePointer fp_j = fp->FFp(j);
              int iaux = fp->FFi(j);
              if (!CheckOrientation(*fp_j, iaux)) {
                IsOriented = false;
                if (!fp_j->IsV()) {
                  SwapEdge(*fp_j, iaux);
                }
                else {
                  IsOrientable = false;
                  break;
                }
              }
              if (!fp_j->IsV()) {
                fp_j->SetV();
                faces.push(fp_j);
              }
            }
          }
        }
      }
      if (!IsOrientable)	break;
    }
    _IsOriented = IsOriented;
    _IsOrientable = IsOrientable;
  }


  void OrientConnectedComponents(MeshType& m)
  {
    // Note: for tri-soup isosurfaces, instead of trying 
    // to orient all faces in mesh (m.face), here we try
    // to orient the faces in each connected component.

    // assumes caller has called FaceFace(m);
    bool IsOrientable = true, IsOriented = true;

    VecPiF CCV;
    int TotalCC = ConnectedComponents(m, CCV);

    FaceClearV(m);

    ConnectedComponentIterator ci;
    for (uint i = 0; i < CCV.size(); ++i)
    {
      // for each connected component, we try to orient 
      // coherently its subset of the faces in the mesh

      // load a vector of pointers to faces in this component
      std::vector<FacePointer> FPV;
      for (ci.start(m, CCV[i].second); !ci.completed(); ++ci) {
        FPV.push_back(*ci);
      }

      std::stack<FacePointer> faces;
      typename std::vector<FacePointer>::iterator fpvi;
      for (fpvi = FPV.begin(); fpvi != FPV.end(); ++fpvi)
      {
        FacePointer fi = (*fpvi);
        if (!fi->IsD() && !fi->IsV()) {

          fi->SetV();     // mark face as Visited
          faces.push(fi); // each face added to stack is "oriented"

          while (!faces.empty()) {
            FacePointer fp = faces.top();
            faces.pop();

            // orient consistently the adjacent faces
            for (int j = 0; j < 3; j++) {
              if (!IsBorder(*fp, j) && IsManifold(*fp, j)) {
                FacePointer fp_j = fp->FFp(j);
                int iaux = fp->FFi(j);
                if (!CheckOrientation(*fp_j, iaux)) {
                  IsOriented = false;
                  if (!fp_j->IsV()) {
                    SwapEdge(*fp_j, iaux);
                  }
                  else {
                    IsOrientable = false;
                    break;
                  }
                }
                if (!fp_j->IsV()) {
                  fp_j->SetV();
                  faces.push(fp_j);
                }
              }
            }
          }
        }
        if (!IsOrientable) {
          break;
        }
      }
    }
  }

  //-------------------------------------
  // Taubin smoothing
  //-------------------------------------
  class SmoothInfo {
  public:
    SmoothInfo(const CoordType& _p, const int _n) : sum(_p), cnt(_n) {}
    SmoothInfo() : sum(CoordType(0, 0, 0)), cnt(0) {}
    CoordType sum;
    ScalarType cnt;
  };

  template <class STL_CONT, class ATTR_TYPE>
  class SmoothData
  {
  public:
    typedef typename STL_CONT::value_type       stl_vtype;
    typedef typename STL_CONT::iterator         stl_it;
    typedef typename STL_CONT::const_iterator   stl_cit;

    const STL_CONT& c;
    std::vector<ATTR_TYPE>  data;

    SmoothData(const STL_CONT& _c) : c(_c) {
      data.reserve(c.capacity()); data.resize(c.size());
    };

    SmoothData(const STL_CONT& _c, const ATTR_TYPE& val) : c(_c) {
      data.reserve(c.capacity()); data.resize(c.size()); Init(val);
    };

    ~SmoothData() { data.clear(); }
    void Init(const ATTR_TYPE& val) { std::fill(data.begin(), data.end(), val); }

    ATTR_TYPE& operator[](const stl_vtype& v) { return data[&v - &*c.begin()]; }
    ATTR_TYPE& operator[](const stl_vtype* v) { return data[v - &*c.begin()]; }
    ATTR_TYPE& operator[](const stl_cit& cont) { return data[&(*cont) - &*c.begin()]; }
    ATTR_TYPE& operator[](const stl_it& cont) { return data[&(*cont) - &*c.begin()]; }
    ATTR_TYPE& operator[](size_t i) { return data[i]; }

    const ATTR_TYPE& operator[](const stl_vtype& v)  const { return data[&v - &*c.begin()]; }
    const ATTR_TYPE& operator[](const stl_vtype* v)  const { return data[v - &*c.begin()]; }
    const ATTR_TYPE& operator[](const stl_cit& cont) const { return data[&(*cont) - &*c.begin()]; }
    const ATTR_TYPE& operator[](const stl_it& cont)  const { return data[&(*cont) - &*c.begin()]; }
    const ATTR_TYPE& operator[](size_t i)            const { return data[i]; }
  };

  // Accumulate over a SmoothData all positions of adjacent vertices
  void AccumulateSmoothInfo(MeshType& m, SmoothData<VertContainer, SmoothInfo>& TD)
  {
    // bool cotangentFlag = false;
    float weight = 1.0f;
    FaceIterator fi;
    for (fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int j = 0; j < 3; ++j) {
          if (!(*fi).IsB(j)) {

            // if (cotangentFlag) {
            //   float angle = Angle(fi->P1(j) - fi->P2(j), fi->P0(j) - fi->P2(j));
            //   weight = tan((M_PI * 0.5) - angle);
            // }

            TD[(*fi).V0(j)].sum += (*fi).P1(j) * weight;
            TD[(*fi).V1(j)].sum += (*fi).P0(j) * weight;
            TD[(*fi).V0(j)].cnt += weight;
            TD[(*fi).V1(j)].cnt += weight;
          }
        }
      }
    }

    // Reset data for border vertices
    for (fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int j = 0; j < 3; ++j) {
          if ((*fi).IsB(j)) {
            TD[(*fi).V0(j)].sum = (*fi).P0(j);
            TD[(*fi).V1(j)].sum = (*fi).P1(j);
            TD[(*fi).V0(j)].cnt = 1;
            TD[(*fi).V1(j)].cnt = 1;
          }
        }
      }
    }

    for (fi = m.face.begin(); fi != m.face.end(); ++fi) {
      if (!(*fi).IsD()) {
        for (int j = 0; j < 3; ++j) {
          if ((*fi).IsB(j)) {
            TD[(*fi).V(j)].sum += (*fi).V1(j)->P();
            TD[(*fi).V1(j)].sum += (*fi).V(j)->P();
            ++TD[(*fi).V(j)].cnt;
            ++TD[(*fi).V1(j)].cnt;
          }
        }
      }
    }
  }


  //-------------------------------------------------------
  // Taubin smoothing: lambda/mu values
  // 
  // Assume:  0 < lambda, 
  // Assume:  mu is neg. scale factor s.t. mu < (-lambda).
  // Assume:  (mu + lambda) < 0   (i.e. |mu| > lambda )
  //
  // Let kpb be the pass-band frequency, Taubin says that:
  //             kpb = (1/lambda + 1/mu) > 0
  //
  // Taubin: "values of kpb from 0.01 to 0.1 produce good results"
  //
  // kpb * mu - mu/lambda = 1
  // mu = 1 / (kpb - 1/lambda)
  //
  // So if
  //    lambda == 0.5, kpb==0.01 -> mu = 1/(0.01 - 2) = -0.502
  //    lambda == 0.5, kpb==0.1  -> mu = 1/(0.1  - 2) = -0.526
  //    lambda == 0.5, kpb==0.2  -> mu = 1/(0.2  - 2) = -0.555 (NBN: isosurfs)
  //-------------------------------------------------------

  void TaubinSmooth(MeshType& m, int steps, float lambda, float mu)
  {
    SmoothInfo lpz(CoordType(0, 0, 0), 0);
    SmoothData<VertContainer, SmoothInfo> TD(m.vert, lpz);
    VertexIterator vi;
    for (int i = 0; i < steps; ++i)
    {
      TD.Init(lpz);
      AccumulateSmoothInfo(m, TD);
      for (vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
        if (!(*vi).IsD() && TD[*vi].cnt > 0) {
          CoordType Delta = TD[*vi].sum / TD[*vi].cnt - (*vi).P();
          (*vi).P() = (*vi).P() + Delta * lambda;
        }
      }

      TD.Init(lpz);
      AccumulateSmoothInfo(m, TD);
      for (vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
        if (!(*vi).IsD() && TD[*vi].cnt > 0) {
          CoordType Delta = TD[*vi].sum / TD[*vi].cnt - (*vi).P();
          (*vi).P() = (*vi).P() + Delta * mu;
        }
      } // for (vi)
    } // for (steps)
  }


  /////////////////////////////////////////////////////////
  // 
  // start: isoSurf/iso_5_quadric.hpp
  // 
  /////////////////////////////////////////////////////////

  template<class T> int IsNAN(T t) { 
    return std::isnan(t) || std::isinf(t);
  }

  // Options for "Quadric Edge Collapse" algorithm
  class TriQuadricParams {
  public:

  //double    BoundaryQuadricWeight = 0.5;
    double    BoundaryQuadricWeight = 1.0;

    bool      PreserveBoundary = false;

  //double    QualityQuadricWeight = 0.001;
  //double    QualityQuadricWeight = 0.010;
    double    QualityQuadricWeight = 0.100;

  //bool      AreaCheck = false;
    bool      AreaCheck = true;

    bool      UseFacePlaneArea = true;
  //bool      UseFacePlaneArea = false;

    bool      OptimalPlacement = true;    // uses Cholesky solve
  //bool      OptimalPlacement = false;   // uses V(1) of edge

    double    QuadricEpsilon = 1e-15;

    // Faces with a "quality" LOWER than QualityThr are penalized. 
    // (i.e. higher value => better quality)
    double    QualityThr = 0.3;
    bool      QualityCheck = true;

    bool      QualityQuadric = true;  // improve mesh quality in planar regions
    bool      QualityWeight = true;   // if true, try to preserve edges in regions of high scalar variance

  //double    QualityWeightFactor = 100.0;  // see also COST_FAC in mesh load.
    double    QualityWeightFactor = 10.0;   // see also COST_FAC in mesh load.

    double    ScaleFactor = 1.0;
    bool      ScaleIndependent = true;    // scales up quadrics for tiny faces
  //bool      ScaleIndependent = false;   // testing...

    TriQuadricParams() { ; }
  };

  // Detach face f from the chain of faces incident on vertex z
  template <class FaceType>
  void VFDetach(FaceType& f, int z)
  {
    if (f.V(z)->VFp() == &f) {
      int fz = f.V(z)->VFi();       // detaching the first face
      f.V(z)->VFp() = f.VFp(fz);
      f.V(z)->VFi() = f.VFi(fz);
    }
    else {
      // scan the list of faces to find face f to be detached
      VFIterator x(f.V(z)->VFp(), f.V(z)->VFi());
      VFIterator y;

      for (;;) {
        y = x;
        ++x;
        assert(x.f != 0);
        if (x.f == &f) {
          y.f->VFp(y.z) = f.VFp(z);   // found target face
          y.f->VFi(y.z) = f.VFi(z);
          break;
        }
      }
    }
  }

  template <class FaceType>
  void VFDetach(FaceType& f) {
    VFDetach(f, 0);
    VFDetach(f, 1);
    VFDetach(f, 2);
  }


  // Heap element: wraps a pointer to a "LocalModification" object

  class MyTriEdgeCollapse;
  typedef MyTriEdgeCollapse LocModType;

  struct HeapElem {
    HeapElem() : locModPtr(nullptr), pri(0) {}
    ~HeapElem() {}
    LocModType* locModPtr;    // pointer to instance of local modifier
    float pri;

    HeapElem(LocModType* _locModPtr);
    bool IsUpToDate() const;

    // Note:
    // In the current case of edge collapse, "priority" means
    // "quadric error". Since we want edges with small errors 
    // to be collapsed first, we need to invert the usual heap 
    // ordering, putting small errors near the top of the heap. 

    bool operator < (const HeapElem& h) const { return (pri > h.pri); }
  };

  typedef typename std::vector<HeapElem>  HeapType;


  class VertexPair {
  public:

    inline VertexPair() { v[0] = v[1] = nullptr; }
    inline VertexPair(VertexType* v0, VertexType* v1) { V(0) = v0; V(1) = v1; }
    void Sort() { if (V(0) < V(0)) std::swap(V(0), V(0)); }
    VertexType*& V(int i) { return v[i]; }
    VertexType* cV(int i) const { return v[i]; }

  private:            // v[0] gets deleted, 
    VertexType* v[2]; // v[1] gets moved to new position
  };


  class EdgeCollapser
  {
    // Utility routines for collapsing an edge in a trimesh

  public:

    typedef typename std::vector<VFIterator>  VFIVec;

  private:
    struct EdgeSet {
      VFIVec av0, av1, av01;
      VFIVec& AV0()  { return av0; }   // Faces incident only on v0
      VFIVec& AV01() { return av01; }  // Faces incident only on both v0 and v1
    };

    static void FindSets(VertexPair& p, EdgeSet& es) {
      VertexType* v0 = p.V(0);
      VertexType* v1 = p.V(1);

      es.AV0().clear();
      es.AV01().clear();

      for (VFIterator x = VFIterator(v0); !x.End(); ++x) {
        bool foundV1 = false;
        for (int j = 0; j < 3; ++j) {
          if (x.f->V(j) == v1) {
            foundV1 = true;
            break;
          }
        }

        // v1 not found ==> face incident only on v0
        if (!foundV1) es.AV0().push_back(x);
        else          es.AV01().push_back(x);
      }
    }

  public:

    // Main collapsing function: collapse edge indicated by VertexPair c. 
    // Note: v[0] gets deleted, v[1] survives (shifted to position p)

    static int Do(MyMesh2& m, VertexPair& c, const Point3<ScalarType>& p)
    {
      EdgeSet es; //, es1;
      FindSets(c, es);

      int n_face_del = 0;

      static int VtoE[3][3] = { -1,  0,  2,
                                 0, -1,  1,
                                 2,  1, -1 };

      std::vector<VertexPointer> topVertices; topVertices.reserve(2);
      std::vector<VertexPointer> fan1V2S; fan1V2S.reserve(2);
      std::vector<VertexPointer> v2s; v2s.reserve(2);
      std::map <VertexPointer, bool> toSel;

      for (auto i = es.AV01().begin(); i != es.AV01().end(); ++i) {
        FaceType& f = *((*i).f);
        assert(f.V((*i).z) == c.V(0));
        iso::VFDetach(f, ((*i).z + 1) % 3);
        iso::VFDetach(f, ((*i).z + 2) % 3);
        iso::DeleteFace(m, f);
        n_face_del++;
      }

      // Low-level update of VF Adjacency:
      // for all the faces incident in v[0]
      // - v[0] will be deleted so we substitute v[0] with v[1]
      // - prepend that face to the list of the faces incident on v[1]

      for (auto i = es.AV0().begin(); i != es.AV0().end(); ++i) {
        FaceType& f = *((*i).f);
        (*i).f->V((*i).z) = c.V(1);	          // For each face in es.AV0, 
        (*i).f->VFp((*i).z) = c.V(1)->VFp();  // substitute v0 with v1
        (*i).f->VFi((*i).z) = c.V(1)->VFi();
        c.V(1)->VFp() = (*i).f;
        c.V(1)->VFi() = (*i).z;
      }

      iso::DeleteVertex(m, *(c.V(0)));
      c.V(1)->P() = p;
      return n_face_del;
    }

  };


  class MyTriEdgeCollapse
  {
    // After a collapse, this class is responsible for 
    // updating the heap of class LocalOptimization 

  protected:

    typedef iso::Quadric                QuadricType;
    typedef TriQuadricParams            QParams;
    typedef typename MyTriEdgeCollapse  MYTYPE;

    VertexPair  m_pos;          // the pair to collapse
    ScalarType  m_priority;     // priority in the heap
    int         localMark;      // mark for up_dating
    CoordType   optimalPos;     // Local storage for computed optimal position

    static int& GlobalMark() {
      static int im = 0;        // mark for up_dating
      return im;
    }

  public:

    MyTriEdgeCollapse() : localMark(0), m_priority(0) {}

    MyTriEdgeCollapse(const VertexPair& p, int mark, TriQuadricParams* pp) {
      localMark = mark;
      m_pos = p;
      m_priority = ComputePriority(pp);
    }

    ~MyTriEdgeCollapse() {}

    // Adjusting the size of the heap vs. num faces in mesh. 
    // The heap is cleared when this ratio is exceeded, so 
    // make it larger than a minimum expected size to avoid 
    // clearing of the heap too frequently.
    // 
    // For symmetric edge collapse, 4 or 5 works well.
    // For non-symmetric edge collapse, 8 or 9 is better.

    static float HeapFaceRatio(TriQuadricParams* _pp) { return IsSymmetric(_pp) ? 4.0f : 8.0f; }
    static bool IsSymmetric(TriQuadricParams* _pp)    { return  ((QParams*)_pp)->OptimalPlacement; }
    static bool IsVertexStable(TriQuadricParams* _pp) { return !((QParams*)_pp)->OptimalPlacement; }

    //-------------------------------------
    // NBN: add miscellaneous routines
    //-------------------------------------

    // Return a shape quality measure (2*AreaTri/(MaxEdge^2)) for 
    // triangle defined by p0,p1,p2. Range is [0.0, 0.866]. E.g.,
    // 0.866 : for equilateral triangles, sqrt(3)/2
    // 0.500 : for a half square
    // 0.000 : for lines
    ScalarType FQuality(Point3<ScalarType> const& p0, Point3<ScalarType> const& p1, Point3<ScalarType> const& p2)
    {
      Point3<ScalarType> d10 = p1 - p0;
      Point3<ScalarType> d20 = p2 - p0;
      Point3<ScalarType> d12 = p1 - p2;
      Point3<ScalarType> x = d10 ^ d20;

      ScalarType a = Norm(x);
      if (a == 0) return 0; // zero area ==> zero quality
      ScalarType b = SquaredNorm(d10);
      if (b == 0) return 0; // zero area ==> zero quality
      ScalarType t = b;
      t = SquaredNorm(d20); if (b < t) b = t;
      t = SquaredNorm(d12); if (b < t) b = t;
      return a / b;
    }

    ScalarType FaceQuality(const FaceType& t) {
      return FQuality(t.cP(0), t.cP(1), t.cP(2));
    }

    // Evaluate the priority (error) for an edge collapse
    // 
    // Simulates the collapse, computing a quadric error 
    // generated by this collapse. Error is weighted by:
    // 
    //  - aspect ratio of the triangles involved,
    //  - other measures removed (NBN:)

    virtual ScalarType ComputePriority(TriQuadricParams* _pp)
    {
      QParams* pp = (QParams*)_pp;

      VertexType* v[2];   // two vertices defining the edge
      v[0] = m_pos.V(0);
      v[1] = m_pos.V(1);

      ScalarType origArea = 0;
      if (pp->AreaCheck) {
        for (VFIterator x(v[0]); !x.End(); ++x)	    // Accumulate areas for...
          origArea += DoubleArea(*x.F());           // all faces in v0

        for (VFIterator x(v[1]); !x.End(); ++x)	    // Accumulate areas for...
          if (x.V1() != v[0] && x.V2() != v[0])     // all faces involving v1
            origArea += DoubleArea(*x.F());
      }

      // Move the two vertexes into new position (storing the old ones)
      CoordType OldPos0 = v[0]->P(), OldPos1 = v[1]->P();
      ComputePosition(_pp);

      //---------------------------------
      // Now Simulate the collapse 
      //---------------------------------
      v[0]->P() = v[1]->P() = this->optimalPos;

      // Rescan faces and compute the worst quality that happens after collapse

      ScalarType newQual = std::numeric_limits<ScalarType>::max();
      if (pp->QualityCheck) {
        for (VFIterator x(v[0]); !x.End(); ++x)
          if (x.V1() != v[1] && x.V2() != v[1])
            newQual = std::min(newQual, FaceQuality(*x.F()));

        for (VFIterator x(v[1]); !x.End(); ++x)
          if (x.V1() != v[0] && x.V2() != v[0])
            newQual = std::min(newQual, FaceQuality(*x.F()));
      }

      ScalarType newArea = 0;
      if (pp->AreaCheck) {
        for (VFIterator x(v[0]); !x.End(); ++x)   // calculate the "change of area" 
          newArea += DoubleArea(*x.F());          // that this collapse would cause

        for (VFIterator x(v[1]); !x.End(); ++x)
          if (x.V1() != v[0] && x.V2() != v[0])
            newArea += DoubleArea(*x.F());
      }

      // add quadrics for both vertices
      QuadricType qq = v[0]->Qd(); qq += v[1]->Qd();

      double QuadErr = pp->ScaleFactor * qq.Apply(Point3d::Construct(v[1]->P()));
      assert(!IsNAN(QuadErr));

      // Add scaled variance to cost of collapsing this edge
      if (pp->QualityWeight) {
        QuadErr *= (1 + v[1]->cQ());  // *= : add 1 in case cost is zero
      //QuadErr += (    v[1]->cQ());  // += : simple addition
      }


      // The lower the newQual measure (below QualityThr), the 
      // larger the scaling of the cost for removing that edge.
      if (newQual > pp->QualityThr) {
        newQual = pp->QualityThr;             // clamp the newQual measure
      }

      QuadErr = std::max(QuadErr, pp->QuadricEpsilon);
      if (QuadErr <= pp->QuadricEpsilon) {
        QuadErr *= Distance(OldPos0, OldPos1);
      }

      ScalarType error = 0;
      if (!pp->QualityCheck) {
        error = (ScalarType)(QuadErr);
      }
      else {
        error = (ScalarType)(QuadErr / newQual);
      }

    //if (pp->AreaCheck && ((fabs(origArea - newArea) / (origArea + newArea)) > 0.01))
      if (pp->AreaCheck && ((fabs(origArea - newArea) / (origArea + newArea)) > 0.10))
      {
        error = std::numeric_limits<ScalarType>::max();
      }

      // Restore old position of v0 and v1
      v[0]->P() = OldPos0;
      v[1]->P() = OldPos1;

      this->m_priority = error;
      return this->m_priority;
    }


    void AddCollapseToHeap(HeapType& h_ret, VertexType* v0, VertexType* v1, TriQuadricParams* _pp)
    {
      QParams* pp = (QParams*)_pp;

      // TODO: consider setting a cutoff value such that,
      // if (edge.pri > value), do not collapse this edge
      constexpr ScalarType maxAdmitErr = std::numeric_limits<ScalarType>::max();

      h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v0, v1), this->GlobalMark(), _pp)));
      if (h_ret.back().pri > maxAdmitErr) {
        delete h_ret.back().locModPtr;        // NBN: protect this edge?
        h_ret.pop_back();
      }
      else {
        std::push_heap(h_ret.begin(), h_ret.end());
      }

      if (!IsSymmetric(pp)) {
        h_ret.push_back(HeapElem(new MYTYPE(VertexPair(v1, v0), this->GlobalMark(), _pp)));
        if (h_ret.back().pri > maxAdmitErr) {
          delete h_ret.back().locModPtr;
          h_ret.pop_back();                   // NBN: protect this edge?
        }
        else {
          std::push_heap(h_ret.begin(), h_ret.end());
        }
      }
    }

    void Execute(MyMesh2& m, TriQuadricParams* /*_pp*/)
    {
      CoordType newPos = this->optimalPos;

      // v0 is deleted and v1 takes the new position
      (m_pos.V(1))->Qd() += (m_pos.V(0))->Qd();

      EdgeCollapser::Do(m, m_pos, newPos);
    }

    void UpdateHeap(HeapType& h_ret, TriQuadricParams* _pp)
    {
      this->GlobalMark()++;
      VertexType* v[2];
      v[0] = m_pos.V(0); v[1] = m_pos.V(1);
      v[1]->IMark() = this->GlobalMark();

      // First loop around the surviving vertex to unmark the Visit flags
      for (VFIterator vfi(v[1]); !vfi.End(); ++vfi) {
        vfi.V1()->ClearV();
        vfi.V2()->ClearV();
        vfi.V1()->IMark() = this->GlobalMark();
        vfi.V2()->IMark() = this->GlobalMark();
      }

      // Second Loop
      for (VFIterator vfi(v[1]); !vfi.End(); ++vfi) {
        if (!(vfi.V1()->IsV()) && vfi.V1()->IsRW()) {
          vfi.V1()->SetV();
          AddCollapseToHeap(h_ret, vfi.V0(), vfi.V1(), _pp);
        }
        if (!(vfi.V2()->IsV()) && vfi.V2()->IsRW()) {
          vfi.V2()->SetV();
          AddCollapseToHeap(h_ret, vfi.V2(), vfi.V0(), _pp);
        }
        if (vfi.V1()->IsRW() && vfi.V2()->IsRW()) {
          AddCollapseToHeap(h_ret, vfi.V1(), vfi.V2(), _pp);
        }
      } // end second loop around surviving vertex.
    }

    // has any data changed?
    bool IsUpToDate() const {
      VertexType* v0 = m_pos.cV(0);
      VertexType* v1 = m_pos.cV(1);
      if (v0->IsD() || v1->IsD() ||
        localMark < v0->IMark() ||
        localMark < v1->IMark())
      {
        return false;
      }
      return true;
    }

    ScalarType Priority() const {
      return m_priority;            // priority to be used in the heap
    }

    static void Init(MyMesh2& m, HeapType& h_ret, TriQuadricParams* _pp)
    {
      QParams* pp = (QParams*)_pp;
      h_ret.clear();
      iso::VertexFace(m);
      iso::FaceBorderFromVF(m);

      if (pp->PreserveBoundary) {

        // Preserve boundary edges by disabling the [W]rite 
        // flag for both vertices of each boundary edge
        for (auto pf = m.face.begin(); pf != m.face.end(); ++pf) {
          if (!(*pf).IsD() && (*pf).IsW()) {
            for (int j = 0; j < 3; ++j) {
              if ((*pf).IsB(j)) {
                (*pf).V(j)->ClearW();
                (*pf).V1(j)->ClearW();
              }
            }
          }
        }
      }

      InitQuadric(m, pp);

      // Initialize the heap with all the possible collapses
      if (IsSymmetric(pp))
      {
        // if the collapse is symmetric (e.g. u->v == v->u)
        for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
          if (!(*vi).IsD() && (*vi).IsRW()) {
            VFIterator x;
            for (x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F() != 0; ++x) {
              x.V1()->ClearV();
              x.V2()->ClearV();
            }
            for (x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F() != 0; ++x) {
              if ((x.V0() < x.V1()) && x.V1()->IsRW() && !x.V1()->IsV()) {
                x.V1()->SetV();
                h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(), x.V1()), GlobalMark(), _pp)));
              }
              if ((x.V0() < x.V2()) && x.V2()->IsRW() && !x.V2()->IsV()) {
                x.V2()->SetV();
                h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V0(), x.V2()), GlobalMark(), _pp)));
              }
            }
          }
        }
      }
      else
      {
        // collapse is not symmetric (e.g. u->v != v->u)
        for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
          if (!(*vi).IsD() && (*vi).IsRW()) {
            VFIterator x;
            UnMarkAll(m);
            for (x.F() = (*vi).VFp(), x.I() = (*vi).VFi(); x.F() != 0; ++x) {
              if (x.V()->IsRW() && x.V1()->IsRW() && !IsMarked(m, x.F()->V1(x.I()))) {
                h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V(), x.V1()), GlobalMark(), _pp)));
              }
              if (x.V()->IsRW() && x.V2()->IsRW() && !IsMarked(m, x.F()->V2(x.I()))) {
                h_ret.push_back(HeapElem(new MYTYPE(VertexPair(x.V(), x.V2()), GlobalMark(), _pp)));
              }
            }
          }
        }
      }
    }


    static void InitQuadric(MyMesh2& m, TriQuadricParams* _pp)
    {
      QParams* pp = (QParams*)_pp;
      for (VertexIterator pv = m.vert.begin(); pv != m.vert.end(); ++pv) {
        if (!(*pv).IsD() && (*pv).IsW()) {
          (pv->Qd()).SetZero();
        }
      }

      for (FaceIterator fi = m.face.begin(); fi != m.face.end(); ++fi) {
        if (!(*fi).IsD() && (*fi).IsR()) {
          if ((*fi).V(0)->IsR() && (*fi).V(1)->IsR() && (*fi).V(2)->IsR()) {
            Plane3<ScalarType, false> facePlane;
            facePlane.SetDirection(((*fi).V(1)->cP() - (*fi).V(0)->cP()) ^ ((*fi).V(2)->cP() - (*fi).V(0)->cP()));
            if (!pp->UseFacePlaneArea) {
              facePlane.Normalize();
            }
            facePlane.SetOffset(facePlane.Direction().dot((*fi).V(0)->cP()));

            QuadricType q;
            q.ByPlane(facePlane);

            // Add "facePlane" quadric to each vertex
            for (int j = 0; j < 3; ++j) {
              if ((*fi).V(j)->IsW()) {
                (fi->V(j))->Qd() += q;
              }
            }

            for (int j = 0; j < 3; ++j) {

              // Handle border faces and cases where QualityQuadric is [ON]

              if (fi->IsB(j) || pp->QualityQuadric) {
                Plane3<ScalarType, false> borderPlane;
                QuadricType bq;
                // Border quadric record the squared distance from the plane 
                // orthogonal to the face and passing through the edge. 
                borderPlane.SetDirection(facePlane.Direction() ^ ((fi->V1(j)->cP() - fi->V(j)->cP()).normalized()));
                if (fi->IsB(j)) {
                  // allow amplification of boundary edges using "BoundaryQuadricWeight"
                  borderPlane.SetDirection(borderPlane.Direction() * (ScalarType)(pp->BoundaryQuadricWeight));
                }
                else {
                  // allow adjustement of effect of Quality using "QualityQuadricWeight"
                  borderPlane.SetDirection(borderPlane.Direction() * (ScalarType)(pp->QualityQuadricWeight));
                }
                borderPlane.SetOffset(borderPlane.Direction().dot(fi->V(j)->cP()));
                bq.ByPlane(borderPlane);
                if (fi->V(j)->IsW())  (fi->V(j))->Qd() += bq;
                if (fi->V1(j)->IsW())	(fi->V1(j))->Qd() += bq;
              }
            }
          }
        }
      }

      if (pp->ScaleIndependent) {

        UpdateBounding(m);
        // Make all quadric factors independent of mesh size
        pp->ScaleFactor = 1e8 * pow(1.0 / m.m_bbox.Diag(), 6); // scaling factor

#if (1)
        nnMSG(1, "          (ScaleIndependent) scale factor: %0.3e\n", pp->ScaleFactor);
#endif
      }

      if (pp->QualityWeight) {
        // Note: "cost" vi->cQ() is calculated when mesh is loaded
        for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
          if (!(*vi).IsD() && (*vi).IsW()) {
            vi->Qd() *= ((1 + vi->cQ()) * pp->QualityWeightFactor);
          }
        }
      }
    }


    static void Finalize(MyMesh2& m, HeapType& /*h_ret*/, TriQuadricParams* _pp)
    {
      // Final clean up after simplification process (if required)
      QParams* pp = (QParams*)_pp;

      // If boundary edges were protected, re-enable [W]rite flags
      // Note: here we just make ALL vertices [W]ritable
      if (pp->PreserveBoundary) {
        for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi) {
          if (!(*vi).IsD()) (*vi).SetW();
        }
      }

#if (0)
      {
        // NBN: check range of final quadric error vaules
        double cQ_err = 0, Qd_err = 0;
        FILE* g_fpQE = fopen("QErr_Check.txt", "w");
        fprintf(g_fpQE, " Quadric   Variance\n");
        for (int i = 0; i < m.vert.size(); ++i) {
          VertexType& vi = m.vert[i];
          if (!vi.IsD()) {
            Qd_err = vi.Qd().m_val;   // final quadric + variance
            cQ_err = vi.cQ();         // color variance
            fprintf(g_fpQE, "%10.3e  %10.3e\n", Qd_err, cQ_err);
          }
        }
        fclose(g_fpQE);
      }
#endif
    }

    void ComputePosition(TriQuadricParams* _pp)
    {
      QParams* pp = (QParams*)_pp;

      CoordType newPos;
      if (pp->OptimalPlacement == false) {
        newPos = m_pos.V(1)->P();
      }
      else
      {
        // initialize with average
        newPos = (m_pos.V(0)->P() + m_pos.V(1)->P()) / 2.0;

        if ((m_pos.V(0)->Qd().Apply(newPos) + m_pos.V(1)->Qd().Apply(newPos)) > 2.0 * pp->QuadricEpsilon)
        {
          // Find optimal vertex location, i.e. that point p which
          // minimizes error ("cost") associated with edge [v1,v2]

          QuadricType q = (m_pos.V(0))->Qd();
          q += (m_pos.V(1))->Qd();

          bool border = ((m_pos.V(0))->IsB() && (m_pos.V(1))->IsB());

          // determinant of symmetric 3x3 matrix
          double det = q.sym_det();
          double absdet = std::abs(det);
          constexpr double det_TOL = 1e-13;

          if ((absdet > det_TOL) && !border)
          {
            // q is PSD/invertible: solve for optimal vertex location
            Point3d p;
            q.Minimum(p);
            newPos = p;
          }
          else {

            // q not invertable (or edge is on border).
            // select the vertex that minimizes error
            CoordType p0 = m_pos.V(0)->P();
            CoordType p1 = m_pos.V(1)->P();
            CoordType p2 = (p0 + p1) / 2;

            double error0 = q.Apply(p0);
            double error1 = q.Apply(p1);
            double error2 = q.Apply(p2);

            double err = std::min(error0, std::min(error1, error2));

#if (1)
            if (err < 1e-15) {
              // NBN: checking for negative error
              // nnTRC(1, "err < 1e-15 : %12.3e\n", err);
              err = 0.0;
            }
#endif

            if (error0 == err) newPos = p0;
            else if (error1 == err) newPos = p1;
            else if (error2 == err) newPos = p2;
          }
        }
      }
      this->optimalPos = newPos;
    }

  };  // end class MyTriEdgeCollapse


  // Define the HeapElem routines declared above
  // (reqires that LocModType has beed defined).
  HeapElem::HeapElem(LocModType* _locModPtr) {
    locModPtr = _locModPtr;
    pri = float(locModPtr->Priority());
  }

  bool HeapElem::IsUpToDate() const {
    return locModPtr->IsUpToDate();
  }



  //-------------------------------------------------------
  // Driver class for edge_collapse modification of a mesh. 
  //-------------------------------------------------------
  class LocalOptimization
  {
  public:

    MeshType& m_m;    // the mesh to optimize
    HeapType m_h;     // the heap of operations

    LocalOptimization(MeshType& msh, TriQuadricParams* _pp)
      : m_m(msh)
    {
      ClearTermination();
      HeapFaceRatio = 5;
      pp = _pp;
    }

    // termination conditions:
    enum LOTermination {
      LOnSimplices = 0x01,  // test number of simplicies	
      LOnVertices  = 0x02,  // test number of verticies
      LOnOps       = 0x04,  // test number of operations
      LOMetric     = 0x08,  // test Metric (error, quality...instance dependent)
      LOTime       = 0x10   // test how much time is passed since the start
    };

    int m_tf;   // Termination Flag

    int nPerformedOps, nTargetOps, nTargetSimplices, nTargetVertices;
    float	timeBudget;
    clock_t	start;
    ScalarType currMetric;
    ScalarType targetMetric;

    TriQuadricParams* pp;

    // Ratio between Heap size and number of faces in current mesh.
    // Use this to determine when to clear the heap.
    float HeapFaceRatio;

    void SetTerminationFlag(int v) { m_tf |= v; }
    void ClearTerminationFlag(int v) { m_tf &= ~v; }
    bool IsTerminationFlag(int v) { return ((m_tf & v) != 0); }

    void SetTargetSimplices(int ts) { nTargetSimplices = ts; SetTerminationFlag(LOnSimplices); }
    void SetTargetVertices(int tv) { nTargetVertices = tv;  SetTerminationFlag(LOnVertices); }
    void SetTargetOperations(int to) { nTargetOps = to;	     SetTerminationFlag(LOnOps); }
    void SetTargetMetric(ScalarType tm) { targetMetric = tm;     SetTerminationFlag(LOMetric); }
    void SetTimeBudget(float tb) { timeBudget = tb;       SetTerminationFlag(LOTime); }

    void ClearTermination() {
      m_tf = 0;
      nTargetSimplices = 0;
      nTargetOps = 0;
      targetMetric = 0;
      timeBudget = 0;
      nTargetVertices = 0;
    }

    ~LocalOptimization() {
      typename HeapType::iterator i;
      for (i = m_h.begin(); i != m_h.end(); i++)
        delete (*i).locModPtr;
    }

    // Main cycle of optimization:
    bool DoOptimization()
    {
      assert(((m_tf & LOnSimplices) == 0) || (nTargetSimplices != -1));
      assert(((m_tf & LOnVertices)  == 0) || (nTargetVertices  != -1));
      assert(((m_tf & LOnOps)       == 0) || (nTargetOps       != -1));
      assert(((m_tf & LOMetric)     == 0) || (targetMetric     != -1));
      assert(((m_tf & LOTime)       == 0) || (timeBudget       != -1));

      start = clock();
      nPerformedOps = 0;
      while (!GoalReached() && !m_h.empty()) {

        if (m_h.size() > m_m.FaceNumber() * HeapFaceRatio)  ClearHeap();
        std::pop_heap(m_h.begin(), m_h.end());
        LocModType* locMod = m_h.back().locModPtr;
        currMetric = m_h.back().pri;
        m_h.pop_back();

        if (locMod->IsUpToDate()) {
          nPerformedOps++;
          locMod->Execute(m_m, this->pp);
          locMod->UpdateHeap(m_h, this->pp);
        }
        delete locMod;
      }
      return !(m_h.empty());
    }

    // Remove all "out-of-date" operations from the heap 
    // (e.g. collapses with recently modified vertices). 
    // Called when heap is larger than a specified ratio.

    void ClearHeap() {
      for (auto hi = m_h.begin(); hi != m_h.end();) {
        if (!(*hi).locModPtr->IsUpToDate()) {
          delete (*hi).locModPtr;
          *hi = m_h.back();
          if (&*hi == &m_h.back()) {
            hi = m_h.end();
            m_h.pop_back();
            break;
          }
          m_h.pop_back();
          continue;
        }
        ++hi;
      }
      std::make_heap(m_h.begin(), m_h.end());
    }


    template <class MeshModifierType>
    void Init() {

      // Initialize flags and quadrics for all vertices in 
      // the mesh. Called only at the start of decimation.

      iso::InitVertexIMark(m_m);
      HeapFaceRatio = MeshModifierType::HeapFaceRatio(pp);
      MeshModifierType::Init(m_m, m_h, pp);
      std::make_heap(m_h.begin(), m_h.end());
      if (!m_h.empty()) currMetric = m_h.front().pri;

#if (1)
      if (pp->OptimalPlacement) {
        nnMSG(1, "          OptimalPlacement: ON, heap size is %d\n", (int)m_h.size());
      }
      else {
        nnMSG(1, "          OptimalPlacement: OFF, heap size is %d\n", (int)m_h.size());
      }
#endif
    }

    template <class MeshModifierType>
    void Finalize() {
      MeshModifierType::Finalize(m_m, m_h, pp);
    }

    // Check if the process is to end or not: the process ends 
    // when any of the termination conditions is verified. 
    bool GoalReached() {
      if (IsTerminationFlag(LOnSimplices) && (m_m.FaceNumber() <= nTargetSimplices)) return true;
      if (IsTerminationFlag(LOnVertices)  && (m_m.VertexNumber() <= nTargetVertices)) return true;
      if (IsTerminationFlag(LOnOps)       && (nPerformedOps == nTargetOps)) return true;
      if (IsTerminationFlag(LOMetric)     && (currMetric > targetMetric))  return true;

      if (IsTerminationFlag(LOTime)) {
        clock_t cur = clock();
        if (cur < start)
          return true; // should not happen
        else if ((cur - start) / (double)CLOCKS_PER_SEC > timeBudget)
          return true;
      }
      return false;
    }
  }; // end class LocalOptimization


  //-------------------------------------------------------
  // Driver routine for Quadric edge collapse
  //-------------------------------------------------------

  void EdgeCollapse(MeshType& m, float fac)
  {
    nnMSG(1, "[proc:%02d] Entering Quadric EdgeCollapse(m, %0.2f)\n", g_procid, fac);

    int Nv = m.VN(), Nf = m.FN();
    int TargetNf = (int)(fac * (float)Nf);
    constexpr double TargetError = std::numeric_limits<float>::max();
    nnMSG(1, "          mesh has %d verts and %d tri faces\n"
             "          target Nfaces (%0.2f * %d) = %d\n", Nv, Nf, fac, Nf, TargetNf);

    TriQuadricParams qparams;

    // decimator initialization
    iso::LocalOptimization DeciSession(m, &qparams);
    DeciSession.Init<MyTriEdgeCollapse>();

    DeciSession.SetTargetSimplices(TargetNf);
    DeciSession.SetTargetOperations(100000);
    DeciSession.SetTargetMetric(TargetError);

#ifndef _DEBUG
    // time for each DoOptimization() iteration
    DeciSession.SetTimeBudget(0.5f);
#endif

    double NfToDel = m.m_fn - TargetNf, progress = 0.0;
    while (DeciSession.DoOptimization() && m.m_fn > TargetNf) {
      progress = 100.0 - (100.0 * (m.m_fn - TargetNf) / (NfToDel));
      nnMSG(1, "[proc:%02d] Simplifying... %5.2lf \n", g_procid, progress);
    };
    DeciSession.Finalize<MyTriEdgeCollapse>();
  }

}   // end namespace iso
