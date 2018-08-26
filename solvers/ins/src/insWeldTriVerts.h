/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

// holmes/solvers/gradient/src/weldTriVerts.h
// weld shared verts in tri mesh, remove degenerate tris
// 2018/06/01
//---------------------------------------------------------
#ifndef HOLMES__WELD_TRI_VERTS__H
#define HOLMES__WELD_TRI_VERTS__H

#include <vector>
#include <set>

//---------------------------------------------------------
class vertexPos
//---------------------------------------------------------
{
public:
  vertexPos() { 
    mPos[0]=0.0;  mPos[1]=0.0;  mPos[2]=0.0; 
  }
  vertexPos(const dfloat *p) {
    mPos[0]=p[0]; mPos[1]=p[1]; mPos[2]=p[2];
  }
  
  void set(const dfloat *p) {
    mPos[0]=p[0]; mPos[1]=p[1]; mPos[2]=p[2]; 
  }

  void set(int index, const dfloat *pos) {
    const dfloat *p = &pos[index*3];
    mPos[0]=p[0]; mPos[1]=p[1]; mPos[2]=p[2];
  }

  dfloat getX() const { return mPos[0]; };
  dfloat getY() const { return mPos[1]; };
  dfloat getZ() const { return mPos[2]; };

  // coords of vertex position
  dfloat mPos[3];
};


//---------------------------------------------------------
template <typename Type> class vertexLess
//---------------------------------------------------------
{
public:

  typedef std::vector<Type> vertexVector;

  bool operator()(int v1,int v2) const {
    const Type& a = get(v1);
    const Type& b = get(v2);

    // convert to scaled int to avoid floating point ambiguity
    // alternative: enable tolerance for merging vertices
    const dfloat scale = 1.0e5;

    int ixA = (int)(a.getX()*scale);
    int ixB = (int)(b.getX()*scale);

    if (ixA < ixB) return true;
    if (ixA > ixB) return false;

    int iyA = (int)(a.getY()*scale);
    int iyB = (int)(b.getY()*scale);

    if (iyA < iyB) return true;
    if (iyA > iyB) return false;

    int izA = (int)(a.getZ()*scale);
    int izB = (int)(b.getZ()*scale);

    if (izA < izB) return true;
    if (izA > izB) return false;

    return false;
  }

  static void setSearch(const Type& match, vertexVector *list) {
    mFind = match;
    mList = list;
  }

private:

  const Type& get(int index) const {
    if ( index == -1 ) return mFind;
    vertexVector &vlist = *mList;
    return vlist[index];
  }

  static Type          mFind;   // vertex to find
  static vertexVector *mList;   // pointer to list of vertices
//static dfloat        s_tol;   // tolerance for merging vertices
};


//---------------------------------------------------------
template <typename Type> class vertexPool
//---------------------------------------------------------
{
public:

  typedef std::set<int, vertexLess<Type> > vertexSet;
  typedef std::vector<Type> vertexVector;

  int getVertex(const Type& vtx){
    vertexLess<Type>::setSearch(vtx,&mVtxs);
    typename vertexSet::iterator found;
    found = mVertSet.find( -1 );
    if (found != mVertSet.end()){
      return (*found);
    }
    int idx = (int)mVtxs.size();
    mVtxs.push_back( vtx );
    mVertSet.insert( idx );
    return idx;
  }

  const dfloat* getPos(int idx) const { return mVtxs[idx].mPos; }
  const Type&   get(int idx) const    { return mVtxs[idx]; }
  int           getSize() const       { return (int)mVtxs.size(); }

  void clear(int reservesize){
    mVertSet.clear();
    mVtxs.clear();
    mVtxs.reserve(reservesize);
  }

  const vertexVector& getVertexList() const { return mVtxs; }
  void  set(const Type& vtx)    { mVtxs.push_back(vtx); }
  int   getVertexCount() const  { return (int)mVtxs.size(); }
  Type* getBuffer()             { return &mVtxs[0]; }

private:
  
  // 1. use a set to accumulate unique vertices, 
  // 2. use a vector to associate an implicit id,
  // 3. use these ids for element connectivity.

  vertexSet     mVertSet;   // set of unique vertices
  vertexVector  mVtxs;      // same verts, with id from array index
};


//---------------------------------------------------------
// define static members
//---------------------------------------------------------
template <typename Type>             Type   vertexLess<Type>::mFind;
template <typename Type> std::vector<Type> *vertexLess<Type>::mList=0;


// simplify syntax
typedef vertexPool<vertexPos> vertexLookup;

#endif  // HOLMES__WELD_TRI_VERTS__H
