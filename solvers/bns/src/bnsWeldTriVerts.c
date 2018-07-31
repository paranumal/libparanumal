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

// holmes/solvers/bns/src/weldTriVerts.cpp
// weld shared verts in tri mesh, remove degenerate tris
// 2018/06/02
//---------------------------------------------------------
#include "bns.h"
#include "bnsWeldTriVerts.h"


//---------------------------------------------------------
int bnsWeldTriVerts(bns_t *bns, int Ntris, double *isoq){
  int Max_N = Ntris*12;                   // upper bound for node data
  int Max_T = Ntris*3;    // upper bound for connectivity ids

  std::vector<int>& refT = bns->iso_tris;
  refT.resize(Max_T);   // storage for 3 node ids for each tri

  std::vector< std::vector<double> > Q;
  Q.resize(bns->isoNfields, vector<double>(Max_T)); // for multiple field variables

  // std::vector<double> Q;
  // Q.resize(Max_T*bns->isoNfields);      // storage for 1 scalar for each node

  mesh_t *mesh = bns->mesh;
  // const int procid = 0;
  // const int procid = mesh->m_procid;
  int Nfields = bns->isoNfields;
  int dim = mesh->dim;
  
  //---------------------------------------------
  // build set of unique verts for tri mesh
  //---------------------------------------------
  vertexLookup* VL = new vertexLookup;

  double p[3][3]={{0.0}};
  double q[bns->isoNfields][3]={{0.0}}; 
  int id=0, v[3]={0}; 
  int nodesk=0;
  vertexPos vpos;

  int sk=0, Nbadfaces=0, n0,n1,n2; //, id0,id1,id2;
  
  int fld=0;  //  3    1        3     0
  int stride = (dim + Nfields)*dim + fld;
  int noff   = (dim + Nfields); 
  int koff   = 0;
  int skT=0; 

  for (int k=0; k<Ntris; ++k) {

    int  koff = nodesk + k*stride;  // offset to kth triangle
    n0 = koff + 0*noff;             // offsets to each node
    n1 = koff + 1*noff;
    n2 = koff + 2*noff;

    // get {x,x,x, y,y,y, z,z,z} for 3 nodes in kth tri:
    for (int i=0; i<3; ++i) {   // 0 1 2
      p[0][i] = isoq[n0+i];     // x,y,z
      p[1][i] = isoq[n1+i];     // x,y,z
      p[2][i] = isoq[n2+i];     // x,y,z
    }

    // get {q, q, q}, for 3 nodes in kth tri:
    for(int fld = 0; fld<bns->isoNfields; fld++){
      q[fld][0] = isoq[n0+(bns->dim + fld)];
      q[fld][1] = isoq[n1+(bns->dim + fld)];
      q[fld][2] = isoq[n2+(bns->dim + fld)];
    }

    for (int i=0; i<3; ++i) {     // for each node in the tri
      vpos.set(p[i]);
      id = VL->getVertex(vpos);   // insert node, get id

      for(int fld = 0; fld<bns->isoNfields; fld++)
        Q[fld][id] = q[fld][i];               // scalar for this node

      refT[skT++] = id;           // id of ith node in triangle k
      v[i] = id;                  // for checking validity of tri
    }

    if (v[0]==v[1] || v[0]==v[2] || v[1]==v[2]) {
      
      // degenerate face: reset vertex ids and roll back insertion marker skT
      ++Nbadfaces;

      refT[skT-1] = refT[skT-2] = refT[skT-3] = 0;
      skT -= 3;
    }
  }
  int Ntris_Final = Ntris;
  if (Nbadfaces>0) {
    printf("\n*** removed %d degenerate faces ***\n\n", Nbadfaces);
    Ntris_Final -= Nbadfaces;

    // resize truncated elem list
    refT.resize(Ntris_Final*3);
  }

  // load node coods and scalar
  int nno = VL->getVertexCount();
  std::vector<double>& refN = bns->iso_nodes;
  refN.resize(nno*(bns->dim + bns->isoNfields));
  int skN=0;
  for (int i=0; i<nno; ++i) {
    refN[skN++] = VL->getPos(i)[0];
    refN[skN++] = VL->getPos(i)[1];
    refN[skN++] = VL->getPos(i)[2];
    for(int fld=0; fld<bns->isoNfields; fld++)
      refN[skN++] = Q[fld][i];
  }

#if (1)
  // printf("\n\n num isotri before: %8d\n", Ntris);
  // printf(    " num isotri  after: %8d\n", Ntris_Final);
  // printf(    " num nodes  before: %8d\n", Max_N);
  // printf(    " num nodes   after: %8d\n\n", nno);
  printf(  "num nodes  before: %8d  and after: %8d\n", Max_N, nno);
#endif

  delete VL;            // clean up helper object
  return Ntris_Final;   // return num good triangles
}

