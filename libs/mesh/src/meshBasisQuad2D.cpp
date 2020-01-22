/*

The MIT License (MIT)

Copyright (c) 2020 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "mesh.hpp"
#include "mesh2D.hpp"

// ------------------------------------------------------------------------
// QUAD 2D NODES
// ------------------------------------------------------------------------
void mesh_t::NodesQuad2D(int _N, dfloat *_r, dfloat *_s){
  int _Nq = _N+1;

  dfloat *r1D = (dfloat*) malloc(_Nq*sizeof(dfloat));
  JacobiGLL(_N, r1D); //Gauss-Legendre-Lobatto nodes

  //Tensor product
  for (int j=0;j<_Nq;j++) {
    for (int i=0;i<_Nq;i++) {
      _r[i+j*_Nq] = r1D[i];
      _s[i+j*_Nq] = r1D[j];
    }
  }

  free(r1D);
}

void mesh_t::FaceNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_faceNodes){
  int _Nq = _N+1;
  int _Nfp = _Nq;
  int _Np = _Nq*_Nq;

  int cnt[4];
  for (int i=0;i<4;i++) cnt[i]=0;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  for (int n=0;n<_Np;n++) {
    if(fabs(_s[n]+1)<NODETOL)
      _faceNodes[0*_Nfp+(cnt[0]++)] = n;
    if(fabs(_r[n]-1)<NODETOL)
      _faceNodes[1*_Nfp+(cnt[1]++)] = n;
    if(fabs(_s[n]-1)<NODETOL)
      _faceNodes[2*_Nfp+(cnt[2]++)] = n;
    if(fabs(_r[n]+1)<NODETOL)
      _faceNodes[3*_Nfp+(cnt[3]++)] = n;
  }
}

void mesh_t::VertexNodesQuad2D(int _N, dfloat *_r, dfloat *_s, int *_vertexNodes){
  int _Nq = _N+1;
  int _Np = _Nq*_Nq;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  for(int n=0;n<_Np;++n){
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]+1)*(_s[n]+1)<NODETOL)
      _vertexNodes[0] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]+1)*(_s[n]+1)<NODETOL)
      _vertexNodes[1] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]-1)*(_s[n]-1)<NODETOL)
      _vertexNodes[2] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]-1)*(_s[n]-1)<NODETOL)
      _vertexNodes[3] = n;
  }
}

void mesh_t::EquispacedNodesQuad2D(int _N, dfloat *_r, dfloat *_s){
  int _Nq = _N+1;

  //Equispaced 1D nodes
  dfloat *r1D = (dfloat*) malloc(_Nq*sizeof(dfloat));
  EquispacedNodes1D(_N, r1D);

  //Tensor product
  for (int j=0;j<_Nq;j++) {
    for (int i=0;i<_Nq;i++) {
      _r[i+j*_Nq] = r1D[i];
      _s[i+j*_Nq] = r1D[j];
    }
  }

  free(r1D);
}

void mesh_t::EquispacedEToVQuad2D(int _N, int *_EToV){
  int _Nq = _N+1;
  int _Nverts = 3;

  //Tensor product
  int cnt=0;
  for (int j=0;j<_N;j++) {
    for (int i=0;i<_N;i++) {
      _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq;
      _EToV[cnt*_Nverts+1] = i+1+(j  )*_Nq;
      _EToV[cnt*_Nverts+2] = i+1+(j+1)*_Nq;
      cnt++;

      _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq;
      _EToV[cnt*_Nverts+1] = i+1+(j+1)*_Nq;
      _EToV[cnt*_Nverts+2] = i  +(j+1)*_Nq;
      cnt++;
    }
  }
}

void mesh_t::SEMFEMEToVQuad2D(int _N, int *_EToV){
  int _Nq = _N+1;
  int _Nverts = 4;

  //Tensor product
  int cnt=0;
  for (int j=0;j<_N;j++) {
    for (int i=0;i<_N;i++) {
      _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq;
      _EToV[cnt*_Nverts+1] = i+1+(j  )*_Nq;
      _EToV[cnt*_Nverts+2] = i+1+(j+1)*_Nq;
      _EToV[cnt*_Nverts+3] = i  +(j+1)*_Nq;
      cnt++;
    }
  }
}

// ------------------------------------------------------------------------
// ORTHONORMAL BASIS POLYNOMIALS
// ------------------------------------------------------------------------
void mesh_t::OrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *P){
  *P = JacobiP(a,0,0,i)*JacobiP(b,0,0,j);
}

void mesh_t::GradOrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *Pr, dfloat *Ps){
  *Pr = GradJacobiP(a,0,0,i)*JacobiP(b,0,0,j);
  *Ps = JacobiP(a,0,0,i)*GradJacobiP(b,0,0,j);
}

// ------------------------------------------------------------------------
// 2D VANDERMONDE MATRICES
// ------------------------------------------------------------------------

void mesh_t::VandermondeQuad2D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *V){

  int _Nq = _N+1;
  int _Np = _Nq*_Nq;

  for(int n=0; n<Npoints; n++){
    for(int j=0; j<_Nq; j++){
      for(int i=0; i<_Nq; i++){
        int id = n*_Np+i+j*_Nq;
        OrthonormalBasisQuad2D(_r[n], _s[n], i, j, V+id);
      }
    }
  }
}

void mesh_t::GradVandermondeQuad2D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *Vr, dfloat *Vs){

  int _Nq = _N+1;
  int _Np = _Nq*_Nq;

  for(int n=0; n<Npoints; n++){
    for(int j=0; j<_Nq; j++){
      for(int i=0; i<_Nq; i++){
        int id = n*_Np+i+j*_Nq;
        GradOrthonormalBasisQuad2D(_r[n], _s[n], i, j, Vr+id, Vs+id);
      }
    }
  }
}

// ------------------------------------------------------------------------
// 2D OPERATOR MATRICES
// ------------------------------------------------------------------------
void mesh_t::MassMatrixQuad2D(int _Np, dfloat *V, dfloat *_MM){

  // masMatrix = inv(V')*inv(V) = inv(V*V')
  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      dfloat res = 0;
      for(int i=0;i<_Np;++i){
        res += V[n*_Np+i]*V[m*_Np+i];
      }
      _MM[n*_Np + m] = res;
    }
  }
  matrixInverse(_Np, _MM);
}

void mesh_t::DmatrixQuad2D(int _N, int Npoints, dfloat *_r, dfloat *_s,
                                                dfloat *_Dr, dfloat *_Ds){

  int _Np = (_N+1)*(_N+1);

  dfloat *V  = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  dfloat *Vr = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  dfloat *Vs = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));

  VandermondeQuad2D(_N, Npoints, _r, _s, V);
  GradVandermondeQuad2D(_N, Npoints, _r, _s, Vr, Vs);

  //Dr = Vr/V, Ds = Vs/V
  matrixRightSolve(_Np, _Np, Vr, _Np, _Np, V, _Dr);
  matrixRightSolve(_Np, _Np, Vs, _Np, _Np, V, _Ds);

  free(V); free(Vr); free(Vs);
}

void mesh_t::InterpolationMatrixQuad2D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut,
                               dfloat *I){

  int _Np = (_N+1)*(_N+1);

  // need NpointsIn = _Np
  if (NpointsIn != _Np)
    LIBP_ABORT(string("Invalid Interplation operator requested."))

  dfloat *VIn = (dfloat*) malloc(NpointsIn*_Np*sizeof(dfloat));
  dfloat *VOut= (dfloat*) malloc(NpointsOut*_Np*sizeof(dfloat));

  VandermondeQuad2D(_N, NpointsIn,   rIn, sIn, VIn);
  VandermondeQuad2D(_N, NpointsOut, rOut, sOut, VOut);

  matrixRightSolve(NpointsOut, _Np, VOut, NpointsIn, _Np, VIn, I);

  free(VIn); free(VOut);
}