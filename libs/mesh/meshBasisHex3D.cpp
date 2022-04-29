/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim WarburtonTim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

// ------------------------------------------------------------------------
// HEX 3D NODES
// ------------------------------------------------------------------------
void mesh_t::NodesHex3D(const int _N,
                        memory<dfloat>& _r,
                        memory<dfloat>& _s,
                        memory<dfloat>& _t){
  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  memory<dfloat> r1D;
  JacobiGLL(_N, r1D); //Gauss-Legendre-Lobatto nodes

  _r.malloc(_Np);
  _s.malloc(_Np);
  _t.malloc(_Np);

  //Tensor product
  for (int k=0;k<_Nq;k++) {
    for (int j=0;j<_Nq;j++) {
      for (int i=0;i<_Nq;i++) {
        _r[i+j*_Nq+k*_Nq*_Nq] = r1D[i];
        _s[i+j*_Nq+k*_Nq*_Nq] = r1D[j];
        _t[i+j*_Nq+k*_Nq*_Nq] = r1D[k];
      }
    }
  }
}

void mesh_t::FaceNodesHex3D(const int _N,
                            const memory<dfloat> _r,
                            const memory<dfloat> _s,
                            const memory<dfloat> _t,
                            memory<int>& _faceNodes){
  int _Nq = _N+1;
  int _Nfp = _Nq*_Nq;
  int _Np = _Nq*_Nq*_Nq;

  int cnt[6];
  for (int i=0;i<6;i++) cnt[i]=0;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  _faceNodes.malloc(6*_Nfp);
  for (int n=0;n<_Np;n++) {
    if(std::abs(_t[n]+1)<NODETOL)
      _faceNodes[0*_Nfp+(cnt[0]++)] = n;
    if(std::abs(_s[n]+1)<NODETOL)
      _faceNodes[1*_Nfp+(cnt[1]++)] = n;
    if(std::abs(_r[n]-1)<NODETOL)
      _faceNodes[2*_Nfp+(cnt[2]++)] = n;
    if(std::abs(_s[n]-1)<NODETOL)
      _faceNodes[3*_Nfp+(cnt[3]++)] = n;
    if(std::abs(_r[n]+1)<NODETOL)
      _faceNodes[4*_Nfp+(cnt[4]++)] = n;
    if(std::abs(_t[n]-1)<NODETOL)
      _faceNodes[5*_Nfp+(cnt[5]++)] = n;
  }
}

void mesh_t::VertexNodesHex3D(const int _N,
                              const memory<dfloat> _r,
                              const memory<dfloat> _s,
                              const memory<dfloat> _t,
                              memory<int>& _vertexNodes){
  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  _vertexNodes.malloc(8);
  for(int n=0;n<_Np;++n){
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]+1)*(_s[n]+1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      _vertexNodes[0] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]+1)*(_s[n]+1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      _vertexNodes[1] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]-1)*(_s[n]-1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      _vertexNodes[2] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]-1)*(_s[n]-1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      _vertexNodes[3] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]+1)*(_s[n]+1)+(_t[n]-1)*(_t[n]-1)<NODETOL)
      _vertexNodes[4] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]+1)*(_s[n]+1)+(_t[n]-1)*(_t[n]-1)<NODETOL)
      _vertexNodes[5] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]-1)*(_s[n]-1)+(_t[n]-1)*(_t[n]-1)<NODETOL)
      _vertexNodes[6] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]-1)*(_s[n]-1)+(_t[n]-1)*(_t[n]-1)<NODETOL)
      _vertexNodes[7] = n;
  }
}

/*Find a matching array between nodes on matching faces */
void mesh_t::FaceNodeMatchingHex3D(const memory<dfloat> _r,
                                   const memory<dfloat> _s,
                                   const memory<dfloat> _t,
                                   const memory<int> _faceNodes,
                                   const memory<int> _faceVertices,
                                   memory<int>& R){

  const int _Nfaces = 6;
  const int _Nverts = 8;
  const int _NfaceVertices = 4;

  const int _Nfp = _faceNodes.length()/_Nfaces;

  const dfloat NODETOL = 1.0e-5;

  dfloat V0[4][2] = {{-1.0,-1.0},{ 1.0,-1.0},{ 1.0, 1.0},{-1.0, 1.0}};
  dfloat V1[4][2] = {{-1.0,-1.0},{-1.0, 1.0},{ 1.0, 1.0},{ 1.0,-1.0}};

  dfloat EX0[_Nverts], EY0[_Nverts];
  dfloat EX1[_Nverts], EY1[_Nverts];

  memory<dfloat> x0(_Nfp);
  memory<dfloat> y0(_Nfp);

  memory<dfloat> x1(_Nfp);
  memory<dfloat> y1(_Nfp);

  R.malloc(_Nfaces*_Nfaces*_NfaceVertices*_Nfp);

  for (int fM=0;fM<_Nfaces;fM++) {

    for (int v=0;v<_Nverts;v++) {
      EX0[v] = 0.0; EY0[v] = 0.0;
    }
    //setup top element with face fM on the bottom
    for (int v=0;v<_NfaceVertices;v++) {
      int fv = _faceVertices[fM*_NfaceVertices + v];
      EX0[fv] = V0[v][0]; EY0[fv] = V0[v][1];
    }

    for(int n=0;n<_Nfp;++n){ /* for each face node */
      const int fn = _faceNodes[fM*_Nfp+n];

      /* (r,s,t) coordinates of interpolation nodes*/
      dfloat rn = _r[fn];
      dfloat sn = _s[fn];
      dfloat tn = _t[fn];

      /* physical coordinate of interpolation node */
      x0[n] =
        +0.125*(1-rn)*(1-sn)*(1-tn)*EX0[0]
        +0.125*(1+rn)*(1-sn)*(1-tn)*EX0[1]
        +0.125*(1+rn)*(1+sn)*(1-tn)*EX0[2]
        +0.125*(1-rn)*(1+sn)*(1-tn)*EX0[3]
        +0.125*(1-rn)*(1-sn)*(1+tn)*EX0[4]
        +0.125*(1+rn)*(1-sn)*(1+tn)*EX0[5]
        +0.125*(1+rn)*(1+sn)*(1+tn)*EX0[6]
        +0.125*(1-rn)*(1+sn)*(1+tn)*EX0[7];

      y0[n] =
        +0.125*(1-rn)*(1-sn)*(1-tn)*EY0[0]
        +0.125*(1+rn)*(1-sn)*(1-tn)*EY0[1]
        +0.125*(1+rn)*(1+sn)*(1-tn)*EY0[2]
        +0.125*(1-rn)*(1+sn)*(1-tn)*EY0[3]
        +0.125*(1-rn)*(1-sn)*(1+tn)*EY0[4]
        +0.125*(1+rn)*(1-sn)*(1+tn)*EY0[5]
        +0.125*(1+rn)*(1+sn)*(1+tn)*EY0[6]
        +0.125*(1-rn)*(1+sn)*(1+tn)*EY0[7];
    }

    for (int fP=0;fP<_Nfaces;fP++) { /*For each neighbor face */
      for (int rot=0;rot<_NfaceVertices;rot++) { /* For each face rotation */
        // Zero vertices
        for (int v=0;v<_Nverts;v++) {
          EX1[v] = 0.0; EY1[v] = 0.0;
        }
        //setup bottom element with face fP on the top
        for (int v=0;v<_NfaceVertices;v++) {
          int fv = _faceVertices[fP*_NfaceVertices + ((v+rot)%_NfaceVertices)];
          EX1[fv] = V1[v][0]; EY1[fv] = V1[v][1];
        }

        for(int n=0;n<_Nfp;++n){ /* for each node */
          const int fn = _faceNodes[fP*_Nfp+n];

          /* (r,s,t) coordinates of interpolation nodes*/
          dfloat rn = _r[fn];
          dfloat sn = _s[fn];
          dfloat tn = _t[fn];

          /* physical coordinate of interpolation node */
          x1[n] =  0.125*(1-rn)*(1-sn)*(1-tn)*EX1[0]
                  +0.125*(1+rn)*(1-sn)*(1-tn)*EX1[1]
                  +0.125*(1+rn)*(1+sn)*(1-tn)*EX1[2]
                  +0.125*(1-rn)*(1+sn)*(1-tn)*EX1[3]
                  +0.125*(1-rn)*(1-sn)*(1+tn)*EX1[4]
                  +0.125*(1+rn)*(1-sn)*(1+tn)*EX1[5]
                  +0.125*(1+rn)*(1+sn)*(1+tn)*EX1[6]
                  +0.125*(1-rn)*(1+sn)*(1+tn)*EX1[7];

          y1[n] =  0.125*(1-rn)*(1-sn)*(1-tn)*EY1[0]
                  +0.125*(1+rn)*(1-sn)*(1-tn)*EY1[1]
                  +0.125*(1+rn)*(1+sn)*(1-tn)*EY1[2]
                  +0.125*(1-rn)*(1+sn)*(1-tn)*EY1[3]
                  +0.125*(1-rn)*(1-sn)*(1+tn)*EY1[4]
                  +0.125*(1+rn)*(1-sn)*(1+tn)*EY1[5]
                  +0.125*(1+rn)*(1+sn)*(1+tn)*EY1[6]
                  +0.125*(1-rn)*(1+sn)*(1+tn)*EY1[7];
        }

        /* for each node on this face find the neighbor node */
        for(int n=0;n<_Nfp;++n){
          const dfloat xM = x0[n];
          const dfloat yM = y0[n];

          int m=0;
          for(;m<_Nfp;++m){ /* for each neighbor node */
            const dfloat xP = x1[m];
            const dfloat yP = y1[m];

            /* distance between target and neighbor node */
            const dfloat dist = pow(xM-xP,2) + pow(yM-yP,2);

            /* if neighbor node is close to target, match */
            if(dist<NODETOL){
              R[fM*_Nfaces*_NfaceVertices*_Nfp
                + fP*_NfaceVertices*_Nfp
                + rot*_Nfp + n] = m;
              break;
            }
          }

          /*Check*/
          const dfloat xP = x1[m];
          const dfloat yP = y1[m];

          /* distance between target and neighbor node */
          const dfloat dist = pow(xM-xP,2) + pow(yM-yP,2);
          //This shouldn't happen
          LIBP_ABORT("Unable to match face node, face: " << fM
                     << ", matching face: " << fP
                     << ", rotation: " << rot
                     << ", node: " << n
                     << ". Is the reference node set not symmetric?",
                     dist>NODETOL);
        }
      }
    }
  }
}

void mesh_t::EquispacedNodesHex3D(const int _N,
                                  memory<dfloat>& _r,
                                  memory<dfloat>& _s,
                                  memory<dfloat>& _t){
  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  //Equispaced 1D nodes
  memory<dfloat> r1D;
  EquispacedNodes1D(_N, r1D);

  //Tensor product
  _r.malloc(_Np);
  _s.malloc(_Np);
  _t.malloc(_Np);
  for (int k=0;k<_Nq;k++) {
    for (int j=0;j<_Nq;j++) {
      for (int i=0;i<_Nq;i++) {
        _r[i+j*_Nq+k*_Nq*_Nq] = r1D[i];
        _s[i+j*_Nq+k*_Nq*_Nq] = r1D[j];
        _t[i+j*_Nq+k*_Nq*_Nq] = r1D[k];
      }
    }
  }
}

void mesh_t::EquispacedEToVHex3D(const int _N, memory<int>& _EToV){
  const int _Nq = _N+1;
  const int _Nelements = 6*_N*_N*_N;
  const int _Nverts = 4;

  _EToV.malloc(_Nelements*_Nverts);

  //Tensor product
  int cnt=0;
  for (int k=0;k<_N;k++) {
    for (int j=0;j<_N;j++) {
      for (int i=0;i<_N;i++) {
        //tet 1 (0,3,2,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i+1+(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i  +(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;

        //tet 2 (0,1,3,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i+1+(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i+1+(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;

        //tet 3 (0,2,6,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i  +(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i  +(j+1)*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;

        //tet 4 (0,6,4,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i  +(j+1)*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i  +(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;

        //tet 5 (0,5,1,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i+1+(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i+1+(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;

        //tet 6 (0,4,5,7)
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i  +(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i+1+(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;
      }
    }
  }
}

void mesh_t::SEMFEMEToVHex3D(const int _N, memory<int>& _EToV){
  const int _Nq = _N+1;
  const int _Nelements = _N*_N*_N;
  const int _Nverts = 8;

  _EToV.malloc(_Nelements*_Nverts);

  //Tensor product
  int cnt=0;
  for (int k=0;k<_N;k++) {
    for (int j=0;j<_N;j++) {
      for (int i=0;i<_N;i++) {
        _EToV[cnt*_Nverts+0] = i  +(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+1] = i+1+(j  )*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+2] = i+1+(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+3] = i  +(j+1)*_Nq+(k  )*_Nq*_Nq;
        _EToV[cnt*_Nverts+4] = i  +(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+5] = i+1+(j  )*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+6] = i+1+(j+1)*_Nq+(k+1)*_Nq*_Nq;
        _EToV[cnt*_Nverts+7] = i  +(j+1)*_Nq+(k+1)*_Nq*_Nq;
        cnt++;
      }
    }
  }
}

// ------------------------------------------------------------------------
// ORTHONORMAL BASIS POLYNOMIALS
// ------------------------------------------------------------------------
void mesh_t::OrthonormalBasisHex3D(const dfloat a, const dfloat b, const dfloat c,
                                   const int i, const int j, const int k,
                                   dfloat& P){
  P = JacobiP(a,0,0,i)*JacobiP(b,0,0,j)*JacobiP(c,0,0,k);
}

void mesh_t::GradOrthonormalBasisHex3D(const dfloat a, const dfloat b, const dfloat c,
                                       const int i, const int j, const int k,
                                       dfloat& Pr, dfloat& Ps, dfloat& Pt){
  Pr = GradJacobiP(a,0,0,i)*JacobiP(b,0,0,j)*JacobiP(c,0,0,k);
  Ps = JacobiP(a,0,0,i)*GradJacobiP(b,0,0,j)*JacobiP(c,0,0,k);
  Pt = JacobiP(a,0,0,i)*JacobiP(b,0,0,j)*GradJacobiP(c,0,0,k);
}

// ------------------------------------------------------------------------
// 2D VANDERMONDE MATRICES
// ------------------------------------------------------------------------

void mesh_t::VandermondeHex3D(const int _N,
                              const memory<dfloat> _r,
                              const memory<dfloat> _s,
                              const memory<dfloat> _t,
                              memory<dfloat>& V){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;
  const int Npoints = _r.length();

  V.malloc(Npoints*_Np);
  for(int n=0; n<Npoints; n++){
    for(int k=0; k<_Nq; k++){
      for(int j=0; j<_Nq; j++){
        for(int i=0; i<_Nq; i++){
          int id = n*_Np+i+j*_Nq+k*_Nq*_Nq;
          OrthonormalBasisHex3D(_r[n], _s[n], _t[n], i, j, k, V[id]);
        }
      }
    }
  }
}

void mesh_t::GradVandermondeHex3D(const int _N,
                                  const memory<dfloat> _r,
                                  const memory<dfloat> _s,
                                  const memory<dfloat> _t,
                                  memory<dfloat>& Vr,
                                  memory<dfloat>& Vs,
                                  memory<dfloat>& Vt){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;
  const int Npoints = _r.length();

  Vr.malloc(Npoints*_Np);
  Vs.malloc(Npoints*_Np);
  Vt.malloc(Npoints*_Np);
  for(int n=0; n<Npoints; n++){
    for(int k=0; k<_Nq; k++){
      for(int j=0; j<_Nq; j++){
        for(int i=0; i<_Nq; i++){
          int id = n*_Np+i+j*_Nq+k*_Nq*_Nq;
          GradOrthonormalBasisHex3D(_r[n], _s[n], _t[n], i, j, k, Vr[id], Vs[id], Vt[id]);
        }
      }
    }
  }
}

// ------------------------------------------------------------------------
// 2D OPERATOR MATRICES
// ------------------------------------------------------------------------
void mesh_t::MassMatrixHex3D(const int _Np,
                             const memory<dfloat> V,
                             memory<dfloat>& _MM){

  // massMatrix = inv(V')*inv(V) = inv(V*V')
  _MM.malloc(_Np*_Np);
  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      dfloat res = 0;
      for(int i=0;i<_Np;++i){
        res += V[n*_Np+i]*V[m*_Np+i];
      }
      _MM[n*_Np + m] = res;
    }
  }
  linAlg_t::matrixInverse(_Np, _MM);
}

void mesh_t::LumpedMassMatrixHex3D(const int _N,
                                   const memory<dfloat> _gllw,
                                   memory<dfloat>& _MM){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  // LumpedMassMatrix = gllw \ctimes gllw \ctimes gllw
  _MM.malloc(_Np*_Np, 0.0);
  for(int k=0;k<_Nq;++k){
    for(int n=0;n<_Nq;++n){
      for(int m=0;m<_Nq;++m){
        int id = n+m*_Nq+k*_Nq*_Nq;
        _MM[id+id*_Np] = _gllw[n]*_gllw[m]*_gllw[k];
      }
    }
  }
}

void mesh_t::invLumpedMassMatrixHex3D(const int _N,
                                      const memory<dfloat> _gllw,
                                      memory<dfloat>& _invMM){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  // invLumpedMassMatrix = invgllw \ctimes invgllw
  _invMM.malloc(_Np*_Np, 0.0);
  for(int k=0;k<_Nq;++k){
    for(int n=0;n<_Nq;++n){
      for(int m=0;m<_Nq;++m){
        int id = n+m*_Nq+k*_Nq*_Nq;
        _invMM[id+id*_Np] = 1.0/(_gllw[n]*_gllw[m]*_gllw[k]);
      }
    }
  }
}

void mesh_t::DmatrixHex3D(const int _N,
                          const memory<dfloat> _r,
                          const memory<dfloat> _s,
                          const memory<dfloat> _t,
                          memory<dfloat>& _D){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  memory<dfloat> V, Vr, Vs, Vt;
  VandermondeHex3D(_N, _r, _s, _t, V);
  GradVandermondeHex3D(_N, _r, _s, _t, Vr, Vs, Vt);

  //Dr = Vr/V, Ds = Vs/V, Dt = Vt/V
  _D.malloc(3*_Np*_Np);
  memory<dfloat> _Dr = _D + 0*_Np*_Np;
  memory<dfloat> _Ds = _D + 1*_Np*_Np;
  memory<dfloat> _Dt = _D + 2*_Np*_Np;
  linAlg_t::matrixRightSolve(_Np, _Np, Vr, _Np, _Np, V, _Dr);
  linAlg_t::matrixRightSolve(_Np, _Np, Vs, _Np, _Np, V, _Ds);
  linAlg_t::matrixRightSolve(_Np, _Np, Vt, _Np, _Np, V, _Dt);
}

void mesh_t::InterpolationMatrixHex3D(const int _N,
                                      const memory<dfloat> rIn,
                                      const memory<dfloat> sIn,
                                      const memory<dfloat> tIn,
                                      const memory<dfloat> rOut,
                                      const memory<dfloat> sOut,
                                      const memory<dfloat> tOut,
                                      memory<dfloat>& I){

  const int _Nq = _N+1;
  const int _Np = _Nq*_Nq*_Nq;

  const int NpointsIn  = rIn.length();
  const int NpointsOut = rOut.length();

  // need NpointsIn = _Np
  LIBP_ABORT("Invalid Interplation operator requested.",
             NpointsIn != _Np);

  memory<dfloat> VIn;
  memory<dfloat> VOut;
  VandermondeHex3D(_N, rIn, sIn, tIn, VIn);
  VandermondeHex3D(_N, rOut, sOut, tOut, VOut);

  I.malloc(NpointsIn*NpointsOut);
  linAlg_t::matrixRightSolve(NpointsOut, _Np, VOut,
                             NpointsIn, _Np, VIn, I);
}

} //namespace libp
