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
#include "mesh3D.hpp"

// ------------------------------------------------------------------------
// TET 3D NODES
// ------------------------------------------------------------------------
void mesh_t::NodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t){

  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  EquispacedNodesTet3D(_N, _r, _s, _t); //make equispaced nodes on reference tet
  WarpBlendTransformTet3D(_N, _Np, _r, _s, _t); //apply warp&blend transform
}

void mesh_t::FaceNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes){
  int _Nfp = (_N+1)*(_N+2)/2;
  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  int cnt[4];
  for (int i=0;i<4;i++) cnt[i]=0;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  for (int n=0;n<_Np;n++) {
    if(fabs(_t[n]+1)<NODETOL)
      _faceNodes[0*_Nfp+(cnt[0]++)] = n;
    if(fabs(_s[n]+1)<NODETOL)
      _faceNodes[1*_Nfp+(cnt[1]++)] = n;
    if(fabs(_r[n]+_s[n]+_t[n]+1.0)<NODETOL)
      _faceNodes[2*_Nfp+(cnt[2]++)] = n;
    if(fabs(_r[n]+1)<NODETOL)
      _faceNodes[3*_Nfp+(cnt[3]++)] = n;
  }
}

void mesh_t::VertexNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t, int *_vertexNodes){
  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  dfloat deps = 1.;
  while((1.+deps)>1.)
    deps *= 0.5;

  const dfloat NODETOL = 1000.*deps;

  for(int n=0;n<_Np;++n){
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]+1)*(_s[n]+1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      vertexNodes[0] = n;
    if( (_r[n]-1)*(_r[n]-1)+(_s[n]+1)*(_s[n]+1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      vertexNodes[1] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]-1)*(_s[n]-1)+(_t[n]+1)*(_t[n]+1)<NODETOL)
      vertexNodes[2] = n;
    if( (_r[n]+1)*(_r[n]+1)+(_s[n]+1)*(_s[n]+1)+(_t[n]-1)*(_t[n]-1)<NODETOL)
      vertexNodes[3] = n;
  }
}

// Create equidistributed nodes on reference tet
void mesh_t::EquispacedNodesTet3D(int _N, dfloat *_r, dfloat *_s, dfloat *_t){

  int sk = 0;
  for (int k=0;k<_N+1;k++) {
    for (int n=0;n<_N+1-k;n++) {
      for (int m=0;m<_N+1-n-k;m++) {
        _r[sk] = -1.0 + 2.0*m/_N;
        _s[sk] = -1.0 + 2.0*n/_N;
        _t[sk] = -1.0 + 2.0*k/_N;
        sk++;
      }
    }
  }
}

void mesh_t::EquispacedEToVTet3D(int _N, int *_EToV){
  int _Nverts = 4;

  int cnt=0;
  int sk=0;
  for (int k=0;k<_N;k++) {
    int layershift = (_N+1-k)*(_N+2-k)/2; //number of nodes in this layer
    for (int j=0;j<_N-k;j++) {

      int shift = _N+1-j-k; //number of nodes in this row

      for (int i=0;i<_N-j-k;i++) {
        _EToV[cnt*_Nverts+0] = sk  ;
        _EToV[cnt*_Nverts+1] = sk+1;
        _EToV[cnt*_Nverts+2] = sk+shift;
        _EToV[cnt*_Nverts+3] = sk+layershift;
        cnt++;

        //middle octohedron becomes 4 tets
        if (i<_N-j-k-1) {
          _EToV[cnt*_Nverts+0] = sk+1;
          _EToV[cnt*_Nverts+1] = sk+shift;
          _EToV[cnt*_Nverts+2] = sk+layershift;
          _EToV[cnt*_Nverts+3] = sk+layershift+1;
          cnt++;

          _EToV[cnt*_Nverts+0] = sk+1;
          _EToV[cnt*_Nverts+1] = sk+shift+1;
          _EToV[cnt*_Nverts+2] = sk+shift;
          _EToV[cnt*_Nverts+3] = sk+layershift+1;
          cnt++;

          _EToV[cnt*_Nverts+0] = sk+shift;
          _EToV[cnt*_Nverts+1] = sk+shift+1;
          _EToV[cnt*_Nverts+2] = sk+shift+layershift-1;
          _EToV[cnt*_Nverts+3] = sk+layershift+1;
          cnt++;

          _EToV[cnt*_Nverts+0] = sk+shift;
          _EToV[cnt*_Nverts+1] = sk+shift+layershift-1;
          _EToV[cnt*_Nverts+2] = sk+layershift;
          _EToV[cnt*_Nverts+3] = sk+layershift+1;
          cnt++;
        }
        //top corner tet
        if (i<_N-j-k-2) {
          _EToV[cnt*_Nverts+0] = sk+shift+1;
          _EToV[cnt*_Nverts+1] = sk+layershift+shift;
          _EToV[cnt*_Nverts+2] = sk+layershift+shift-1;
          _EToV[cnt*_Nverts+3] = sk+layershift+1;
          cnt++;
        }
        sk++;
      }
      sk++;
      layershift--;
    }
    sk++;
  }
}

void mesh_t::SEMFEMEToVTet3D(int _N, int *_EToV){
  EquispacedEToVTet3D(_N, _EToV);
}

// ------------------------------------------------------------------------
// ORTHONORMAL BASIS POLYNOMIALS
// ------------------------------------------------------------------------
void mesh_t::OrthonormalBasisTet3D(dfloat _r, dfloat _s, dfloat _t, int i, int j, int k, dfloat *P){
  // First convert to abc coordinates
  dfloat a, b, c;
  if(fabs(_s+_t)>1e-8)
    a = 2.0*(1.+_r)/(-_s-_t)-1.0;
  else
    a = -1.0;

  if(fabs(_t-1)>1e-8)
    b = 2.0*(1+_s)/(1.-_t)-1.0;
  else
    b = -1.;

  c = _t;

  dfloat p1 = JacobiP(a,0,0,i);
  dfloat p2 = JacobiP(b,2*i+1,0,j);
  dfloat p3 = JacobiP(c,2*(i+j)+2,0,k);

  *P = 2.*sqrt(2.0)*p1*p2*p3*pow(1.0-b,i)*pow(1.0-c,i+j);
}

void mesh_t::GradOrthonormalBasisTet3D(dfloat _r, dfloat _s, dfloat _t,
                                       int i, int j, int k, dfloat *Pr, dfloat *Ps, dfloat *Pt){
  // First convert to abc coordinates
  dfloat a, b, c;
  if(fabs(_s+_t)>1e-8)
    a = 2.0*(1.+_r)/(-_s-_t)-1.0;
  else
    a = -1.0;

  if(fabs(_t-1)>1e-8)
    b = 2.0*(1+_s)/(1.-_t)-1.0;
  else
    b = -1.;

  c = _t;

  dfloat p1 = JacobiP(a,0,0,i);
  dfloat p2 = JacobiP(b,2*i+1,0,j);
  dfloat p3 = JacobiP(c,2*(i+j)+2,0,k);

  dfloat p1a = GradJacobiP(a,0,0,i);
  dfloat p2b = GradJacobiP(b,2*i+1,0,j);
  dfloat p3c = GradJacobiP(c,2*(i+j)+2,0,k);

  *Pr = p1a*p2*p3;
  if(i>0)
    *Pr *= pow(0.5*(1.0-b), i-1);
  if(i+j>0)
    *Pr *= pow(0.5*(1.0-c), i+j-1);

  *Ps = 0.5*(1.0+a)*(*Pr);
  dfloat tmp = p2b*pow(0.5*(1.0-b), i);
  if(i>0)
    tmp += -0.5*i*p2*pow(0.5*(1.0-b), i-1);
  if(i+j>0)
    tmp *= pow(0.5*(1.0-c), i+j-1);
  tmp *= p1*p3;
  *Ps += tmp;

  *Pt = 0.5*(1.0+a)*(*Pr) + 0.5*(1.0+b)*tmp;
  tmp = p3c*pow(0.5*(1-c), i+j);
  if(i+j>0)
    tmp -= 0.5*(i+j)*(p3*pow(0.5*(1.0-c), i+j-1));
  tmp *= p1*p2*pow(0.5*(1-b), i);
  *Pt += tmp;

  *Pr *= pow(2, 2*i+j+1.5);
  *Ps *= pow(2, 2*i+j+1.5);
  *Pt *= pow(2, 2*i+j+1.5);
}

// ------------------------------------------------------------------------
// 3D VANDERMONDE MATRICES
// ------------------------------------------------------------------------

void mesh_t::VandermondeTet3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t, dfloat *V){

  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  for(int n=0; n<Npoints; n++){
    int sk=0;
    for(int i=0; i<_N+1; i++){
      for(int j=0; j<_N+1-i; j++){
        for(int k=0; k<_N+1-i-j; k++){
          int id = n*_Np+sk;
          OrthonormalBasisTet3D(_r[n], _s[n], _t[n], i, j, k, V+id);
          sk++;
        }
      }
    }
  }
}

void mesh_t::GradVandermondeTet3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                                  dfloat *Vr, dfloat *Vs, dfloat *Vt){

  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  for(int n=0; n<Npoints; n++){
    int sk=0;
    for(int i=0; i<_N+1; i++){
      for(int j=0; j<_N+1-i; j++){
        for(int k=0; k<_N+1-i-j; k++){
          int id = n*_Np+sk;
          GradOrthonormalBasisTet3D(_r[n], _s[n], _t[n], i, j, k, Vr+id, Vs+id, Vt+id);
          sk++;
        }
      }
    }
  }
}

// ------------------------------------------------------------------------
// 3D OPERATOR MATRICES
// ------------------------------------------------------------------------
void mesh_t::MassMatrixTet3D(int _Np, dfloat *V, dfloat *_MM){

  // massMatrix = inv(V')*inv(V) = inv(V*V')
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

void mesh_t::DmatrixTet3D(int _N, int Npoints, dfloat *_r, dfloat *_s, dfloat *_t,
                                               dfloat *_Dr, dfloat *_Ds, dfloat *_Dt){

  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  dfloat *V  = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  dfloat *Vr = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  dfloat *Vs = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  dfloat *Vt = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));

  VandermondeTet3D(_N, Npoints, _r, _s, _t, V);
  GradVandermondeTet3D(_N, Npoints, _r, _s, _t, Vr, Vs, Vt);

  //Dr = Vr/V, Ds = Vs/V
  matrixRightSolve(_Np, _Np, Vr, _Np, _Np, V, _Dr);
  matrixRightSolve(_Np, _Np, Vs, _Np, _Np, V, _Ds);
  matrixRightSolve(_Np, _Np, Vt, _Np, _Np, V, _Dt);

  free(V); free(Vr); free(Vs); free(Vt);
}

void mesh_t::LIFTmatrixTet3D(int _N, int *_faceNodes,
                             dfloat *_r, dfloat *_s, dfloat *_t, dfloat *_LIFT){

  int _Nfp = (_N+1)*(_N+2)/2;
  int _Np = (_N+1)*(_N+2)*(_N+3)/6;
  int _Nfaces = 4;

  dfloat *E = (dfloat *) calloc(_Np*_Nfaces*_Nfp, sizeof(dfloat));

  dfloat *r2D = (dfloat *) malloc(_Nfp*sizeof(dfloat));
  dfloat *s2D = (dfloat *) malloc(_Nfp*sizeof(dfloat));
  dfloat *V2D = (dfloat *) malloc(_Nfp*_Nfp*sizeof(dfloat));
  dfloat *MM2D = (dfloat *) malloc(_Nfp*_Nfp*sizeof(dfloat));

  for (int f=0;f<_Nfaces;f++) {
    dfloat *rFace, *sFace;
    if (f==0) {rFace = _r; sFace = _s;}
    if (f==1) {rFace = _r; sFace = _t;}
    if (f==2) {rFace = _s; sFace = _t;}
    if (f==3) {rFace = _s; sFace = _t;}

    for (int i=0;i<_Nfp;i++) {
      r2D[i] = rFace[_faceNodes[f*_Nfp+i]];
      s2D[i] = sFace[_faceNodes[f*_Nfp+i]];
    }

    VandermondeTri2D(_N, _Nfp, r2D, s2D, V2D);
    MassMatrixTri2D(_Nfp, V2D, MM2D);

    for (int j=0;j<_Nfp;j++) {
      int fid = _faceNodes[f*_Nfp+j];
      for (int i=0;i<_Nfp;i++) {
        E[fid*_Nfaces*_Nfp + i + f*_Nfp] = MM2D[j*_Nfp+i];
      }
    }
  }

  dfloat *V = (dfloat *) malloc(_Np*_Np*sizeof(dfloat));
  VandermondeTet3D(_N, _Np, _r, _s, _t, V);

  for (int n=0;n<_Np;n++) {
    for (int m=0;m<_Nfaces*_Nfp;m++) {

      _LIFT[n*_Nfaces*_Nfp+m]=0.0;

      // LIFT = V*(V'*E);
      for (int j=0;j<_Np;j++) {
        for (int i=0;i<_Np;i++) {
          _LIFT[n*_Nfaces*_Nfp+m] += V[n*_Np+i]*V[j*_Np+i]*E[j*_Nfaces*_Nfp+m];
        }
      }
    }
  }

  free(V); free(r2D); free(s2D); free(V2D); free(MM2D); free(E);
}

void mesh_t::InterpolationMatrixTet3D(int _N,
                               int NpointsIn, dfloat *rIn, dfloat *sIn, dfloat *tIn,
                               int NpointsOut, dfloat *rOut, dfloat *sOut, dfloat *tOut,
                               dfloat *I){

  int _Np = (_N+1)*(_N+2)*(_N+3)/6;

  // need NpointsIn = _Np
  if (NpointsIn != _Np)
    LIBP_ABORT(string("Invalid Interplation operator requested."))

  dfloat *VIn = (dfloat*) malloc(NpointsIn*_Np*sizeof(dfloat));
  dfloat *VOut= (dfloat*) malloc(NpointsOut*_Np*sizeof(dfloat));

  VandermondeTet3D(_N, NpointsIn,   rIn, sIn, tIn, VIn);
  VandermondeTet3D(_N, NpointsOut, rOut, sOut, tOut, VOut);

  matrixRightSolve(NpointsOut, _Np, VOut, NpointsIn, _Np, VIn, I);

  free(VIn); free(VOut);
}

void mesh_t::DegreeRaiseMatrixTet3D(int Nc, int Nf, dfloat *P){

  int Npc = (Nc+1)*(Nc+2)*(Nc+3)/6;
  int Npf = (Nf+1)*(Nf+2)*(Nf+3)/6;

  dfloat *rc = (dfloat *) malloc(Npc*sizeof(dfloat));
  dfloat *sc = (dfloat *) malloc(Npc*sizeof(dfloat));
  dfloat *tc = (dfloat *) malloc(Npc*sizeof(dfloat));
  dfloat *rf = (dfloat *) malloc(Npf*sizeof(dfloat));
  dfloat *sf = (dfloat *) malloc(Npf*sizeof(dfloat));
  dfloat *tf = (dfloat *) malloc(Npf*sizeof(dfloat));

  NodesTet3D(Nc, rc, sc, tc);
  NodesTet3D(Nf, rf, sf, tf);

  InterpolationMatrixTet3D(Nc, Npc, rc, sc, tc, Npf, rf, sf, tf, P);

  free(rc); free(sc); free(tc); free(rf); free(sf); free(tf);
}

void mesh_t::CubatureWeakDmatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                    int _cubNp, dfloat *_cubr, dfloat *_cubs, dfloat *_cubt, dfloat *_cubw,
                                    dfloat *_cubDrW, dfloat *_cubDsW, dfloat *_cubDtW, dfloat *_cubProject){

  dfloat *V = (dfloat*) malloc(_Np*_Np*sizeof(dfloat));
  VandermondeTet3D(_N, _Np, _r, _s, _t, V);

  dfloat *cubV  = (dfloat*) malloc(_cubNp*_Np*sizeof(dfloat));
  dfloat *cubVr = (dfloat*) malloc(_cubNp*_Np*sizeof(dfloat));
  dfloat *cubVs = (dfloat*) malloc(_cubNp*_Np*sizeof(dfloat));
  dfloat *cubVt = (dfloat*) malloc(_cubNp*_Np*sizeof(dfloat));
  VandermondeTet3D(N, _cubNp, _cubr, _cubs, _cubt, cubV);
  GradVandermondeTet3D(_N, _cubNp, _cubr, _cubs, _cubt, cubVr, cubVs, cubVt);

  // cubDrW = V*transpose(cVr)*diag(cubw);
  // cubDsW = V*transpose(cVs)*diag(cubw);
  // cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)
  for(int n=0;n<_cubNp;++n){
    for(int m=0;m<_Np;++m){
      // scale by cubw
      cubVr[n*Np+m] *= _cubw[n];
      cubVs[n*Np+m] *= _cubw[n];
      cubVt[n*Np+m] *= _cubw[n];
      cubV[n*Np+m]  *= _cubw[n];
    }
  }

  for(int n=0;n<_Np;++n){
    for(int m=0;m<_cubNp;++m){
      dfloat resP = 0, resDrW = 0, resDsW = 0, resDtW = 0;

      for(int i=0;i<Np;++i){
        dfloat Vni = V[n*Np+i];
        resDrW += Vni*cubVr[m*Np+i];
        resDsW += Vni*cubVs[m*Np+i];
        resDtW += Vni*cubVt[m*Np+i];
        resP   += Vni*cubV[m*Np+i];
      }

      cubDrW[n*cubNp+m] = resDrW;
      cubDsW[n*cubNp+m] = resDsW;
      cubDtW[n*cubNp+m] = resDtW;
      cubProject[n*cubNp+m] = resP;
    }
  }

  free(V); free(cubV); free(cubVr); free(cubVs); free(cubVt);
}

void mesh_t::CubatureSurfaceMatricesTet3D(int _N, int _Np, dfloat *_r, dfloat *_s, dfloat *_t, int *_faceNodes,
                                    int _intNfp, dfloat *_intr, dfloat *_ints, dfloat *_intw,
                                    dfloat *_intInterp, dfloat *_intLIFT){

  int _Nfaces = 4;
  int _Nfp = (_N+1)*(_N+2)/2;

  dfloat *V = (dfloat*) malloc(_Np*_Np*sizeof(dfloat));
  VandermondeTet3D(_N, _Np, _r, _s, _t, V);

  dfloat *ir = (dfloat*) calloc(_intNfp*_Nfaces, sizeof(dfloat));
  dfloat *is = (dfloat*) calloc(_intNfp*_Nfaces, sizeof(dfloat));
  dfloat *it = (dfloat*) calloc(_intNfp*_Nfaces, sizeof(dfloat));
  dfloat *iw = (dfloat*) calloc(_intNfp*_Nfaces, sizeof(dfloat));

  for(int n=0;n<_intNfp;++n){
    ir[0*_intNfp + n] =  _intr[n];
    ir[1*_intNfp + n] =  _intr[n];
    ir[2*_intNfp + n] =  _intr[n];
    ir[3*_intNfp + n] = -1.0;

    is[0*_intNfp + n] =  _ints[n];
    is[1*_intNfp + n] = -1.0;
    is[2*_intNfp + n] =  _ints[n];
    is[3*_intNfp + n] =  _intr[n];

    it[0*_intNfp + n] = -1.0;
    it[1*_intNfp + n] =  _ints[n];
    it[2*_intNfp + n] = -(1.0 + _intr[n] + _ints[n]);
    it[3*_intNfp + n] =  _ints[n];

    iw[0*_intNfp + n] =  _intw[n];
    iw[1*_intNfp + n] =  _intw[n];
    iw[2*_intNfp + n] =  _intw[n];
    iw[3*_intNfp + n] =  _intw[n];
  }

  dfloat *sInterp = (dfloat*) malloc(_intNfp*_Nfaces*_Np*sizeof(dfloat));
  InterpolationMatrixTet3D(_N, _Np, _r, _s, _t, _Nfaces*_intNfp, ir, is, it, sInterp);

  for(int n=0;n<_intNfp;++n){
    for(int m=0;m<_Nfp;++m){
      intInterp[0*_intNfp*_Nfp + n*_Nfp + m] = sInterp[(n+0*_intNfp)*_Np+_faceNodes[0*_Nfp+m]];
      intInterp[1*_intNfp*_Nfp + n*_Nfp + m] = sInterp[(n+1*_intNfp)*_Np+_faceNodes[1*_Nfp+m]];
      intInterp[2*_intNfp*_Nfp + n*_Nfp + m] = sInterp[(n+2*_intNfp)*_Np+_faceNodes[2*_Nfp+m]];
      intInterp[3*_intNfp*_Nfp + n*_Nfp + m] = sInterp[(n+3*_intNfp)*_Np+_faceNodes[3*_Nfp+m]];
    }
  }

  // integration node lift matrix
  //iLIFT = V*V'*sInterp'*diag(iw(:));
  for(int n=0;n<_Nfaces*_intNfp;++n){
    for(int m=0;m<_Np;++m){
      intLIFT[m*_Nfaces*_intNfp+n] = 0.0;
      for(int j=0;j<_Np;++j){
        for(int i=0;i<_Np;++i){
          intLIFT[m*_Nfaces*_intNfp+n] += V[m*Np+i]*V[j*_Np+i]*sInterp[n*_Np+j]*iw[n];
        }
      }
    }
  }

  free(V); free(ir);  free(is); free(it); free(iw);  free(sInterp);
}

void mesh_t::SEMFEMInterpMatrixTet3D(int _N,
                                    int _Np, dfloat *_r, dfloat *_s, dfloat *_t,
                                    int _NpFEM, dfloat *_rFEM, dfloat *_sFEM, dfloat *_tFEM,
                                    dfloat *I){

  dfloat *IQN = (dfloat*) malloc(_NpFEM*_Np*sizeof(dfloat));
  InterpolationMatrixTet3D(_N, _Np, _r, _s, _t, _NpFEM, _rFEM, _sFEM, _tFEM, IQN);

  dfloat *IQTIQ = (dfloat*) malloc(_Np*_Np*sizeof(dfloat));
  // IQTIQ = IQN'*IQN
  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      IQTIQ[n*_Np+m] = 0.0;
      for(int i=0;i<_NpFEM;++i){
        IQTIQ[n*_Np+m] += IQN[i*_Np+n]*IQN[i*_Np+m];
      }
    }
  }

  // I = IQN/(IQN'*IQN)  - pseudo inverse
  matrixRightSolve(_NpFEM, _Np, IQN, _Np, _Np, IQTIQ, I);

  free(IQN); free(IQTIQ);
}

// ------------------------------------------------------------------------
// Warp & Blend routines
//  Warburton, T. (2006). An explicit construction of interpolation nodes on the simplex.
//                       Journal of engineering mathematics, 56(3), 247-262.
// ------------------------------------------------------------------------

static void xyztorst(int Npoints, dfloat *x, dfloat *y, dfloat *z, dfloat *r, dfloat *s, dfloat *t) {
  // vertices of tetrahedron
  dfloat v1[3] = {-1.0, -1./sqrt(3.), -1./sqrt(6.)};
  dfloat v2[3] = { 1.0, -1./sqrt(3.), -1./sqrt(6.)};
  dfloat v3[3] = { 0.0,  2./sqrt(3.), -1./sqrt(6.)};
  dfloat v4[3] = { 0.0,  0.,           3./sqrt(6.)};

  dfloat *XYZ = (dfloat *) malloc(3*Npoints*sizeof(dfloat));
  dfloat *RST = (dfloat *) malloc(3*Npoints*sizeof(dfloat));
  dfloat *A = (dfloat *) malloc(3*3*sizeof(dfloat));

  for (int i=0;i<3;i++) {
    A[0*3+i] = 0.5*(v2[i]-v1[i]);
    A[1*3+i] = 0.5*(v3[i]-v1[i]);
    A[2*3+i] = 0.5*(v4[i]-v1[i]);
  }

  for (int n=0;n<Npoints;n++) {
    XYZ[3*n+0] = x[n]-0.5*(v2[0]+v3[0]+v4[0]-v1[0]);
    XYZ[3*n+1] = y[n]-0.5*(v2[1]+v3[1]+v4[1]-v1[1]);
    XYZ[3*n+2] = z[n]-0.5*(v2[2]+v3[2]+v4[2]-v1[2]);
  }

  matrixRightSolve(Npoints, 3, XYZ, 3, 3, A, RST);

  for (int n=0;n<Npoints;n++) {
    r[n] = RST[3*n+0];
    s[n] = RST[3*n+1];
    t[n] = RST[3*n+2];
  }

  free(XYZ); free(RST); free(A);
}

void mesh_t::WarpShiftFace3D(int _N, int Npoints, dfloat alpha,
                             dfloat *L1, dfloat *L2, dfloat *L3,
                             dfloat *w1, dfloat *w2) {
  // Compute scaled warp function at order N
  // based on rout interpolation nodes

  dfloat *dL32 = (dfloat*) malloc(Npoints*sizeof(dfloat));
  dfloat *dL13 = (dfloat*) malloc(Npoints*sizeof(dfloat));
  dfloat *dL21 = (dfloat*) malloc(Npoints*sizeof(dfloat));

  dfloat *warpf1 = (dfloat*) malloc(Npoints*sizeof(dfloat));
  dfloat *warpf2 = (dfloat*) malloc(Npoints*sizeof(dfloat));
  dfloat *warpf3 = (dfloat*) malloc(Npoints*sizeof(dfloat));

  for (int n=0;n<Npoints;n++) {
    dL32[n] = L3[n]-L2[n];
    dL13[n] = L1[n]-L3[n];
    dL21[n] = L2[n]-L1[n];
  }

  Warpfactor(_N, Npoints, dL32, warpf1);
  Warpfactor(_N, Npoints, dL13, warpf2);
  Warpfactor(_N, Npoints, dL21, warpf3);

  for (int n=0;n<Npoints;n++) {
    dfloat blend1 = 4.0*L2[n]*L3[n];
    dfloat blend2 = 4.0*L3[n]*L1[n];
    dfloat blend3 = 4.0*L1[n]*L2[n];

    dfloat warp1 = blend1*warpf1[n]*(1.0+alpha*alpha*L1[n]*L1[n]);
    dfloat warp2 = blend2*warpf2[n]*(1.0+alpha*alpha*L2[n]*L2[n]);
    dfloat warp3 = blend3*warpf3[n]*(1.0+alpha*alpha*L3[n]*L3[n]);

    w1[n] = 1.*warp1 + cos(2.*M_PI/3.)*warp2 + cos(4.*M_PI/3.)*warp3;
    w2[n] = 0.*warp1 + sin(2.*M_PI/3.)*warp2 + sin(4.*M_PI/3.)*warp3;
  }

  free(dL32); free(dL21); free(dL13);
  free(warpf1); free(warpf2); free(warpf3);
}

void mesh_t::WarpBlendTransformTet3D(int _N, int _Npoints, dfloat *_r, dfloat *_s, dfloat *_t, dfloat alphaIn){

  const dfloat alpopt[15] = {0.0000,0.0000,0.00000,0.1002,1.1332,1.5608,1.3413,
                             1.2577,1.1603,1.10153,0.6080,0.4523,0.8856,0.8717,0.9655};

  dfloat alpha;
  if (alphaIn==-1) {
    if (_N<16) {
      alpha = alpopt[_N-1];
    } else {
      alpha = 1.;
    }
  } else {
    alpha = alphaIn;
  }

  // vertices of tetrahedron
  dfloat v1[3] = {-1.0, -1./sqrt(3.), -1./sqrt(6.)};
  dfloat v2[3] = { 1.0, -1./sqrt(3.), -1./sqrt(6.)};
  dfloat v3[3] = { 0.0,  2./sqrt(3.), -1./sqrt(6.)};
  dfloat v4[3] = { 0.0,  0.,           3./sqrt(6.)};

  // orthogonal axis tangents on faces 1-4
  dfloat t1[4][4], t2[4][4];
  for (int v=0;v<3;v++) {
    t1[0][v] = v2[v]-v1[v];             t1[1][v] = v2[v]-v1[v];
    t1[2][v] = v3[v]-v2[v];             t1[3][v] = v3[v]-v1[v];
    t2[0][v] = v3[v]-0.5*(v1[v]+v2[v]); t2[1][v] = v4[v]-0.5*(v1[v]+v2[v]);
    t2[2][v] = v4[v]-0.5*(v2[v]+v3[v]); t2[3][v] = v4[v]-0.5*(v1[v]+v3[v]);
  }
  // normalize tangents
  for (int v=0;v<4;v++) {
    dfloat normt1 = sqrt(t1[v][0]*t1[v][0]+t1[v][1]*t1[v][1]+t1[v][2]*t1[v][2]);
    dfloat normt2 = sqrt(t2[v][0]*t2[v][0]+t2[v][1]*t2[v][1]+t2[v][2]*t2[v][2]);
    for (int i=0;i<4;i++) {
      t1[v][i] /= normt1;
      t2[v][i] /= normt2;
    }
  }

  // Convert r s coordinates to points in equilateral triangle
  dfloat *L1 = (dfloat*) malloc(_Npoints*sizeof(dfloat));
  dfloat *L2 = (dfloat*) malloc(_Npoints*sizeof(dfloat));
  dfloat *L3 = (dfloat*) malloc(_Npoints*sizeof(dfloat));
  dfloat *L4 = (dfloat*) malloc(_Npoints*sizeof(dfloat));

  dfloat *_x = (dfloat*) malloc(_Npoints*sizeof(dfloat));
  dfloat *_y = (dfloat*) malloc(_Npoints*sizeof(dfloat));
  dfloat *_z = (dfloat*) malloc(_Npoints*sizeof(dfloat));

  dfloat *shiftx = (dfloat*) calloc(_Npoints,sizeof(dfloat));
  dfloat *shifty = (dfloat*) calloc(_Npoints,sizeof(dfloat));
  dfloat *shiftz = (dfloat*) calloc(_Npoints,sizeof(dfloat));

  for (int n=0;n<_Npoints;n++) {
    L1[n] =  0.5*(1.+_t[n]);
    L2[n] =  0.5*(1.+_s[n]);
    L3[n] = -0.5*(1.0+_r[n]+_s[n]+_t[n]);
    L4[n] =  0.5*(1.+_r[n]);

    _x[n] =  L3[n]*v1[0]+L4[n]*v2[0]+L2[n]*v3[0]+L1[n]*v4[0];
    _y[n] =  L3[n]*v1[1]+L4[n]*v2[1]+L2[n]*v3[1]+L1[n]*v4[1];
    _z[n] =  L3[n]*v1[2]+L4[n]*v2[2]+L2[n]*v3[2]+L1[n]*v4[2];
  }

  dfloat *warp1 = (dfloat*) calloc(_Npoints,sizeof(dfloat));
  dfloat *warp2 = (dfloat*) calloc(_Npoints,sizeof(dfloat));

  for (int f=0;f<4;f++) {
    dfloat *La, *Lb, *Lc, *Ld;
    if(f==0) {La = L1; Lb = L2; Lc = L3; Ld = L4;}
    if(f==1) {La = L2; Lb = L1; Lc = L3; Ld = L4;}
    if(f==2) {La = L3; Lb = L1; Lc = L4; Ld = L2;}
    if(f==3) {La = L4; Lb = L1; Lc = L3; Ld = L2;}

    // compute warp tangential to face
    WarpShiftFace3D(N, _Npoints, alpha, Lb, Lc, Ld, warp1, warp2);

    for (int n=0;n<_Npoints;n++) {
      dfloat blend = Lb[n]*Lc[n]*Ld[n];

      // modify linear blend
      dfloat denom = (Lb[n]+.5*La[n])*(Lc[n]+.5*La[n])*(Ld[n]+.5*La[n]);
      if (denom>1.e-10) {
        blend = (1+(alpha*La[n])*(alpha*La[n]))*blend/denom;
      }

      // compute warp & blend
      shiftx[n] += (blend*warp1[n])*t1[f][0] + (blend*warp2[n])*t2[f][0];
      shifty[n] += (blend*warp1[n])*t1[f][1] + (blend*warp2[n])*t2[f][1];
      shiftz[n] += (blend*warp1[n])*t1[f][2] + (blend*warp2[n])*t2[f][2];

      // fix face warp
      if ((La[n]<1.e-10) && ((Lb[n]<1.e-10)||(Lc[n]<1.e-10)||(Ld[n]<1.e-10))) {
        shiftx[n] = warp1[n]*t1[f][0] + warp2[n]*t2[f][0];
        shifty[n] = warp1[n]*t1[f][1] + warp2[n]*t2[f][1];
        shiftz[n] = warp1[n]*t1[f][2] + warp2[n]*t2[f][2];
      }
    }
  }
  for (int n=0;n<_Npoints;n++) {
    _x[n] += shiftx[n];
    _y[n] += shifty[n];
    _z[n] += shiftz[n];
  }

  xyztorst(_Npoints, _x, _y, _z, _r, _s, _t);

  free(L1); free(L2); free(L3); free(L4);
  free(warp1); free(warp2);
  free(shiftx); free(shifty); free(shiftz);
  free(_x); free(_y); free(_z);
}

// ------------------------------------------------------------------------
// Gimbutas Xiao Cubature
// ------------------------------------------------------------------------

static const int cubTetNps[15] = {1, 4, 6, 11, 14, 23, 31, 44, 57, 74, 95, 122, 146, 177, 214};

static const dfloat cubR1[1] = {-5.000000000000000e-01};
static const dfloat cubS1[1] = {-5.000000000000000e-01};
static const dfloat cubT1[1] = {-4.999999999999999e-01};
static const dfloat cubW1[1] = { 1.333333333333333e+00};

static const dfloat cubR2[4] = {-6.399406129792692e-01,-6.881337590016283e-01,-5.678471416303044e-01, 6.431450819352392e-01};
static const dfloat cubS2[4] = {-2.693709623707308e-01,-8.507682582880892e-02,-9.992489699425411e-01,-7.526663993430834e-01};
static const dfloat cubT2[4] = {-9.861535288527451e-01,-2.364692878613065e-01,-1.385965858443282e-01,-9.201339027170031e-01};
static const dfloat cubW2[4] = { 4.006945857826349e-01, 3.717034355820901e-01, 4.254585806686327e-01, 1.354767312999743e-01};

static const dfloat cubR3[6] = {-6.197659951214324e-01,-6.826296735451189e-01,-9.781895755776219e-01, 1.424521042982300e-01,-6.583661496700222e-01,-6.759970166029514e-01};
static const dfloat cubS3[6] = {-1.202821047014497e-01,-7.503902756695054e-01,-3.091116885605387e-01,-7.170344960609916e-01,-9.242567364352863e-01, 2.828595829913932e-01};
static const dfloat cubT3[6] = {-9.771933411108853e-01, 1.713256113104316e-01,-4.369523957529076e-01,-7.061632198256604e-01,-6.943637138181457e-01,-6.322992990158042e-01};
static const dfloat cubW3[6] = { 1.767573695259264e-01, 2.992202225567431e-01, 1.707552169476622e-01, 3.360053077460037e-01, 1.874992128805373e-01, 1.630960036764600e-01};

static const dfloat cubR4[11] = {-8.055907108248335e-01,-9.408610095870402e-01,-1.345795219044632e-01,-5.194466701438546e-01,-7.411772524221795e-01,-7.569160173321442e-01,-9.846824781744662e-02,-1.614673722409738e-01,-8.655534102132334e-01, 5.050170140193101e-01,-9.190189865448195e-01};
static const dfloat cubS4[11] = {-7.867916548760125e-01,-3.415340805147064e-01,-7.923117671780141e-01,-3.911031951310060e-01, 7.601440783237155e-02,-9.820174798133283e-01,-1.340930190372888e-01,-8.933175209285095e-01, 4.824577641872456e-01,-8.371901631942812e-01,-6.506118826055384e-01};
static const dfloat cubT4[11] = { 3.687808309060797e-01,-3.641928795732112e-01,-2.923535215814058e-01,-7.463965481692164e-01,-3.396191703250704e-01,-3.870120231406193e-01,-8.810867674011320e-01,-9.044371288818271e-01,-9.296321404528026e-01,-8.638012581235871e-01,-9.728785962403947e-01};
static const dfloat cubW4[11] = { 1.419573788739868e-01, 1.469789765713298e-01, 2.066348213549953e-01, 2.578810827328453e-01, 1.015502869940778e-01, 1.059022400907004e-01, 9.262662125018049e-02, 7.991091352874596e-02, 7.385839849543516e-02, 7.369782554124912e-02, 5.233478789978663e-02};

static const dfloat cubR5[14] = {-8.145294993782177e-01,-9.100740825129933e-02,-9.089925917487004e-01,-9.089925917487010e-01,-3.782281614733988e-01,-8.653155155798035e-01,-3.782281614733992e-01,-3.782281614733983e-01, 4.435884981346523e-01,-8.145294993782176e-01,-8.145294993782173e-01,-9.100740825129877e-02,-9.100740825129949e-02,-9.089925917487007e-01};
static const dfloat cubS5[14] = {-8.145294993782174e-01,-9.089925917487008e-01,-9.100740825129933e-02,-9.089925917487003e-01,-3.782281614733988e-01,-3.782281614733984e-01,-8.653155155798030e-01,-3.782281614733984e-01,-8.145294993782171e-01,-8.145294993782171e-01, 4.435884981346519e-01,-9.100740825129915e-02,-9.089925917487007e-01,-9.100740825129915e-02};
static const dfloat cubT5[14] = { 4.435884981346527e-01,-9.100740825129856e-02,-9.100740825129978e-02,-9.100740825129978e-02,-8.653155155798035e-01,-3.782281614733988e-01,-3.782281614733987e-01,-3.782281614733988e-01,-8.145294993782176e-01,-8.145294993782176e-01,-8.145294993782176e-01,-9.089925917487014e-01,-9.089925917487002e-01,-9.089925917487014e-01};
static const dfloat cubW5[14] = { 9.799072415514931e-02, 5.672802770277523e-02, 5.672802770277523e-02, 5.672802770277538e-02, 1.502505676240207e-01, 1.502505676240207e-01, 1.502505676240207e-01, 1.502505676240207e-01, 9.799072415514914e-02, 9.799072415514931e-02, 9.799072415514914e-02, 5.672802770277523e-02, 5.672802770277523e-02, 5.672802770277523e-02};

static const dfloat cubR6[23] = {-9.223278313102310e-01,-8.704611261398938e-01,-8.704496791057892e-01,-4.441926613398445e-01,-8.678026751706386e-01,-3.497606828459493e-01,-3.616114393021378e-01,-3.432236575375570e-01,-8.898019550185483e-01,-7.507000727250277e-01,-8.681501536799804e-01,-9.852909523238611e-01, 2.349114402945376e-01,-4.411598941080236e-01,-4.245498103470715e-01, 1.894346037515911e-01,-8.664208004365244e-01, 2.530804034177644e-01,-8.799788339594619e-01,-4.484273990602982e-01,-8.973495876695942e-01,-9.188478978663642e-01, 8.075400026643635e-01};
static const dfloat cubS6[23] = {-9.513620515037141e-01,-4.643116036328488e-01,-9.530644088538914e-01,-8.725342094100046e-01,-8.326423718798897e-01,-3.412405629016029e-01,-3.916614693004364e-01,-9.234226585235092e-01,-2.961216053305913e-01,-6.957923773801387e-01, 2.486427271068594e-01,-5.774046828368271e-01,-8.736000381148603e-01,-4.883584314700274e-01, 1.546915627794538e-01,-8.696440144732590e-01, 6.012655096203365e-02,-5.031009197622103e-01,-5.739176335276283e-01,-8.920077183281708e-01, 6.822779033246367e-01,-9.824360844449622e-01,-9.542683523719533e-01};
static const dfloat cubT6[23] = { 8.058575980272226e-01, 2.735350171170275e-01,-2.182758986579757e-01, 1.898193780435916e-01, 2.601091102219790e-01,-3.463329907619084e-01,-9.112333112855833e-01,-3.594251326046153e-01,-2.378313821873799e-01,-5.975308652711581e-01,-4.928126505135995e-01,-4.976310094449407e-01,-4.831017020321489e-01,-4.608601407334560e-01,-8.707587238532627e-01,-8.667934039847933e-01,-8.460145657980651e-01,-8.757689336328026e-01,-9.483146274785937e-01,-8.799677016676628e-01,-9.254704957232888e-01,-8.227992990621793e-01,-9.413285578336420e-01};
static const dfloat cubW6[23] = { 9.461059802212705e-03, 4.201254651027517e-02, 3.230838250325908e-02, 6.518676786992283e-02, 7.071109854422583e-02, 5.765239559396446e-02, 8.951442161674263e-02, 7.976179688190549e-02, 8.348596704174839e-02, 8.577869596411661e-02, 5.812053074750553e-02, 3.008755637085715e-02, 6.215084550210757e-02, 1.501325393254235e-01, 7.191334750441594e-02, 6.635817345535235e-02, 8.408848251402726e-02, 6.286404062968159e-02, 3.400576568938985e-02, 5.295213019877631e-02, 2.123397224667169e-02, 1.389778096492792e-02, 9.655035855822626e-03};

static const dfloat cubR7[31] = {-9.615840130228295e-01,-3.530863164208045e-01,-8.750519549536997e-01,-9.719338310533894e-01,-8.751111745741820e-01,-4.769824616468345e-01,-4.791588240034766e-01,-4.834389052651127e-01,-2.056511940101652e-01, 1.009118832497191e-01,-9.483336645365394e-01,-7.850751384336935e-01,-4.957389387412098e-01, 2.036471552636233e-01,-9.336686003379285e-01,-7.802080473459577e-01,-1.792838080101212e-01,-6.471235539701039e-01,-6.017788558942327e-01, 5.742261282730587e-02,-9.921406970258241e-01,-5.530595949397141e-01,-4.044295312329518e-01,-4.124599269525862e-01,-8.942447412422904e-01,-8.939997233309063e-01, 1.982018872400515e-01, 6.414906014831900e-01,-8.802116098600259e-01,-8.759778161267122e-01, 2.710430211675223e-01};
static const dfloat cubS7[31] = {-9.960063483633999e-01,-8.781556308290985e-01,-8.720611393484018e-01,-5.384301995246593e-01,-3.521038580216352e-01,-2.997693909856581e-01,-9.549340330123351e-01, 7.573396996314106e-04,-4.202774361826186e-01,-6.855162880279936e-01, 2.630969478360095e-01,-7.774977331461141e-01,-5.648911115672932e-01,-9.107394641867267e-01,-8.475170594320620e-01,-1.255038204709026e-01,-8.519825957217687e-01,-5.623747549049908e-01,-9.967729347600200e-01,-5.312070053280753e-01,-2.989431863254339e-01, 1.017779562844258e-01,-8.639777902014777e-01, 2.559664587171953e-01, 1.980384797597949e-01,-4.493803786257354e-01,-3.606832478059763e-01,-8.743698543095261e-01, 6.490881333907904e-01,-8.957266218839776e-01,-9.989991331114559e-01};
static const dfloat cubT7[31] = { 3.026697916964751e-01, 2.303419766237413e-01, 6.257451111419999e-01, 1.323261490613942e-01, 1.726425716602438e-01,-9.760482462216850e-01, 1.293002049353824e-01,-8.346618087852156e-01,-6.579022443624429e-01,-9.125130596785365e-01,-4.391895237961879e-01,-4.830167652136928e-02,-1.025214892328492e-01,-4.200119312689738e-01,-4.503969302040128e-01,-4.468566345222848e-01,-5.104766981613168e-01,-6.484677794670978e-01,-4.725052493228704e-01,-5.660034013082684e-01,-5.163107738828342e-01,-5.489429246054712e-01,-9.047011137982147e-01,-9.450352436143311e-01,-8.900408342089723e-01,-9.100446737462401e-01,-9.295724094050908e-01,-8.821791743487821e-01,-8.864482666201061e-01,-8.856356909664419e-01,-8.802110555361827e-01};
static const dfloat cubW7[31] = { 1.027757488266732e-02, 1.600505625095977e-02, 3.095574577048298e-02, 2.988182052718841e-02, 3.105464927642637e-02, 3.643988575686199e-02, 3.839504479364609e-02, 5.732154928530760e-02, 7.037555420473829e-02, 5.578470881013387e-02, 2.761610955168380e-02, 8.040202078827963e-02, 8.894827622097576e-02, 4.188150327399958e-02, 3.775038348257681e-02, 8.898598589450760e-02, 8.558511858061323e-02, 1.097842250407599e-01, 2.783027354905359e-02, 3.274990651463104e-02, 2.979822864344559e-02, 2.957375890028240e-02, 4.924774548191803e-02, 2.226933475152646e-02, 4.506115373700136e-02, 3.741490441857740e-02, 2.812368778917140e-02, 2.745379817637410e-02, 2.560374906157176e-02, 2.407610331221044e-02, 1.668547660576086e-02};

static const dfloat cubR8[44] = {-9.822888114425899e-01,-9.275906561429631e-01, 2.226730673837682e-01,-9.444464840128696e-01,-6.870266299198677e-01,-6.089078361547958e-01,-9.170116819634433e-01,-5.580450700023645e-01,-1.698468314600223e-01,-4.118669329426861e-02,-8.755112010779593e-01,-5.749615477649417e-01,-2.054490850771143e-01,-9.091052514742826e-01,-5.677206147398416e-01,-9.094185002228330e-01,-1.431844585279626e-01,-7.149707401822105e-01,-6.374756176581090e-01,-5.009584678400951e-01,-2.754950826768038e-02,-5.450759000264900e-01,-6.834818133082194e-01,-9.423454821981293e-01,-5.409187698728505e-01,-2.937288472386135e-02,-9.547781807266793e-01, 4.410223314357435e-01,-6.936967917392711e-02,-9.212395626309886e-01,-9.091309392957807e-01,-1.506998918408216e-01,-5.674717022186972e-01,-5.910163725055539e-01, 4.391466601433898e-01,-9.075988357026712e-01, 6.593440659526384e-04,-9.178994037090353e-01,-9.202866441531541e-01, 5.247357699089055e-01,-3.490076435591544e-03,-9.026219907219561e-01,-9.398955044481340e-01, 9.147632521166054e-01};
static const dfloat cubS8[44] = {-9.647856877462442e-01,-6.218510622599334e-01,-5.385130064107014e-01,-9.483863022549348e-01,-9.070277411778146e-01,-5.007258338842491e-01,-9.018293113450074e-01,-4.841666536913819e-01,-9.137552332442402e-01,-6.001865419341569e-01,-7.520430213443628e-01,-9.323265238055517e-01,-2.342389976685463e-01,-5.573356090065297e-02,-1.281983833195988e-02,-4.995451074139229e-01,-6.355079015532782e-01,-5.438396318898419e-01,-9.383550558001381e-01,-6.706089881993432e-01,-9.151393977087543e-01,-5.993385862443937e-01,-1.762075927031755e-01, 2.221214382715361e-02,-3.296343941797995e-02,-5.807153942587078e-01,-4.990561561653583e-01,-9.264387005305038e-01,-7.637654966347461e-02,-8.990899861365506e-01, 4.971491214915261e-01,-9.380989415677337e-01,-8.912339687875325e-01, 4.033135521287796e-01,-9.104395690079359e-01, 4.474189900834977e-01,-9.284748582425598e-01,-1.024439844950352e-01,-6.712261380303836e-01,-6.987245764371494e-01,-1.721574331223598e-01,-9.908576312509686e-01, 8.816164056632010e-01,-9.980562432619695e-01};
static const dfloat cubT8[44] = { 8.971278595980133e-01, 4.627347200721783e-01,-9.888484539391230e-01,-6.691882821816188e-02, 5.029930723036760e-01, 8.358636492361879e-02, 4.343547141878999e-01,-9.752081532796095e-01,-1.040166136001711e-02,-3.598087864992038e-01,-1.926904499798824e-01, 1.128390196074010e-02,-8.365327435904173e-01,-1.286660109086743e-01,-9.216188536049289e-01,-1.371046408187091e-02,-8.761060769063105e-01,-7.437125204667648e-01,-4.813072216807797e-01,-5.001888985484900e-01,-4.495279454001869e-01,-1.020909426726044e-01,-5.852165297187530e-01,-5.409185738546890e-01,-5.159648853650391e-01,-6.027600446615375e-01,-5.551307150144428e-01,-5.912621521853492e-01,-9.891556333856432e-01,-6.536087826332979e-01,-6.734013134012919e-01,-7.275875175496033e-01,-8.974924697322526e-01,-9.062626934382535e-01,-8.880187106357958e-01,-9.130023710960993e-01,-9.837758984885883e-01,-9.086292685077639e-01,-9.310358490421146e-01,-8.951232201425936e-01,-8.404842580661891e-01,-9.556758578066753e-01,-9.865606423556936e-01,-9.679840315374033e-01};
static const dfloat cubW8[44] = { 2.326911380407538e-03, 2.008502598013748e-02, 1.316868905984896e-02, 1.350522833433149e-02, 2.348338158017796e-02, 2.082185227170751e-02, 2.398518653047993e-02, 2.203915573406974e-02, 3.000608795736148e-02, 1.706584449682616e-02, 3.633991457929743e-02, 4.305733418779512e-02, 3.758091453334270e-02, 3.221112523708480e-02, 4.298833779190321e-02, 4.873345427684407e-02, 5.542332055396455e-02, 5.208361326526754e-02, 3.351878142745155e-02, 7.048390494136669e-02, 4.363636633805530e-02, 7.438998373030328e-02, 7.113840326380720e-02, 3.615671798064746e-02, 5.346507352177339e-02, 6.903949858534457e-02, 3.032496084600489e-02, 1.962654290056966e-02, 9.196571088500877e-03, 2.263363901068635e-02, 2.318560341668196e-02, 3.007898500466981e-02, 3.133189805159075e-02, 2.686228112905257e-02, 2.739338179275024e-02, 2.386507174305719e-02, 1.002820046914775e-02, 2.923266480869545e-02, 1.831436341104725e-02, 2.051792721654784e-02, 1.582056092943098e-02, 3.395867359123967e-03, 2.772112756025174e-03, 2.018593860153416e-03};

static const dfloat cubR9[57] = {-5.533046643727986e-01,-9.158670291802449e-01,-9.701523520430464e-01,-9.206457589681702e-01,-9.337366660891633e-01,-6.380957350966594e-01,-7.135713192043529e-01,-6.107566394690580e-01,-1.298676464011230e-01,-1.664358333494954e-01,-7.514287025389632e-01,-1.113095043790728e-01,-7.342331339208183e-01,-9.689597680156343e-01,-9.266490702197298e-01,-5.373715022569363e-01,-8.943748201422934e-01,-9.451345695722863e-01,-6.401909968058997e-01,-5.456802613699343e-01,-2.510060418460261e-01,-2.796865735769462e-01,-5.420319818586138e-01,-3.327616188639549e-01, 1.459127371604473e-01,-7.046258637970575e-01,-1.555427966561573e-01,-1.704761554019259e-01,-2.925992831309229e-01,-6.698743181853890e-01,-7.130341551605860e-01,-7.386686248489396e-01,-7.349912291618884e-01,-3.211354118094386e-01,-8.899252123406390e-01, 2.908989268740044e-01, 4.407614406718393e-01,-9.103378816977798e-01,-5.683089329465612e-01,-9.648263734898368e-01,-9.539504662043419e-01, 2.670695115549040e-01,-9.581114473978323e-01, 2.887166196779363e-01, 4.757968987501355e-01, 1.732388127357021e-01,-6.257054157709013e-01,-3.253520660505721e-02,-9.365768769292344e-01,-9.395805700596319e-01, 5.660476761353139e-01,-9.510355234372726e-01,-1.189571903053532e-01,-5.633928989032095e-01,-9.082033404150098e-01, 7.479685689385656e-01,-9.250656883579180e-01};
static const dfloat cubS9[57] = {-9.192284574478844e-01,-9.160905359846262e-01,-5.732924390755566e-01,-9.766416244652930e-01,-9.197772435191071e-01,-9.490813120195858e-01, 2.246816043641724e-01,-2.385259301169455e-01,-9.883852374828955e-01,-5.463253600713984e-01,-6.152564094597048e-01,-8.778467993652843e-01,-2.864957262025621e-01,-2.760982372512223e-01,-8.351066097292752e-02,-7.939769723872107e-01,-6.320292181652613e-01,-7.515287540235962e-01,-9.217593295662161e-01,-7.165676802509563e-01,-5.555048525674710e-01,-8.250864840889946e-02,-3.271309457569846e-01,-4.580675737837868e-01,-6.662500542586419e-01,-6.950633656529526e-01,-8.331524229587274e-01,-2.732866397752712e-01,-7.807750833428199e-01,-6.360455530450386e-01,-2.425112117576626e-01,-9.083825243710147e-01, 1.143206439498955e-01,-9.692208003786650e-01,-2.967816456232522e-01,-7.936857063434514e-01,-9.862172663214093e-01, 4.196054391172780e-01, 1.310327041650577e-01, 2.088525999212051e-01,-9.372243266147463e-01,-4.478451203530572e-01,-6.322900043647489e-01,-9.606988791583196e-01,-8.972526359227828e-01,-4.524900616954256e-01, 4.609979885240922e-01,-9.333455935924870e-01,-5.384537758565796e-01, 5.037786889766072e-01,-6.888421833945486e-01, 1.540226645667018e-03, 3.057012311455353e-02,-9.832534620977386e-01,-9.082819610479546e-01,-9.305195511631422e-01, 7.695471905246293e-01};
static const dfloat cubT9[57] = { 4.474481896141336e-01, 7.503920530568159e-01, 4.552224596612802e-01,-5.365671713622035e-02, 4.273822481233936e-01, 3.641172259369879e-01,-9.986675641825126e-01,-1.589896515942665e-01,-1.056316993155485e-01,-9.844032143628079e-01, 3.122764470087541e-01,-5.028623426553014e-02,-9.589708571964229e-01,-2.151326389458592e-01,-7.282854719917632e-02,-9.942708440232070e-01, 1.582881449664418e-01,-1.749137514146293e-01,-4.698119799191238e-02, 3.519795905700290e-02,-2.721492944747828e-01,-9.262003377957502e-01,-8.345329893053417e-01,-6.003539397970914e-01,-8.445114632975066e-01,-8.118584756178400e-01,-3.937162708510041e-01,-6.749472724987970e-01,-7.641311876161534e-01,-3.815071188418233e-01,-3.061183874904804e-01,-5.290457907046182e-01,-7.551928096236605e-01,-5.162609760038315e-01,-6.061178761247058e-01,-5.948965692091915e-01,-5.118426729186494e-01,-5.773640152816579e-01,-6.276181223461234e-01,-5.925662730794791e-01,-5.979697800078242e-01,-9.594102826767850e-01,-6.591221643854689e-01,-7.318109291016270e-01,-9.789502598685401e-01,-7.336525734932146e-01,-9.258602682046155e-01,-9.255429370952352e-01,-9.417017379403116e-01,-9.213734346055054e-01,-9.255453718785110e-01,-9.013357205929660e-01,-9.338512861491433e-01,-8.723247837732444e-01,-9.213986974859054e-01,-9.039405322254990e-01,-9.137203010999236e-01};
static const dfloat cubW9[57] = { 7.800212910650462e-03, 9.458250846547485e-03, 9.474373461649973e-03, 9.599714961394884e-03, 1.568338205232868e-02, 1.676387616329851e-02, 1.050685463763321e-02, 1.310048147306960e-02, 1.303920343219988e-02, 1.902192811379036e-02, 1.986073113428314e-02, 1.549993751121349e-02, 1.771500052129017e-02, 2.015837198880603e-02, 2.095650303125563e-02, 1.093106704497207e-02, 3.504065691166473e-02, 2.651841414209441e-02, 3.516745122395425e-02, 4.351961177069496e-02, 3.372350242246109e-02, 3.385668187512971e-02, 4.696048351121531e-02, 5.583379695105220e-02, 4.246040614681167e-02, 3.750295819431537e-02, 5.096144788527118e-02, 3.420250930904183e-02, 4.703915299077518e-02, 6.546949941566481e-02, 6.028184410562054e-02, 3.067931896236044e-02, 4.762204359764656e-02, 2.196718770411703e-02, 4.416143256303730e-02, 2.984094802302178e-02, 5.729612835236330e-03, 1.731670818374603e-02, 2.847239961322439e-02, 1.730034429779180e-02, 9.395696925810014e-03, 1.374864781813647e-02, 1.823631387807064e-02, 2.051164372655097e-02, 9.520052954932823e-03, 9.375535059554770e-03, 1.744782605744033e-02, 1.769855472899909e-02, 1.222025919372541e-02, 1.219542757986468e-02, 7.262216915703982e-03, 1.764147054086365e-02, 8.801520740585440e-03, 1.106164768380024e-02, 1.069494738034133e-02, 8.988533895500545e-03, 7.334736333120193e-03};

static const dfloat cubR10[74] = { 6.517784730010368e-02,-7.619105379214477e-01,-7.757512270699543e-01,-1.820930647800567e-02,-2.607747102875660e-01,-7.709240273679512e-01,-9.385960782586074e-01,-9.376373624955783e-01,-7.467210182428361e-01,-3.400316248731436e-01,-9.652919051084328e-01, 2.405166025801018e-01,-4.027298204212719e-01,-7.613548130161666e-01,-9.302815128201264e-02,-1.650474962583800e-01,-9.663808170607030e-01,-1.912459710232673e-01,-9.296208354624638e-01,-6.976708558367206e-01, 5.611444916483015e-01,-9.538054824488602e-01,-6.607723048792105e-01, 2.678560687806029e-01,-6.459632979621646e-01,-9.513179433594532e-01, 2.671751306285224e-01,-2.979381466761734e-01, 1.720291110545033e-01,-9.394121244845574e-01,-3.433089570523171e-01,-8.667471154851962e-01,-3.959961853465120e-01,-9.746471367028977e-01,-6.709055682201694e-01,-1.725977431011307e-01,-3.209738747577238e-01,-7.426376635037392e-01,-7.585737767282890e-01,-8.339501062356643e-01, 4.564189619948546e-01,-6.994624194964894e-01,-7.350659785104570e-01,-5.116605778881920e-01, 3.183363322888239e-02,-1.747983276416606e-01,-6.455712538545869e-01,-9.754391428751119e-01, 6.298045658810059e-01,-9.329687289891602e-01,-9.355894271645430e-01,-9.412897662791638e-01,-6.519924690908426e-01,-6.484007855015702e-01, 1.985476022784536e-01,-9.340318800246465e-01,-9.377107750198432e-01,-9.285218311769171e-01,-9.324417488177632e-01,-3.241021603478104e-01,-6.384795008595457e-01,-4.957182683133580e-01, 1.193201777114047e-01,-9.866658391033001e-01,-9.229088686052509e-01, 8.934396965837035e-01,-7.928126795048904e-01,-9.923462161339119e-01, 7.480502655303353e-01,-4.858560384748986e-01,-9.903102209448552e-01,-6.660804613648394e-01,-3.635875983588155e-01,-7.326801009347781e-01};
static const dfloat cubS10[74] = {-7.757512270700414e-01, 6.517784729987100e-02,-7.619105379214348e-01,-7.709240273683173e-01,-1.820930647889407e-02,-2.607747102877043e-01,-7.467210182422573e-01,-9.385960782585098e-01,-9.376373624955957e-01, 2.405166025797735e-01,-3.400316248725770e-01,-9.652919051086308e-01,-9.302815128242126e-02,-4.027298204213297e-01,-7.613548130165168e-01,-1.912459710232630e-01,-1.650474962581840e-01,-9.663808170608019e-01, 5.611444916487768e-01,-9.296208354625831e-01,-6.976708558369646e-01, 2.678560687807144e-01,-9.538054824488867e-01,-6.607723048791950e-01, 2.671751306287957e-01,-6.459632979621062e-01,-9.513179433594716e-01,-9.394121244845912e-01,-2.979381466764564e-01, 1.720291110549545e-01,-3.959961853467096e-01,-3.433089570523432e-01,-8.667471154853744e-01,-1.725977431009272e-01,-9.746471367029086e-01,-6.709055682199777e-01,-7.585737767283822e-01,-3.209738747581960e-01,-7.426376635036127e-01,-6.994624194965924e-01,-8.339501062360731e-01, 4.564189619941376e-01, 3.183363322854889e-02,-7.350659785107161e-01,-5.116605778886708e-01,-9.754391428751308e-01,-1.747983276415959e-01,-6.455712538545152e-01,-9.355894271645714e-01, 6.298045658815959e-01,-9.329687289891651e-01,-6.484007855013953e-01,-9.412897662792561e-01,-6.519924690907442e-01,-9.377107750198261e-01, 1.985476022790324e-01,-9.340318800246411e-01,-3.241021603471450e-01,-9.285218311769089e-01,-9.324417488177208e-01, 1.193201777109193e-01,-6.384795008598869e-01,-4.957182683139326e-01, 8.934396965847128e-01,-9.866658391026820e-01,-9.229088686050392e-01, 7.480502655298643e-01,-7.928126795039943e-01,-9.923462161339409e-01,-4.858560384751444e-01,-9.903102209444652e-01,-6.660804613655873e-01,-3.635875983593522e-01,-7.326801009352152e-01};
static const dfloat cubT10[74] = {-5.275160823084778e-01,-5.275160823086784e-01,-5.275160823086025e-01,-9.500919558655234e-01,-9.500919558655199e-01,-9.500919558656692e-01, 6.229544589963890e-01, 6.229544589959419e-01, 6.229544589969793e-01,-9.351930725981743e-01,-9.351930725981473e-01,-9.351930725981669e-01,-7.428872152800948e-01,-7.428872152804732e-01,-7.428872152800838e-01,-6.773257156576789e-01,-6.773257156578382e-01,-6.773257156578504e-01,-9.338528003489132e-01,-9.338528003488458e-01,-9.338528003488543e-01,-6.532782814526801e-01,-6.532782814525367e-01,-6.532782814524829e-01,-6.698938893072007e-01,-6.698938893072497e-01,-6.698938893069852e-01,-9.346788398934703e-01,-9.346788398934202e-01,-9.346788398934556e-01,-3.939477421158102e-01,-3.939477421162904e-01,-3.939477421158435e-01,-1.818495519762574e-01,-1.818495519761044e-01,-1.818495519759868e-01,-1.778146850102663e-01,-1.778146850100838e-01,-1.778146850099773e-01,-9.230064362620708e-01,-9.230064362620121e-01,-9.230064362620096e-01,-7.851070768298830e-01,-7.851070768296919e-01,-7.851070768296160e-01,-2.041912756285660e-01,-2.041912756288685e-01,-2.041912756288133e-01,-7.612464097273302e-01,-7.612464097278605e-01,-7.612464097282768e-01, 2.416830208712831e-01, 2.416830208714436e-01, 2.416830208715526e-01,-3.268049472339450e-01,-3.268049472345868e-01,-3.268049472349763e-01, 1.850657403417951e-01, 1.850657403414275e-01, 1.850657403424453e-01,-9.851224085376650e-01,-9.851224085376012e-01,-9.851224085375707e-01,-9.838649888761907e-01,-9.838649888764406e-01,-9.838649888756519e-01,-9.628913698913086e-01,-9.628913698913258e-01,-9.628913698912855e-01,-5.424318845750749e-01, 9.709306628334797e-01,-1.758615903775734e-03,-9.092372049227418e-01, 1.980403028050551e-01};
static const dfloat cubW10[74] = { 2.756079111468018e-02, 2.756079111473660e-02, 2.756079111473208e-02, 1.651728622129342e-02, 1.651728622129794e-02, 1.651728622126273e-02, 8.827127869984224e-03, 8.827127869990035e-03, 8.827127869980037e-03, 8.185538986460163e-03, 8.185538986463684e-03, 8.185538986433970e-03, 3.387610113488516e-02, 3.387610113487611e-02, 3.387610113486918e-02, 1.675748787500219e-02, 1.675748787498494e-02, 1.675748787495920e-02, 9.313281339996723e-03, 9.313281340003173e-03, 9.313281340004842e-03, 1.753135755552536e-02, 1.753135755551546e-02, 1.753135755551263e-02, 1.827560624231835e-02, 1.827560624231354e-02, 1.827560624230717e-02, 1.218397324826652e-02, 1.218397324827248e-02, 1.218397324827924e-02, 3.576925978360336e-02, 3.576925978361185e-02, 3.576925978361609e-02, 1.498958846787545e-02, 1.498958846785876e-02, 1.498958846786880e-02, 3.745674993988912e-02, 3.745674993992178e-02, 3.745674993992857e-02, 1.667130533373333e-02, 1.667130533371325e-02, 1.667130533372230e-02, 3.643106559328947e-02, 3.643106559329400e-02, 3.643106559330701e-02, 1.575595091854842e-02, 1.575595091860258e-02, 1.575595091855973e-02, 9.531695272862311e-03, 9.531695272852877e-03, 9.531695272855819e-03, 2.221923927343592e-02, 2.221923927340425e-02, 2.221923927340948e-02, 1.356664812179151e-02, 1.356664812180096e-02, 1.356664812179412e-02, 1.518962478003038e-02, 1.518962478002415e-02, 1.518962478003236e-02, 1.094939643308055e-02, 1.094939643309527e-02, 1.094939643309776e-02, 1.053465486327635e-03, 1.053465486335802e-03, 1.053465486342203e-03, 2.256754452207803e-03, 2.256754452196206e-03, 2.256754452192911e-03, 6.210469597922555e-02, 5.092217550654402e-04, 1.567347452275818e-02, 3.721542358429886e-02, 1.522263115666367e-02};

static const dfloat cubR11[95] = {-7.679538362911502e-01,-9.609773918290878e-01, 4.577464989682887e-01,-9.847579184704128e-01,-9.757022456214157e-01,-5.321454743397008e-02,-2.733049486061983e-01,-3.220780592357142e-01,-9.662227196205052e-01,-9.653902083247552e-01,-7.118900754072411e-01, 3.617214902690057e-01,-7.140107196669039e-01,-5.230512706308128e-01,-5.163574270884314e-02,-9.644032204152739e-01,-9.485865642419784e-01, 6.223979628844538e-01,-9.409645149934525e-01,-6.420356285802742e-01,-1.099587574693295e-02,-4.620133965830218e-01,-6.950687687402711e-01, 1.127037388627261e-01,-9.455231122055394e-01,-9.414276020716339e-01,-7.573236015256012e-01, 8.935000911028185e-02,-6.043209058679191e-01,-9.407945482659268e-01,-9.762651325032106e-01,-7.474011870194273e-01, 6.566822622786648e-01,-9.354023733633352e-01,-3.256618598446160e-01,-3.035899334160957e-03, 1.133388069533189e-01,-3.069892546981613e-01,-9.611523296968980e-01,-8.134322794944013e-01,-6.853508675135351e-01, 4.750316521804119e-01,-9.322209866584998e-01,-9.431610059527277e-01,-4.396671109028148e-01,-3.784511526616307e-01,-3.986203922912018e-01,-7.449087452771068e-01,-4.270101772449363e-01,-7.768065929337287e-01,-1.782389476407426e-02,-7.489277736663206e-01,-7.702356064364334e-01,-7.977739065333606e-02,-4.998935175251013e-01,-4.984767801597670e-01,-9.187915156179884e-01,-9.531077681358991e-01,-3.828531583576022e-01, 2.706389544135504e-01,-8.194129287473002e-01,-8.147665099539618e-01,-3.818126395293394e-01,-7.995290508026287e-01,-7.584974927219390e-01, 3.556702729058751e-01,-6.901428251973104e-01,-6.840899148041231e-01,-9.756798215296782e-01,-7.237861161653755e-01,-3.290784063369139e-01,-1.260874584214489e-02,-7.662036066132658e-01,-9.799751554857663e-01,-1.547579082025432e-01,-7.889561391572769e-01,-1.882833225713522e-01,-9.878548432217089e-01,-9.223142707500790e-01,-9.467509756719994e-01, 8.132400421776769e-01,-6.296772266165893e-01,-9.365909943466163e-01, 5.141113529511696e-01,-9.258450401605933e-01,-9.303159525936109e-01, 2.267111719712813e-01, 5.610706757272649e-02,-1.326212642699352e-01,-9.317357523187260e-01,-9.565128767575304e-01,-7.720609303767829e-01,-3.444577745519064e-01,-4.079830601730307e-01,-6.063462709910203e-01};
static const dfloat cubS11[95] = { 4.577464989693064e-01,-7.679538362913773e-01,-9.609773918292007e-01,-5.321454743520816e-02,-9.847579184706120e-01,-9.757022456215160e-01,-9.662227196204634e-01,-2.733049486047525e-01,-3.220780592358554e-01, 3.617214902691905e-01,-9.653902083250243e-01,-7.118900754074279e-01,-5.163574270829095e-02,-7.140107196671319e-01,-5.230512706315071e-01, 6.223979628860419e-01,-9.644032204153778e-01,-9.485865642419543e-01,-1.099587574519130e-02,-9.409645149929431e-01,-6.420356285814927e-01, 1.127037388631749e-01,-4.620133965824464e-01,-6.950687687402429e-01,-7.573236015326951e-01,-9.455231122053777e-01,-9.414276020709444e-01,-9.407945482658525e-01, 8.935000911195126e-02,-6.043209058673978e-01, 6.566822622788856e-01,-9.762651325028601e-01,-7.474011870194225e-01,-3.035899335218925e-03,-9.354023733634328e-01,-3.256618598456624e-01,-9.611523296970399e-01, 1.133388069531991e-01,-3.069892546977568e-01, 4.750316521818557e-01,-8.134322794935307e-01,-6.853508675140709e-01,-4.396671109093365e-01,-9.322209866583300e-01,-9.431610059534287e-01,-7.449087452768705e-01,-3.784511526607067e-01,-3.986203922900602e-01,-1.782389476426723e-02,-4.270101772455154e-01,-7.768065929339733e-01,-7.977739064995708e-02,-7.489277736654356e-01,-7.702356064370101e-01,-9.187915156177521e-01,-4.998935175247330e-01,-4.984767801581871e-01, 2.706389544134378e-01,-9.531077681357768e-01,-3.828531583575322e-01,-3.818126395269364e-01,-8.194129287468309e-01,-8.147665099537936e-01, 3.556702729076395e-01,-7.995290508029310e-01,-7.584974927226493e-01,-9.756798215293145e-01,-6.901428251965293e-01,-6.840899148027154e-01,-1.260874584188131e-02,-7.237861161651740e-01,-3.290784063371327e-01,-1.547579082019079e-01,-7.662036066132224e-01,-9.799751554861351e-01,-9.878548432216758e-01,-7.889561391574725e-01,-1.882833225679202e-01, 8.132400421778441e-01,-9.223142707500882e-01,-9.467509756719902e-01, 5.141113529510306e-01,-6.296772266166528e-01,-9.365909943466179e-01, 2.267111719745459e-01,-9.258450401605582e-01,-9.303159525935650e-01,-9.317357523188697e-01, 5.610706757214612e-02,-1.326212642690617e-01,-9.565128767605066e-01,-7.720609303752631e-01,-3.444577745515867e-01,-4.079830601730581e-01,-6.063462709910926e-01};
static const dfloat cubT11[95] = {-7.288152708490122e-01,-7.288152708490550e-01,-7.288152708477885e-01, 1.367471152724070e-02, 1.367471152812739e-02, 1.367471152733010e-02,-4.383942725383680e-01,-4.383942725392109e-01,-4.383942725394879e-01,-6.844412065367120e-01,-6.844412065362393e-01,-6.844412065365075e-01,-7.113022669934016e-01,-7.113022669941305e-01,-7.113022669923007e-01,-7.094081782286773e-01,-7.094081782288734e-01,-7.094081782271893e-01,-4.060039806819337e-01,-4.060039806810906e-01,-4.060039806785093e-01,-9.556215735400153e-01,-9.556215735399612e-01,-9.556215735397813e-01, 6.442743158095920e-01, 6.442743158098897e-01, 6.442743158014782e-01,-5.442345549775491e-01,-5.442345549782510e-01,-5.442345549787743e-01,-9.330159427561004e-01,-9.330159427561836e-01,-9.330159427562670e-01,-7.358998674567107e-01,-7.358998674568270e-01,-7.358998674567119e-01,-8.451972225576958e-01,-8.451972225579445e-01,-8.451972225585827e-01,-9.762485051739834e-01,-9.762485051739271e-01,-9.762485051726448e-01, 3.150491035205434e-01, 3.150491035209673e-01, 3.150491035146120e-01,-4.780197097722490e-01,-4.780197097710318e-01,-4.780197097723581e-01,-7.783593350567001e-01,-7.783593350569291e-01,-7.783593350564662e-01,-4.010592292476484e-01,-4.010592292482036e-01,-4.010592292446700e-01,-8.283818669742238e-02,-8.283818669749098e-02,-8.283818669921417e-02,-9.346780279198313e-01,-9.346780279200211e-01,-9.346780279201987e-01, 1.599207822769532e-02, 1.599207822741243e-02, 1.599207823167453e-02,-7.976437293829682e-01,-7.976437293826510e-01,-7.976437293804048e-01, 3.499125615299393e-01, 3.499125615300912e-01, 3.499125615286900e-01,-9.345267316555559e-01,-9.345267316557261e-01,-9.345267316555975e-01,-9.906332969910832e-02,-9.906332969919528e-02,-9.906332969842246e-02,-3.490569505225247e-02,-3.490569504842889e-02,-3.490569505272526e-02,-9.441747957558063e-01,-9.441747957558688e-01,-9.441747957555076e-01,-9.478431319878280e-01,-9.478431319878011e-01,-9.478431319876015e-01,-3.705501792204114e-01,-3.705501792203661e-01,-3.705501792172956e-01,-9.917500509834563e-01,-9.917500509834758e-01,-9.917500509840943e-01, 8.695386302787249e-01, 3.161827911273822e-01,-9.666266763448393e-01,-7.760508194806258e-01,-1.809611870271467e-01};
static const dfloat cubW11[95] = { 8.497233737804958e-03, 8.497233737812155e-03, 8.497233737828392e-03, 1.951372627176398e-03, 1.951372627182381e-03, 1.951372627096651e-03, 1.301998872398560e-02, 1.301998872404748e-02, 1.301998872395198e-02, 9.358078434992697e-03, 9.358078434974141e-03, 9.358078434968938e-03, 2.625573745688143e-02, 2.625573745682995e-02, 2.625573745703841e-02, 4.223936341047088e-03, 4.223936341046989e-03, 4.223936341074736e-03, 1.559706111217679e-02, 1.559706111224566e-02, 1.559706111238764e-02, 1.498965770388662e-02, 1.498965770387997e-02, 1.498965770393060e-02, 5.786148186164216e-03, 5.786148186153031e-03, 5.786148186293616e-03, 1.657203299779843e-02, 1.657203299780267e-02, 1.657203299778386e-02, 3.541235839082939e-03, 3.541235839111025e-03, 3.541235839096940e-03, 1.719845284625162e-02, 1.719845284624525e-02, 1.719845284627750e-02, 1.107318891220779e-02, 1.107318891220495e-02, 1.107318891224347e-02, 6.856663583113830e-03, 6.856663583126162e-03, 6.856663583292700e-03, 8.445303975913391e-03, 8.445303975919216e-03, 8.445303975732527e-03, 3.677904306744947e-02, 3.677904306745668e-02, 3.677904306746672e-02, 2.673380190572628e-02, 2.673380190566716e-02, 2.673380190568145e-02, 2.886549696168238e-02, 2.886549696166923e-02, 2.886549696173047e-02, 2.314310172287898e-02, 2.314310172278861e-02, 2.314310172283090e-02, 8.102945026016736e-03, 8.102945026019225e-03, 8.102945025991761e-03, 2.321444526942162e-02, 2.321444526944185e-02, 2.321444526939263e-02, 2.027279380665854e-02, 2.027279380666717e-02, 2.027279380665430e-02, 9.602352588021107e-03, 9.602352587982557e-03, 9.602352587981283e-03, 2.077696376975419e-02, 2.077696376971826e-02, 2.077696376977017e-02, 9.634057010055147e-03, 9.634057010051936e-03, 9.634057009965303e-03, 6.818068346919792e-03, 6.818068346716189e-03, 6.818068346905608e-03, 3.725019917228112e-03, 3.725019917221903e-03, 3.725019917241533e-03, 8.102738004743620e-03, 8.102738004744398e-03, 8.102738004771255e-03, 1.204124614555485e-02, 1.204124614556417e-02, 1.204124614561331e-02, 4.360434875712713e-03, 4.360434875708343e-03, 4.360434875655932e-03, 1.481640621184332e-03, 2.399497001991400e-02, 1.551667731231699e-02, 3.968612441682814e-02, 3.603811827413330e-02};

static const dfloat cubR12[122] = {-1.374503437769854e-01,-9.291591635595258e-01, 4.494310271081625e-02,-9.422439932245119e-01, 2.430423677648110e-01,-7.956608618114338e-01,-8.208949817859721e-01,-7.300121447810803e-01,-1.786133718196966e-01,-4.798880252744967e-01, 2.028841008979334e-01,-8.468644315814380e-01,-4.801936918283287e-01,-2.162492056004726e-01,-7.704649664738020e-01,-3.898780668455788e-01,-9.732218982337848e-01, 3.101592620596687e-01,-8.401130268256528e-01,-7.587042371149260e-01, 2.042970932723433e-01,-9.597313293543706e-01,-7.633141983579799e-01, 9.190092199783506e-02,-5.324156538842341e-01,-9.692334606205502e-01,-3.791473546039409e-02,-5.583164813732576e-01,-5.451029263909702e-01,-7.995004912557284e-01,-7.397458104759506e-01,-3.672764754873818e-01,-1.477440514886856e-01,-7.439273446422778e-01, 4.341532062019434e-01,-8.507478973455769e-01,-8.816060112431570e-01,-8.576780070370036e-01, 7.039700978984957e-01,-8.286368529395126e-01,-9.533416421732388e-01,-7.167787839002487e-02,-5.351640020872266e-01,-9.823447747059300e-01,-8.775290251116272e-01,-6.281991826324685e-01,-9.614945926988339e-01,-4.155034713787327e-01,-9.512801798601274e-01,-8.029791970793696e-01,-9.507695663607698e-01,-1.909952628734674e-01,-5.821770063190647e-01,-8.104747619419612e-01,-8.759635613946177e-01,-8.399099562227985e-01,-4.075544098241169e-01,-9.816558398598512e-01,-6.751330775192570e-01,-3.149360040996975e-01,-9.527447173711141e-01,-9.642973907015797e-01, 3.882799450590653e-01,-4.054196893657805e-01,-2.393278735475145e-01,-9.501368719194545e-01, 1.213317675400288e-01,-9.634372956535556e-01,-4.582195921420686e-01, 1.497455178272146e-01,-9.649408756348019e-01,-2.457872011647745e-01,-3.331460493589239e-02,-9.795878526898599e-01,-9.400848917408163e-01, 5.351913398839814e-02,-5.771247308970763e-01,-5.272439333234028e-01,-9.908408885142427e-01,-9.542078132340630e-01, 9.092369697334675e-01,-9.487762812978140e-01, 3.912186294449625e-01,-6.733380907456182e-01,-1.942367507430980e-01,-7.507934072992886e-01,-1.073689634700975e-01,-5.367046095654924e-01,-7.316400271438055e-01, 1.509721024051824e-02,-9.312283446813046e-01,-9.530099739999303e-01, 6.860650921566884e-01,-9.451599502645341e-01,-6.947721730584699e-01,-7.633656849146985e-01, 3.796790207398015e-01,-5.940915240375572e-01,-8.337496877757204e-01,-2.659802125250012e-01,-9.424883056920581e-01,-2.178037054639842e-02,-9.913843356785719e-01,-8.017801206802969e-01, 4.884770153263237e-01,-7.647093103679242e-01,-9.927319239202493e-01, 7.178088783690493e-01,-9.928858092454618e-01,-6.612014489340987e-01, 5.879776998244531e-01,-4.187164755261538e-01,-9.301990868248645e-01,-9.956181752033222e-01,-8.543417457662290e-01,-5.317102318519026e-01, 3.850028351761021e-01,-7.646313205708825e-01,-9.893278459606969e-01,-3.350060800467011e-01,-3.893408313278781e-01,-5.035752568875519e-01};
static const dfloat cubS12[122] = { 4.494310271080339e-02,-1.374503437769800e-01,-9.291591635595342e-01,-7.956608618114207e-01,-9.422439932245108e-01, 2.430423677648021e-01,-1.786133718196977e-01,-8.208949817859824e-01,-7.300121447810836e-01,-8.468644315814318e-01,-4.798880252744660e-01, 2.028841008979130e-01,-7.704649664737874e-01,-4.801936918283312e-01,-2.162492056004864e-01, 3.101592620596832e-01,-3.898780668455887e-01,-9.732218982337789e-01, 2.042970932723419e-01,-8.401130268256584e-01,-7.587042371149224e-01, 9.190092199784240e-02,-9.597313293543740e-01,-7.633141983579822e-01,-3.791473546039571e-02,-5.324156538842368e-01,-9.692334606205503e-01,-7.995004912557280e-01,-5.583164813732471e-01,-5.451029263909755e-01,-1.477440514886944e-01,-7.397458104759530e-01,-3.672764754873972e-01,-8.507478973455825e-01,-7.439273446422684e-01, 4.341532062019415e-01, 7.039700978984880e-01,-8.816060112431580e-01,-8.576780070370065e-01,-7.167787839004998e-02,-8.286368529395153e-01,-9.533416421732391e-01,-8.775290251116095e-01,-5.351640020872259e-01,-9.823447747059338e-01,-4.155034713787540e-01,-6.281991826324663e-01,-9.614945926988320e-01,-9.507695663607681e-01,-9.512801798601279e-01,-8.029791970793906e-01,-8.104747619419593e-01,-1.909952628734482e-01,-5.821770063190658e-01,-4.075544098241037e-01,-8.759635613946225e-01,-8.399099562227862e-01,-3.149360040996973e-01,-9.816558398598504e-01,-6.751330775192623e-01, 3.882799450590572e-01,-9.527447173711093e-01,-9.642973907015800e-01,-9.501368719194521e-01,-4.054196893657869e-01,-2.393278735475080e-01,-4.582195921420612e-01, 1.213317675400408e-01,-9.634372956535532e-01,-2.457872011647652e-01, 1.497455178272235e-01,-9.649408756347968e-01,-9.400848917408172e-01,-3.331460493590697e-02,-9.795878526898629e-01,-5.272439333234166e-01, 5.351913398837858e-02,-5.771247308970848e-01, 9.092369697334570e-01,-9.908408885142483e-01,-9.542078132340662e-01,-6.733380907456206e-01,-9.487762812978144e-01, 3.912186294449516e-01,-1.073689634701154e-01,-1.942367507431101e-01,-7.507934072993036e-01, 1.509721024050988e-02,-5.367046095654978e-01,-7.316400271438089e-01, 6.860650921566785e-01,-9.312283446813061e-01,-9.530099739999284e-01,-7.633656849146908e-01,-9.451599502645367e-01,-6.947721730584682e-01,-8.337496877757223e-01, 3.796790207397898e-01,-5.940915240375597e-01,-2.178037054640109e-02,-2.659802125249980e-01,-9.424883056920595e-01, 4.884770153263270e-01,-9.913843356785754e-01,-8.017801206802900e-01, 7.178088783690437e-01,-7.647093103679297e-01,-9.927319239202509e-01, 5.879776998244501e-01,-9.928858092454629e-01,-6.612014489341050e-01,-9.956181752033096e-01,-4.187164755261620e-01,-9.301990868248643e-01, 3.850028351760949e-01,-8.543417457662337e-01,-5.317102318518995e-01,-7.646313205708875e-01,-9.893278459607124e-01,-3.350060800467284e-01,-3.893408313278872e-01,-5.035752568875422e-01};
static const dfloat cubT12[122] = {-9.783335953742984e-01,-9.783335953742984e-01,-9.783335953743082e-01,-5.051375127288734e-01,-5.051375127288813e-01,-5.051375127288567e-01,-2.704795016132610e-01,-2.704795016132427e-01,-2.704795016132353e-01,-8.761316440420067e-01,-8.761316440420128e-01,-8.761316440420128e-01,-5.330921360973955e-01,-5.330921360973915e-01,-5.330921360973797e-01,-9.470592969803197e-01,-9.470592969803197e-01,-9.470592969803185e-01,-6.054798293317695e-01,-6.054798293317559e-01,-6.054798293317541e-01,-3.688553942854864e-01,-3.688553942854803e-01,-3.688553942854779e-01,-4.604361500348142e-01,-4.604361500348209e-01,-4.604361500348155e-01,-9.708010098004605e-02,-9.708010098005836e-02,-9.708010098004360e-02,-7.452336625479586e-01,-7.452336625479734e-01,-7.452336625479695e-01,-8.394779642140980e-01,-8.394779642141005e-01,-8.394779642141005e-01,-9.646860796183304e-01,-9.646860796183329e-01,-9.646860796183329e-01,-1.463436264971988e-01,-1.463436264972110e-01,-1.463436264972110e-01, 3.950378019047747e-01, 3.950378019047796e-01, 3.950378019047723e-01, 5.197246710055239e-03, 5.197246710044225e-03, 5.197246710045449e-03, 7.050289433002711e-01, 7.050289433002650e-01, 7.050289433002906e-01,-4.163529688654979e-01,-4.163529688655126e-01,-4.163529688654997e-01, 1.234279274415084e-01, 1.234279274415292e-01, 1.234279274415279e-01,-2.827507852118435e-02,-2.827507852117946e-02,-2.827507852118680e-02,-4.712378369863621e-01,-4.712378369863733e-01,-4.712378369863746e-01,-4.051155651672539e-01,-4.051155651672461e-01,-4.051155651672552e-01,-6.996748797444104e-01,-6.996748797444178e-01,-6.996748797444153e-01,-9.390174410276441e-01,-9.390174410276441e-01,-9.390174410276404e-01,-4.701265063342596e-02,-4.701265063341494e-02,-4.701265063342840e-02,-9.491504697678957e-01,-9.491504697678982e-01,-9.491504697678996e-01,-9.641882679851499e-01,-9.641882679851536e-01,-9.641882679851548e-01,-7.691042574015277e-01,-7.691042574015287e-01,-7.691042574015213e-01,-9.476008784874947e-01,-9.476008784874947e-01,-9.476008784874957e-01,-7.467525735312169e-01,-7.467525735312133e-01,-7.467525735312096e-01,-8.018267734754447e-01,-8.018267734754533e-01,-8.018267734754545e-01, 4.032978082376928e-01, 4.032978082376916e-01, 4.032978082376990e-01,-9.518378089265165e-01,-9.518378089265154e-01,-9.518378089265154e-01,-7.697511112365405e-01,-7.697511112365454e-01,-7.697511112365416e-01,-6.953125589674611e-01,-6.953125589674550e-01,-6.953125589674575e-01,-9.603676440808753e-01,-9.603676440808777e-01,-9.603676440808777e-01,-9.338904416448865e-01,-9.338904416448853e-01,-9.338904416448853e-01, 3.445337375543313e-01, 3.445337375543459e-01, 3.445337375543349e-01,-9.989508575579646e-01,-9.989508575579671e-01,-9.989508575579706e-01, 2.938939617126564e-01, 9.679835378821102e-01,-9.949817598598508e-01,-8.319775060163606e-01,-4.892742293373437e-01};
static const dfloat cubW12[122] = { 2.781735057653978e-03, 2.781735057654092e-03, 2.781735057652762e-03, 8.398769485703079e-03, 8.398769485703389e-03, 8.398769485703702e-03, 1.751926923808389e-02, 1.751926923808319e-02, 1.751926923808276e-02, 1.144500673832101e-02, 1.144500673832079e-02, 1.144500673832059e-02, 2.070223832089109e-02, 2.070223832088912e-02, 2.070223832088700e-02, 3.755144648977120e-03, 3.755144648976951e-03, 3.755144648977461e-03, 1.585064692214937e-02, 1.585064692214994e-02, 1.585064692214965e-02, 9.837382530578264e-03, 9.837382530577457e-03, 9.837382530577046e-03, 1.040938129049506e-02, 1.040938129049636e-02, 1.040938129049624e-02, 2.601141444106337e-02, 2.601141444106351e-02, 2.601141444106238e-02, 2.254462784978858e-02, 2.254462784978603e-02, 2.254462784978717e-02, 1.158500980094870e-02, 1.158500980094966e-02, 1.158500980094936e-02, 4.168974300566039e-03, 4.168974300565897e-03, 4.168974300565954e-03, 1.028842110893551e-02, 1.028842110893561e-02, 1.028842110893553e-02, 4.205298136714887e-03, 4.205298136715934e-03, 4.205298136715269e-03, 1.216104814495457e-02, 1.216104814495443e-02, 1.216104814495468e-02, 4.217582910822220e-03, 4.217582910822276e-03, 4.217582910821895e-03, 2.504450458257267e-02, 2.504450458257196e-02, 2.504450458257253e-02, 1.540237379264268e-02, 1.540237379264239e-02, 1.540237379264282e-02, 7.109722872505770e-03, 7.109722872505939e-03, 7.109722872505064e-03, 4.452836223595802e-03, 4.452836223596298e-03, 4.452836223596240e-03, 1.560511141396748e-02, 1.560511141396692e-02, 1.560511141396706e-02, 1.122469707341839e-02, 1.122469707341878e-02, 1.122469707341904e-02, 5.606952565186992e-03, 5.606952565187232e-03, 5.606952565188124e-03, 4.295440580382450e-03, 4.295440580382917e-03, 4.295440580382775e-03, 1.485317436346403e-02, 1.485317436346361e-02, 1.485317436346333e-02, 6.743238395204040e-04, 6.743238395202980e-04, 6.743238395202797e-04, 1.108774447134587e-02, 1.108774447134572e-02, 1.108774447134648e-02, 1.527946705889213e-02, 1.527946705889241e-02, 1.527946705889255e-02, 2.926965483792584e-02, 2.926965483792640e-02, 2.926965483792612e-02, 4.903655570817210e-03, 4.903655570817012e-03, 4.903655570817195e-03, 1.066746136540665e-02, 1.066746136540648e-02, 1.066746136540661e-02, 1.046201140582620e-02, 1.046201140582678e-02, 1.046201140582673e-02, 1.620898349121864e-02, 1.620898349121822e-02, 1.620898349121807e-02, 3.781331883907675e-03, 3.781331883907420e-03, 3.781331883907420e-03, 1.719274762893588e-03, 1.719274762893348e-03, 1.719274762893291e-03, 2.795111755340936e-03, 2.795111755340936e-03, 2.795111755340964e-03, 2.956085825040871e-03, 2.956085825040178e-03, 2.956085825039980e-03, 3.486457943478327e-03, 3.486457943478129e-03, 3.486457943477903e-03, 1.480798009773252e-02, 2.669332563147237e-04, 6.590607255576522e-03, 3.402287040173104e-02, 3.933995650399713e-02};

static const dfloat cubR13[146] = {-9.732571421128717e-01,-8.777094417633484e-01, 6.360523243885940e-01,-8.341709591339117e-01,-8.942891180332886e-01,-1.998041407397467e-01,-9.386582648478321e-01,-8.680617765856911e-01,-2.616421416492865e-02,-8.963403716926903e-01, 5.343670444851365e-01,-6.635447485924715e-01,-3.045068568243497e-02,-3.787172459437460e-02,-9.683387312981971e-01,-2.591148880836474e-01,-1.871316765763450e-01,-5.751164814667691e-01,-4.250587282142239e-01,-3.818279810172570e-01,-9.643478781984359e-01,-9.807820953032897e-01, 1.932519690055645e-03,-4.558451169144688e-01,-8.422980060971172e-01,-5.570417863326089e-01, 6.665546184094190e-02, 5.600514449454901e-01,-9.636995078165799e-01,-6.792382718083578e-01,-4.466624483634102e-01,-5.622004059640262e-01,-7.951382421541775e-01,-9.648753815273471e-01,-6.819315530681420e-01,-3.605709513818484e-02,-9.754549036654413e-01,-6.958032696001026e-01, 3.239668349489298e-01, 8.762706778079965e-02,-6.344367164934905e-01,-8.453627023133832e-01,-8.272048561870711e-01,-8.763851254913780e-01, 5.592441857943067e-01,-9.821739752750656e-01,-9.608620724515218e-01, 8.825622987044712e-01,-9.421912729503297e-01,-1.191341671759996e-01,-6.528395283098720e-01,-9.545982997764215e-01,-9.577406806012169e-01,-8.192963331930044e-01,-9.570233745958590e-01,-3.587854022699037e-01,-1.096610415336468e-01, 2.349732966499787e-01,-5.393659894285496e-01,-8.395229976753453e-01,-8.109423756048937e-01,-6.496446495980878e-01, 3.192749999106354e-01, 3.825698293210655e-01,-9.585509414775775e-01,-4.358122649164298e-01,-9.648449161180356e-01,-3.658619438859911e-01, 1.691044838712373e-01,-2.821139978456132e-01,-3.723169573141468e-01,-8.067307637292832e-01,-8.219486660185765e-01,-3.227915674738000e-01, 9.423053270234544e-02,-3.004482562658674e-01,-8.419492219637048e-01, 1.223709283249110e-01,-9.327737093395267e-01,-8.416881016492157e-01, 7.576374936718195e-01,-9.370963799046415e-01,-9.776898685073001e-01, 6.217453302896785e-01,-6.876169316746629e-01,-9.759681789639764e-01, 4.071202836332542e-01,-5.275153882760935e-01,-5.980787986510294e-01, 6.427195366309985e-02,-2.181407851211924e-01,-2.024851137176574e-01,-7.545018742660632e-01,-8.447241584640430e-01,-9.808572706033634e-01, 2.137909346010534e-01,-6.708617375695154e-01,-4.771603842649455e-01,-9.290809718499888e-01,-8.691646621839529e-01,-8.491305015738382e-01, 2.769494379663832e-01,-7.993602927334967e-01,-8.686786982101000e-01,-5.131728755180972e-01,-8.370884347828356e-01,-9.827937094740478e-01, 7.689841364496782e-01, 1.086475398267162e-02,-2.578131788008660e-01,-9.468553183080501e-01,-7.657159432089945e-01,-7.387495169446457e-01,-9.696856835049740e-01,-7.783009954008024e-01,-6.747640155947847e-01,-1.995744100753914e-01,-9.504514061880363e-01,-9.850355556452155e-01,-1.071706000748585e-01,-5.861552747394648e-01,-5.783256911891492e-01,-1.163888222392147e-01,-9.851710780454848e-01, 3.910577496982700e-01,-4.623142964636475e-01,-9.457501779316387e-01,-9.556721239442023e-01,-5.281837494500825e-01,-9.948116694696253e-01,-8.483115481646346e-01,-3.164138515358303e-01,-7.583425904617824e-01,-6.647725606680912e-01, 4.172862178469927e-01,-9.915229809889643e-01,-9.415148665544931e-01, 3.212250088359117e-01,-7.250992104065092e-01,-9.979391240452457e-01,-3.987004611045780e-01,-8.097659985419874e-01,-9.793097854869544e-01,-3.808899633138076e-01,-6.662734259168809e-01,-4.988822261209018e-01};
static const dfloat cubS13[146] = { 6.360523243888241e-01,-9.732571421129482e-01,-8.777094417634068e-01,-1.998041407402172e-01,-8.341709591339200e-01,-8.942891180333804e-01,-2.616421416513242e-02,-9.386582648478748e-01,-8.680617765856307e-01,-6.635447485925710e-01,-8.963403716927343e-01, 5.343670444850736e-01,-9.683387312982029e-01,-3.045068568251808e-02,-3.787172459441498e-02,-5.751164814668782e-01,-2.591148880838861e-01,-1.871316765765378e-01,-9.643478781984429e-01,-4.250587282142386e-01,-3.818279810172600e-01,-4.558451169144027e-01,-9.807820953033244e-01, 1.932519689933552e-03, 6.665546184083668e-02,-8.422980060971939e-01,-5.570417863327513e-01,-6.792382718085200e-01, 5.600514449454066e-01,-9.636995078165774e-01,-7.951382421540705e-01,-4.466624483632516e-01,-5.622004059640536e-01,-3.605709513812018e-02,-9.648753815273329e-01,-6.819315530681272e-01, 3.239668349490811e-01,-9.754549036654717e-01,-6.958032696002108e-01,-8.453627023133928e-01, 8.762706778065467e-02,-6.344367164936563e-01, 5.592441857940212e-01,-8.272048561872356e-01,-8.763851254914318e-01, 8.825622987046170e-01,-9.821739752750327e-01,-9.608620724514576e-01,-6.528395283098161e-01,-9.421912729502248e-01,-1.191341671758371e-01,-8.192963331928110e-01,-9.545982997764149e-01,-9.577406806012398e-01,-1.096610415335899e-01,-9.570233745958655e-01,-3.587854022699748e-01,-8.395229976754708e-01, 2.349732966498296e-01,-5.393659894287345e-01, 3.192749999104236e-01,-8.109423756049685e-01,-6.496446495981585e-01,-4.358122649165939e-01, 3.825698293210474e-01,-9.585509414775819e-01, 1.691044838711478e-01,-9.648449161180522e-01,-3.658619438861107e-01,-8.067307637292715e-01,-2.821139978457123e-01,-3.723169573141921e-01, 9.423053270227139e-02,-8.219486660186259e-01,-3.227915674737588e-01, 1.223709283248278e-01,-3.004482562660676e-01,-8.419492219637827e-01, 7.576374936719125e-01,-9.327737093395108e-01,-8.416881016491287e-01, 6.217453302896253e-01,-9.370963799046809e-01,-9.776898685073103e-01, 4.071202836331183e-01,-6.876169316746668e-01,-9.759681789639969e-01, 6.427195366303383e-02,-5.275153882762045e-01,-5.980787986511875e-01,-7.545018742660826e-01,-2.181407851212919e-01,-2.024851137177368e-01, 2.137909346010763e-01,-8.447241584640276e-01,-9.808572706034018e-01,-9.290809718499947e-01,-6.708617375695477e-01,-4.771603842648413e-01, 2.769494379662987e-01,-8.691646621840267e-01,-8.491305015739224e-01,-5.131728755184237e-01,-7.993602927334694e-01,-8.686786982101975e-01, 7.689841364495766e-01,-8.370884347827919e-01,-9.827937094740744e-01,-9.468553183080493e-01, 1.086475398262476e-02,-2.578131788009341e-01,-9.696856835050209e-01,-7.657159432090259e-01,-7.387495169446614e-01,-1.995744100753992e-01,-7.783009954008309e-01,-6.747640155948592e-01,-1.071706000747267e-01,-9.504514061880432e-01,-9.850355556451971e-01,-1.163888222392223e-01,-5.861552747395741e-01,-5.783256911892162e-01,-4.623142964636497e-01,-9.851710780455284e-01, 3.910577496982430e-01,-5.281837494499044e-01,-9.457501779316365e-01,-9.556721239441894e-01,-3.164138515359573e-01,-9.948116694695791e-01,-8.483115481646086e-01, 4.172862178471399e-01,-7.583425904617823e-01,-6.647725606680003e-01, 3.212250088360542e-01,-9.915229809889501e-01,-9.415148665544906e-01,-3.987004611046259e-01,-7.250992104065455e-01,-9.979391240453227e-01,-8.097659985422143e-01,-9.793097854868257e-01,-3.808899633138974e-01,-6.662734259166581e-01,-4.988822261208352e-01};
static const dfloat cubT13[146] = {-7.850857405124128e-01,-7.850857405123810e-01,-7.850857405123014e-01,-7.173578209259969e-02,-7.173578209280664e-02,-7.173578209275028e-02,-1.671157444014210e-01,-1.671157444013193e-01,-1.671157444016121e-01,-9.744819241998242e-01,-9.744819241998829e-01,-9.744819241998707e-01,-9.633388584248760e-01,-9.633388584248747e-01,-9.633388584248969e-01,-9.786369538728992e-01,-9.786369538728601e-01,-9.786369538728601e-01,-2.287654125700825e-01,-2.287654125700593e-01,-2.287654125700678e-01,-5.653053074721941e-01,-5.653053074721549e-01,-5.653053074721600e-01,-6.673156694110401e-01,-6.673156694110401e-01,-6.673156694110806e-01,-9.171136653203762e-01,-9.171136653203908e-01,-9.171136653204471e-01,-1.959989035185533e-01,-1.959989035185520e-01,-1.959989035186035e-01,-3.171359702664373e-01,-3.171359702663871e-01,-3.171359702663613e-01,-6.527086616834269e-01,-6.527086616833829e-01,-6.527086616832848e-01,-6.078276489738966e-01,-6.078276489738356e-01,-6.078276489736553e-01,-8.556542041156293e-01,-8.556542041155865e-01,-8.556542041157358e-01,-9.395262509780005e-01,-9.395262509779857e-01,-9.395262509779869e-01,-2.858350315639968e-01,-2.858350315639539e-01,-2.858350315640151e-01, 7.316353135704496e-01, 7.316353135703885e-01, 7.316353135706825e-01,-5.745301816005585e-01,-5.745301816005292e-01,-5.745301816005276e-01,-8.560843095458567e-01,-8.560843095458724e-01,-8.560843095457389e-01,-8.586879747074871e-01,-8.586879747074100e-01,-8.586879747075643e-01,-9.882066229268871e-01,-9.882066229268982e-01,-9.882066229269264e-01,-8.383976238671185e-01,-8.383976238671356e-01,-8.383976238670720e-01,-5.388382811109242e-01,-5.388382811108720e-01,-5.388382811109336e-01,-9.494902992099282e-01,-9.494902992099221e-01,-9.494902992098999e-01,-9.799734500951530e-01,-9.799734500951247e-01,-9.799734500951433e-01,-9.831756826831526e-01,-9.831756826831698e-01,-9.831756826832224e-01,-7.069590818776603e-01,-7.069590818777044e-01,-7.069590818776689e-01,-7.435351729944846e-01,-7.435351729945029e-01,-7.435351729945116e-01,-9.386777667357825e-01,-9.386777667357629e-01,-9.386777667357223e-01,-8.248722268950576e-01,-8.248722268949975e-01,-8.248722268949988e-01,-3.882095055336683e-01,-3.882095055337357e-01,-3.882095055336431e-01, 7.710309368438679e-02, 7.710309368446273e-02, 7.710309368439903e-02,-5.586542742084717e-01,-5.586542742084810e-01,-5.586542742085601e-01, 1.812118664621056e-01, 1.812118664620359e-01, 1.812118664616538e-01,-9.491019921927196e-01,-9.491019921927479e-01,-9.491019921927601e-01,-8.061962568736880e-01,-8.061962568737149e-01,-8.061962568737394e-01, 4.741511436586612e-01, 4.741511436587114e-01, 4.741511436586735e-01,-3.473605789290306e-01,-3.473605789290257e-01,-3.473605789289559e-01, 4.265756190794634e-02, 4.265756190798427e-02, 4.265756190813004e-02,-7.191302118321772e-01,-7.191302118320853e-01,-7.191302118320804e-01,-9.435723751890991e-01,-9.435723751890734e-01,-9.435723751890759e-01, 4.296060513257585e-01, 4.296060513257070e-01, 4.296060513259066e-01, 1.595370691702228e-01, 1.595370691701504e-01, 1.595370691700034e-01,-9.941710667172281e-01,-9.941710667171791e-01,-9.941710667172869e-01,-3.881871612926001e-01,-3.881871612925659e-01,-3.881871612924909e-01, 1.217387955564122e-01, 1.217387955564991e-01, 1.217387955564306e-01, 4.292979956264287e-01, 9.379293564605812e-01,-8.573301100584970e-01,-1.179722249622360e-03,-5.033533216373864e-01};
static const dfloat cubW13[146] = { 2.840980725841578e-03, 2.840980725838184e-03, 2.840980725840815e-03, 9.422995106205527e-03, 9.422995106213770e-03, 9.422995106190889e-03, 6.515806867772963e-03, 6.515806867763020e-03, 6.515806867761803e-03, 3.624190270983565e-03, 3.624190270978883e-03, 3.624190270979492e-03, 3.323481837533718e-03, 3.323481837532035e-03, 3.323481837533619e-03, 7.453651529601763e-03, 7.453651529609273e-03, 7.453651529608312e-03, 1.011108061748424e-02, 1.011108061748533e-02, 1.011108061748594e-02, 6.347638435625232e-03, 6.347638435621045e-03, 6.347638435624043e-03, 1.528366691775590e-02, 1.528366691775081e-02, 1.528366691776114e-02, 4.146858859962316e-03, 4.146858859963236e-03, 4.146858859962613e-03, 2.056909001662421e-02, 2.056909001661827e-02, 2.056909001662252e-02, 9.022163389018846e-03, 9.022163389022340e-03, 9.022163389022085e-03, 6.184839575665225e-03, 6.184839575659540e-03, 6.184839575665112e-03, 1.701957325506264e-02, 1.701957325507297e-02, 1.701957325506590e-02, 7.370232899285110e-03, 7.370232899281885e-03, 7.370232899279707e-03, 8.673226162974020e-04, 8.673226162988997e-04, 8.673226163001441e-04, 1.218846751209406e-02, 1.218846751211666e-02, 1.218846751210588e-02, 2.969389984180547e-03, 2.969389984181013e-03, 2.969389984178114e-03, 1.105996179111814e-02, 1.105996179111763e-02, 1.105996179111948e-02, 1.242710471565888e-02, 1.242710471566187e-02, 1.242710471567059e-02, 1.262849562762590e-02, 1.262849562762402e-02, 1.262849562762580e-02, 1.908277532574823e-03, 1.908277532574413e-03, 1.908277532572475e-03, 7.510407203198469e-03, 7.510407203195513e-03, 7.510407203195993e-03, 2.369350844363481e-02, 2.369350844363312e-02, 2.369350844363128e-02, 9.246862858535350e-03, 9.246862858532974e-03, 9.246862858533895e-03, 5.523697568553409e-03, 5.523697568558994e-03, 5.523697568554398e-03, 1.787802063828094e-03, 1.787802063827245e-03, 1.787802063825718e-03, 2.784520816196351e-03, 2.784520816198076e-03, 2.784520816196492e-03, 6.013899523999122e-03, 6.013899523998826e-03, 6.013899523993691e-03, 1.344172911438606e-02, 1.344172911438769e-02, 1.344172911438968e-02, 2.048092950331958e-02, 2.048092950332282e-02, 2.048092950332268e-02, 5.240726869874686e-03, 5.240726869875294e-03, 5.240726869871207e-03, 1.446990890238765e-02, 1.446990890239090e-02, 1.446990890238935e-02, 1.287823094631799e-02, 1.287823094631513e-02, 1.287823094631803e-02, 1.421162330323071e-02, 1.421162330324019e-02, 1.421162330324924e-02, 1.721299942608860e-03, 1.721299942606993e-03, 1.721299942605169e-03, 1.047888097554055e-02, 1.047888097554181e-02, 1.047888097554007e-02, 6.402096370957237e-03, 6.402096370953929e-03, 6.402096370963122e-03, 2.642693254567321e-02, 2.642693254566501e-02, 2.642693254567590e-02, 2.953156526667655e-03, 2.953156526664727e-03, 2.953156526665066e-03, 2.493408796090700e-02, 2.493408796091125e-02, 2.493408796090715e-02, 2.848671606431755e-03, 2.848671606429393e-03, 2.848671606429973e-03, 5.295388238559288e-03, 5.295388238560391e-03, 5.295388238563516e-03, 3.217049829531768e-03, 3.217049829535106e-03, 3.217049829538189e-03, 3.413588783402041e-03, 3.413588783406071e-03, 3.413588783397247e-03, 2.689053313984700e-03, 2.689053313985096e-03, 2.689053313986850e-03, 4.273381649551712e-03, 4.273381649555487e-03, 4.273381649546677e-03, 1.150852063876846e-02, 3.307459566891739e-04, 2.267424349796184e-02, 1.872018224494073e-02, 3.443152615917819e-02};

static const dfloat cubR14[177] = {-8.567523903984056e-01, 2.937321767696985e-01,-5.692613549698415e-01,-5.193473480512563e-01, 1.412680487794378e-01,-8.319937625402691e-01,-9.708892474738398e-01,-9.593489098771848e-01, 4.440149956802136e-01,-9.749418803839179e-01,-3.684585382998631e-01, 9.951912643713334e-02,-2.308758062699288e-01,-8.300439114686674e-01, 3.502705143586354e-02,-2.829109456728378e-01,-7.457066474604837e-01,-6.644950594335102e-01,-5.505346489626760e-01,-7.902266429023022e-01, 3.067760377339515e-01,-8.215286306248032e-01, 2.308477306581685e-02,-3.871904504425462e-01,-9.854112892373914e-01,-5.863008773397973e-01, 8.095953613083051e-02,-8.052559327293417e-01,-7.985311971756309e-01,-3.820306143645372e-01,-2.019721875456377e-01,-9.823300269139608e-01,-9.815052727635920e-01,-9.569583489308782e-01,-9.678392999746168e-01, 1.665560319874374e-01,-8.623336520822803e-01,-9.645317122756909e-01, 7.979301210773273e-01,-9.945047520567687e-01,-3.276980727669719e-01,-2.307085927903063e-01, 6.079556793213509e-02,-5.326462754985429e-01,-5.452872262463172e-01,-5.383521084278374e-02,-8.432515189485366e-01,-2.476363726136155e-01, 1.044364208690578e-01,-8.160190642126887e-01,-2.982800071098369e-01, 3.225964073564471e-01,-5.414443953322645e-01,-9.731776707058488e-01,-9.554207271113230e-01,-9.619077322874070e-01, 6.687677811518801e-01,-9.846649841209297e-01,-3.987370345109075e-01,-4.774806610095652e-01,-5.679869113087561e-01, 1.435556152935313e-02,-8.753919064423815e-01,-2.114016862036898e-01,-9.655172118841630e-01,-5.214829017076167e-02,-9.623380148719576e-01,-9.608693011846281e-01,-8.406253388749148e-01,-9.605000973076959e-01, 6.116223916452550e-01,-8.028459184018204e-01,-9.609541747446682e-01,-6.013382886926650e-01, 5.419800458524395e-01,-6.433130187729085e-01,-6.902954537605610e-01,-9.470225601038479e-02,-9.079006849812619e-01,-2.062635049164693e-01,-9.184813465187585e-01,-8.917774773820670e-01, 6.822193673149513e-02,-5.886734938918742e-01, 3.004201958278954e-01,-8.003296766514031e-01,-5.529337075291529e-01,-8.457093705404592e-01,-6.135265713434932e-01,-8.572698174976345e-01,-9.611864986280678e-01,-6.019932601376445e-01, 4.123065091794484e-01,-8.118099987139180e-01,-8.023473755211947e-01, 5.755357972483064e-01,-5.921536175986672e-01,-9.597575687639389e-01, 5.231101115221591e-01,-9.072484447659781e-01,-2.757518889867298e-01,-2.994770969636995e-01,-9.233890657799098e-01,-1.713576967924220e-01,-6.779908687722540e-01,-4.296978649684835e-01,-4.923548227148208e-01,-8.059174914812960e-01,-7.971189807552446e-01, 3.909606769611357e-01,-9.624516025049893e-01, 3.513898553375933e-01,-8.084228945618219e-01,-9.635359245664379e-01,-4.129022988444946e-01,-8.047600052292334e-01,-9.835803988600644e-01, 7.089094525142836e-03,-5.806062924564268e-01,-5.846458394348144e-01,-2.353981070162650e-01,-9.664022120806307e-01, 1.566894708545493e-01,-9.613552545036474e-01,-8.472233495894415e-01, 6.958108952006254e-01,-7.770896833229748e-01,-9.793952947663495e-01,-7.690278404598587e-01,-6.449042333727806e-01,-2.200647816871925e-01,-9.410792865416311e-01,-9.235556209423426e-01,-5.931991810169928e-01,-6.089556124846054e-01,-9.650342420780642e-01,-2.414487369555118e-01, 1.603915248028213e-01,-9.917797371775030e-01, 1.138438287833665e-02,-8.175856314288508e-01, 9.053183821168509e-01,-9.814952923071467e-01,-9.904223104721384e-01,-2.160132996322185e-01,-6.171636437984327e-01,-2.280396216261566e-01,-8.641058019347291e-01,-8.389926351662794e-01, 9.587125236990479e-02,-3.354078332290266e-01,-3.010696584863675e-01,-6.923879278393141e-01,-7.969656816595116e-01,-7.958731300763120e-01, 3.602938355273282e-01,-7.919304509816861e-01,-4.251550483406301e-01,-9.887971646106728e-01,-5.741149257080890e-01,-9.574385047232588e-01,-9.558614116822212e-01,-9.734489361183414e-01,-8.771711114944296e-01, 8.465999303170888e-01,-8.319822467713663e-01, 5.291032756854286e-02,-9.935203813401260e-01,-5.301187750231947e-01,-9.961517088434864e-01, 3.128792709078354e-02,-3.941637613412473e-01,-9.862845926370784e-01,-6.682249594899097e-01,-8.443681536636525e-01,-3.365520947720209e-01,-5.135579361166598e-01};
static const dfloat cubS14[177] = {-5.692613549714406e-01,-8.567523903996869e-01, 2.937321767690869e-01,-8.319937625396620e-01,-5.193473480562776e-01, 1.412680487798826e-01, 4.440149956816773e-01,-9.708892474731111e-01,-9.593489098770824e-01, 9.951912644081932e-02,-9.749418803837925e-01,-3.684585382974073e-01, 3.502705143446270e-02,-2.308758062704557e-01,-8.300439114681392e-01,-6.644950594364012e-01,-2.829109456709421e-01,-7.457066474618167e-01, 3.067760377326236e-01,-5.505346489628263e-01,-7.902266429009204e-01,-3.871904504420789e-01,-8.215286306227646e-01, 2.308477306517287e-02, 8.095953613145350e-02,-9.854112892368967e-01,-5.863008773374992e-01,-3.820306143652253e-01,-8.052559327273110e-01,-7.985311971787540e-01,-9.815052727643349e-01,-2.019721875491939e-01,-9.823300269138270e-01, 1.665560319844528e-01,-9.569583489309840e-01,-9.678392999737244e-01, 7.979301210771497e-01,-8.623336520833290e-01,-9.645317122753843e-01,-2.307085927852067e-01,-9.945047520572895e-01,-3.276980727638241e-01,-5.452872262449433e-01, 6.079556793294856e-02,-5.326462754990520e-01,-2.476363726174116e-01,-5.383521084314242e-02,-8.432515189493912e-01,-2.982800071080640e-01, 1.044364208685466e-01,-8.160190642139122e-01,-9.731776707060575e-01, 3.225964073566463e-01,-5.414443953300053e-01, 6.687677811551758e-01,-9.554207271112243e-01,-9.619077322873879e-01,-4.774806610081035e-01,-9.846649841206343e-01,-3.987370345104631e-01,-8.753919064418382e-01,-5.679869113075616e-01, 1.435556153035758e-02,-5.214829016948568e-02,-2.114016861997780e-01,-9.655172118841225e-01,-8.406253388768438e-01,-9.623380148726018e-01,-9.608693011840216e-01,-8.028459184008299e-01,-9.605000973077215e-01, 6.116223916439572e-01, 5.419800458525613e-01,-9.609541747446029e-01,-6.013382886931526e-01,-9.470225600927738e-02,-6.433130187737914e-01,-6.902954537588261e-01,-9.184813465183521e-01,-9.079006849819555e-01,-2.062635049177409e-01,-5.886734938916006e-01,-8.917774773814240e-01, 6.822193673246273e-02,-5.529337075266429e-01, 3.004201958305422e-01,-8.003296766508911e-01,-8.572698174958018e-01,-8.457093705411737e-01,-6.135265713416332e-01, 4.123065091803993e-01,-9.611864986280264e-01,-6.019932601368023e-01, 5.755357972487538e-01,-8.118099987142103e-01,-8.023473755205793e-01, 5.231101115219017e-01,-5.921536175995283e-01,-9.597575687639091e-01,-2.994770969640125e-01,-9.072484447659134e-01,-2.757518889843976e-01,-6.779908687730015e-01,-9.233890657804160e-01,-1.713576967945722e-01,-8.059174914815999e-01,-4.296978649692322e-01,-4.923548227152426e-01,-9.624516025052400e-01,-7.971189807552889e-01, 3.909606769611551e-01,-9.635359245662354e-01, 3.513898553388829e-01,-8.084228945607991e-01,-9.835803988600607e-01,-4.129022988459327e-01,-8.047600052299576e-01,-5.846458394339499e-01, 7.089094525170311e-03,-5.806062924561544e-01, 1.566894708543361e-01,-2.353981070165389e-01,-9.664022120805861e-01, 6.958108951988492e-01,-9.613552545036335e-01,-8.472233495900784e-01,-7.690278404605405e-01,-7.770896833226055e-01,-9.793952947658595e-01,-9.410792865416145e-01,-6.449042333713295e-01,-2.200647816864868e-01,-6.089556124842277e-01,-9.235556209417903e-01,-5.931991810168001e-01, 1.603915248026716e-01,-9.650342420782030e-01,-2.414487369559702e-01,-8.175856314282485e-01,-9.917797371769650e-01, 1.138438287667271e-02,-9.904223104714033e-01, 9.053183821180540e-01,-9.814952923078711e-01,-2.280396216249042e-01,-2.160132996332856e-01,-6.171636437980160e-01, 9.587125237040388e-02,-8.641058019354901e-01,-8.389926351650465e-01,-6.923879278393541e-01,-3.354078332285020e-01,-3.010696584868394e-01, 3.602938355270157e-01,-7.969656816600444e-01,-7.958731300766900e-01,-9.887971646102336e-01,-7.919304509810771e-01,-4.251550483392493e-01,-9.558614116826611e-01,-5.741149257106475e-01,-9.574385047228751e-01, 8.465999303149507e-01,-9.734489361185846e-01,-8.771711114950046e-01,-9.935203813408606e-01,-8.319822467700391e-01, 5.291032756969321e-02, 3.128792709372000e-02,-5.301187750200479e-01,-9.961517088435393e-01,-3.941637613413972e-01,-9.862845926380797e-01,-6.682249594898159e-01,-8.443681536623480e-01,-3.365520947699970e-01,-5.135579361169699e-01};
static const dfloat cubT14[177] = {-8.677184314001082e-01,-8.677184314013502e-01,-8.677184314014397e-01,-7.899269381888595e-01,-7.899269381834939e-01,-7.899269381852465e-01,-5.137768383309290e-01,-5.137768383287411e-01,-5.137768383292588e-01,-7.561187077554274e-01,-7.561187077568812e-01,-7.561187077559455e-01,-9.741073336963032e-01,-9.741073336962174e-01,-9.741073336967880e-01,-3.068873474294131e-01,-3.068873474319239e-01,-3.068873474323990e-01,-9.660147458686635e-01,-9.660147458680035e-01,-9.660147458683426e-01,-8.143656919989729e-01,-8.143656920003850e-01,-8.143656919991664e-01,-5.092473695541364e-01,-5.092473695566246e-01,-5.092473695560799e-01,-1.418225572800657e-02,-1.418225573432502e-02,-1.418225572608982e-02, 1.658074872242221e-01, 1.658074872273402e-01, 1.658074872337371e-01,-2.417583830790724e-01,-2.417583830712365e-01,-2.417583830824184e-01,-9.710647567190727e-01,-9.710647567187028e-01,-9.710647567195663e-01,-4.470885823916320e-01,-4.470885823924011e-01,-4.470885823880192e-01,-9.828620661897888e-01,-9.828620661897289e-01,-9.828620661896076e-01,-8.552768975909628e-01,-8.552768975906884e-01,-8.552768975901996e-01,-9.901373495478656e-01,-9.901373495471454e-01,-9.901373495460957e-01,-8.079743413184145e-01,-8.079743413186852e-01,-8.079743413191628e-01,-7.514393217564428e-01,-7.514393217582320e-01,-7.514393217532622e-01,-1.391173203618309e-01,-1.391173203605498e-01,-1.391173203594304e-01,-5.709767437822296e-01,-5.709767437793426e-01,-5.709767437801685e-01,-7.709328117426426e-01,-7.709328117429819e-01,-7.709328117421687e-01, 7.638326549330181e-01, 7.638326549376111e-01, 7.638326549308504e-01,-8.482763759340863e-01,-8.482763759353786e-01,-8.482763759345421e-01,-9.796875824148232e-01,-9.796875824146811e-01,-9.796875824145904e-01,-5.716892714592209e-01,-5.716892714564207e-01,-5.716892714580704e-01, 3.264553641769556e-02, 3.264553641816582e-02, 3.264553641853452e-02,-5.877709654590555e-01,-5.877709654604684e-01,-5.877709654585925e-01,-9.471568116500158e-01,-9.471568116501859e-01,-9.471568116504013e-01, 3.165057593730721e-01, 3.165057593814518e-01, 3.165057593799808e-01,-8.491267504154898e-01,-8.491267504155167e-01,-8.491267504146484e-01,-9.613784230132565e-01,-9.613784230132503e-01,-9.613784230132515e-01,-9.711989251592792e-01,-9.711989251593954e-01,-9.711989251594053e-01,-5.175225692848496e-01,-5.175225692843031e-01,-5.175225692851814e-01,-2.272623686536590e-01,-2.272623686539199e-01,-2.272623686533577e-01,-2.720298208326884e-01,-2.720298208350962e-01,-2.720298208332739e-01,-6.313900937014347e-01,-6.313900937008885e-01,-6.313900937012632e-01,-5.794310362105732e-01,-5.794310362107408e-01,-5.794310362108394e-01, 2.012427029337758e-01, 2.012427029353985e-01, 2.012427029338565e-01,-8.418369626352739e-01,-8.418369626355140e-01,-8.418369626340908e-01,-9.548891517575462e-01,-9.548891517576674e-01,-9.548891517574666e-01,-8.872322911059731e-01,-8.872322911052898e-01,-8.872322911067606e-01, 5.255128185486104e-01, 5.255128185486962e-01, 5.255128185480410e-01,-1.939516984007722e-01,-1.939516983998096e-01,-1.939516983993527e-01, 1.257104144430931e-01, 1.257104144429155e-01, 1.257104144429070e-01,-9.539085457686953e-01,-9.539085457687393e-01,-9.539085457687344e-01,-2.020190142730135e-01,-2.020190142699284e-01,-2.020190142710576e-01,-9.334007793387326e-01,-9.334007793395176e-01,-9.334007793399144e-01,-9.387834349447086e-01,-9.387834349445531e-01,-9.387834349447185e-01,-3.927728152701174e-01,-3.927728152699886e-01,-3.927728152709291e-01,-6.711345804459720e-01,-6.711345804455909e-01,-6.711345804457723e-01,-7.674550237911505e-01,-7.674550237908994e-01,-7.674550237914602e-01, 2.058826639314879e-01, 2.058826639322951e-01, 2.058826639318492e-01, 4.874148421139309e-01, 4.874148421164918e-01, 4.874148421202836e-01,-9.959798827017677e-01,-9.959798827003800e-01,-9.959798827032251e-01,-2.274076994581263e-01,-2.274076994588807e-01,-2.274076994575580e-01,-5.050174432265367e-01,-5.050174432283195e-01,-5.050174432250477e-01,-8.175087159753148e-01, 9.588537779156814e-01, 4.674878470335996e-03, 5.331044609837571e-01,-9.903437156880256e-01,-4.593261916505825e-01};
static const dfloat cubW14[177] = { 5.760248903305241e-03, 5.760248903325903e-03, 5.760248903365261e-03, 8.932868689521610e-03, 8.932868689511656e-03, 8.932868689476301e-03, 2.076049498090074e-03, 2.076049498174885e-03, 2.076049498042698e-03, 4.911437360394771e-03, 4.911437360421683e-03, 4.911437360432558e-03, 4.632630603019636e-03, 4.632630603023865e-03, 4.632630602996896e-03, 1.478616415369785e-02, 1.478616415372260e-02, 1.478616415358613e-02, 5.161692502570879e-03, 5.161692502629611e-03, 5.161692502599545e-03, 1.140073505154521e-02, 1.140073505147236e-02, 1.140073505157061e-02, 4.040956960491586e-03, 4.040956960533503e-03, 4.040956960535977e-03, 1.214382348234666e-02, 1.214382348228025e-02, 1.214382348235584e-02, 1.340088420032186e-03, 1.340088420047175e-03, 1.340088420071598e-03, 2.752950520631146e-03, 2.752950520621246e-03, 2.752950520637114e-03, 1.331401247764424e-03, 1.331401247780917e-03, 1.331401247756442e-03, 2.684423028380161e-03, 2.684423028363544e-03, 2.684423028306820e-03, 5.069562772386248e-03, 5.069562772383122e-03, 5.069562772402412e-03, 1.094533055932680e-02, 1.094533055925776e-02, 1.094533055915119e-02, 2.717780254308358e-03, 2.717780254386791e-03, 2.717780254453796e-03, 5.079032159261387e-03, 5.079032159297902e-03, 5.079032159282812e-03, 2.316987739291757e-03, 2.316987739296989e-03, 2.316987739314384e-03, 5.266486904775362e-03, 5.266486904824209e-03, 5.266486904782348e-03, 1.395506591225304e-02, 1.395506591215987e-02, 1.395506591221747e-02, 7.232545731546133e-03, 7.232545731564914e-03, 7.232545731587993e-03, 2.098806854227512e-03, 2.098806854174692e-03, 2.098806854244964e-03, 4.095418260649108e-03, 4.095418260645926e-03, 4.095418260643578e-03, 2.013844392669503e-03, 2.013844392679473e-03, 2.013844392681934e-03, 2.006292561902968e-02, 2.006292561903576e-02, 2.006292561909558e-02, 7.636349045050322e-03, 7.636349044960237e-03, 7.636349045015122e-03, 1.339120189395926e-02, 1.339120189406734e-02, 1.339120189406358e-02, 7.281319055054554e-03, 7.281319054993121e-03, 7.281319054909172e-03, 9.799335256329167e-03, 9.799335256279062e-03, 9.799335256278751e-03, 5.441700720819044e-03, 5.441700720816442e-03, 5.441700720874283e-03, 4.832027392359751e-03, 4.832027392367500e-03, 4.832027392366722e-03, 2.623461512808826e-03, 2.623461512803495e-03, 2.623461512803665e-03, 1.530232186139421e-02, 1.530232186146322e-02, 1.530232186140284e-02, 1.282421163940995e-02, 1.282421163935870e-02, 1.282421163937584e-02, 1.960238059385917e-02, 1.960238059391475e-02, 1.960238059397032e-02, 6.461243759462251e-03, 6.461243759546609e-03, 6.461243759527885e-03, 6.574681086137257e-03, 6.574681086085157e-03, 6.574681086080772e-03, 4.619145926356672e-03, 4.619145926334002e-03, 4.619145926347734e-03, 1.763943232158073e-02, 1.763943232157564e-02, 1.763943232166799e-02, 3.780109952877583e-03, 3.780109952860302e-03, 3.780109952877427e-03, 3.494590792272494e-03, 3.494590792242060e-03, 3.494590792284953e-03, 4.191790076637947e-03, 4.191790076513666e-03, 4.191790076570659e-03, 1.217415884392229e-02, 1.217415884388091e-02, 1.217415884391409e-02, 1.239994032767911e-02, 1.239994032767526e-02, 1.239994032768864e-02, 4.107419260796993e-03, 4.107419260774153e-03, 4.107419260786400e-03, 3.036323024481491e-03, 3.036323024512674e-03, 3.036323024549628e-03, 4.083230290834423e-04, 4.083230290736291e-04, 4.083230290554777e-04, 1.259565293606578e-02, 1.259565293606549e-02, 1.259565293608774e-02, 1.343412553142801e-02, 1.343412553139515e-02, 1.343412553145282e-02, 2.386756526414242e-02, 2.386756526415388e-02, 2.386756526415473e-02, 1.341142697054699e-02, 1.341142697054656e-02, 1.341142697055976e-02, 4.126697319202634e-03, 4.126697319174223e-03, 4.126697319148965e-03, 3.824606804156526e-03, 3.824606804148607e-03, 3.824606804194654e-03, 5.154031158903054e-04, 5.154031159091442e-04, 5.154031158631511e-04, 3.008884311995277e-03, 3.008884312101371e-03, 3.008884312053670e-03, 3.475550074368887e-03, 3.475550074426672e-03, 3.475550074439499e-03, 1.985270760196320e-02, 1.768406287405760e-04, 2.162569384673413e-02, 6.346297860541651e-03, 3.834603804842544e-03, 2.542527774527081e-02};

static const dfloat cubR15[214] = {-4.643883053659482e-01,-9.806925443236449e-01,-1.956932081750188e-01,-4.653464282631503e-01,-8.149337537832231e-01,-4.517628670350377e-02,-4.884081161872981e-01, 4.882356091005899e-02,-5.830844011837657e-01,-8.537213192645825e-01,-8.644349360961447e-01,-6.912514924034796e-01,-9.046811078786205e-01, 7.242451190904076e-01,-9.697631291508340e-01,-5.294284984115514e-01,-9.793001867986513e-01, 4.853831669745099e-01,-9.823417544632196e-01, 1.123723224050577e-01,-8.562808212160420e-01,-9.030284534793107e-01,-8.554526246331162e-01,-4.522022018995108e-01,-9.788306410188258e-02, 4.995437228294657e-02,-9.726554719749202e-01,-9.811851305179435e-01,-3.875842243781824e-01,-6.286608571589433e-01,-9.818041165754183e-01,-9.733219200647203e-01,-2.562832364766001e-01,-9.419117214263670e-01,-9.503164941544398e-01,-9.823039070211337e-01,-9.702381351695488e-01,-8.624671096921959e-01, 5.412036808840410e-01,-8.623412224408559e-01,-9.352617785668350e-01, 3.334698653996355e-01,-3.139388520443420e-01,-3.391089213808151e-01,-8.178199640831652e-01,-2.343477353567101e-01,-2.420509859922935e-01,-5.488898624223500e-01,-3.324304545273839e-01,-3.254598517074725e-01,-6.494234734897651e-01,-1.526783886694892e-01,-2.587710909624500e-01,-9.679538407320751e-01,-7.427837756926913e-01,-8.758623466224407e-01,-3.952261469503152e-01,-9.725790978827638e-01,-6.391093116702669e-01,-9.635979517435672e-02,-7.808489049037033e-01,-9.707195613500970e-01,-9.744570392779863e-01,-4.035227593725479e-01,-3.014031249231439e-01,-9.516274946739665e-01,-6.548972678664049e-01,-6.132923621309210e-01,-3.138664551037936e-01,-9.546022538653359e-02,-6.798164731284529e-01,-9.560289916888933e-01,-6.088046555182679e-01,-4.331865389582874e-01,-9.677867047429870e-01,-6.557745244765812e-01,-6.168827753406187e-01,-5.228043384933917e-02,-9.773491934323073e-01,-8.310805670240774e-01, 2.162461842795751e-01,-8.327094277904781e-01,-8.020801730068206e-02,-6.544711989541597e-01,-7.062305631011258e-01,-6.302726376643720e-01, 2.999032710991442e-01,-9.669781313762058e-01,-9.724792423631858e-01, 4.830320594159364e-01, 2.850680325148247e-01,-9.395331951982926e-01,-7.010033354706134e-01,-9.012961985134473e-01,-3.721507544955594e-01,-2.600352041533796e-02,-8.275151966459791e-01,-6.948959476233765e-01,-3.422511127474528e-01,-9.728630503483611e-01,-9.644622821619812e-01, 1.192942880143547e-01,-4.202945993671708e-01,-5.281878639376225e-01,-8.297194275561455e-01,-6.103564448832565e-01,-8.139907302831556e-01, 2.648172082209886e-01,-6.082110736395772e-01,-8.573137339002826e-01, 4.337977760414746e-01,-9.754970037176441e-01,-7.783784126016599e-01, 7.162639199129448e-01,-9.731209597556876e-01,-6.220433556990723e-01, 4.417289788115173e-01,-8.424518463630635e-01,-8.299165880486181e-01,-9.625609846688696e-01,-9.843203335838651e-01,-5.796441119303599e-01,-7.985509492426148e-01,-9.700181365782471e-01,-8.405481417066596e-01,-1.936469840431067e-01,-8.602202507982447e-01,-6.136151417127925e-01,-4.075088764786616e-02,-5.522186453777402e-01,-9.501484282830269e-01, 1.226968646946504e-01,-7.832941949220155e-01,-8.268740275259179e-01, 2.583717741907658e-01,-7.194071680947965e-01,-3.165025561289078e-01,-1.003014508800699e-01,-6.694990509019798e-01,-6.461076941893704e-01,-8.981825530809362e-01,-4.788819486190393e-01,-6.154755688380665e-01,-3.602959727174426e-02,-7.747128044924225e-01, 8.492219612122176e-02,-3.371830044757640e-01, 4.118361381028008e-01,-9.038252917125962e-01,-5.420754838488659e-01,-9.682347712499504e-01,-2.273086239226234e-01,-8.696174120649703e-01,-8.168862545018932e-01,-5.782851509401810e-01, 2.215129797696634e-01,-2.466935687505188e-01, 5.421054226985865e-02,-8.378958752522317e-01,-5.691582757189227e-01,-9.643242866696536e-01, 3.775943212697339e-01,-8.659878997763408e-01,-8.510616249988688e-01,-8.698563628160738e-03,-9.569288671106377e-01,-9.610278135124122e-01,-5.557071327165761e-01,-8.210674250927418e-01,-9.641091452966988e-01, 5.663803816093957e-01, 1.733733570675947e-01,-2.599944640090449e-01,-9.517880688934234e-01,-1.041698906825117e-02,-2.663877770650462e-01,-8.827477102198451e-01,-8.521120631804805e-01,-8.205558156903667e-01, 5.508934755698189e-01,-9.601757867817463e-01,-9.587545875355591e-01, 8.867584170146786e-01,-8.508177150690071e-01,-8.131586061133455e-01, 6.557556977303685e-01,-9.618045345702616e-01,-2.867610747004070e-01, 1.215870279553673e-01, 2.716752503710635e-01,-7.276837572332995e-01,-9.961559530711328e-01, 7.206156381491000e-01,-7.986661217522745e-01,-9.655962539865030e-01,-9.936730692955971e-01, 7.703324161118421e-01,-9.611662894911912e-01, 4.196958364782455e-01,-9.938526304411064e-01,-4.414964783738483e-01,-9.939816946579979e-01,-8.102525057291776e-01,-5.725005009030300e-01, 8.328958951698762e-02,-9.977370129880159e-01,-4.735718614332430e-01, 5.744698714694374e-02,-2.341648955521002e-01,-9.999971439746738e-01,-3.773645077780259e-01,-6.745431359676690e-01,-5.593856381682771e-01,-4.562031817459618e-01};
static const dfloat cubS15[214] = {-1.956932081750901e-01,-4.643883053659112e-01,-9.806925443236437e-01,-4.517628670354546e-02,-4.653464282631740e-01,-8.149337537831614e-01,-5.830844011837819e-01,-4.884081161873131e-01, 4.882356091000834e-02,-6.912514924034424e-01,-8.537213192645864e-01,-8.644349360961118e-01,-9.697631291508356e-01,-9.046811078786193e-01, 7.242451190904010e-01, 4.853831669744998e-01,-5.294284984115346e-01,-9.793001867986568e-01,-8.562808212160142e-01,-9.823417544632069e-01, 1.123723224050740e-01,-4.522022018994907e-01,-9.030284534793108e-01,-8.554526246331079e-01,-9.726554719749126e-01,-9.788306410190260e-02, 4.995437228295948e-02,-6.286608571589389e-01,-9.811851305179350e-01,-3.875842243781447e-01,-2.562832364765835e-01,-9.818041165754111e-01,-9.733219200647211e-01,-9.823039070211251e-01,-9.419117214263756e-01,-9.503164941544310e-01, 5.412036808840583e-01,-9.702381351695494e-01,-8.624671096921911e-01, 3.334698653996048e-01,-8.623412224408524e-01,-9.352617785668299e-01,-8.178199640831895e-01,-3.139388520443826e-01,-3.391089213807830e-01,-5.488898624223250e-01,-2.343477353567271e-01,-2.420509859923017e-01,-6.494234734897472e-01,-3.324304545273779e-01,-3.254598517075228e-01,-9.679538407320714e-01,-1.526783886695001e-01,-2.587710909624583e-01,-3.952261469503215e-01,-7.427837756926992e-01,-8.758623466224110e-01,-9.635979517436784e-02,-9.725790978827612e-01,-6.391093116702623e-01,-9.744570392779863e-01,-7.808489049037070e-01,-9.707195613500912e-01,-9.516274946739771e-01,-4.035227593725730e-01,-3.014031249231665e-01,-3.138664551038582e-01,-6.548972678663051e-01,-6.132923621310278e-01,-9.560289916889030e-01,-9.546022538645867e-02,-6.798164731283882e-01,-9.677867047429876e-01,-6.088046555182602e-01,-4.331865389582947e-01,-5.228043384936472e-02,-6.557745244765563e-01,-6.168827753406025e-01, 2.162461842795877e-01,-9.773491934323004e-01,-8.310805670240743e-01,-6.544711989541513e-01,-8.327094277905337e-01,-8.020801730069999e-02, 2.999032710991041e-01,-7.062305631011467e-01,-6.302726376643931e-01, 4.830320594159287e-01,-9.669781313762033e-01,-9.724792423631912e-01,-7.010033354705972e-01, 2.850680325148471e-01,-9.395331951982968e-01,-2.600352041530780e-02,-9.012961985134290e-01,-3.721507544955477e-01,-3.422511127474248e-01,-8.275151966459768e-01,-6.948959476233220e-01, 1.192942880143550e-01,-9.728630503483602e-01,-9.644622821619796e-01,-8.297194275561447e-01,-4.202945993671151e-01,-5.281878639376197e-01, 2.648172082209769e-01,-6.103564448832770e-01,-8.139907302831462e-01, 4.337977760414659e-01,-6.082110736395807e-01,-8.573137339002923e-01, 7.162639199129395e-01,-9.754970037176485e-01,-7.783784126016626e-01, 4.417289788115011e-01,-9.731209597556929e-01,-6.220433556990796e-01,-9.625609846688702e-01,-8.424518463630619e-01,-8.299165880486141e-01,-7.985509492426254e-01,-9.843203335838525e-01,-5.796441119303628e-01,-1.936469840431012e-01,-9.700181365782451e-01,-8.405481417066629e-01,-4.075088764784493e-02,-8.602202507982122e-01,-6.136151417127935e-01, 1.226968646946293e-01,-5.522186453777678e-01,-9.501484282830229e-01, 2.583717741907725e-01,-7.832941949220492e-01,-8.268740275259213e-01,-1.003014508800369e-01,-7.194071680947886e-01,-3.165025561288871e-01,-8.981825530809275e-01,-6.694990509019569e-01,-6.461076941893930e-01,-3.602959727173487e-02,-4.788819486190639e-01,-6.154755688380403e-01,-3.371830044757918e-01,-7.747128044924128e-01, 8.492219612122189e-02,-5.420754838488806e-01, 4.118361381027958e-01,-9.038252917126203e-01,-8.696174120649623e-01,-9.682347712499534e-01,-2.273086239226363e-01, 2.215129797696696e-01,-8.168862545018911e-01,-5.782851509401727e-01,-8.378958752522282e-01,-2.466935687505224e-01, 5.421054226984534e-02, 3.775943212697254e-01,-5.691582757189002e-01,-9.643242866696480e-01,-8.698563628129792e-03,-8.659878997763302e-01,-8.510616249988442e-01,-5.557071327165682e-01,-9.569288671106436e-01,-9.610278135124100e-01, 5.663803816094106e-01,-8.210674250927283e-01,-9.641091452967021e-01,-9.517880688934220e-01, 1.733733570675994e-01,-2.599944640090250e-01,-8.827477102198569e-01,-1.041698906823088e-02,-2.663877770650632e-01, 5.508934755698234e-01,-8.521120631804924e-01,-8.205558156903746e-01, 8.867584170146727e-01,-9.601757867817466e-01,-9.587545875355588e-01, 6.557556977303718e-01,-8.508177150690146e-01,-8.131586061133439e-01, 1.215870279553617e-01,-9.618045345702638e-01,-2.867610747004197e-01,-9.961559530711124e-01, 2.716752503711084e-01,-7.276837572333040e-01,-9.655962539865013e-01, 7.206156381491010e-01,-7.986661217522719e-01,-9.611662894911879e-01,-9.936730692955931e-01, 7.703324161118531e-01,-4.414964783738899e-01, 4.196958364782518e-01,-9.938526304411267e-01,-5.725005009030275e-01,-9.939816946579966e-01,-8.102525057291788e-01,-4.735718614332526e-01, 8.328958951699263e-02,-9.977370129880007e-01,-9.999971439746830e-01, 5.744698714694076e-02,-2.341648955520927e-01,-3.773645077780107e-01,-6.745431359676521e-01,-5.593856381682558e-01,-4.562031817458799e-01};
static const dfloat cubT15[214] = {-3.592259421353274e-01,-3.592259421353629e-01,-3.592259421353948e-01,-6.745435312501461e-01,-6.745435312500934e-01,-6.745435312501350e-01,-9.773310435389496e-01,-9.773310435389447e-01,-9.773310435389521e-01, 4.094077477641385e-01, 4.094077477641581e-01, 4.094077477641520e-01,-8.498008820609469e-01,-8.498008820609518e-01,-8.498008820609360e-01,-9.766544817643010e-01,-9.766544817643134e-01,-9.766544817643071e-01,-2.737497467258065e-01,-2.737497467258163e-01,-2.737497467258102e-01, 2.106832800118864e-01, 2.106832800119096e-01, 2.106832800119195e-01,-9.794158362061375e-01,-9.794158362061326e-01,-9.794158362061410e-01,-2.569787944988005e-03,-2.569787944940278e-03,-2.569787944970872e-03, 2.114092731167252e-01, 2.114092731167276e-01, 2.114092731167472e-01, 8.745321226019229e-01, 8.745321226019106e-01, 8.745321226019229e-01,-7.084984360223180e-01,-7.084984360223070e-01,-7.084984360223009e-01,-5.358668643919273e-01,-5.358668643919274e-01,-5.358668643919191e-01,-5.291322624916689e-01,-5.291322624916613e-01,-5.291322624916905e-01,-9.747114162286691e-01,-9.747114162286677e-01,-9.747114162286642e-01,-6.926862202753598e-01,-6.926862202753991e-01,-6.926862202753991e-01,-6.205966796360076e-01,-6.205966796359800e-01,-6.205966796359994e-01, 1.387226926542950e-02, 1.387226926545037e-02, 1.387226926545527e-02,-2.919517952726106e-01,-2.919517952726388e-01,-2.919517952726143e-01, 7.260255055317858e-01, 7.260255055317980e-01, 7.260255055317980e-01,-3.434466210302873e-01,-3.434466210303130e-01,-3.434466210303326e-01,-4.179439148988016e-01,-4.179439148987066e-01,-4.179439148988169e-01,-2.686943097961418e-01,-2.686943097961932e-01,-2.686943097962190e-01, 9.777899219516327e-03, 9.777899219537199e-03, 9.777899219524893e-03,-6.750622663334673e-01,-6.750622663334280e-01,-6.750622663334610e-01,-4.078164238232052e-01,-4.078164238231988e-01,-4.078164238231999e-01,-4.326113559546607e-01,-4.326113559546422e-01,-4.326113559546595e-01,-9.634000703336254e-01,-9.634000703336254e-01,-9.634000703336254e-01,-5.435746856765372e-01,-5.435746856765189e-01,-5.435746856765354e-01,-6.445315018459322e-01,-6.445315018459200e-01,-6.445315018459200e-01,-7.005495265756838e-01,-7.005495265756863e-01,-7.005495265756765e-01,-1.353377429832427e-01,-1.353377429832280e-01,-1.353377429832464e-01,-1.819689555040141e-01,-1.819689555040055e-01,-1.819689555040116e-01,-2.217981091390599e-01,-2.217981091390978e-01,-2.217981091390794e-01,-8.404700330545780e-01,-8.404700330545438e-01,-8.404700330545770e-01,-9.682729685016156e-01,-9.682729685016083e-01,-9.682729685016169e-01,-9.623885035936355e-01,-9.623885035936379e-01,-9.623885035936391e-01,-8.465646633567417e-01,-8.465646633567429e-01,-8.465646633567491e-01, 6.349294190805568e-01, 6.349294190805519e-01, 6.349294190805556e-01, 3.625153947568439e-01, 3.625153947568402e-01, 3.625153947568280e-01, 4.213262328010296e-03, 4.213262328004177e-03, 4.213262328010296e-03,-4.854137198411078e-01,-4.854137198410846e-01,-4.854137198410982e-01,-6.203297910338710e-01,-6.203297910338782e-01,-6.203297910338957e-01,-6.482035517428409e-01,-6.482035517428127e-01,-6.482035517428103e-01,-8.637888248962520e-01,-8.637888248962520e-01,-8.637888248962471e-01, 2.137892981722631e-01, 2.137892981722398e-01, 2.137892981722814e-01,-8.696128852711745e-01,-8.696128852711710e-01,-8.696128852711672e-01,-9.730263871530412e-01,-9.730263871530387e-01,-9.730263871530399e-01,-9.659353625413271e-01,-9.659353625413256e-01,-9.659353625413256e-01, 6.516080723754117e-02, 6.516080723756197e-02, 6.516080723755341e-02,-8.263415743275953e-01,-8.263415743275868e-01,-8.263415743276001e-01,-9.696210982671161e-01,-9.696210982671148e-01,-9.696210982671161e-01,-8.441117588811560e-01,-8.441117588811744e-01,-8.441117588811657e-01,-2.742519115966724e-01,-2.742519115966700e-01,-2.742519115966724e-01, 4.736638133396158e-01, 4.736638133396402e-01, 4.736638133396244e-01,-7.812038112199684e-01,-7.812038112199684e-01,-7.812038112199574e-01,-9.615908241651252e-01,-9.615908241651276e-01,-9.615908241651241e-01,-8.404475236468387e-01,-8.404475236468447e-01,-8.404475236468387e-01,-8.782255966989736e-01,-8.782255966989724e-01,-8.782255966989613e-01,-9.678280426973700e-01,-9.678280426973690e-01,-9.678280426973724e-01,-9.917793765480273e-01,-9.917793765480284e-01,-9.917793765480211e-01,-8.730214186846984e-01,-8.730214186846921e-01,-8.730214186846847e-01,-5.478355400666467e-01,-5.478355400666755e-01,-5.478355400666706e-01,-9.563532624103227e-01,-9.563532624103251e-01,-9.563532624103240e-01,-8.154930573250420e-01,-8.154930573250592e-01,-8.154930573250556e-01,-9.843467276632510e-01,-9.843467276632510e-01,-9.843467276632608e-01, 3.767347012902068e-01, 3.767347012902080e-01, 3.767347012901908e-01,-6.119807150957207e-01,-6.119807150957191e-01,-6.119807150957195e-01,-8.232849476201566e-01,-8.232849476201444e-01,-8.232849476201517e-01,-8.679064766659740e-01, 2.362940790295578e-02,-3.218430854948759e-01,-6.313904547620893e-01};
static const dfloat cubW15[214] = { 3.522723551354820e-03, 3.522723551352486e-03, 3.522723551352938e-03, 9.232955535331875e-03, 9.232955535327506e-03, 9.232955535330659e-03, 4.237026901463632e-03, 4.237026901464212e-03, 4.237026901463377e-03, 6.106499343749692e-03, 6.106499343748844e-03, 6.106499343750059e-03, 1.627360858046573e-03, 1.627360858046587e-03, 1.627360858046813e-03, 1.148548912222280e-03, 1.148548912221887e-03, 1.148548912221821e-03, 2.564399625003663e-03, 2.564399625004229e-03, 2.564399625002914e-03, 5.856076670044469e-03, 5.856076670043464e-03, 5.856076670043295e-03, 1.423643710751563e-03, 1.423643710751888e-03, 1.423643710751167e-03, 3.978120910726397e-03, 3.978120910728038e-03, 3.978120910727048e-03, 1.363809538497426e-03, 1.363809538497801e-03, 1.363809538497208e-03, 4.075959956903368e-04, 4.075959956909096e-04, 4.075959956900526e-04, 2.860089500953389e-03, 2.860089500953686e-03, 2.860089500953841e-03, 4.850559515934299e-03, 4.850559515933083e-03, 4.850559515932247e-03, 1.352971840215884e-02, 1.352971840216170e-02, 1.352971840215909e-02, 6.132094336152466e-03, 6.132094336152919e-03, 6.132094336153004e-03, 1.461675478856591e-02, 1.461675478855672e-02, 1.461675478856139e-02, 5.990023547122631e-03, 5.990023547123226e-03, 5.990023547122334e-03, 9.631452064974138e-03, 9.631452064972767e-03, 9.631452064974280e-03, 5.614990995964819e-03, 5.614990995964790e-03, 5.614990995963984e-03, 1.144090371383849e-03, 1.144090371383791e-03, 1.144090371383924e-03, 7.816866183700298e-03, 7.816866183699294e-03, 7.816866183700680e-03, 1.639566856148552e-02, 1.639566856149825e-02, 1.639566856148552e-02, 6.721523996645333e-03, 6.721523996646124e-03, 6.721523996646986e-03, 6.261479412226208e-03, 6.261479412226463e-03, 6.261479412227057e-03, 1.527689367197913e-02, 1.527689367198168e-02, 1.527689367197715e-02, 4.042885330071795e-03, 4.042885330072389e-03, 4.042885330072658e-03, 1.229463521022889e-02, 1.229463521022621e-02, 1.229463521022872e-02, 5.427737795462986e-03, 5.427737795463169e-03, 5.427737795462180e-03, 1.976199256754921e-03, 1.976199256754936e-03, 1.976199256754398e-03, 6.919752984901465e-03, 6.919752984901423e-03, 6.919752984901507e-03, 1.046027736959284e-02, 1.046027736959218e-02, 1.046027736959269e-02, 1.447481953654798e-02, 1.447481953654770e-02, 1.447481953654855e-02, 2.391550565412188e-03, 2.391550565412258e-03, 2.391550565412131e-03, 1.511375252393240e-02, 1.511375252392830e-02, 1.511375252393155e-02, 9.511095825311733e-03, 9.511095825312935e-03, 9.511095825312142e-03, 4.036967691998535e-03, 4.036967691999214e-03, 4.036967691998252e-03, 1.508260348095485e-03, 1.508260348095230e-03, 1.508260348095372e-03, 3.937346811175521e-03, 3.937346811175041e-03, 3.937346811175210e-03, 3.608140871956763e-03, 3.608140871956792e-03, 3.608140871956891e-03, 3.365526474882623e-03, 3.365526474883641e-03, 3.365526474884108e-03, 5.192349870771271e-03, 5.192349870771441e-03, 5.192349870772318e-03, 1.334407891882406e-02, 1.334407891882760e-02, 1.334407891882411e-02, 8.025809123555732e-03, 8.025809123554318e-03, 8.025809123554217e-03, 1.070887983557521e-02, 1.070887983557467e-02, 1.070887983557649e-02, 1.309736406162602e-02, 1.309736406162684e-02, 1.309736406162661e-02, 8.885655025694786e-03, 8.885655025694404e-03, 8.885655025693740e-03, 1.344113214741205e-02, 1.344113214741205e-02, 1.344113214741270e-02, 5.698568685004529e-03, 5.698568685005180e-03, 5.698568685005053e-03, 3.847402885188293e-03, 3.847402885188307e-03, 3.847402885188477e-03, 4.769112898361405e-03, 4.769112898361265e-03, 4.769112898362141e-03, 1.118686556204748e-02, 1.118686556204881e-02, 1.118686556204810e-02, 5.635998890136026e-03, 5.635998890136154e-03, 5.635998890136069e-03, 5.044085510105154e-03, 5.044085510104333e-03, 5.044085510104730e-03, 1.014203848072567e-02, 1.014203848072546e-02, 1.014203848072654e-02, 2.681590717335358e-03, 2.681590717334821e-03, 2.681590717335202e-03, 4.059907280598535e-03, 4.059907280598903e-03, 4.059907280598535e-03, 3.336448036565914e-03, 3.336448036565589e-03, 3.336448036566042e-03, 1.076375181138509e-02, 1.076375181138425e-02, 1.076375181138517e-02, 6.573170212491274e-03, 6.573170212491260e-03, 6.573170212491358e-03, 8.877430748418416e-04, 8.877430748418500e-04, 8.877430748417468e-04, 1.503300080640668e-03, 1.503300080640583e-03, 1.503300080640908e-03, 5.744672404119040e-03, 5.744672404119181e-03, 5.744672404119252e-03, 1.906496203051725e-03, 1.906496203050283e-03, 1.906496203049802e-03, 2.124665847459818e-03, 2.124665847459747e-03, 2.124665847459747e-03, 7.203086774524117e-04, 7.203086774524287e-04, 7.203086774520808e-04, 6.473262420273394e-04, 6.473262420272856e-04, 6.473262420267072e-04, 2.381499975344257e-03, 2.381499975344469e-03, 2.381499975345459e-03, 3.013274913110685e-03, 3.013274913110755e-03, 3.013274913111859e-03, 2.067897521108355e-03, 2.067897521107776e-03, 2.067897521108949e-03, 1.499773420159665e-02, 1.925991642486046e-02, 1.254797412755794e-02, 1.209377740627947e-02};


void mesh_t::CubatureNodesTet3D(int cubTetN, int *_cubNp, dfloat **_cubr, dfloat **_cubs, dfloat **_cubt, dfloat **_cubw){

  if (cubTetN>15)
    LIBP_ABORT(string("Requested Cubature order unavailable."))

  int cubTetNp = cubTetNps[cubTetN-1];

  *_cubNp = cubTetNp;

  *_cubr = (dfloat*) calloc(cubTetNp, sizeof(dfloat));
  *_cubs = (dfloat*) calloc(cubTetNp, sizeof(dfloat));
  *_cubt = (dfloat*) calloc(cubTetNp, sizeof(dfloat));
  *_cubw = (dfloat*) calloc(cubTetNp, sizeof(dfloat));

  const dfloat *cubTetR=NULL, *cubTetS=NULL, *cubTetT=NULL, *cubTetW=NULL;
  switch(cubTetN){
    case 1:  cubTetR = cubR1;  cubTetS = cubS1;  cubTetT = cubT1;  cubTetW = cubW1; break;
    case 2:  cubTetR = cubR2;  cubTetS = cubS2;  cubTetT = cubT2;  cubTetW = cubW2; break;
    case 3:  cubTetR = cubR3;  cubTetS = cubS3;  cubTetT = cubT3;  cubTetW = cubW3; break;
    case 4:  cubTetR = cubR4;  cubTetS = cubS4;  cubTetT = cubT4;  cubTetW = cubW4; break;
    case 5:  cubTetR = cubR5;  cubTetS = cubS5;  cubTetT = cubT5;  cubTetW = cubW5; break;
    case 6:  cubTetR = cubR6;  cubTetS = cubS6;  cubTetT = cubT6;  cubTetW = cubW6; break;
    case 7:  cubTetR = cubR7;  cubTetS = cubS7;  cubTetT = cubT7;  cubTetW = cubW7; break;
    case 8:  cubTetR = cubR8;  cubTetS = cubS8;  cubTetT = cubT8;  cubTetW = cubW8; break;
    case 9:  cubTetR = cubR9;  cubTetS = cubS9;  cubTetT = cubT9;  cubTetW = cubW9; break;
    case 10: cubTetR = cubR10; cubTetS = cubS10; cubTetT = cubT10; cubTetW = cubW10; break;
    case 11: cubTetR = cubR11; cubTetS = cubS11; cubTetT = cubT11; cubTetW = cubW11; break;
    case 12: cubTetR = cubR12; cubTetS = cubS12; cubTetT = cubT12; cubTetW = cubW12; break;
    case 13: cubTetR = cubR13; cubTetS = cubS13; cubTetT = cubT13; cubTetW = cubW13; break;
    case 14: cubTetR = cubR14; cubTetS = cubS14; cubTetT = cubT14; cubTetW = cubW14; break;
    case 15: cubTetR = cubR15; cubTetS = cubS15; cubTetT = cubT15; cubTetW = cubW15; break;
    default:
      LIBP_ABORT(string("Requested Cubature order unavailable."))
  }

  for(int n=0;n<cubTetNp;++n){
    _cubr[0][n] = cubTetR[n];
    _cubs[0][n] = cubTetS[n];
    _cubt[0][n] = cubTetT[n];
    _cubw[0][n] = cubTetW[n];
  }
}