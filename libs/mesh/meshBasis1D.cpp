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
// 1D NODES
// ------------------------------------------------------------------------
void mesh_t::Nodes1D(const int _N, memory<dfloat>& _r){
  JacobiGLL(_N, _r); //Gauss-Legendre-Lobatto nodes
}

void mesh_t::EquispacedNodes1D(const int _N, memory<dfloat>& _r){
  int _Nq = _N+1;

  dfloat dr = 2.0/_N;

  _r.malloc(_Nq);
  for (int i=0;i<_Nq;i++) _r[i] = -1.0 + i*dr;
}

// ------------------------------------------------------------------------
// ORTHONORMAL BASIS POLYNOMIALS
// ------------------------------------------------------------------------
void mesh_t::OrthonormalBasis1D(const dfloat a, const int i, dfloat& P){
  P = JacobiP(a,0,0,i); //Legendre Polynomials
}

void mesh_t::GradOrthonormalBasis1D(const dfloat a, const int i, dfloat& Pr){
  Pr = GradJacobiP(a,0,0,i);
}

// ------------------------------------------------------------------------
// 1D VANDERMONDE MATRICES
// ------------------------------------------------------------------------
void mesh_t::Vandermonde1D(const int _N,
                           const memory<dfloat> _r,
                           memory<dfloat>& V){

  const int _Np = (_N+1);
  const int Npoints = _r.length();

  V.malloc(Npoints*_Np);
  for(int n=0; n<Npoints; n++){
    for(int i=0; i<_Np; i++){
      int id = n*_Np+i;
      OrthonormalBasis1D(_r[n], i, V[id]);
    }
  }
}

void mesh_t::GradVandermonde1D(const int _N,
                               const memory<dfloat> _r,
                               memory<dfloat>& Vr){

  const int _Np = (_N+1);
  const int Npoints = _r.length();

  Vr.malloc(Npoints*_Np);
  for(int n=0; n<Npoints; n++){
    for(int i=0; i<_Np; i++){
      int id = n*_Np+i;
      GradOrthonormalBasis1D(_r[n], i, Vr[id]);
    }
  }
}

// ------------------------------------------------------------------------
// 1D OPERATOR MATRICES
// ------------------------------------------------------------------------
void mesh_t::MassMatrix1D(const int _Np,
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

void mesh_t::Dmatrix1D(const int _N,
                       const memory<dfloat> _rIn,
                       const memory<dfloat> _rOut,
                       memory<dfloat>& _Dr){


  const int _Np = _N+1;
  const int NpointsIn  = _rIn.length();
  const int NpointsOut = _rOut.length();

  // need NpointsIn = (_N+1)
  LIBP_ABORT("Invalid Differentiation operator requested.",
             NpointsIn != _N+1);

  memory<dfloat> V;
  memory<dfloat> Vr;
  Vandermonde1D(_N, _rIn, V);
  GradVandermonde1D(_N, _rOut, Vr);

  //D = Vr/V
  _Dr.malloc(NpointsOut*_Np);
  linAlg_t::matrixRightSolve(NpointsOut, _Np, Vr, _Np, _Np, V, _Dr);
}

void mesh_t::InterpolationMatrix1D(const int _N,
                                   const memory<dfloat> _rIn,
                                   const memory<dfloat> _rOut,
                                   memory<dfloat>& I){

  const int _Np = _N+1;
  const int NpointsIn  = _rIn.length();
  const int NpointsOut = _rOut.length();

  // need NpointsIn = (_N+1)
  LIBP_ABORT("Invalid Interplation operator requested.",
             NpointsIn != _N+1);

  memory<dfloat> VIn;
  memory<dfloat> VOut;
  Vandermonde1D(_N, _rIn, VIn);
  Vandermonde1D(_N, _rOut, VOut);

  I.malloc(NpointsIn*NpointsOut);
  linAlg_t::matrixRightSolve(NpointsOut, _Np, VOut,
                             NpointsIn, _Np, VIn, I);
}

void mesh_t::DegreeRaiseMatrix1D(const int Nc, const int Nf,
                                 memory<dfloat>& P){

  memory<dfloat> rc;
  memory<dfloat> rf;
  Nodes1D(Nc, rc);
  Nodes1D(Nf, rf);

  InterpolationMatrix1D(Nc, rc, rf, P);
}

void mesh_t::CubatureWeakDmatrix1D(const int _Nq, const int _cubNq,
                                   const memory<dfloat> _cubProject,
                                   const memory<dfloat> _cubD,
                                   memory<dfloat>& _cubPDT){

  // cubPDT = cubProject*cubD';
  _cubPDT.malloc(_Nq*_cubNq);
  for(int n=0;n<_Nq;++n){
    for(int m=0;m<_cubNq;++m){
      _cubPDT[n*_cubNq+m] = 0.0;
      for(int k=0;k<_cubNq;++k){
        _cubPDT[n*_cubNq+m] += _cubProject[n*_cubNq+k]*_cubD[m*_cubNq+k];
      }
    }
  }
}

// ------------------------------------------------------------------------
// 1D JACOBI POLYNOMIALS
// ------------------------------------------------------------------------
static dfloat mygamma(const dfloat x){
  dfloat lgam = lgamma(x);
  dfloat gam  = signgam*exp(lgam);
  return gam;
}

dfloat mesh_t::JacobiP(const dfloat a, const dfloat alpha,
                       const dfloat beta, const int _N){

  const dfloat ax = a;

  memory<dfloat> P(_N+1);

  // Zero order
  const dfloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*mygamma(1+alpha)*mygamma(1+beta)/mygamma(1+alpha+beta);
  const dfloat p0     = 1.0/sqrt(gamma0);

  if (_N==0){ return p0;}
  P[0] = p0;

  // first order
  const dfloat gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  const dfloat p1     = ((alpha+beta+2)*ax/2 + (alpha-beta)/2)/sqrt(gamma1);
  if (_N==1){ return p1;}

  P[1] = p1;

  /// Repeat value in recurrence.
  dfloat aold = 2/(2+alpha+beta)*sqrt((alpha+1.)*(beta+1.)/(alpha+beta+3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for(int i=1;i<=_N-1;++i){
    dfloat h1 = 2.*i+alpha+beta;
    dfloat anew = 2./(h1+2.)*sqrt( (i+1.)*(i+1.+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    dfloat bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    P[i+1] = 1./anew*( -aold*P[i-1] + (ax-bnew)*P[i]);
    aold =anew;
  }
  return P[_N];
}

dfloat mesh_t::GradJacobiP(const dfloat a, const dfloat alpha,
                           const dfloat beta, const int _N){

  dfloat PNr = 0;

  if(_N>0)
    PNr = sqrt(_N*(_N+alpha+beta+1.))*JacobiP(a, alpha+1.0, beta+1.0, _N-1);

  return PNr;
}

// ------------------------------------------------------------------------
// 1D GAUSS-LEGENDRE-LOBATTO QUADRATURE
// ------------------------------------------------------------------------
void mesh_t::JacobiGLL(const int _N, memory<dfloat>& _x){

  _x.malloc(_N+1);

  _x[0] = -1.;
  _x[_N] =  1.;

  if(_N>1){
    memory<dfloat> wtmp;
    memory<dfloat> xp1 = _x + 1;
    JacobiGQ(1,1, _N-2, xp1, wtmp);
  }
}

void mesh_t::JacobiGLL(const int _N,
                       memory<dfloat>& _x,
                       memory<dfloat>& _w){

  _x.malloc(_N+1);
  _w.malloc(_N+1);

  _x[0] = -1.;
  _x[_N] =  1.;

  if(_N>1){
    memory<dfloat> wtmp;
    memory<dfloat> xp1 = _x + 1;
    JacobiGQ(1,1, _N-2, xp1, wtmp);
  }

  memory<dfloat> V;
  memory<dfloat> _MM;
  Vandermonde1D(_N, _x, V);
  MassMatrix1D(_N+1, V, _MM);

  // use weights from mass lumping
  for(int n=0;n<=_N;++n){
    dfloat res = 0;
    for(int m=0;m<=_N;++m){
      res += _MM[n*(_N+1)+m];
    }
    _w[n] = res;
  }
}

// ------------------------------------------------------------------------
// 1D GAUSS QUADRATURE
// ------------------------------------------------------------------------
void mesh_t::JacobiGQ(const dfloat alpha, const dfloat beta,
                      const int _N,
                      memory<dfloat>& _x,
                      memory<dfloat>& _w){

  // function NGQ = JacobiGQ(alpha,beta,_N, _x, _w)
  // Purpose: Compute the _N'th order Gauss quadrature points, _x,
  //          and weights, _w, associated with the Jacobi
  //          polynomial, of type (alpha,beta) > -1 ( <> -0.5).
  if (_x.length()==0) _x.malloc(_N+1);
  if (_w.length()==0) _w.malloc(_N+1);

  if (_N==0){
    _x[0] = (alpha-beta)/(alpha+beta+2);
    _w[0] = 2;
  }

  // Form symmetric matrix from recurrence.
  memory<dfloat> J((_N+1)*(_N+1), 0.0);
  memory<dfloat> h1(_N+1);

  for(int n=0;n<=_N;++n){
    h1[n] = 2*n+alpha+beta;
  }

  // J = J + J';
  for(int n=0;n<=_N;++n){
    // J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
    J[n*(_N+1)+n]+= -0.5*(alpha*alpha-beta*beta)/((h1[n]+2)*h1[n])*2; // *2 for symm

    //    diag(2./(h1(1:_N)+2).*sqrt((1:_N).*((1:_N)+alpha+beta).*((1:_N)+alpha).*((1:_N)+beta)./(h1(1:_N)+1)./(h1(1:_N)+3)),1);
    if(n<_N){
      J[n*(_N+1)+n+1]   += (2./(h1[n]+2.))*sqrt((n+1)*(n+1+alpha+beta)*(n+1+alpha)*(n+1+beta)/((h1[n]+1)*(h1[n]+3)));
      J[(n+1)*(_N+1)+n] += (2./(h1[n]+2.))*sqrt((n+1)*(n+1+alpha+beta)*(n+1+alpha)*(n+1+beta)/((h1[n]+1)*(h1[n]+3)));
    }
  }

  dfloat eps = 1;
  while(1+eps>1){
    eps = eps/2.;
  }
  // printf("MACHINE PRECISION %e\n", eps);

  if (alpha+beta<10*eps) J[0] = 0;

  // Compute quadrature by eigenvalue solve

  //  [V,D] = eig(J);
  memory<dfloat> WI(_N+1);
  memory<dfloat> VR((_N+1)*(_N+1));

  // _x = diag(D);
  linAlg_t::matrixEigenVectors(_N+1, J, VR, _x, WI);

  //_w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*.gamma(beta+1)/gamma(alpha+beta+1);
  for(int n=0;n<=_N;++n){
    _w[n] = pow(VR[0*(_N+1)+n],2)*(pow(2,alpha+beta+1)/(alpha+beta+1))*mygamma(alpha+1)*mygamma(beta+1)/mygamma(alpha+beta+1);
  }

  // sloppy sort
  for(int n=0;n<=_N;++n){
    for(int m=n+1;m<=_N;++m){
      if(_x[n]>_x[m]){
        dfloat tmpx = _x[m];
        dfloat tmpw = _w[m];
        _x[m] = _x[n];
        _w[m] = _w[n];
        _x[n] = tmpx;
        _w[n] = tmpw;
      }
    }
  }

#if 0
  for(int n=0;n<=_N;++n){
    printf("zgl[%d] = % e, wgl[%d] = % e\n", n, _x[0][n], n, _w[0][n]);
  }
#endif
}

/*
// C0 basis
int meshContinuousVandermonde1D(int _N, int Npoints, dfloat *_r, dfloat **V, dfloat **Vr){

  int _Np = (_N+1);

  *V  = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*_Np, sizeof(dfloat));

  for(int n=0; n<Npoints; n++){

    int sk = 0;
    for(int i=0; i<=_N; i++){
      int id = n*_Np+sk;
      if(i==0){
        V[0][id] = 0.5*(1-_r[n]);
        Vr[0][id] = -0.5;
      }
      else  if(i==1){
        V[0][id] = 0.5*(1+_r[n]);
        Vr[0][id] = +0.5;
      }
      else{
        // 0.25*(1+_r)*(1-_r)*P^{0,0}_{i-2}(_r)
        dfloat P =  meshJacobiP(_r[n], 0, 0, i-2);
        dfloat Pr = meshGradJacobiP(_r[n], 0, 0, i-2);
        V[0][id]  = 0.25*(1+_r[n])*(1-_r[n])*P;
        Vr[0][id] = 0.25*( (-2*_r[n])*P + (1+_r[n])*(1-_r[n])*Pr);
      }

      sk++;
    }
  }

  return _Np;
}
*/

/*
void meshContinuousFilterMatrix1D(int _N, int Nlow, dfloat *_r, dfloat **F){

  dfloat *VC0, *VrC0;
  dfloat *L = (dfloat*) calloc((_N+1)*(_N+1), sizeof(dfloat));
  dfloat *LinvF = (dfloat*) calloc((_N+1)*(_N+1), sizeof(dfloat));

  int _Np = meshContinuousVandermonde1D(_N, _N+1, _r, &VC0, &VrC0);
  //  int _Np = meshVandermonde1D(_N, _N+1, _r, &VC0, &VrC0); use
  printf("CONTINUOUS VANDERMONDE MATRIX: [\n");
  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      printf("% e ", VC0[n*_Np+m]);
    }
    printf("\n");
  }
  printf("\n");

  *F = (dfloat *) calloc(_Np*_Np, sizeof(dfloat));

  for(int n=0;n<=Nlow;++n){
    L[n*(_N+1)+n] = 1;
  }

  matrixRightSolve(_Np, _Np, L, _Np, _Np, VC0, LinvF);

  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      dfloat res = 0;
      printf("% e ", LinvF[n*_Np+m]);
    }
    printf("\n");
  }
  printf("\n");

  printf("FILTER MATRIX: [\n");
  for(int n=0;n<_Np;++n){
    for(int m=0;m<_Np;++m){
      dfloat res = 0;
      for(int i=0;i<_Np;++i){
        res += VC0[n*_Np+i]*LinvF[i*_Np+m];
      }
      F[0][n*_Np+m] = res;
      printf("% e ", res);
    }
    printf("\n");
  }
  printf("\n");

  free(VC0);
  free(VrC0);
  free(L);
  free(LinvF);
}
*/

// ------------------------------------------------------------------------
// 1D INTERPOLATION MATRICES
// ------------------------------------------------------------------------

/*

*/

/*

void meshCubatureWeakDmatrices1D(int _N, int _Np, dfloat *V,
                                 int cubNp, dfloat *cubr, dfloat *cubw,
                                 dfloat **cubDrT, dfloat **cubProject){

  dfloat *cubV, *cubVr;

  meshVandermonde1D(_N, cubNp, cubr, &cubV, &cubVr);

  // cubDrT = V*transpose(cVr)*diag(cubw);
  // cubProject = V*cV'*diag(cubw); %% relies on (transpose(cV)*diag(cubw)*cV being the identity)

  for(int n=0;n<cubNp;++n){
    for(int m=0;m<_Np;++m){
      // scale by cubw
      cubVr[n*_Np+m] *= cubw[n];
      cubV[n*_Np+m]  *= cubw[n];
    }
  }

  *cubDrT = (dfloat*) calloc(cubNp*_Np, sizeof(dfloat));
  *cubProject = (dfloat*) calloc(cubNp*_Np, sizeof(dfloat));

  for(int n=0;n<_Np;++n){
    for(int m=0;m<cubNp;++m){
      dfloat resP = 0, resDrT = 0;

      for(int i=0;i<_Np;++i){
        dfloat Vni = V[n*_Np+i];
        resDrT += Vni*cubVr[m*_Np+i];
        resP   += Vni*cubV[m*_Np+i];
      }

      cubDrT[0][n*cubNp+m] = resDrT;
      cubProject[0][n*cubNp+m] = resP;
    }
  }

  free(cubV);
  free(cubVr);
}
*/

} //namespace libp
