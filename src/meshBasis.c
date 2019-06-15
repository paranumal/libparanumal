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

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>

void matrixInverse(int N, dfloat *A);

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);

dfloat mygamma(dfloat x){

  dfloat lgam = lgamma(x);
  dfloat gam  = signgam*exp(lgam);
  return gam;
}

dfloat meshJacobiP(dfloat a, dfloat alpha, dfloat beta, int N){
  
  dfloat ax = a; 

  dfloat *P = (dfloat *) calloc((N+1), sizeof(dfloat));

  // Zero order
  dfloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*mygamma(1+alpha)*mygamma(1+beta)/mygamma(1+alpha+beta);
  dfloat p0     = 1.0/sqrt(gamma0);

  if (N==0){ free(P); return p0;}
  P[0] = p0; 

  // first order
  dfloat gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  dfloat p1     = ((alpha+beta+2)*ax/2 + (alpha-beta)/2)/sqrt(gamma1);
  if (N==1){free(P); return p1;} 

  P[1] = p1;

  /// Repeat value in recurrence.
  dfloat aold = 2/(2+alpha+beta)*sqrt((alpha+1.)*(beta+1.)/(alpha+beta+3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for(int i=1;i<=N-1;++i){
    dfloat h1 = 2.*i+alpha+beta;
    dfloat anew = 2./(h1+2.)*sqrt( (i+1.)*(i+1.+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    dfloat bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    P[i+1] = 1./anew*( -aold*P[i-1] + (ax-bnew)*P[i]);
    aold =anew;
  }
  
  dfloat pN = P[N]; 
  free(P);
  return pN;

}

dfloat meshGradJacobiP(dfloat a, dfloat alpha, dfloat beta, int N){

  dfloat PNr = 0;

  if(N>0)
    PNr = sqrt(N*(N+alpha+beta+1.))*meshJacobiP(a, alpha+1.0, beta+1.0, N-1);

  return PNr;
}


void meshOrthonormalBasis1D(dfloat a, int i, dfloat *P, dfloat *Pr){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p1a = meshGradJacobiP(a,0,0,i);

  *P = p1;
  *Pr = p1a;
}


void meshOrthonormalBasisTri2D(dfloat a, dfloat b, int i, int j, dfloat *P, dfloat *Pr, dfloat *Ps){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p2 = meshJacobiP(b,2*i+1,0,j);

  dfloat p1a = meshGradJacobiP(a,0,0,i);
  dfloat p2b = meshGradJacobiP(b,2*i+1,0,j);

  *P = sqrt(2.0)*p1*p2*pow(1.0-b,i);

  *Pr = p1a*p2;
  if(i>0)
    *Pr *= pow(0.5*(1.0-b),i-1);

  *Ps = p1a*p2*0.5*(1.0+a);
  if(i>0)
    *Ps *= pow(0.5*(1.0-b),i-1);

  dfloat tmp = p2b*pow(0.5*(1.0-b), i);
  if(i>0)
    tmp -= 0.5*i*p2*pow(0.5*(1.0-b), i-1);

  *Ps += p1*tmp;

  // normalize
  *Pr *= pow(2.0, i+0.5);
  *Ps *= pow(2.0, i+0.5);
  
}

void meshOrthonormalBasisQuad2D(dfloat a, dfloat b, int i, int j, dfloat *P, dfloat *Pr, dfloat *Ps){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p2 = meshJacobiP(b,0,0,j);
  dfloat p1a = meshGradJacobiP(a,0,0,i);
  dfloat p2b = meshGradJacobiP(b,0,0,j);

  *P = p1*p2;
  *Pr = p1a*p2;
  *Ps = p1*p2b;
  
}

// TW: check normalization
void meshOrthonormalBasisTet3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P, dfloat *Pr, dfloat *Ps, dfloat *Pt){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p2 = meshJacobiP(b,2*i+1,0,j);
  dfloat p3 = meshJacobiP(c,2*(i+j)+2,0,k);

  dfloat p1a = meshGradJacobiP(a,0,0,i);
  dfloat p2b = meshGradJacobiP(b,2*i+1,0,j);
  dfloat p3c = meshGradJacobiP(c,2*(i+j)+2,0,k);

  *P = 2.*sqrt(2.0)*p1*p2*p3*pow(1.0-b,i)*pow(1.0-c,i+j);

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

void meshOrthonormalBasisHex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k, dfloat *P, dfloat *Pr, dfloat *Ps, dfloat *Pt){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p2 = meshJacobiP(b,0,0,j);
  dfloat p3 = meshJacobiP(c,0,0,k);
  dfloat p1a = meshGradJacobiP(a,0,0,i);
  dfloat p2b = meshGradJacobiP(b,0,0,j);
  dfloat p3c = meshGradJacobiP(c,0,0,k);

  *P = p1*p2*p3;
  *Pr = p1a*p2*p3;
  *Ps = p1*p2b*p3;
  *Pt = p1*p2*p3c;
  
}


int meshVandermonde1D(int N, int Npoints, dfloat *r, dfloat **V, dfloat **Vr){

  int Np = (N+1);
  
  *V  = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  for(int n=0; n<Npoints; n++){

    int sk = 0;
    for(int i=0; i<=N; i++){
      int id = n*Np+sk;
      meshOrthonormalBasis1D(r[n], i, V[0]+id, Vr[0]+id);
      sk++;
    }
  }
  
  return Np;
}



int meshVandermondeTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat **V, dfloat **Vr, dfloat **Vs){

  int Np = (N+1)*(N+2)/2;
  
  *V  = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vs = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  for(int n=0; n<Npoints; n++){

    dfloat a, b;
    
    // First convert to ab coordinates
    if(fabs(s[n]-1.0)>1e-8)
      a = 2.0*(1.+r[n])/(1.0-s[n])-1.0;
    else
      a = -1.0; 
    
    b = s[n];
    
    int sk=0;
    for(int i=0; i<=N; i++){
      for(int j=0; j<=N-i; j++){
	int id = n*Np+sk;
	meshOrthonormalBasisTri2D(a, b, i, j, V[0]+id,Vr[0]+id, Vs[0]+id);
	sk++;
      }
    }
  }

  return Np;
}


int meshVandermondeQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat **V, dfloat **Vr, dfloat **Vs){

  int Np = (N+1)*(N+1);
  
  *V  = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vs = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  for(int n=0; n<Npoints; n++){

    int sk = 0;
    for(int i=0; i<=N; i++){
      for(int j=0; j<=N; j++){
	int id = n*Np+sk;
	meshOrthonormalBasisQuad2D(r[n], s[n], i, j, V[0]+id, Vr[0]+id, Vs[0]+id);
	sk++;
      }
    }
  }

  return Np;
}



int meshVandermondeTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t,
			  dfloat **V, dfloat **Vr, dfloat **Vs, dfloat **Vt){

  int Np = (N+1)*(N+2)*(N+3)/6; 

  *V  = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vs = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vt = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  
  for(int n=0; n<Npoints; n++){
    // First convert to abc coordinates
    dfloat a, b, c;
    
    if(fabs(s[n]+t[n])>1e-8)
      a = 2.0*(1.+r[n])/(-s[n]-t[n])-1.0;
    else
      a = -1.0; 

    if(fabs(t[n]-1)>1e-8)
      b = 2.0*(1+s[n])/(1.-t[n])-1.0;
    else
      b = -1.;
    
    c = t[n];

    int sk=0;
    for(int i=0; i<=N; i++){
      for(int j=0; j<=N-i; j++){
	for(int k=0; k<=N-i-j; k++){
	  int id = n*Np+sk;	
	  meshOrthonormalBasisTet3D(a, b, c, i, j, k, V[0]+id, Vr[0]+id, Vs[0]+id, Vt[0]+id);
	  sk++;
	}
      }
    }
  }

  return Np;
}


int meshVandermondeHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat **V, dfloat **Vr, dfloat **Vs, dfloat **Vt){

  int Np = (N+1)*(N+1)*(N+1);
  
  *V  = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vs = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Vt = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  for(int n=0; n<Npoints; n++){

    int sk = 0;
    for(int i=0; i<=N; i++){
      for(int j=0; j<=N; j++){
	for(int k=0; k<=N; k++){
	  int id = n*Np+sk;
	  meshOrthonormalBasisHex3D(r[n], s[n], t[n], i, j, k, V[0]+id, Vr[0]+id, Vs[0]+id, Vt[0]+id);
	  sk++;
	}
      }
    }
  }

  return Np;
}

extern "C"
{
  void dgesv_ (int *NrowsA, int *NcolsA, double *A, int *LDA, int *ipiv,  double *B, int *LDB, int *info);
}

// C = A/B  = trans(trans(B)\trans(A))
// assume row major
void matrixRightSolve(int NrowsA, int NcolsA, dfloat *A, int NrowsB, int NcolsB, dfloat *B, dfloat *C){

  int info;

  int NrowsX = NcolsB;
  int NcolsX = NrowsB;
  
  int NrowsY = NcolsA;
  int NcolsY = NrowsA;
  
  int lwork = NrowsX*NcolsX;
  
  // compute inverse mass matrix
  double *tmpX = (double*) calloc(NrowsX*NcolsX, sizeof(double));
  double *tmpY = (double*) calloc(NrowsY*NcolsY, sizeof(double));

  int    *ipiv = (int*) calloc(NrowsX, sizeof(int));
  double *work = (double*) calloc(lwork, sizeof(double));

  for(int n=0;n<NrowsX*NcolsX;++n){
    tmpX[n] = B[n];
  }

  for(int n=0;n<NrowsY*NcolsY;++n){
    tmpY[n] =A[n];
  }
  
  dgesv_(&NrowsX, &NcolsY, tmpX, &NrowsX, ipiv, tmpY, &NrowsY, &info); // ?

  for(int n=0;n<NrowsY*NcolsY;++n){
    C[n] = tmpY[n];
  }

  if(info)
    printf("matrixRightSolve: dgesv reports info = %d when inverting matrix\n", info);

  free(work);
  free(ipiv);
  free(tmpX);
  free(tmpY);
}

void matrixPrint(FILE *fp, const char *mess, int Nrows, int Ncols, dfloat *A){
#if 0
  fprintf(fp, "%s=[\n", mess);
  for(int n=0;n<Nrows;++n){
    for(int m=0;m<Ncols;++m){
      fprintf(fp," % e ", A[n*Ncols+m]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"];\n");
#endif
}


void matrixCompare(FILE *fp, const char *mess, int Nrows, int Ncols, dfloat *A, dfloat *B){

  dfloat maxd = 0;
  for(int n=0;n<Nrows;++n){
    for(int m=0;m<Ncols;++m){
      dfloat newd = fabs(A[n*Ncols+m]-B[n*Ncols+m]);
      if(newd>maxd)
	maxd = newd;
    }
  }
  fprintf(fp,"%s=% e\n", mess, maxd);
}


void meshDmatrix1D(int N, int Npoints, dfloat *r, dfloat **Dr){

  dfloat *V, *Vr;
  int Np = meshVandermonde1D(N, Npoints, r, &V, &Vr);
  
  *Dr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  
  matrixRightSolve(Np, Np, Vr, Np, Np, V, Dr[0]);
  
  free(V);
  free(Vr);
}




void meshDmatricesTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat **Dr, dfloat **Ds){
  
  dfloat *V, *Vr, *Vs;
  int Np = meshVandermondeTri2D(N, Npoints, r, s, &V, &Vr, &Vs);
  
  *Dr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Ds = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  
  matrixRightSolve(Np, Np, Vr, Np, Np, V, Dr[0]);
  matrixRightSolve(Np, Np, Vs, Np, Np, V, Ds[0]);

  free(V);
  free(Vr);
  free(Vs);
}

void meshDmatricesQuad2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat **Dr, dfloat **Ds){

  dfloat *V, *Vr, *Vs;
  int Np = meshVandermondeQuad2D(N, Npoints, r, s, &V, &Vr, &Vs);

  *Dr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Ds = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  matrixRightSolve(Np, Np, Vr, Np, Np, V, Dr[0]);
  matrixRightSolve(Np, Np, Vs, Np, Np, V, Ds[0]);
  
  free(V);
  free(Vr);
  free(Vs);
}

void meshDmatricesTet3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat **Dr, dfloat **Ds, dfloat **Dt){

  dfloat *V, *Vr, *Vs, *Vt;
  int Np = meshVandermondeTet3D(N, Npoints, r, s, t, &V, &Vr, &Vs, &Vt);

  *Dr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Ds = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Dt = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  matrixRightSolve(Np, Np, Vr, Np, Np, V, Dr[0]);
  matrixRightSolve(Np, Np, Vs, Np, Np, V, Ds[0]);
  matrixRightSolve(Np, Np, Vt, Np, Np, V, Dt[0]);

  free(V);
  free(Vr);
  free(Vs);
  free(Vt);
}


void meshDmatricesHex3D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *t, dfloat **Dr, dfloat **Ds, dfloat **Dt){

  dfloat *V, *Vr, *Vs, *Vt;
  int Np = meshVandermondeHex3D(N, Npoints, r, s, t, &V, &Vr, &Vs, &Vt);

  *Dr = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Ds = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));
  *Dt = (dfloat *) calloc(Npoints*Np, sizeof(dfloat));

  matrixRightSolve(Np, Np, Vr, Np, Np, V, Dr[0]);
  matrixRightSolve(Np, Np, Vs, Np, Np, V, Ds[0]);
  matrixRightSolve(Np, Np, Vt, Np, Np, V, Dt[0]);

  free(V);
  free(Vr);
  free(Vs);
  free(Vt);

}

// assumes NpointsIn = (N+1)*(N+2)/2
void meshInterpolationMatrixTri2D(int N,
				  int NpointsIn, dfloat *rIn, dfloat *sIn,
				  int NpointsOut, dfloat *rOut, dfloat *sOut,
				  dfloat **I){

  dfloat *VIn, *VrIn, *VsIn;
  dfloat *VOut, *VrOut, *VsOut;
  
  int NpIn  = meshVandermondeTri2D(N, NpointsIn, rIn, sIn, &VIn, &VrIn, &VsIn);
  int NpOut = meshVandermondeTri2D(N, NpointsOut, rOut, sOut, &VOut, &VrOut, &VsOut);
  
  *I = (dfloat *) calloc(NpointsIn*NpointsOut, sizeof(dfloat));
  
  matrixRightSolve(NpointsOut, NpOut, VOut, NpointsIn, NpIn, VIn, *I);

  free(VIn);
  free(VrIn);
  free(VsIn);

  free(VOut);
  free(VrOut);
  free(VsOut);
}


// masMatrix = inv(V')*inv(V) = inv(V*V')
void meshMassMatrix1D(int N, int Npoints, dfloat *r, dfloat **MM){

  dfloat *V, *Vr;
  int Np = meshVandermonde1D(N, Npoints, r, &V, &Vr);

  *MM = (dfloat *) calloc(Np*Np, sizeof(dfloat));
  
  for(int n=0;n<Npoints;++n){
    for(int m=0;m<Npoints;++m){
      dfloat res = 0;
      for(int i=0;i<Np;++i){
	res += V[n*Np+i]*V[m*Np+i];
      }
      MM[0][n*Npoints + m] = res;
    }
  }
  
  matrixInverse(Np, MM[0]);

  free(V);
  free(Vr);
}


// masMatrix = inv(V')*inv(V) = inv(V*V')
void meshMassMatrixTri2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat **MM){

  dfloat *V, *Vr, *Vs;
  
  int Np = meshVandermondeTri2D(N, Npoints, r, s, &V, &Vr, &Vs);

  *MM = (dfloat *) calloc(Np*Np, sizeof(dfloat));
  
  for(int n=0;n<Npoints;++n){
    for(int m=0;m<Npoints;++m){
      dfloat res = 0;
      for(int i=0;i<Np;++i){
	res += V[n*Np+i]*V[m*Np+i];
      }
      MM[0][n*Npoints + m] = res;
    }
  }
  
  matrixInverse(Np, MM[0]);

  free(V);
  free(Vr);
  free(Vs);
}




#if TEST_MESH_BASIS==1
// mpic++ -I../libs/gatherScatter/ -I../../occa/include  -I../include -o meshBasis meshBasis.c matrixInverse.c readArray.c -Ddfloat=double -llapack -lblas -lm -DDHOLMES='"../"' -DdfloatFormat='"%lf"' -DTEST_MESH_BASIS=1 

// to run with degree 2:
// ./meshBasis 2
int main(int argc, char **argv){

  dfloat *r, *s, *t;
  dfloat *Dr, *Ds, *Dt;
  dfloat *Vr, *Vs, *Vt, *V;
  dfloat *MM;

  dfloat *cubr, *cubs, *cubt;
  dfloat *cubInterp;
  
  dfloat *fileDr, *fileDs, *fileDt, *fileCubInterp, *fileMM;
  
  int Nrows, Ncols;

  int N = atoi(argv[1]);
  int cubNp, Np;

  char fname[BUFSIZ];

  { // TRIANGLE TEST
    sprintf(fname, DHOLMES "/nodes/triangleN%02d.dat", N);
    
    FILE *fp = fopen(fname, "r");
    
    readDfloatArray(fp, "Nodal r-coordinates", &r,&Np,&Ncols);
    readDfloatArray(fp, "Nodal s-coordinates", &s,&Np,&Ncols);

    readDfloatArray(fp, "Nodal Dr differentiation matrix", &(fileDr), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Ds differentiation matrix", &(fileDs), &Np, &Ncols);

    readDfloatArray(fp, "Cubature r-coordinates", &cubr,&cubNp,&Ncols);
    readDfloatArray(fp, "Cubature s-coordinates", &cubs,&cubNp,&Ncols);

    printf("cubNp = %d\n", cubNp);
    
    readDfloatArray(fp, "Cubature Interpolation Matrix", &fileCubInterp,&cubNp,&Ncols);

    readDfloatArray(fp, "Nodal Mass Matrix", &fileMM,&Np,&Ncols);
    
    fclose(fp);
    
    meshDmatricesTri2D(N, Np, r, s, &Dr, &Ds);
    meshVandermondeTri2D(N, Np, r, s, &V, &Vr, &Vs);
    meshInterpolationMatrixTri2D(N, Np, r, s, cubNp, cubr, cubs, &cubInterp);
    meshMassMatrixTri2D(N, Np, r, s, &MM);
    
    matrixCompare(stdout, "TRI2D: |Dr-fileDr|", Np, Np, Dr, fileDr);
    matrixCompare(stdout, "TRI2D: |Ds-fileDs|", Np, Np, Ds, fileDs);
    matrixCompare(stdout, "TRI2D: |MM-fileMM|", Np, Np, MM, fileMM);
    matrixCompare(stdout, "TRI2D: |cubInterp-fileCubInterp|", cubNp, Np, cubInterp, fileCubInterp);


    matrixPrint(stdout, "TRI2D: cubInterp", cubNp, Np, cubInterp);
    matrixPrint(stdout, "TRI2D: fileCubInterp", cubNp, Np, fileCubInterp);
    
  }

  { // QUAD TEST
    sprintf(fname, DHOLMES "/nodes/quadrilateralN%02d.dat", N);
    
    FILE *fp = fopen(fname, "r");
    
    readDfloatArray(fp, "Nodal r-coordinates", &r,&Np,&Ncols);
    readDfloatArray(fp, "Nodal s-coordinates", &s,&Np,&Ncols);

    readDfloatArray(fp, "Nodal Dr differentiation matrix", &(fileDr), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Ds differentiation matrix", &(fileDs), &Np, &Ncols);
    
    fclose(fp);
    
    meshDmatricesQuad2D(N, Np, r, s, &Dr, &Ds);
    meshVandermondeQuad2D(N, Np, r, s, &V, &Vr, &Vs);

    matrixCompare(stdout, "QUAD2D: |Dr-fileDr|", Np, Np, Dr, fileDr);
    matrixCompare(stdout, "QUAD2D: |Ds-fileDs|", Np, Np, Ds, fileDs);
    
  }

  
  { // TETRAHEDRON TEST
    sprintf(fname, DHOLMES "/nodes/tetN%02d.dat", N);
    
    FILE *fp = fopen(fname, "r");
    
    readDfloatArray(fp, "Nodal r-coordinates", &r,&Np,&Ncols);
    readDfloatArray(fp, "Nodal s-coordinates", &s,&Np,&Ncols);
    readDfloatArray(fp, "Nodal t-coordinates", &t,&Np,&Ncols);

    readDfloatArray(fp, "Nodal Dr differentiation matrix", &(fileDr), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Ds differentiation matrix", &(fileDs), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Dt differentiation matrix", &(fileDt), &Np, &Ncols);
    
    fclose(fp);
    
    meshDmatricesTet3D(N, Np, r, s, t, &Dr, &Ds, &Dt);
    meshVandermondeTet3D(N, Np, r, s, t, &V, &Vr, &Vs, &Vt);
    
    matrixCompare(stdout, "TET3D: |Dr-fileDr|", Np, Np, Dr, fileDr);
    matrixCompare(stdout, "TET3D: |Ds-fileDs|", Np, Np, Ds, fileDs);
    matrixCompare(stdout, "TET3D: |Dt-fileDt|", Np, Np, Dt, fileDt);
  }

  if(0){ // HEX TEST [ does not have Dr,Ds,Dt in node files ]
    sprintf(fname, DHOLMES "/nodes/hexN%02d.dat", N);
    
    FILE *fp = fopen(fname, "r");
    
    readDfloatArray(fp, "Nodal r-coordinates", &r,&Np,&Ncols);
    readDfloatArray(fp, "Nodal s-coordinates", &s,&Np,&Ncols);
    readDfloatArray(fp, "Nodal t-coordinates", &t,&Np,&Ncols);

    readDfloatArray(fp, "Nodal Dr differentiation matrix", &(fileDr), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Ds differentiation matrix", &(fileDs), &Np, &Ncols);
    readDfloatArray(fp, "Nodal Dt differentiation matrix", &(fileDt), &Np, &Ncols);
    
    fclose(fp);
    
    meshDmatricesHex3D(N, Np, r, s, t, &Dr, &Ds, &Dt);
    meshVandermondeHex3D(N, Np, r, s, t, &V, &Vr, &Vs, &Vt);
    
    matrixCompare(stdout, "HEX3D: |Dr-fileDr|", Np, Np, Dr, fileDr);
    matrixCompare(stdout, "HEX3D: |Ds-fileDs|", Np, Np, Ds, fileDs);
    matrixCompare(stdout, "HEX3D: |Dt-fileDt|", Np, Np, Dt, fileDt);
  }

  
  return 0;
}


#endif






