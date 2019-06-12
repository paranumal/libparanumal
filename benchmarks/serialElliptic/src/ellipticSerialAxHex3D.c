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

#define p_Np (p_Nq*p_Nq*p_Nq)

#include "elliptic.h"

extern "C"
{
  void ax_e_(dfloat *w, const dfloat *u, const dfloat *g, dfloat *ur, dfloat *us, dfloat *ut,dfloat *wk, const dfloat *DT, const dfloat *D);
  void ax_(const int *Nel, dfloat *w, const dfloat *u, const dfloat *g, dfloat *ur, dfloat *us, dfloat *ut,dfloat *wk, const int *, const dfloat *DT, const dfloat *D);
  
  void local_grad3_ (dfloat * __restrict__ qr,dfloat * __restrict__ qs, dfloat * __restrict__ qt, 
		     const dfloat * __restrict__ q, const int *N, const dfloat * __restrict__ DT, const dfloat * __restrict__ D);
  
  void local_grad3_t_ (dfloat * __restrict__ q,
		       const dfloat * __restrict__ qr, 
		       const dfloat * __restrict__ qs, 
		       const dfloat * __restrict__ qt,
		       const int *N, const dfloat * __restrict__ DT, const dfloat * __restrict__ D, dfloat * __restrict__ wk);

  void libxsmm_dgemm_ (char *, char *, int *, int *, int *,
		       const dfloat *, const dfloat * __restrict, int *,
		       const dfloat * __restrict, int *,
		       const dfloat *, dfloat * __restrict, int *);

  void dgemm_ (char *, char *, int *, int *, int *,
	       const dfloat *, const dfloat * __restrict, int *,
	       const dfloat * __restrict, int *,
	       const dfloat *, dfloat * __restrict, int *);
}

#define USE_BLAS 0
#define USE_XSMM 0

// hack
#define p_Nggeo 7
//#define p_Nggeo 6

template < const int rowsA, const int rowsB, const int colsC >
  static void mxm(const dfloat * __restrict__ A,
		    const dfloat * __restrict__ B,
		    const dfloat BETA, 
		    dfloat * __restrict__ C){

#if USE_XSMM || USE_BLAS
  // dgemm (TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  // C = beta*C + A*B
  char TRANSA = 'N';
  char TRANSB = 'N';
  int M = rowsA;
  int N = colsC;
  int K = rowsB;
  dfloat ALPHA = 1;
  int LDA = rowsA;
  int LDB = rowsB;
  int LDC = rowsA;

#if USE_XSMM
  libxsmm_dgemm_
#else
    dgemm_
#endif
    (&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);


#else
  if(BETA)
    for(int j=0;j<colsC;++j){
      for(int i=0;i<rowsA;++i){
	dfloat res = BETA*C[i+j*rowsA];
	for(int k=0;k<rowsB;++k){
	  res += A[i+k*rowsA]*B[k+j*rowsB];
	}
	C[i+j*rowsA] = res;
      }
    }
  else
    for(int j=0;j<colsC;++j){
      for(int i=0;i<rowsA;++i){
	dfloat res = 0;
	for(int k=0;k<rowsB;++k){
	  res += A[i+k*rowsA]*B[k+j*rowsB];
	}
	C[i+j*rowsA] = res;
      }
    }
#endif
}


template < const int p_Nq >
void ellipticSerialElementAxHexKernel3D(const dfloat * __restrict__ ggeo,
					const dfloat * __restrict__ D,
					const dfloat * __restrict__ S,
					const dfloat * __restrict__ MM,
					const dfloat lambda,
					const dfloat * __restrict__ q,
					dfloat * __restrict__ qr,
					dfloat * __restrict__ qs,
					dfloat * __restrict__ qt,
					dfloat * __restrict__ Aq,
					dfloat * __restrict__ wk){

  D    = (dfloat*)__builtin_assume_aligned(D, USE_OCCA_MEM_BYTE_ALIGN) ;
  S    = (dfloat*)__builtin_assume_aligned(S, USE_OCCA_MEM_BYTE_ALIGN) ;
  MM   = (dfloat*)__builtin_assume_aligned(MM, USE_OCCA_MEM_BYTE_ALIGN) ;
  q    = (dfloat*)__builtin_assume_aligned(q, USE_OCCA_MEM_BYTE_ALIGN) ;
  qr    = (dfloat*)__builtin_assume_aligned(qr, USE_OCCA_MEM_BYTE_ALIGN) ;
  qs    = (dfloat*)__builtin_assume_aligned(qs, USE_OCCA_MEM_BYTE_ALIGN) ;
  qt    = (dfloat*)__builtin_assume_aligned(qt, USE_OCCA_MEM_BYTE_ALIGN) ;
  Aq   = (dfloat*)__builtin_assume_aligned(Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
  ggeo = (dfloat*)__builtin_assume_aligned(ggeo, USE_OCCA_MEM_BYTE_ALIGN) ;

dfloat zero = 0, one = 1.0;
const int N = p_Nq-1;

  // grad
#if 0
  mxm<p_Nq,p_Nq,p_Nq*p_Nq>(S, q, zero, qr); // D(:,:)*q(:,:+::)
  for(int k=0;k<p_Nq;++k){
    mxm<p_Nq,p_Nq,p_Nq>(q+k*p_Nq*p_Nq, D, zero, qs+k*p_Nq*p_Nq);
  }
  mxm<p_Nq*p_Nq,p_Nq,p_Nq>(q, D, zero, qt);
#else
  local_grad3_ (qr, qs, qt, q, &N, S, D);
#endif
  
  for(int n=0;n<p_Np;++n){
#if 1
    const dfloat G00 = ggeo[n+G00ID*p_Np], G01 = ggeo[n+G01ID*p_Np], G11 = ggeo[n+G11ID*p_Np];
    const dfloat GWJ = ggeo[n+GWJID*p_Np];
    const dfloat G12 = ggeo[n+G12ID*p_Np], G02 = ggeo[n+G02ID*p_Np], G22 = ggeo[n+G22ID*p_Np];
#else
    const dfloat G00 = ggeo[0+n*6], G01 = ggeo[1+n*6], G02 = ggeo[2+n*6];
    const dfloat G11 = ggeo[3+n*6], G12 = ggeo[4+n*6], G22 = ggeo[5+n*6];
#endif
    dfloat qrn = G00*qr[n] + G01*qs[n] + G02*qt[n];
    dfloat qsn = G01*qr[n] + G11*qs[n] + G12*qt[n];
    dfloat qtn = G02*qr[n] + G12*qs[n] + G22*qt[n];

#if 0
    dfloat qrn = ggeo[n+G00ID*p_Np]*qr[n] + ggeo[n+G01ID*p_Np]*qs[n] + ggeo[n+G02ID*p_Np]*qt[n];
    dfloat qsn = ggeo[n+G01ID*p_Np]*qr[n] + ggeo[n+G11ID*p_Np]*qs[n] + ggeo[n+G12ID*p_Np]*qt[n];
    dfloat qtn = ggeo[n+G02ID*p_Np]*qr[n] + ggeo[n+G12ID*p_Np]*qs[n] + ggeo[n+G22ID*p_Np]*qt[n];
#endif
#if 0
    dfloat qrn = ggeo[n*p_Nggeo+G00ID]*qr[n] + ggeo[n*p_Nggeo+G01ID]*qs[n] + ggeo[n*p_Nggeo+G02ID]*qt[n];
    dfloat qsn = ggeo[n*p_Nggeo+G01ID]*qr[n] + ggeo[n*p_Nggeo+G11ID]*qs[n] + ggeo[n*p_Nggeo+G12ID]*qt[n];
    dfloat qtn = ggeo[n*p_Nggeo+G02ID]*qr[n] + ggeo[n*p_Nggeo+G12ID]*qs[n] + ggeo[n*p_Nggeo+G22ID]*qt[n];
#endif
    qr[n] = qrn;
    qs[n] = qsn;
    qt[n] = qtn;
  }

  // gradT
#if 0
  mxm<p_Nq,p_Nq,p_Nq*p_Nq>(D, qr, zero, Aq); // D(:,:)*q(:,:+::)
  for(int k=0;k<p_Nq;++k){
    mxm<p_Nq,p_Nq,p_Nq>(qs+k*p_Nq*p_Nq, S, one, Aq+k*p_Nq*p_Nq);
  }
  mxm<p_Nq*p_Nq,p_Nq,p_Nq>(qt, S, one, Aq);
#else
  //  const int N = p_Nq-1;
  local_grad3_t_ (Aq, qr, qs, qt, &N, S, D, wk);
#endif
}


template < const int p_Nq >
void ellipticSerialAxHexKernel3D (const hlong Nelements,
				  const dfloat * __restrict__ ggeo ,
				  const dfloat * __restrict__ D ,
				  const dfloat * __restrict__ S ,
				  const dfloat * __restrict__ MM ,
				  const dfloat lambda,
				  const dfloat * __restrict__ q ,
				  dfloat * __restrict__ Aq ){
  
  D    = (dfloat*)__builtin_assume_aligned(D, USE_OCCA_MEM_BYTE_ALIGN) ;
  S    = (dfloat*)__builtin_assume_aligned(S, USE_OCCA_MEM_BYTE_ALIGN) ;
  MM   = (dfloat*)__builtin_assume_aligned(MM, USE_OCCA_MEM_BYTE_ALIGN) ;
  q    = (dfloat*)__builtin_assume_aligned(q, USE_OCCA_MEM_BYTE_ALIGN) ;
  Aq   = (dfloat*)__builtin_assume_aligned(Aq, USE_OCCA_MEM_BYTE_ALIGN) ;
  ggeo = (dfloat*)__builtin_assume_aligned(ggeo, USE_OCCA_MEM_BYTE_ALIGN) ;
  
  dfloat s_q  [p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));

#if 0
  dfloat s_tmp[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  
  dfloat s_qr[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_qs[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_qt[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_wk[p_Np] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
#endif

  // ok 
  dfloat s_D[p_Nq][p_Nq]  __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_S[p_Nq][p_Nq]  __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));

  for(int j=0;j<p_Nq;++j){
    for(int i=0;i<p_Nq;++i){
      s_D[j][i] = D[j*p_Nq+i];
      s_S[j][i] = S[j*p_Nq+i];
    }
  }

  const int c_Np = p_Np;
  const int p_N = p_Nq-1;

  //  ax_ (&Nelements, Aq, q, ggeo, s_qr, s_qs, s_qt, s_wk, &p_N, S, D);

  if(1)
  for(dlong e=0; e<Nelements; ++e){
    
    const dlong element = e;

#if 0
    // MIXED STEFAN + TW VERSION
    ellipticSerialElementAxHexKernel3D<p_Nq>(ggeo+element*p_Np*p_Nggeo,
					     D, S, MM, lambda, q + element*p_Np,
					     s_qr, s_qs, s_qt, Aq+element*p_Np, s_wk);
#endif
   
#if 0
    // STEFAN version
    ax_e_(Aq+element*p_Np, 
	  q+element*p_Np,
	  ggeo + element*p_Nggeo*p_Np, // note layout is wrong
	  s_qr,
	  s_qs,
	  s_qt,
	  s_wk,
	  S, 
	  D);
#endif
    
#if 1
    // TW version
    for(int k = 0; k < p_Nq; k++) {
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
          const dlong base = i + j*p_Nq + k*p_Nq*p_Nq + element*c_Np;
          const dfloat qbase = q[base];
          s_q[k][j][i] = qbase;
        }
      }
    }

    for(int k=0;k<p_Nq;++k){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
#if 1
          const dlong gbase = element*p_Nggeo*c_Np + k*p_Nq*p_Nq + j*p_Nq + i;
          const dfloat r_G00 = ggeo[gbase+G00ID*p_Np];
          const dfloat r_G01 = ggeo[gbase+G01ID*p_Np];
          const dfloat r_G11 = ggeo[gbase+G11ID*p_Np];
          const dfloat r_G12 = ggeo[gbase+G12ID*p_Np];
          const dfloat r_G02 = ggeo[gbase+G02ID*p_Np];
          const dfloat r_G22 = ggeo[gbase+G22ID*p_Np];
#endif

#if 0
          const dlong gbase = element*p_Nggeo*c_Np + (k*p_Nq*p_Nq + j*p_Nq + i)*p_Nggeo;
	  const dfloat * __restrict__ ggeobase = ggeo+gbase;
          const dfloat r_G00 = ggeobase[0];
	  const dfloat r_G01 = ggeobase[1];
          const dfloat r_G02 = ggeobase[2];
          const dfloat r_G11 = ggeobase[3];
          const dfloat r_G12 = ggeobase[4];
          const dfloat r_G22 = ggeobase[5];
#endif

#if 0
          const dlong gbase = element*p_Nggeo*c_Np + (k*p_Nq*p_Nq + j*p_Nq + i);
	  const dfloat * __restrict__ ggeobase = ggeo+gbase;
          const dfloat r_G00 = ggeobase[0*p_Np];
          const dfloat r_G01 = ggeobase[1*p_Np];
          const dfloat r_G02 = ggeobase[2*p_Np];
          const dfloat r_G11 = ggeobase[3*p_Np];
          const dfloat r_G12 = ggeobase[4*p_Np];
          const dfloat r_G22 = ggeobase[5*p_Np];
#endif

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

          for(int m = 0; m < p_Nq; m++) {
            qr += s_S[m][i]*s_q[k][j][m];  
            qs += s_S[m][j]*s_q[k][m][i];           
            qt += s_S[m][k]*s_q[m][j][i]; 
          }

          dfloat Gqr = r_G00*qr;
          Gqr += r_G01*qs;
          Gqr += r_G02*qt;
          
          dfloat Gqs = r_G01*qr;
          Gqs += r_G11*qs;
          Gqs += r_G12*qt;

          dfloat Gqt = r_G02*qr;
          Gqt += r_G12*qs;
          Gqt += r_G22*qt;
          
          s_Gqr[k][j][i] = Gqr;
          s_Gqs[k][j][i] = Gqs;
          s_Gqt[k][j][i] = Gqt;
        }
      }
    }

    for(int k = 0;k < p_Nq; k++){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
#if 1
          const dlong gbase = element*p_Nggeo*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase+GWJID*p_Np];

          dfloat r_Aq = r_GwJ*lambda*s_q[k][j][i];
#endif
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;


          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_D[m][i]*s_Gqr[k][j][m];
          for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_D[m][j]*s_Gqs[k][m][i];
          for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_D[m][k]*s_Gqt[m][j][i];

          const dlong id = element*p_Np +k*p_Nq*p_Nq+ j*p_Nq + i;
          Aq[id] = r_Aqr + r_Aqs + r_Aqt +r_Aq;
        }
      }
    }
#endif

    // ok
#if 0
    for(int k=0;k<p_Nq;++k){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){

#if 0
          const dlong gbase = element*p_Nggeo*c_Np + k*p_Nq*p_Nq + j*p_Nq + i;
          const dfloat r_G00 = ggeo[gbase+G00ID*p_Np];
          const dfloat r_G01 = ggeo[gbase+G01ID*p_Np];
          const dfloat r_G11 = ggeo[gbase+G11ID*p_Np];
          const dfloat r_G02 = ggeo[gbase+G02ID*p_Np];
          const dfloat r_G12 = ggeo[gbase+G12ID*p_Np];
          const dfloat r_G22 = ggeo[gbase+G22ID*p_Np];
#endif

#if 1
          const dlong gbase = element*p_Nggeo*c_Np + (k*p_Nq*p_Nq + j*p_Nq + i);
	  const dfloat * __restrict__ ggeobase = ggeo+gbase;
          const dfloat r_G00 = ggeobase[0*p_Np];
          const dfloat r_G01 = ggeobase[1*p_Np];
          const dfloat r_G02 = ggeobase[2*p_Np];
          const dfloat r_G11 = ggeobase[3*p_Np];
          const dfloat r_G12 = ggeobase[4*p_Np];
          const dfloat r_G22 = ggeobase[5*p_Np];
#endif

          dfloat qr = 0.f;
          dfloat qs = 0.f;
          dfloat qt = 0.f;

	  const dfloat * __restrict__ qei = q + element*p_Np + k*p_Nq*p_Nq + j*p_Nq;
	  const dfloat * __restrict__ qej = q + element*p_Np + k*p_Nq*p_Nq + i;
	  const dfloat * __restrict__ qek = q + element*p_Np + j*p_Nq + i;

          for(int m = 0; m < p_Nq; m++)
	    qr += s_S[m][i]*qei[m]; 
	  for(int m = 0; m < p_Nq; m++) 
	    //	    qs += s_S[m][j]*s_q[k][m][i]; 	    
	    qs += s_S[m][j]*qej[m*p_Nq];
          for(int m = 0; m < p_Nq; m++) {
	    //	    qt += s_S[m][k]*s_q[m][j][i]; 
	    qt += s_S[m][k]*qek[m*p_Nq*p_Nq];
	  }
	  dfloat Gqr = r_G00*qr;
	  Gqr += r_G01*qs;
	  Gqr += r_G02*qt;
	  
	  dfloat Gqs = r_G01*qr;
	  Gqs += r_G11*qs;
	  Gqs += r_G12*qt;

	  dfloat Gqt = r_G02*qr;
	  Gqt += r_G12*qs;
	  Gqt += r_G22*qt;
	  
          s_Gqr[k][j][i] = Gqr;
          s_Gqs[k][j][i] = Gqs;
          s_Gqt[k][j][i] = Gqt;
        }
      }
    }

    for(int k = 0;k < p_Nq; k++){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
#if 0
	  const dlong gbase = element*p_Nggeo*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;
	  const dfloat r_GwJ = ggeo[gbase+p_GWJID*p_Np];

	  dfloat r_Aq = r_GwJ*lambda*s_q[k][j][i];
#endif
          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;


          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_S[i][m]*s_Gqr[k][j][m];
	  for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_S[j][m]*s_Gqs[k][m][i];
	  for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_S[k][m]*s_Gqt[m][j][i];

          const dlong id = element*p_Np +k*p_Nq*p_Nq+ j*p_Nq + i;
          Aq[id] = r_Aqr + r_Aqs + r_Aqt; // +r_Aq;
        }
      }
    }
#endif

#if 0
    // sophisticated 9.48e7
    dfloat s_qr[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
    dfloat s_qs[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
    dfloat s_qt[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
    
    // loop over slices 
    for(int k=0;k<p_Nq;++k){
      for(int j=0;j<p_Nq;++j){
	const dfloat * __restrict__ qei = q + element*p_Np + k*p_Nq*p_Nq + j*p_Nq;
	
        for(int i=0;i<p_Nq;++i){
          dfloat qr = 0.f;
          for(int m = 0; m < p_Nq; m++)
	    qr += s_S[m][i]*qei[m]; 
	  s_qr[k][j][i] = qr;
	}
      }
     
      for(int i=0;i<p_Nq;++i){
	const dfloat * __restrict__ qej = q + element*p_Np + k*p_Nq*p_Nq + i;
	  
	for(int j=0;j<p_Nq;++j){
	  
	  dfloat qs = 0.f;
	  for(int m = 0; m < p_Nq; m++)
	    qs += s_S[m][j]*qej[m*p_Nq]; 
	  
	  s_qs[k][j][i] = qs;
	}
      }
    }
    
    for(int j=0;j<p_Nq;++j){
      for(int i=0;i<p_Nq;++i){
	const dfloat * __restrict__ qek = q + element*p_Np + j*p_Nq + i;
	
	for(int k=0;k<p_Nq;++k){
	  dfloat qt = 0;
          for(int m = 0; m < p_Nq; m++) 
	    qt += s_S[m][k]*qek[m*p_Nq*p_Nq];
	  
	  s_qt[k][j][i] = qt;
	}
      }
    }
    
    for(int k=0;k<p_Nq;++k){
      for(int j=0;j<p_Nq;++j){
	for(int i=0;i<p_Nq;++i){
	  const dlong gbase = element*p_Nggeo*c_Np + k*p_Nq*p_Nq + j*p_Nq + i;
	  const dfloat r_G00 = ggeo[gbase+G00ID*p_Np];
	  const dfloat r_G01 = ggeo[gbase+G01ID*p_Np];
	  const dfloat r_G11 = ggeo[gbase+G11ID*p_Np];
	  const dfloat r_G02 = ggeo[gbase+G02ID*p_Np];
	  const dfloat r_G12 = ggeo[gbase+G12ID*p_Np];
	  const dfloat r_G22 = ggeo[gbase+G22ID*p_Np];
	  
	  const dfloat qr = s_qr[k][j][i];
	  const dfloat qs = s_qs[k][j][i];
	  const dfloat qt = s_qt[k][j][i];
	  
	  dfloat Gqr = r_G00*qr;
	  Gqr += r_G01*qs;
	  Gqr += r_G02*qt;
	
	  dfloat Gqs = r_G01*qr;
	  Gqs += r_G11*qs;
	  Gqs += r_G12*qt;
	  
	  dfloat Gqt = r_G02*qr;
	  Gqt += r_G12*qs;
	  Gqt += r_G22*qt;
	  
	  s_Gqr[k][j][i] = Gqr;
	  s_Gqs[k][j][i] = Gqs;
	  s_Gqt[k][j][i] = Gqt;
	}
      }
    }
    
    for(int k = 0;k < p_Nq; k++){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){

          dfloat r_Aqr = 0, r_Aqs = 0, r_Aqt = 0;

          for(int m = 0; m < p_Nq; m++)
            r_Aqr += s_S[i][m]*s_Gqr[k][j][m];
	  for(int m = 0; m < p_Nq; m++)
            r_Aqs += s_S[j][m]*s_Gqs[k][m][i];
	  for(int m = 0; m < p_Nq; m++)
            r_Aqt += s_S[k][m]*s_Gqt[m][j][i];

          const dlong id = element*p_Np +k*p_Nq*p_Nq+ j*p_Nq + i;
          Aq[id] = r_Aqr + r_Aqs + r_Aqt; // +r_Aq;
        }
      }
    }
#endif

  }
}

void ellipticSerialAxHexKernel3D(const int Nq,
				 const hlong Nelements,
				 const occa::memory &o_ggeo,
				 const occa::memory &o_Dmatrices,
				 const occa::memory &o_Smatrices,
				 const occa::memory &o_MM,
				 const dfloat lambda,
				 const occa::memory &o_q,
				 occa::memory &o_Aq,
				 const occa::memory &o_ggeoNoJW
				 ){
  

  const dfloat *D    = (dfloat*) o_Dmatrices.ptr();
  const dfloat *S    = (dfloat*) o_Smatrices.ptr();
  const dfloat *MM   = (dfloat*) o_MM.ptr();

  const dfloat *q    = (dfloat*)o_q.ptr();
  const dfloat *ggeo = (dfloat*)o_ggeo.ptr();
  //  const dfloat *ggeo = (dfloat*)o_ggeoNoJW.ptr();
  dfloat *Aq  = (dfloat*)o_Aq.ptr();
  
  switch(Nq){
  case  2: ellipticSerialAxHexKernel3D <  2 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  3: ellipticSerialAxHexKernel3D <  3 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  4: ellipticSerialAxHexKernel3D <  4 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  5: ellipticSerialAxHexKernel3D <  5 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  6: ellipticSerialAxHexKernel3D <  6 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  7: ellipticSerialAxHexKernel3D <  7 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  8: ellipticSerialAxHexKernel3D <  8 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case  9: ellipticSerialAxHexKernel3D <  9 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case 10: ellipticSerialAxHexKernel3D < 10 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case 11: ellipticSerialAxHexKernel3D < 11 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;
  case 12: ellipticSerialAxHexKernel3D < 12 > (Nelements, ggeo, D, S, MM, lambda, q, Aq); break;

  }
      
}
