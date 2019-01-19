
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

#include "elliptic.h"

#if USE_STEFAN_MXM==1
extern "C"
{
  void ax_e_(dfloat *w, const dfloat *u, const dfloat *g, dfloat *ur, dfloat *us, dfloat *ut,dfloat *wk);
}
#endif
   
// hack
#define p_Nggeo 7 

template < const int p_Nq >
void ellipticSerialPartialAxHexKernel3D (const hlong Nelements,
					 const dlong  * __restrict__ elementList , // type ?
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
  elementList = (dlong*)__builtin_assume_aligned(elementList, USE_OCCA_MEM_BYTE_ALIGN) ; // type ?

#define p_Np (p_Nq*p_Nq*p_Nq)
  
  dfloat s_q  [p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));

  dfloat s_tmp[p_Nq][p_Nq][p_Nq] __attribute__((aligned(USE_OCCA_MEM_BYTE_ALIGN)));
  
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

  for(dlong e=0; e<Nelements; ++e){
    
    const dlong element = elementList[e];

#if USE_STEFAN_MXM==1
    ax_e_(Aq+element*p_Np, 
      q+element*p_Np,
      ggeo + element*p_Nggeo*p_Np, // note layout is wrong
      s_Gqr[0][0],
      s_Gqs[0][0],
      s_Gqt[0][0],
      s_tmp[0][0]);
#endif
    
#if 1
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

          const dlong gbase = element*p_Nggeo*c_Np + k*p_Nq*p_Nq + j*p_Nq + i;
          const dfloat r_G00 = ggeo[gbase+G00ID*p_Np];
          const dfloat r_G01 = ggeo[gbase+G01ID*p_Np];
          const dfloat r_G11 = ggeo[gbase+G11ID*p_Np];
          const dfloat r_G02 = ggeo[gbase+G02ID*p_Np];
          const dfloat r_G12 = ggeo[gbase+G12ID*p_Np];
          const dfloat r_G22 = ggeo[gbase+G22ID*p_Np];

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
#if 0
          const dlong gbase = element*p_Nggeo*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;
          const dfloat r_GwJ = ggeo[gbase+p_GWJID*p_Np];

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
          Aq[id] = r_Aqr + r_Aqs + r_Aqt; // +r_Aq;
        }
      }
    }
#endif

    // ok
#if 0
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

#undef p_Np


void ellipticSerialPartialAxHexKernel3D(const int Nq,
					const hlong Nelements,
					const occa::memory &o_elementList,
					const occa::memory &o_ggeo,
					const occa::memory &o_Dmatrices,
					const occa::memory &o_Smatrices,
					const occa::memory &o_MM,
					const dfloat lambda,
					const occa::memory &o_q,
					occa::memory &o_Aq){


  const dlong  *elementList = (dlong*) o_elementList.ptr();
  const dfloat *D    = (dfloat*) o_Dmatrices.ptr();
  const dfloat *S    = (dfloat*) o_Smatrices.ptr();
  const dfloat *MM   = (dfloat*) o_MM.ptr();

  const dfloat *q    = (dfloat*)o_q.ptr();
  const dfloat *ggeo = (dfloat*)o_ggeo.ptr();
  dfloat *Aq  = (dfloat*)o_Aq.ptr();
  
  switch(Nq){
  case  2: ellipticSerialPartialAxHexKernel3D <  2 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  3: ellipticSerialPartialAxHexKernel3D <  3 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  4: ellipticSerialPartialAxHexKernel3D <  4 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  5: ellipticSerialPartialAxHexKernel3D <  5 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  6: ellipticSerialPartialAxHexKernel3D <  6 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  7: ellipticSerialPartialAxHexKernel3D <  7 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  8: ellipticSerialPartialAxHexKernel3D <  8 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case  9: ellipticSerialPartialAxHexKernel3D <  9 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case 10: ellipticSerialPartialAxHexKernel3D < 10 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case 11: ellipticSerialPartialAxHexKernel3D < 11 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;
  case 12: ellipticSerialPartialAxHexKernel3D < 12 > (Nelements, elementList, ggeo, D, S, MM, lambda, q, Aq); break;

  }
      
}
