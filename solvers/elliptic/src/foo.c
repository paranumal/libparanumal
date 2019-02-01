
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

// hack
#define p_Nggeo 7 

template < int p_Nq >
void ellipticSerialPartialAxHexKernel3D (const hlong Nelements,
					 const occa::memory &o_elementList,
					 const occa::memory &o_ggeo,
					 const occa::memory &o_Dmatrices,
					 const occa::memory &o_Smatrices,
					 const occa::memory &o_MM,
					 const dfloat lambda,
					 const occa::memory &o_q,
					 occa::memory &o_Aq){
  
#define p_Np (p_Nq*p_Nq*p_Nq)
  
  dfloat *Aq = (dfloat*) o_Aq.ptr();

  const dlong  *elementList = (dlong*) o_elementList.ptr();
  
  const dfloat *ggeo  = (dfloat*) o_ggeo.ptr();
  const dfloat *D     = (dfloat*) o_Dmatrices.ptr();
  const dfloat *S     = (dfloat*) o_Smatrices.ptr();
  const dfloat *q  = (dfloat*) o_q.ptr();

  dfloat s_q[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqr[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqs[p_Nq][p_Nq][p_Nq];
  dfloat s_Gqt[p_Nq][p_Nq][p_Nq];
  
  for(dlong e=0; e<Nelements; ++e){
    
    const dlong element = elementList[e];

    for(int k = 0; k < p_Nq; k++) {
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){
          const dlong base = i + j*p_Nq + element*p_Np;
          const dfloat qbase = q[base + k*p_Nq*p_Nq];
          s_q[k][j][i] = qbase;
        }
      }
    }

    // Layer by layer                                                                                                                                                                          
    for(int k = 0;k < p_Nq; k++){
      for(int j=0;j<p_Nq;++j){
        for(int i=0;i<p_Nq;++i){

          const dlong gbase = element*p_Nggeo*p_Np + k*p_Nq*p_Nq + j*p_Nq + i;

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
            qr += S[m*p_Nq+i]*s_q[k][j][m];
            qs += S[m*p_Nq+j]*s_q[k][m][i];	    
            qt += S[m*p_Nq+k]*s_q[m][j][i];
          }

          s_Gqr[k][j][i] = (r_G00*qr + r_G01*qs + r_G02*qt);
          s_Gqs[k][j][i] = (r_G01*qr + r_G11*qs + r_G12*qt);
          s_Gqt[k][j][i] = (r_G02*qr + r_G12*qs + r_G22*qt);
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
	  
          for(int m = 0; m < p_Nq; m++){
            r_Aqr += D[m*p_Nq+i]*s_Gqr[k][j][m];
            r_Aqs += D[m*p_Nq+j]*s_Gqs[k][m][i];
            r_Aqt += D[m*p_Nq+k]*s_Gqt[m][j][i];
          }

          const dlong id = element*p_Np +k*p_Nq*p_Nq+ j*p_Nq + i;
          Aq[id] = r_Aqr + r_Aqs + r_Aqt; // +r_Aq;
        }
      }
    }
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


  switch(Nq){
  case  2: ellipticSerialPartialAxHexKernel3D <  2 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  3: ellipticSerialPartialAxHexKernel3D <  3 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  4: ellipticSerialPartialAxHexKernel3D <  4 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  5: ellipticSerialPartialAxHexKernel3D <  5 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  6: ellipticSerialPartialAxHexKernel3D <  6 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  7: ellipticSerialPartialAxHexKernel3D <  7 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  8: ellipticSerialPartialAxHexKernel3D <  8 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case  9: ellipticSerialPartialAxHexKernel3D <  9 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case 10: ellipticSerialPartialAxHexKernel3D < 10 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case 11: ellipticSerialPartialAxHexKernel3D < 11 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  case 12: ellipticSerialPartialAxHexKernel3D < 12 > (Nelements, o_elementList, o_ggeo, o_Dmatrices, o_Smatrices, o_MM, lambda, o_q, o_Aq); break;
  }
      
}
