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

#include "adaptive.h"

template < int p_Nq >
void adaptiveSerialScaledAddKernel(const hlong Nelements, const dfloat alpha, const dfloat * __restrict__ cpu_a, const dfloat beta, dfloat * __restrict__ cpu_b){

  cpu_a = (dfloat*)__builtin_assume_aligned(cpu_a, USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_b = (dfloat*)__builtin_assume_aligned(cpu_b, USE_OCCA_MEM_BYTE_ALIGN) ;
  
#define p_Np (p_Nq*p_Nq*p_Nq)

  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const dfloat ai = cpu_a[e*p_Np+i];
      const dfloat bi = cpu_b[e*p_Np+i];
      cpu_b[e*p_Np+i] = alpha*ai + beta*bi;
    }
  }

#undef p_Np
}

void adaptiveSerialScaledAdd(const int Nq, const hlong Nelements, const dfloat alpha, occa::memory &o_a, const dfloat beta, occa::memory &o_b){

  const dfloat * __restrict__ cpu_a = (dfloat*)__builtin_assume_aligned(o_a.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  dfloat * __restrict__ cpu_b = (dfloat*)__builtin_assume_aligned(o_b.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  
  switch(Nq){
  case  2: adaptiveSerialScaledAddKernel <  2 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  3: adaptiveSerialScaledAddKernel <  3 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  4: adaptiveSerialScaledAddKernel <  4 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  5: adaptiveSerialScaledAddKernel <  5 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  6: adaptiveSerialScaledAddKernel <  6 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  7: adaptiveSerialScaledAddKernel <  7 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  8: adaptiveSerialScaledAddKernel <  8 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case  9: adaptiveSerialScaledAddKernel <  9 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  case 10: adaptiveSerialScaledAddKernel < 10 > (Nelements, alpha, cpu_a, beta, cpu_b); break;
  }
}

void adaptiveScaledAdd(adaptive_t *adaptive, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  mesh_t *mesh = adaptive->mesh;

  setupAide &options = adaptive->options;
  
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  
  dlong Ntotal = mesh->Nelements*mesh->Np;

  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  if(serial == 1 && continuous==1){
    adaptiveSerialScaledAdd(mesh->Nq, mesh->Nelements, alpha, o_a, beta, o_b);
    return;
  }
  
  // if not Serial
  adaptive->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);
}
