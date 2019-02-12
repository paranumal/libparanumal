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
dfloat adaptiveSerialWeightedInnerProductKernel(const hlong Nelements,
						const dfloat * __restrict__ cpu_w,
						const dfloat * __restrict__ cpu_a,
						const dfloat * __restrict__ cpu_b){


  cpu_a = (dfloat*)__builtin_assume_aligned(cpu_a, USE_OCCA_MEM_BYTE_ALIGN);
  cpu_b = (dfloat*)__builtin_assume_aligned(cpu_b, USE_OCCA_MEM_BYTE_ALIGN);
  cpu_w = (dfloat*)__builtin_assume_aligned(cpu_w, USE_OCCA_MEM_BYTE_ALIGN);

#define p_Np (p_Nq*p_Nq*p_Nq)

  dfloat wab = 0;

  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong id = e*p_Np+i;
      const dfloat ai = cpu_a[id];
      const dfloat bi = cpu_b[id];
      const dfloat wi = cpu_w[id];
      wab += ai*bi*wi;
    }
  }

  return wab;

#undef p_Np
}

dfloat adaptiveSerialWeightedInnerProduct(const int Nq, const hlong Nelements, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  const dfloat * __restrict__ cpu_a = (dfloat*)__builtin_assume_aligned(o_a.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_b = (dfloat*)__builtin_assume_aligned(o_b.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_w = (dfloat*)__builtin_assume_aligned(o_w.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  switch(Nq){
  case  2: return adaptiveSerialWeightedInnerProductKernel <  2 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  3: return adaptiveSerialWeightedInnerProductKernel <  3 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  4: return adaptiveSerialWeightedInnerProductKernel <  4 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  5: return adaptiveSerialWeightedInnerProductKernel <  5 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  6: return adaptiveSerialWeightedInnerProductKernel <  6 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  7: return adaptiveSerialWeightedInnerProductKernel <  7 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  8: return adaptiveSerialWeightedInnerProductKernel <  8 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  9: return adaptiveSerialWeightedInnerProductKernel <  9 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 10: return adaptiveSerialWeightedInnerProductKernel < 10 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 11: return adaptiveSerialWeightedInnerProductKernel < 11 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 12: return adaptiveSerialWeightedInnerProductKernel < 12 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  }

  return -99;
}

dfloat adaptiveWeightedInnerProduct(adaptive_t *adaptive, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  setupAide &options = adaptive->options;

  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int enableReductions = 1;
  options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);

  
  mesh_t *mesh = adaptive->mesh;
  dfloat *tmp = adaptive->tmp;
  dlong Nblock = adaptive->Nblock;
  dlong Nblock2 = adaptive->Nblock2;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  if(serial==1 && continuous==1){
    
    dfloat wab = adaptiveSerialWeightedInnerProduct(mesh->Nq, mesh->Nelements, o_w, o_a, o_b);
    dfloat globalwab = 0;
    
    MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
    
    return globalwab;
  }
  
  occa::memory &o_tmp = adaptive->o_tmp;
  occa::memory &o_tmp2 = adaptive->o_tmp2;

  if(continuous==1)
    adaptive->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    adaptive->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

  /* add a second sweep if Nblock>Ncutoff */
  dlong Ncutoff = 100;
  dlong Nfinal;
  if(Nblock>Ncutoff){

    mesh->sumKernel(Nblock, o_tmp, o_tmp2);

    o_tmp2.copyTo(tmp);

    Nfinal = Nblock2;
	
  }
  else{
    o_tmp.copyTo(tmp);
    
    Nfinal = Nblock;

  }    

  dfloat wab = 0;
  for(dlong n=0;n<Nfinal;++n){
    wab += tmp[n];
  }

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  return globalwab;
}


