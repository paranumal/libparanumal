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

template < int p_Nq >
dfloat ellipticSerialWeightedInnerProductKernel(const hlong Nelements,
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

dfloat ellipticSerialWeightedInnerProduct(const int Nq, const hlong Nelements, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  const dfloat * __restrict__ cpu_a = (dfloat*)__builtin_assume_aligned(o_a.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_b = (dfloat*)__builtin_assume_aligned(o_b.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_w = (dfloat*)__builtin_assume_aligned(o_w.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  switch(Nq){
  case  2: return ellipticSerialWeightedInnerProductKernel <  2 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  3: return ellipticSerialWeightedInnerProductKernel <  3 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  4: return ellipticSerialWeightedInnerProductKernel <  4 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  5: return ellipticSerialWeightedInnerProductKernel <  5 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  6: return ellipticSerialWeightedInnerProductKernel <  6 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  7: return ellipticSerialWeightedInnerProductKernel <  7 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  8: return ellipticSerialWeightedInnerProductKernel <  8 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case  9: return ellipticSerialWeightedInnerProductKernel <  9 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 10: return ellipticSerialWeightedInnerProductKernel < 10 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 11: return ellipticSerialWeightedInnerProductKernel < 11 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  case 12: return ellipticSerialWeightedInnerProductKernel < 12 > (Nelements, cpu_w, cpu_a, cpu_b); break;
  }

  return -99;
}

dfloat ellipticWeightedInnerProduct(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b){

  setupAide &options = elliptic->options;

  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int enableReductions = 1;
  // options.getArgs("DEBUG ENABLE REDUCTIONS", enableReductions);


  mesh_t *mesh = elliptic->mesh;
  dfloat *tmp = elliptic->tmp;
  dlong Nblock = elliptic->Nblock;
  dlong Nblock2 = elliptic->Nblock2;
  dlong Ntotal = mesh->Nelements*mesh->Np;

  if(serial==1 && continuous==1){

    dfloat wab = ellipticSerialWeightedInnerProduct(mesh->Nq, mesh->Nelements, o_w, o_a, o_b);
    dfloat globalwab = 0;

    MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    return globalwab;
  }

  occa::memory &o_tmp = elliptic->o_tmp;
  occa::memory &o_tmp2 = elliptic->o_tmp2;

  if(continuous==1)
    elliptic->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    elliptic->innerProductKernel(Ntotal, o_a, o_b, o_tmp);

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


