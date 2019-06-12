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
dfloat ellipticSerialUpdatePCGKernel(const hlong Nelements,
				     const dfloat * __restrict__ cpu_invDegree,
				     const dfloat * __restrict__ cpu_p,
				     const dfloat * __restrict__ cpu_Ap,
				     const dfloat alpha,
				     dfloat * __restrict__ cpu_x,
				     dfloat * __restrict__ cpu_r){

#define p_Np (p_Nq*p_Nq*p_Nq)

  cpu_p  = (dfloat*)__builtin_assume_aligned(cpu_p,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_Ap = (dfloat*)__builtin_assume_aligned(cpu_Ap, USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_x  = (dfloat*)__builtin_assume_aligned(cpu_x,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN) ;

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN) ;
  
  dfloat rdotr = 0;
  
  cpu_p = (dfloat*)__builtin_assume_aligned(cpu_p, USE_OCCA_MEM_BYTE_ALIGN) ;

  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong n = e*p_Np+i;
      cpu_x[n] += alpha*cpu_p[n];

      const dfloat rn = cpu_r[n] - alpha*cpu_Ap[n];
      rdotr += rn*rn*cpu_invDegree[n];
      cpu_r[n] = rn;
    }
  }

#undef p_Np
  
  return rdotr;
}
				     
dfloat ellipticSerialUpdatePCG(const int Nq, const hlong Nelements,
			       occa::memory &o_invDegree, occa::memory &o_p, occa::memory &o_Ap, const dfloat alpha,
			       occa::memory &o_x, occa::memory &o_r){

  const dfloat * __restrict__ cpu_p  = (dfloat*)__builtin_assume_aligned(o_p.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_Ap = (dfloat*)__builtin_assume_aligned(o_Ap.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(o_invDegree.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat * __restrict__ cpu_x  = (dfloat*)__builtin_assume_aligned(o_x.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  dfloat * __restrict__ cpu_r  = (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat rdotr = 0;
  
  switch(Nq){
  case  2: rdotr = ellipticSerialUpdatePCGKernel <  2 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break; 
  case  3: rdotr = ellipticSerialUpdatePCGKernel <  3 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  4: rdotr = ellipticSerialUpdatePCGKernel <  4 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  5: rdotr = ellipticSerialUpdatePCGKernel <  5 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  6: rdotr = ellipticSerialUpdatePCGKernel <  6 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  7: rdotr = ellipticSerialUpdatePCGKernel <  7 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  8: rdotr = ellipticSerialUpdatePCGKernel <  8 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case  9: rdotr = ellipticSerialUpdatePCGKernel <  9 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case 10: rdotr = ellipticSerialUpdatePCGKernel < 10 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case 11: rdotr = ellipticSerialUpdatePCGKernel < 11 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  case 12: rdotr = ellipticSerialUpdatePCGKernel < 12 > (Nelements, cpu_invDegree, cpu_p, cpu_Ap, alpha, cpu_x, cpu_r); break;
  }

  return rdotr;
}

dfloat ellipticUpdatePCG(elliptic_t *elliptic,
			 occa::memory &o_p, occa::memory &o_Ap, const dfloat alpha,
			 occa::memory &o_x, occa::memory &o_r){

  setupAide &options = elliptic->options;
  
  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  mesh_t *mesh = elliptic->mesh;

  if(serial==1 && continuous==1){
    
    dfloat rdotr1 = ellipticSerialUpdatePCG(mesh->Nq, mesh->Nelements, 
					    elliptic->o_invDegree,
					    o_p, o_Ap, alpha, o_x, o_r);

    dfloat globalrdotr1 = 0;
    if(enableReductions)      
      MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
    else
      globalrdotr1 = 1;
    
    return globalrdotr1;
  }
  
  dfloat rdotr1 = 0;
  
  if(!continuous){
    
    // x <= x + alpha*p
    ellipticScaledAdd(elliptic,  alpha, o_p,  1.f, o_x);
    
    // [
    // r <= r - alpha*A*p
    ellipticScaledAdd(elliptic, -alpha, o_Ap, 1.f, o_r);
    
    // dot(r,r)
    if(enableReductions)
      rdotr1 = ellipticWeightedNorm2(elliptic, elliptic->o_invDegree, o_r);
    else
      rdotr1 = 1;
  }else{
    
    // x <= x + alpha*p
    // r <= r - alpha*A*p
    // dot(r,r)
    elliptic->updatePCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
			      elliptic->o_invDegree, o_p, o_Ap, alpha, o_x, o_r, elliptic->o_tmpNormr);

    elliptic->o_tmpNormr.copyTo(elliptic->tmpNormr);

    rdotr1 = 0;
    for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
      rdotr1 += elliptic->tmpNormr[n];
    }
    
    dfloat globalrdotr1 = 0;
    MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    rdotr1 = globalrdotr1;
    
  }

  return rdotr1;
}

