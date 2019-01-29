
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
void ellipticSerialUpdate0NBFPCGKernel(const hlong Nelements,
				      const dfloat * __restrict__ cpu_invDegree,
				      const dfloat * __restrict__ cpu_u,
				      const dfloat * __restrict__ cpu_r,
				      const dfloat * __restrict__ cpu_w,
				      dfloat * __restrict__ localdots){
  

#define p_Np (p_Nq*p_Nq*p_Nq)

  cpu_u  = (dfloat*)__builtin_assume_aligned(cpu_u,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_w  = (dfloat*)__builtin_assume_aligned(cpu_w,  USE_OCCA_MEM_BYTE_ALIGN) ;

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN) ;

  dfloat udotr = 0;
  dfloat udotw = 0;
  
  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong n = e*p_Np+i;

      dfloat un = cpu_u[n];
      dfloat rn = cpu_r[n];
      dfloat wn = cpu_w[n];

      dfloat invDeg = cpu_invDegree[n];
      
      udotr += un*rn*invDeg;
      udotw += un*wn*invDeg;
    }
  }

  localdots[0] = udotr;
  localdots[1] = udotw;

#undef p_Np
}
				     
void ellipticSerialUpdate0NBFPCG(const int Nq, const hlong Nelements,
				 occa::memory &o_invDegree,
				 occa::memory &o_u, occa::memory &o_r, occa::memory &o_w,
				 dfloat * __restrict__ localdots){
  
  const dfloat * __restrict__ cpu_u  = (dfloat*)__builtin_assume_aligned(o_u.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_r  = (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_w  = (dfloat*)__builtin_assume_aligned(o_w.ptr(), USE_OCCA_MEM_BYTE_ALIGN); 
  const dfloat * __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(o_invDegree.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  
  switch(Nq){
  case  2: ellipticSerialUpdate0NBFPCGKernel <  2 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break; 
  case  3: ellipticSerialUpdate0NBFPCGKernel <  3 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  4: ellipticSerialUpdate0NBFPCGKernel <  4 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  5: ellipticSerialUpdate0NBFPCGKernel <  5 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  6: ellipticSerialUpdate0NBFPCGKernel <  6 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  7: ellipticSerialUpdate0NBFPCGKernel <  7 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  8: ellipticSerialUpdate0NBFPCGKernel <  8 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case  9: ellipticSerialUpdate0NBFPCGKernel <  9 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case 10: ellipticSerialUpdate0NBFPCGKernel < 10 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case 11: ellipticSerialUpdate0NBFPCGKernel < 11 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  case 12: ellipticSerialUpdate0NBFPCGKernel < 12 > (Nelements, cpu_invDegree, cpu_u, cpu_r, cpu_w, localdots); break;
  }

}

void ellipticNonBlockingUpdate0NBFPCG(elliptic_t *elliptic,
				      occa::memory &o_u, occa::memory &o_r, occa::memory &o_w,
				      dfloat *localdots, dfloat *globaldots, MPI_Request *request){
  
  const cgOptions_t cgOptions = elliptic->cgOptions;
  
  mesh_t *mesh = elliptic->mesh;

  localdots[0] = 0;
  localdots[1] = 0;
  
  if(cgOptions.serial==1 && cgOptions.continuous==1){
    
    ellipticSerialUpdate0NBFPCG(mesh->Nq, mesh->Nelements, elliptic->o_invDegree, o_u, o_r, o_w, localdots);
    
  }
  else{
    if(!cgOptions.continuous){ // e.g. IPDG
      printf("EXITING - NBPCG NOT IMPLEMENTED FOR IPDG\n");
      MPI_Finalize();
      exit(-1);
    }else{
    
      // (u.r)
      // (u.w)
      elliptic->update0NBFPCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
				    elliptic->o_invDegree, o_u, o_r, o_w,
				    elliptic->o_tmpudotr, elliptic->o_tmpudotw);
      
      elliptic->o_tmpudotr.copyTo(elliptic->tmpudotr);
      elliptic->o_tmpudotw.copyTo(elliptic->tmpudotw);
      
      for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
	localdots[0] += elliptic->tmpudotr[n];
	localdots[1] += elliptic->tmpudotw[n];
      }
    }
  }
  
  globaldots[0] = 1;
  globaldots[1] = 1;
  if(cgOptions.enableReductions)      
    MPI_Iallreduce(localdots, globaldots, 2, MPI_DFLOAT, MPI_SUM, mesh->comm, request);

}


// PART 1

template < int p_Nq >
void ellipticSerialUpdate1NBFPCGKernel(const hlong Nelements,
				       const dfloat * __restrict__ cpu_invDegree,
				       const dfloat * __restrict__ cpu_p,
				       const dfloat * __restrict__ cpu_s,
				       const dfloat * __restrict__ cpu_q,
				       const dfloat * __restrict__ cpu_z,
				       const dfloat alpha,
				       dfloat * __restrict__ cpu_x,
				       dfloat * __restrict__ cpu_r,
				       dfloat * __restrict__ cpu_u,
				       dfloat * __restrict__ cpu_w,
				       dfloat * __restrict__ localdots){
  
#define p_Np (p_Nq*p_Nq*p_Nq)

  cpu_p  = (dfloat*)__builtin_assume_aligned(cpu_p,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_s  = (dfloat*)__builtin_assume_aligned(cpu_s,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_q  = (dfloat*)__builtin_assume_aligned(cpu_q,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_z  = (dfloat*)__builtin_assume_aligned(cpu_z,  USE_OCCA_MEM_BYTE_ALIGN) ;

  cpu_x  = (dfloat*)__builtin_assume_aligned(cpu_x,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_r  = (dfloat*)__builtin_assume_aligned(cpu_r,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_u  = (dfloat*)__builtin_assume_aligned(cpu_u,  USE_OCCA_MEM_BYTE_ALIGN) ;
  cpu_w  = (dfloat*)__builtin_assume_aligned(cpu_w,  USE_OCCA_MEM_BYTE_ALIGN) ;

  cpu_invDegree = (dfloat*)__builtin_assume_aligned(cpu_invDegree,  USE_OCCA_MEM_BYTE_ALIGN) ;
  
  dfloat udotr = 0, udots = 0, udotw = 0;
  
  for(hlong e=0;e<Nelements;++e){
    for(int i=0;i<p_Np;++i){
      const hlong n = e*p_Np+i;

      dfloat xn = cpu_x[n];
      dfloat rn = cpu_r[n];
      dfloat un = cpu_u[n];
      dfloat wn = cpu_w[n];

      dfloat pn = cpu_p[n];
      dfloat sn = cpu_s[n];
      dfloat qn = cpu_q[n];
      dfloat zn = cpu_z[n];

      xn = xn + alpha*pn;
      rn = rn - alpha*sn;
      un = un - alpha*qn;
      wn = wn - alpha*zn;
      
      dfloat invDeg = cpu_invDegree[n];

      udotr += un*rn*invDeg;
      udots += un*sn*invDeg;
      udotw += un*wn*invDeg;

      cpu_x[n] = xn;
      cpu_r[n] = rn;
      cpu_u[n] = un;
      cpu_w[n] = wn;
    }
  }

  localdots[0] = udotr;
  localdots[1] = udots;
  localdots[2] = udotw;

#undef p_Np
}

void ellipticSerialUpdate1NBFPCG(const int Nq, const hlong Nelements,
				 occa::memory &o_invDegree,
				 occa::memory &o_p, occa::memory &o_s, occa::memory &o_q, occa::memory &o_z, 
				 const dfloat alpha,
				 occa::memory &o_x, occa::memory &o_r, occa::memory &o_u, occa::memory &o_w, 
				 dfloat *localdots){

  const dfloat * __restrict__ cpu_p  = (dfloat*)__builtin_assume_aligned(o_p.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_s  = (dfloat*)__builtin_assume_aligned(o_s.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_q  = (dfloat*)__builtin_assume_aligned(o_q.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  const dfloat * __restrict__ cpu_z  = (dfloat*)__builtin_assume_aligned(o_z.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  
  dfloat * __restrict__ cpu_x = (dfloat*)__builtin_assume_aligned(o_x.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat * __restrict__ cpu_r = (dfloat*)__builtin_assume_aligned(o_r.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat * __restrict__ cpu_u = (dfloat*)__builtin_assume_aligned(o_u.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  dfloat * __restrict__ cpu_w = (dfloat*)__builtin_assume_aligned(o_w.ptr(), USE_OCCA_MEM_BYTE_ALIGN);
  
  const dfloat * __restrict__ cpu_invDegree = (dfloat*)__builtin_assume_aligned(o_invDegree.ptr(), USE_OCCA_MEM_BYTE_ALIGN) ;
  
  switch(Nq){
  case  2: ellipticSerialUpdate1NBFPCGKernel <  2 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break; 
  case  3: ellipticSerialUpdate1NBFPCGKernel <  3 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  4: ellipticSerialUpdate1NBFPCGKernel <  4 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  5: ellipticSerialUpdate1NBFPCGKernel <  5 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  6: ellipticSerialUpdate1NBFPCGKernel <  6 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  7: ellipticSerialUpdate1NBFPCGKernel <  7 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  8: ellipticSerialUpdate1NBFPCGKernel <  8 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case  9: ellipticSerialUpdate1NBFPCGKernel <  9 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case 10: ellipticSerialUpdate1NBFPCGKernel < 10 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case 11: ellipticSerialUpdate1NBFPCGKernel < 11 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  case 12: ellipticSerialUpdate1NBFPCGKernel < 12 > (Nelements, cpu_invDegree, cpu_p, cpu_s, cpu_q, cpu_z, alpha, cpu_x, cpu_r, cpu_u, cpu_w, localdots); break;
  }
}

void ellipticNonBlockingUpdate1NBFPCG(elliptic_t *elliptic,
				      occa::memory &o_p, occa::memory &o_s, occa::memory &o_q, occa::memory &o_z,
				      const dfloat alpha,
				      occa::memory &o_x, occa::memory &o_r, occa::memory &o_u, occa::memory &o_w,
				      dfloat *localdots, dfloat *globaldots, MPI_Request *request){
  
  const cgOptions_t cgOptions = elliptic->cgOptions;
  
  mesh_t *mesh = elliptic->mesh;

  if(cgOptions.serial==1 && cgOptions.continuous==1){
    
    ellipticSerialUpdate1NBFPCG(mesh->Nq, mesh->Nelements, 
				elliptic->o_invDegree,
				o_p, o_s, o_q, o_z,
				alpha,
				o_x, o_r, o_u, o_w,
				localdots);
  }
  else{
    if(!cgOptions.continuous){ // e.g. IPDG
      printf("EXITING - NBPCG NOT IMPLEMENTED FOR IPDG\n");
      MPI_Finalize();
      exit(-1);
    }else{
    
      // p <= z + beta*p
      // s <= Z + beta*s
      // dot(p,s)
      elliptic->update1NBFPCGKernel(mesh->Nelements*mesh->Np, elliptic->NblocksUpdatePCG,
				    elliptic->o_invDegree, o_p, o_s, o_q, o_z, alpha, o_x, o_r, o_u, o_w,
				    elliptic->o_tmpudotr, elliptic->o_tmpudots, elliptic->o_tmpudotw);
      
      elliptic->o_tmpudotr.copyTo(elliptic->tmpudotr);
      elliptic->o_tmpudots.copyTo(elliptic->tmpudots);
      elliptic->o_tmpudotw.copyTo(elliptic->tmpudotw);
      
      localdots[0] = 0;
      localdots[1] = 0;
      localdots[2] = 0;
      for(int n=0;n<elliptic->NblocksUpdatePCG;++n){
	localdots[0] += elliptic->tmpudotr[n];
	localdots[1] += elliptic->tmpudots[n];
	localdots[2] += elliptic->tmpudotw[n];
      }
    }
  }
  
  if(cgOptions.enableReductions)      
    MPI_Iallreduce(localdots, globaldots, 3, MPI_DFLOAT, MPI_SUM, mesh->comm, request);
  else{
    globaldots[0] = 1;
    globaldots[1] = 1;
  }
}
