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

void BuildLocalContinuousDiagHex3D (adaptive_t* adaptive, level_t *level, dfloat *geo, dfloat *D,
				    dfloat lambda, dlong eM, dfloat *B, dfloat *Br, dfloat *Bs, dfloat *Bt, dfloat *A);

void adaptiveBuildJacobi(adaptive_t* adaptive, level_t *level, dfloat lambda, dfloat **invDiagA){

  setupAide options = adaptive->options;

  dfloat *ggeo = (dfloat*) calloc(NGGEO*level->Np*level->Klocal, sizeof(dfloat));
  level->o_ggeo.copyTo(ggeo, NGGEO*level->Np*level->Klocal*sizeof(dfloat), 0);

  dfloat *D = (dfloat*) calloc(level->Nq*level->Nq, sizeof(dfloat));
  level->o_D.copyTo(D, level->Nq*level->Nq*sizeof(dfloat));
  
  dlong diagNnum = level->Np*level->Klocal;
  
  dfloat *diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  
  if(adaptive->rank==0) printf("Building diagonal...");fflush(stdout);
  
#pragma omp parallel for 
  for(dlong eM=0;eM<level->Klocal;++eM)
    BuildLocalContinuousDiagHex3D(adaptive, level, ggeo, D, lambda, eM, B, Br, Bs, Bt, diagA + eM*level->Np);

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) 
    ogsGatherScatter(diagA, ogsDfloat, ogsAdd, level->ogs);
  
  *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
  for (dlong n=0;n<level->Klocal*level->Np;n++) {
    (*invDiagA)[n] = 1/diagA[n];
  }

  if(adaptive->rank==0) printf("done.\n");

  free(diagA);
  free(B); free(Br); free(Bs); free(Bt);
}

void BuildLocalContinuousDiagHex3D(adaptive_t* adaptive, level_t *level, dfloat *ggeo, dfloat *D,
				   dfloat lambda, dlong eM, dfloat *A) {

  for (int nz=0;nz<level->Nq;nz++) {
    for (int ny=0;ny<level->Nq;ny++) {
      for (int nx=0;nx<level->Nq;nx++) {

	int idn = nx+ny*level->Nq+nz*level->Nq*level->Nq;

	if(1){ // TW: fix later
	  A[idn] = 0;

	  int id = nx+ny*level->Nq+nz*level->Nq*level->Nq;
	  dlong base = eM*level->Np*level->Nggeo;
    
	  dfloat Grs = level->ggeo[base + id + GGEO_RS*level->Np];
	  A[idn] += 2*Grs*D[nx+nx*level->Nq]*D[ny+ny*level->Nq];

	  dfloat Grt = level->ggeo[base + id + GGEO_RT*level->Np];
	  A[idn] += 2*Grt*D[nx+nx*level->Nq]*D[nz+nz*level->Nq];
    
	  dfloat Gst = level->ggeo[base + id + GGEO_ST*level->Np];
	  A[idn] += 2*Gst*D[ny+ny*level->Nq]*D[nz+nz*level->Nq];

	  for (int k=0;k<level->Nq;k++) {
	    int iid = k+ny*level->Nq+nz*level->Nq*level->Nq;
	    dfloat Grr = level->ggeo[base + iid + GGEO_RR*level->Np];
	    A[idn] += Grr*D[nx+k*level->Nq]*D[nx+k*level->Nq];
	  }

	  for (int k=0;k<level->Nq;k++) {
	    int iid = nx+k*level->Nq+nz*level->Nq*level->Nq;
	    dfloat Gss = level->ggeo[base + iid + GGEO_SS*level->Np];
	    A[idn] += Gss*D[ny+k*level->Nq]*D[ny+k*level->Nq];
	  }
    
	  for (int k=0;k<level->Nq;k++) {
	    int iid = nx+ny*level->Nq+k*level->Nq*level->Nq;
	    dfloat Gtt = level->ggeo[base + iid + GGEO_TT*level->Np];
	    A[idn] += Gtt*D[nz+k*level->Nq]*D[nz+k*level->Nq];
	  }
      
	  dfloat JW = level->ggeo[base + id + GGEO_JW*level->Np];
	  A[idn] += JW*lambda;
	} else {
	  A[idn] = 1; //just put a 1 so A is invertable
	}
      }
    }
  }

#if 0
  //add the rank boost for the allNeumann Poisson problem
  if (level->allNeumann) {
    for(int n=0;n<level->Np;++n){
      if (adaptive->mapB[n+eM*level->Np]!=1) { //dont fill rows for masked nodes
        A[n] += adaptive->allNeumannPenalty*adaptive->allNeumannScale*adaptive->allNeumannScale;
      } 
    }
  }
#endif
}
