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

#include "cds.h"

void cdsSolveStep(cds_t *cds, dfloat time, dfloat dt, occa::memory o_U, occa::memory o_S){

  mesh_t *mesh = cds->mesh;
  
  hlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
   
  if(cds->Nsubsteps) {
    cdsSubCycle(cds, time, cds->Nstages, o_U, o_S,  cds->o_NS);
  } else {
    cdsAdvection(cds, time, o_U, o_S, cds->o_NS);
  }    

  cdsHelmholtzRhs(cds, time+dt, cds->Nstages, cds->o_rhsS);

  cdsHelmholtzSolve(cds, time+dt, cds->Nstages, cds->o_rhsS, cds->o_rkS);
   
  //cycle history
  for (int s=cds->Nstages;s>1;s--) {
    o_S.copyFrom( o_S, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
		  (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
		  (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));

    cds->o_NS.copyFrom(cds->o_NS, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
		       (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
		       (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));
  }
  //copy updated scalar
  o_S.copyFrom(cds->o_rkS, cds->NSfields*cds->Ntotal*sizeof(dfloat)); 
}
