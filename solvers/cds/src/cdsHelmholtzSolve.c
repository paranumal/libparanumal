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

// solve -lambda*S + A*S = rhsS
void cdsHelmholtzSolve(cds_t *cds, dfloat time, int stage,occa::memory o_rhsS,occa::memory o_Shat){
  
  mesh_t     *mesh   = cds->mesh; 
  elliptic_t *solver = cds->solver;
  
  int quad3D = (cds->dim==3 && cds->elementType==QUADRILATERALS) ? 1 : 0;  
  
  if (cds->options.compareArgs("DISCRETIZATION","CONTINUOUS")){

    if(!quad3D) 
      cds->helmholtzRhsBCKernel(mesh->Nelements,
                                mesh->o_ggeo,
                                mesh->o_sgeo,
                                mesh->o_Dmatrices,
                                mesh->o_Smatrices,
                                mesh->o_MM,
                                mesh->o_vmapM,
    	                          mesh->o_EToB,
                                mesh->o_sMT,
                                cds->lambda,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                cds->o_mapB,
                                o_rhsS);

    // gather-scatter
    ogsGatherScatter(o_rhsS, ogsDfloat, ogsAdd, mesh->ogs);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_rhsS);

  } else if (cds->options.compareArgs("DISCRETIZATION","IPDG") && !quad3D) {

    occaTimerTic(mesh->device,"velocityRhsIpdg");    
    cds->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  solver->tau,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  o_rhsS);
    occaTimerToc(mesh->device,"velocityRhsIpdg");   
  }

  //copy current solution fields as initial guess? (could use Shat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  o_Shat.copyFrom(cds->o_S,Ntotal*sizeof(dfloat),0,0*cds->sOffset*sizeof(dfloat)); 
 
  if (cds->options.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Shat);
  }
  
  occaTimerTic(mesh->device,"S-Solve");
  cds->Niter = ellipticSolve(solver, cds->lambda, cds->TOL, o_rhsS, o_Shat);
  occaTimerToc(mesh->device,"S-Solve"); 

  if (cds->options.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    cds->helmholtzAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            cds->o_mapB,
                            o_Shat);
  }
}
