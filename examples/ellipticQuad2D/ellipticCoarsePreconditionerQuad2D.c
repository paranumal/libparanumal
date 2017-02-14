#include "ellipticQuad2D.h"

void ellipticCoarsePreconditionerQuad2D(mesh_t *mesh, precon_t *precon,
					occa::memory &o_r, occa::memory &o_z){

  precon->coarsenKernel(mesh->Nelements, precon->o_V1, o_r, precon->o_r1);

  precon->o_r1.copyTo(precon->r1);
  
  xxtSolve(precon->z1, precon->xxt, precon->r1);  

  precon->o_z1.copyTo(precon->z1);

  // adds in coarse contribution
  precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, o_z);

  
}
