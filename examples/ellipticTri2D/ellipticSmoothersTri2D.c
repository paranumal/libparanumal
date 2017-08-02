#include "ellipticTri2D.h"

void overlappingPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t*) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = precon->o_zP;

  occaTimerTic(mesh->device,"OverlappingPatchKernel");
  precon->overlappingPatchKernel(mesh->Nelements,
                                 mesh->o_vmapP,
                                 precon->o_oasForwardDgT,
                                 precon->o_oasDiagInvOpDg,
                                 precon->o_oasBackDgT,
                                 o_r,
                                 o_zP);
  ellipticParallelGatherScatterTri2D(mesh, precon->ogsDg, o_zP, o_zP, solver->type, "add");
  solver->dotMultiplyKernel(mesh->NpP*mesh->Nelements,precon->o_invDegreeDGP,o_zP,o_zP);
  occaTimerToc(mesh->device,"OverlappingPatchKernel");

  // extract block interiors on DEVICE
  occaTimerTic(mesh->device,"restrictKernel");
  precon->restrictKernel(mesh->Nelements, o_zP, o_Sr);
  occaTimerToc(mesh->device,"restrictKernel");
}

void exactFullPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t*) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = precon->o_zP;

  occaTimerTic(mesh->device,"exactFullPatchSolveKernel");
  precon->exactPatchSolverKernel(mesh->Nelements,
                            precon->o_invAP,
                            mesh->o_EToE,
                            precon->o_invDegreeAP,
                            o_r,
                            o_zP);
  solver->precon->patchGatherKernel(mesh->Nelements, mesh->o_EToE, mesh->o_EToF, o_zP, o_Sr);
  occaTimerToc(mesh->device,"exactFullPatchSolveKernel");
}

void approxFullPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t*) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = precon->o_zP;

  occaTimerTic(mesh->device,"approxFullPatchSolveKernel");
  precon->approxPatchSolverKernel(mesh->Nelements,
                            precon->o_patchesIndex,
                            precon->o_invAP,
                            mesh->o_EToE,
                            precon->o_invDegreeAP,
                            o_r,
                            o_zP);
  solver->precon->patchGatherKernel(mesh->Nelements, mesh->o_EToE, mesh->o_EToF, o_zP, o_Sr);
  occaTimerToc(mesh->device,"approxFullPatchSolveKernel");
}

void dampedJacobi(void **args, occa::memory o_r, occa::memory o_x, bool x_is_zero) {

  solver_t *solver = (solver_t *) args[0];
  occa::memory *o_invDiagA = (occa::memory *) args[1];
  dfloat *lambda = (dfloat *) args[2];
  char* options = (char*) args[3];

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  if (x_is_zero) {
    solver->dotMultiplyKernel(mesh->Np*mesh->Nelements,*o_invDiagA,o_r,o_x);
    return;
  }

  dfloat one = 1.; dfloat mone = -1.;

  //res = r-Ax
  ellipticOperator2D(solver, *lambda, o_x, solver->o_res, options);
  ellipticScaledAdd(solver, one, o_r, mone, solver->o_res);

  //smooth the fine problem x = x + S(r-Ax)
  solver->dotMultiplyKernel(mesh->Np*mesh->Nelements,*o_invDiagA,solver->o_res,solver->o_res);
  ellipticScaledAdd(solver, one, solver->o_res, one, o_x);
}