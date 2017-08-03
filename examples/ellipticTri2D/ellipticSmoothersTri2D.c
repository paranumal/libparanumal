#include "ellipticTri2D.h"

void smoothTri2D(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  solver_t *solver = (solver_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  occa::memory o_res = level->o_smootherResidual;

  if (xIsZero) {
    level->device_smoother(level->smootherArgs, o_r, o_x);
    return;
  }

  dfloat one = 1.; dfloat mone = -1.;

  //res = r-Ax
  level->device_Ax(level->AxArgs,o_x,o_res);
  solver->scaledAddKernel(level->Nrows,one, o_r, mone, o_res);

  //smooth the fine problem x = x + S(r-Ax)
  level->device_smoother(level->smootherArgs, o_res, o_res);
  solver->scaledAddKernel(level->Nrows,one, o_res, one, o_x);
}

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

void dampedJacobi(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t *) args[0];
  mesh_t *mesh = solver->mesh;
  occa::memory o_invDiagA = solver->precon->o_invDiagA;

  solver->dotMultiplyKernel(mesh->Np*mesh->Nelements,o_invDiagA,o_r,o_Sr);
}