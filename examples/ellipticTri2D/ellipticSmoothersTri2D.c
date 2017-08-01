#include "ellipticTri2D.h"

void overlappingPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = solver->o_zP;

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

  solver_t *solver = args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = solver->o_zP;

  occaTimerTic(mesh->device,"exactFullPatchSolveKernel");
  precon->patchSolverKernel(mesh->Nelements,
                            precon->o_patchesIndex,
                            precon->o_invAP,
                            mesh->o_EToE,
                            precon->o_invDegreeAP,
                            o_r,
                            solver->o_zP);
#if 0
  meshParallelGather(mesh, precon->hgsDg, solver->o_zP, o_Sr);
#else
  solver->precon->patchGatherKernel(mesh->Nelements, mesh->o_EToE, mesh->o_EToF, solver->o_zP, o_Sr);
#endif
  occaTimerToc(mesh->device,"exactFullPatchSolveKernel");
}

void approxFullPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = solver->o_zP;

  occaTimerTic(mesh->device,"approxFullPatchSolveKernel");
  //precon->approxPatchSolverKernel(mesh->Nelements,
  //                          precon->o_invAP,
  //                          mesh->o_EToE,
  //                          mesh->o_EToF,
  //                          mesh->o_rmapP,
  //                          precon->o_invDegreeAP,
  //                          o_r,
  //                          o_zP);
  precon->patchSolverKernel(mesh->Nelements,
                            precon->o_patchesIndex,
                            precon->o_invAP,
                            mesh->o_EToE,
                            precon->o_invDegreeAP,
                            o_r,
                            solver->o_zP);
#if 0
  meshParallelGather(mesh, precon->hgsDg, solver->o_zP, o_Sr);
#else
  solver->precon->patchGatherKernel(mesh->Nelements, mesh->o_EToE, mesh->o_EToF, solver->o_zP, o_Sr);
#endif
  occaTimerToc(mesh->device,"approxFullPatchSolveKernel");
}

void localFullPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = solver->o_zP;

  occaTimerTic(mesh->device,"LocalFullPatchSolveKernel");
  precon->localPatchSolverKernel(mesh->Nelements,
                                  precon->o_invAP,
                                  mesh->o_EToE,
                                  o_r,
                                  o_Sr);
  occaTimerToc(mesh->device,"LocalFullPatchSolveKernel");
}