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

void smoothChebyshevTri2D(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero) {

  solver_t *solver = (solver_t *) args[0];
  agmgLevel *level = (agmgLevel *) args[1];

  dfloat lambdaN = level->smoother_params[0];
  dfloat lambda1 = level->smoother_params[1];

  dfloat theta = 0.5*(lambdaN+lambda1);
  dfloat delta = 0.5*(lambdaN-lambda1);
  dfloat invTheta = 1.0/theta;
  dfloat sigma = theta/delta;
  dfloat rho_n = 1./sigma;
  dfloat rho_np1;

  dfloat one = 1., mone = -1., zero = 0.0;

  occa::memory o_res = level->o_smootherResidual;
  occa::memory o_Ad  = level->o_smootherResidual2;
  occa::memory o_d   = level->o_smootherUpdate;

  if(xIsZero){ //skip the Ax if x is zero
    //res = Sr
    level->device_smoother(level->smootherArgs, o_r, o_res);

    //d = invTheta*res
    solver->scaledAddKernel(level->Nrows, invTheta, o_res, zero, o_d);
  } else {
    //res = S(r-Ax)
    level->device_Ax(level->AxArgs,o_x,o_res);
    solver->scaledAddKernel(level->Nrows, one, o_r, mone, o_res);
    level->device_smoother(level->smootherArgs, o_res, o_res);

    //d = invTheta*res
    solver->scaledAddKernel(level->Nrows, invTheta, o_res, zero, o_d);
  }

  for (int k=0;k<level->ChebyshevIterations;k++) {
    //x_k+1 = x_k + d_k
    if (xIsZero&&(k==0))
      solver->scaledAddKernel(level->Nrows, one, o_d, zero, o_x);
    else
      solver->scaledAddKernel(level->Nrows, one, o_d, one, o_x);

    //r_k+1 = r_k - SAd_k
    level->device_Ax(level->AxArgs,o_d,o_Ad);
    level->device_smoother(level->smootherArgs, o_Ad, o_Ad);
    solver->scaledAddKernel(level->Nrows, mone, o_Ad, one, o_res);

    rho_np1 = 1.0/(2.*sigma-rho_n);
    dfloat rhoDivDelta = 2.0*rho_np1/delta;

    //d_k+1 = rho_k+1*rho_k*d_k  + 2*rho_k+1*r_k+1/delta
    solver->scaledAddKernel(level->Nrows, rhoDivDelta, o_res, rho_np1*rho_n, o_d);

    rho_n = rho_np1;
  }
  //x_k+1 = x_k + d_k
  solver->scaledAddKernel(level->Nrows, one, o_d, one, o_x);

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
  ellipticParallelGatherScatterTri2D(mesh, precon->ogsDg, o_zP, solver->type, "add");
  solver->dotMultiplyKernel(mesh->NpP*mesh->Nelements,precon->o_invDegreeDGP,o_zP,o_zP);
  occaTimerToc(mesh->device,"OverlappingPatchKernel");

  // extract block interiors on DEVICE
  occaTimerTic(mesh->device,"restrictKernel");
  precon->restrictKernel(mesh->Nelements, o_zP, o_Sr);
  occaTimerToc(mesh->device,"restrictKernel");
}

void FullPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

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

void FacePatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t*) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory o_zP = precon->o_zP;

  occaTimerTic(mesh->device,"approxFacePatchSolveKernel");
  precon->approxFacePatchSolverKernel(mesh->NfacePairs,
                            precon->o_patchesIndex,
                            precon->o_invAP,
                            mesh->o_FPairsToE,
                            precon->o_invDegreeAP,
                            o_r,
                            o_zP);
  solver->precon->facePatchGatherKernel(mesh->Nelements, mesh->o_EToFPairs, o_zP, o_Sr);
  occaTimerToc(mesh->device,"approxFacePatchSolveKernel");
}

void LocalPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t*) args[0];
  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;

  occaTimerTic(mesh->device,"approxBlockJacobiSolveKernel");
  // precon->approxBlockJacobiSolverKernel(mesh->Nelements,
  //                           precon->o_patchesIndex,
  //                           precon->o_invAP,
  //                           precon->o_invDegreeAP,
  //                           o_r,
  //                           o_Sr);
  dfloat tol = 1E-12;
  dfloat *lambda = (dfloat *) args[2];

  precon->CGLocalPatchKernel(mesh->Nelements,
                             mesh->o_vmapM,
                             *lambda,
                             solver->tau,
                             mesh->o_vgeo,
                             mesh->o_sgeo,
                             mesh->o_EToB,
                             mesh->o_DrT,
                             mesh->o_DsT,
                             mesh->o_LIFTT,
                             mesh->o_MM,
                             precon->o_invDegreeAP,
                             o_r,
                             o_Sr,
                             tol);
  occaTimerToc(mesh->device,"approxBlockJacobiSolveKernel");
}

void dampedJacobi(void **args, occa::memory &o_r, occa::memory &o_Sr) {

  solver_t *solver = (solver_t *) args[0];
  char *options = (char *) args[1];
  mesh_t *mesh = solver->mesh;

  occa::memory o_invDiagA = solver->precon->o_invDiagA;

  solver->dotMultiplyKernel(mesh->Np*mesh->Nelements,o_invDiagA,o_r,o_Sr);
}