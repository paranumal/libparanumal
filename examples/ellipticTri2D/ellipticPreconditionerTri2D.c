#include "ellipticTri2D.h"

void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer);
void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);
dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);


void ellipticPatchSmootherTri2D(solver_t *solver,
                                occa::memory &o_r,
                                occa::memory &o_Sr,
                                const char *options) {

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  occa::memory &o_zP  = solver->o_zP;

  if (strstr(options,"PATCHSOLVE")) {
    occaTimerTic(mesh->device,"PatchSolveKernel");
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
    occaTimerToc(mesh->device,"PatchSolveKernel");
  } else if (strstr(options,"APPROXPATCH")) {
    occaTimerTic(mesh->device,"PatchSolveKernel");
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
    meshParallelGather(mesh, precon->hgsDg, solver->o_zP, o_Sr);
    occaTimerToc(mesh->device,"PatchSolveKernel");
  } else if (strstr(options,"LOCALPATCH")) {
    occaTimerTic(mesh->device,"PatchSolveKernel");
    precon->localPatchSolverKernel(mesh->Nelements,
                              precon->o_invAP,
                              mesh->o_EToE,
                              o_r,
                              o_Sr);
    occaTimerToc(mesh->device,"PatchSolveKernel");
  } else {
    occaTimerTic(mesh->device,"PatchSmoothKernel");
    precon->preconKernel(mesh->Nelements,
       mesh->o_vmapP,
       precon->o_oasForwardDgT,
       precon->o_oasDiagInvOpDg,
       precon->o_oasBackDgT,
       o_r,
       o_zP);
    ellipticParallelGatherScatterTri2D(mesh, precon->ogsDg, o_zP, o_zP, solver->type, "add");
    solver->dotMultiplyKernel(mesh->NpP*mesh->Nelements,precon->o_invDegreeDGP,o_zP,o_zP);
    occaTimerToc(mesh->device,"PatchSmoothKernel");

    // extract block interiors on DEVICE
    occaTimerTic(mesh->device,"restrictKernel");
    precon->restrictKernel(mesh->Nelements, o_zP, o_Sr);
    occaTimerToc(mesh->device,"restrictKernel");
  }
}

void ellipticPreconditioner2D(solver_t *solver,
            dfloat lambda,
            occa::memory &o_r,
            occa::memory &o_z,
            const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  ogs_t    *ogs = solver->ogs; // C0 Gather ScatterTri info

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;


  if(strstr(options, "OAS")){
    //L2 project weighting
    //if(strstr(options,"CONTINUOUS")||strstr(options,"PROJECT")) {
    //  ellipticParallelGatherScatterTri2D(mesh,ogs,o_r,o_r,dfloatString,"add");
    //  mesh->dotMultiplyKernel(mesh->Np*mesh->Nelements,mesh->o_projectL2,o_r,o_r);
    //}

    ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
    ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);

    //patch solve
    ellipticPatchSmootherTri2D(solver,o_r,o_z,options);

    if(strstr(options, "COARSEGRID")){ // should split into two parts
      occaTimerTic(mesh->device,"coarseGrid");

      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      occaTimerTic(mesh->device,"coarsenKernel");
      precon->coarsenKernel(mesh->Nelements, precon->o_V1, o_r, precon->o_r1);
      occaTimerToc(mesh->device,"coarsenKernel");

      // solve coarse problem using xxt
      if(strstr(options, "XXT")){
        precon->o_r1.copyTo(precon->r1);
        occaTimerTic(mesh->device,"xxtSolve");
        xxtSolve(precon->z1, precon->xxt,precon->r1);
        occaTimerToc(mesh->device,"xxtSolve");
        precon->o_z1.copyFrom(precon->z1);
      }

      if(strstr(options,"ALMOND")){
        occaTimerTic(mesh->device,"ALMOND");
        parAlmondPrecon(precon->o_z1, precon->parAlmond, precon->o_r1);
        occaTimerToc(mesh->device,"ALMOND");
      }

      // prolongate from P1 to PN
      occaTimerTic(mesh->device,"prolongateKernel");
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, solver->o_res);
      occaTimerToc(mesh->device,"prolongateKernel");

      // add patch and coarse solves together
      dfloat one = 1.;
      ellipticScaledAdd(solver, one, solver->o_res, one, o_z);
      occaTimerToc(mesh->device,"coarseGrid");
    }

    //project weighting
    //if(strstr(options,"CONTINUOUS")||strstr(options,"PROJECT")) {
    //  mesh->dotMultiplyKernel(mesh->Np*mesh->Nelements,mesh->o_projectL2,o_z,o_z);
    //  ellipticParallelGatherScatterTri2D(mesh, ogs, o_z, o_z, dfloatString, "add");
    //}
#if 0
  } else if(){
    ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
    ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);

    //smooth the fine problem z = Sr
    ellipticPatchSmootherTri2D(solver,o_r,o_z,options);

    dfloat one = 1.; dfloat mone = -1.;
    if(strstr(options, "COARSEGRID")){
      occaTimerTic(mesh->device,"coarseGrid");
      //res = r-Az
      ellipticOperator2D(solver, lambda, o_z, solver->o_res, options);
      ellipticScaledAdd(solver, one, o_r, mone, solver->o_res);

      // restrict to coarse problem
      occaTimerTic(mesh->device,"coarsenKernel");
      precon->coarsenKernel(mesh->Nelements, precon->o_V1, solver->o_res, precon->o_r1);
      occaTimerToc(mesh->device,"coarsenKernel");

      occaTimerTic(mesh->device,"ALMOND");
      parAlmondPrecon(precon->o_z1, precon->parAlmond, precon->o_r1);
      occaTimerToc(mesh->device,"ALMOND");

      // prolongate back to fine problem
      occaTimerTic(mesh->device,"prolongateKernel");
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, solver->o_res);
      occaTimerToc(mesh->device,"prolongateKernel");
      ellipticScaledAdd(solver, one, solver->o_res, one, o_z);
    }

    //do another fine smoothing
    //res = r-Az
    ellipticOperator2D(solver, lambda, o_z, solver->o_res, options);
    ellipticScaledAdd(solver, one, o_r, mone, solver->o_res);

    ellipticStartHaloExchange2D(mesh, solver->o_res, sendBuffer, recvBuffer);
    ellipticEndHaloExchange2D(mesh, solver->o_res, recvBuffer);

    //smooth the fine problem z = z + S(r-Az)
    ellipticPatchSmootherTri2D(solver,solver->o_res,solver->o_Sres,options);
    ellipticScaledAdd(solver, one, solver->o_Sres, one, o_z);
#endif
  } else if (strstr(options, "FULLALMOND")||strstr(options, "OMS")) {

    occaTimerTic(mesh->device,"parALMOND");
    parAlmondPrecon(o_z, precon->parAlmond, o_r);
    occaTimerToc(mesh->device,"parALMOND");

  } else if(strstr(options, "BLOCKJACOBI")){

    dfloat invLambda = 1./lambda;

    occaTimerTic(mesh->device,"blockJacobiKernel");
    precon->blockJacobiKernel(mesh->Nelements, invLambda, mesh->o_vgeo, precon->o_invMM, o_r, o_z);
    occaTimerToc(mesh->device,"blockJacobiKernel");

  } else if(strstr(options, "JACOBI")){

    iint Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    solver->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
}
