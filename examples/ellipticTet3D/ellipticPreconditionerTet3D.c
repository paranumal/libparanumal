#include "ellipticTet3D.h"

void ellipticStartHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticEndHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *recvBuffer);
void ellipticParallelGatherScatterTet3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op);
dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
void ellipticOperator3D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);


void ellipticPreconditioner3D(solver_t *solver,
            dfloat lambda,
            occa::memory &o_r,
            occa::memory &o_z,
            const char *options){

  mesh_t *mesh = solver->mesh;
  precon_t *precon = solver->precon;
  ogs_t    *ogs = solver->ogs; // C0 Gather ScatterTri info

  dfloat *sendBuffer = solver->sendBuffer;
  dfloat *recvBuffer = solver->recvBuffer;

  if (strstr(options, "FULLALMOND")||strstr(options, "MULTIGRID")) {

    occaTimerTic(mesh->device,"parALMOND");
    parAlmondPrecon(precon->parAlmond, o_z, o_r);
    occaTimerToc(mesh->device,"parALMOND");


  } else if(strstr(options, "MASSMATRIX")){

    dfloat invLambda = 1./lambda;

    occaTimerTic(mesh->device,"blockJacobiKernel");
    precon->blockJacobiKernel(mesh->Nelements, invLambda, mesh->o_vgeo, precon->o_invMM, o_r, o_z);
    occaTimerToc(mesh->device,"blockJacobiKernel");

  } else if(strstr(options, "JACOBI")){

    iint Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    occaTimerTic(mesh->device,"dotDivideKernel");
    solver->dotMultiplyKernel(Ntotal, o_r, precon->o_invDiagA, o_z);
    occaTimerToc(mesh->device,"dotDivideKernel");
  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
}
