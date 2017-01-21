#include "ellipticHex3D.h"

void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type);

void ellipticOperator(mesh3D *mesh, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq){
  
  mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq); // store A*q in rhsq

  // do parallel gather scatter
  meshParallelGatherScatter3D(mesh, o_Aq, o_Aq, "float");
}

dfloat ellipticScaledAdd(mesh3D *mesh, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){
  
  mesh->ellipticScaleAddKernel(mesh->Ntotal, alpha, o_a, beta, o_b);
  
}

void ellipticWeightedInnerProduct(iint Ntotal, int Nblock, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b, occa::memory &o_tmp, dfloat *tmp){

  mesh->weightedInnerProduct(Ntotal, o_w, o_a, o_b, o_tmp);

  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(iint n=0;n<Nblock;++n)
    wab += o_tmp[n];

  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalwab;
}


int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);
  
  // set up mesh stuff
  mesh3D *meshSetupHex3D(char *, iint);
  mesh3D *mesh = meshSetupHex3D(argv[1], N);

  // set up elliptic stuff
  int B = 256; // block size for reduction

  void ellipticSetupHex3D(mesh3D *mesh, iint B);
  ellipticSetupHex3D(mesh, B);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint Nblock = (Ntotal+B-1)/B;
  
  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    mesh->rhsq[n] = 1;

  mesh->o_rhsq.copyFrom(mesh->rhsq);
  
  meshParallelGatherScatter3D(mesh, mesh->o_rhsq, mesh->o_rhsq, "float");
  
  mesh->o_rhsq.copyTo(mesh->rhsq);

  dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  for(iint n=0;n<Ntotals;++n)
    invDegree[n] = 1./mesh->rhsq[n];

  occa::memory o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);

  dfloat *p   = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *r   = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *Ap  = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  
  occa::memory o_p   = mesh->device.malloc(Ntotal*sizeof(dfloat), p);
  occa::memory o_r   = mesh->device.malloc(Ntotal*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Ntotal*sizeof(dfloat), x);
  occa::memory o_Ap  = mesh->device.malloc(Ntotal*sizeof(dfloat), Ap);
  occa::memory o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), tmp);
  
  // at this point gather-scatter is available
  dfloat lambda = 10;

  dfloat normr = 0;
  dfloat rdotr0, rdotr1;

  const dfloat tol = 1e-4;
  
  do{

    // placeholder conjugate gradient:
    // https://en.wikipedia.org/wiki/Conjugate_gradient_method
    
    // need to look at the local storage version of CG
    ellipticOperator(mesh, lambda, o_p, o_Ap); // eventually add reduction in scatterKernel

    // need local reduction operation
    iint Ntotal = mesh->Np*mesh->Nelements;
    
    dfloat pAp = ellipticWeightedInnerProduct(Ntotal, o_invDegree, o_p, o_Ap);

    dfloat alpha = rdotr0/pAp;

    // should merge --- -
    ellipticScaledAdd(mesh,  alpha, o_p,  1.f, o_x);
    ellipticScaledAdd(mesh, -alpha, o_Ap, 1.f, o_r);

    rdotr1 = ellipticWeightedNorm(Ntotal, o_invDegree, o_r);
    // ----

    dfloat beta = rdotr1/rdotr0;

    ellipticScaledAdd(mesh, 1.f, o_r, beta, o_p);
    
    
  }while(normr>tol);
  
  mesh->o_rhsq.copyTo(mesh->rhsq);

  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
