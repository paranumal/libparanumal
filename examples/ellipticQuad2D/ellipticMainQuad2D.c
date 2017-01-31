#include "ellipticQuad2D.h"

const int B = 256; // block size for reduction (hard coded)


// SOME WRITE OR READ RACE CONDITION ??? non-deterministic 

void ellipticParallelGatherScatter2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type){

  mesh->device.finish();
  occa::tic("meshParallelGatherScatter2D");
  
  // use gather map for gather and scatter
  meshParallelGatherScatter2D(mesh, ogs, o_q, o_gsq, type);

  mesh->device.finish();
  occa::toc("meshParallelGatherScatter2D");
  
}

void ellipticComputeDegreeVector(mesh2D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatter2D(mesh, ogs, o_deg, o_deg, dfloatString);
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

void ellipticOperator2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda, occa::memory &o_q, occa::memory &o_gradq, occa::memory &o_Aq, const char *method){

  mesh->device.finish();
  occa::tic("AxKernel");

  if(!strcmp(method, "H0")){
    // compute local element operations and store result in o_Aq
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);
    
    // parallel gather scatter
    ellipticParallelGatherScatter2D(mesh, ogs, o_Aq, o_Aq, dfloatString);
    
  }
  else{

    printf("IPDG\n");

    mesh->gradientKernel(mesh->Nelements, mesh->o_vgeo, mesh->o_D, o_q, o_gradq);

    dfloat tau = 1000.f;
    mesh->ipdgKernel(mesh->Nelements, mesh->o_vmapM, mesh->o_vmapP, lambda, tau,
		     mesh->o_vgeo, mesh->o_sgeo, mesh->o_D, o_gradq, o_Aq);

  }
  mesh->device.finish();
  occa::toc("AxKernel");

}

dfloat ellipticScaledAdd(mesh2D *mesh, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->device.finish();
  occa::tic("scaledAddKernel");
  
  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  mesh->device.finish();
  occa::toc("scaledAddKernel");
  
}

dfloat ellipticWeightedInnerProduct(mesh2D *mesh,
				    iint Nblock,
				    occa::memory &o_w,
				    occa::memory &o_a,
				    occa::memory &o_b,
				    occa::memory &o_tmp,
				    dfloat *tmp,
				    const char *methodType){

  mesh->device.finish();
  occa::tic("weighted inner product2");

  iint Ntotal = mesh->Nelements*mesh->Np;

  if(!strcmp(methodType,"H0"))
    mesh->weightedInnerProduct2Kernel(Ntotal, o_w, o_a, o_b, o_tmp);
  else
    mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  
  mesh->device.finish();
  occa::toc("weighted inner product2");
  
  o_tmp.copyTo(tmp);

  dfloat wab = 0;
  for(iint n=0;n<Nblock;++n){
    wab += tmp[n];
  }
      
  dfloat globalwab = 0;
  MPI_Allreduce(&wab, &globalwab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalwab;
}

#if 0
dfloat ellipticWeightedInnerProduct(mesh2D *mesh,
				    iint Nblock,
				    occa::memory &o_w,
				    occa::memory &o_a,
				    occa::memory &o_tmp,
				    dfloat *tmp){

  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->device.finish();
  occa::tic("weighted inner product1");
  
  mesh->weightedInnerProduct1Kernel(Ntotal, o_w, o_a, o_tmp);

  mesh->device.finish();
  occa::toc("weighted inner product1");
  
  o_tmp.copyTo(tmp);
  
  dfloat wa2 = 0;
  for(iint n=0;n<Nblock;++n){
    wa2 += tmp[n];
  }
  
  dfloat globalwa2 = 0;
  MPI_Allreduce(&wa2, &globalwa2, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalwa2;
}

dfloat diagnostics(mesh2D *mesh, iint Nblock, occa::memory &o_w, occa::memory &o_a,
		   occa::memory &o_tmp, dfloat *tmp, const char *message){

  dfloat wa2 = ellipticWeightedInnerProduct(mesh, Nblock, o_w, o_a, o_tmp, tmp);

  printf("%s: L2 norm = %g\n", message, sqrt(wa2));
}
#endif

void ellipticProject(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_Pv){


  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->dotMultiplyKernel(Ntotal, mesh->o_projectL2, o_v, o_Pv);
  
  ellipticParallelGatherScatter2D(mesh, ogs, o_Pv, o_Pv, dfloatString);
}

void ellipticPreconditioner2D(mesh2D *mesh,
			      precon_t *precon,
			      ogs_t *ogs,
			      dfloat *sendBuffer,
			      dfloat *recvBuffer,
			      occa::memory &o_r,
			      occa::memory &o_zP,
			      occa::memory &o_z,
			      const char *type,
			      char *preconType,
			      dfloat *invDegree,
			      dfloat *invDegreeP){

  if(!strcmp(preconType, "OAS")){

    // the non-repeatability happens in here somewhere
    
    // count size of halo for this process
    iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
    iint haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
    
    // extract halo on DEVICE
    if(haloBytes){
      
      // WARNING: uses dfloats
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      mesh->Np,
			      mesh->o_haloElementList,
			      o_r,
			      mesh->o_haloBuffer);
    
      // copy extracted halo to HOST 
      mesh->o_haloBuffer.copyTo(sendBuffer);
      
      // start halo exchange HOST<>HOST
      meshHaloExchangeStart2D(mesh,
			      mesh->Np*sizeof(dfloat),
			      sendBuffer,
			      recvBuffer);
      
      // finalize recv on HOST
      meshHaloExchangeFinish2D(mesh);
      
      // copy into halo zone of o_r  HOST>DEVICE
      o_r.copyFrom(recvBuffer, haloBytes, haloOffset);
    }
    
    mesh->device.finish();
    occa::tic("preconKernel");
    
    // compute local precon on DEVICE
    precon->preconKernel(mesh->Nelements,
			 precon->o_vmapPP,
			 precon->o_faceNodesP,
			 precon->o_oasForward,
			 precon->o_oasDiagInvOp,
			 precon->o_oasBack,
			 o_r,
			 o_zP);

    mesh->device.finish();
    occa::toc("preconKernel");
    
    // gather-scatter precon blocks on DEVICE [ sum up all contributions ]
    // something goes wrong here:
    ellipticParallelGatherScatter2D(mesh, precon->ogsP, o_zP, o_zP, type);

    // extract block interiors on DEVICE
    mesh->device.finish();
    occa::tic("restrictKernel");
    
    precon->restrictKernel(mesh->Nelements, o_zP, o_z);
    
    
    mesh->device.finish();
    occa::toc("restrictKernel");   
  }
  else if(!strcmp(preconType, "JACOBI")){

    occa::tic("dotDivideKernel");   
    mesh->device.finish();
    
    // Jacobi preconditioner
    iint Ntotal = mesh->Np*mesh->Nelements;
    mesh->dotDivideKernel(Ntotal, o_r, precon->o_diagA, o_z);

    occa::toc("dotDivideKernel");   
    mesh->device.finish();
  }
  else // turn off preconditioner
    o_z.copyFrom(o_r);
  
}



int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // int specify polynomial degree 
  int N = atoi(argv[2]);
  
  // set up mesh stuff
  mesh2D *meshSetupQuad2D(char *, iint);
  mesh2D *mesh = meshSetupQuad2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // set up elliptic stuff

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 10;
  
  // set up
  ellipticSetupQuad2D(mesh, &ogs, &precon, lambda);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;
  iint Nblock = (Ntotal+B-1)/B;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;

  dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *degree = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  occa::memory o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);
  
  ellipticComputeDegreeVector(mesh, Ntotal, ogs, degree);

  for(iint n=0;n<Ntotal;++n){
    if(degree[n]==0) printf("degree[%d]=%d\n", n, degree[n]);
    invDegree[n] = 1./degree[n];
  }
  o_invDegree.copyFrom(invDegree);


  dfloat *invDegreeP = (dfloat*) calloc(NtotalP, sizeof(dfloat));
  dfloat *degreeP = (dfloat*) calloc(NtotalP, sizeof(dfloat));

  occa::memory o_invDegreeP = mesh->device.malloc(NtotalP*sizeof(dfloat), invDegreeP);
  
  ellipticComputeDegreeVector(mesh, NtotalP, precon->ogsP, degreeP);

  for(iint n=0;n<NtotalP;++n){
    if(degreeP[n]==0) printf("degreeP[%d]=%d\n", n, degreeP[n]);
    invDegreeP[n] = 1./degreeP[n];
  }
  o_invDegreeP.copyFrom(invDegreeP);
  
  dfloat *p   = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *r   = (dfloat*) calloc(Ntotal+Nhalo, sizeof(dfloat));
  dfloat *z   = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *zP  = (dfloat*) calloc(NtotalP, sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Ntotal+Nhalo, sizeof(dfloat));
  dfloat *Ap  = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  
  // at this point gather-scatter is available

  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-10;

  // load rhs into r
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      iint ggid = e*mesh->Np*mesh->Nggeo + n;
      dfloat wJ = mesh->ggeo[ggid+mesh->Np*GWJID];

      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];

      dfloat f = -(2*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn);
      
      r[id] = -wJ*f;

      x[id] = 0; // initial guess
    }
  }

  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  
  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  // need to rename o_r, o_x to avoid confusion
  occa::memory o_p   = mesh->device.malloc(Ntotal*sizeof(dfloat), p);
  occa::memory o_r   = mesh->device.malloc((Ntotal+Nhalo)*sizeof(dfloat), r);
  occa::memory o_z   = mesh->device.malloc(Ntotal*sizeof(dfloat), z);
  occa::memory o_zP  = mesh->device.malloc(NtotalP*sizeof(dfloat), zP); // CAUTION
  occa::memory o_x   = mesh->device.malloc((Ntotal+Nhalo)*sizeof(dfloat), x);
  occa::memory o_Ax  = mesh->device.malloc(Ntotal*sizeof(dfloat), x);
  occa::memory o_Ap  = mesh->device.malloc(Ntotal*sizeof(dfloat), Ap);
  occa::memory o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), tmp);

  dfloat *x4   = (dfloat*) calloc(4*(Ntotal+Nhalo), sizeof(dfloat));
  occa::memory o_gradx  = mesh->device.malloc((Ntotal+Nhalo)*4*sizeof(dfloat), x4);
  occa::memory o_gradp  = mesh->device.malloc((Ntotal+Nhalo)*4*sizeof(dfloat), x4);
  

  char *iterationType = strdup("PCG"); // can be CG or PCG
  char *preconType = strdup("NONE");    // can be JACOBI, OAS, NONE
  char *methodType = strdup("IPDG");   // can be IPDG or H0
  
  // use this for OAS precon pairwise halo exchange
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);

  dfloat *sendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));
  dfloat *recvBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));

  // copy initial guess for x to DEVICE
  o_x.copyFrom(x);

  // compute A*x
  ellipticOperator2D(mesh, ogs, lambda, o_x, o_gradx, o_Ax, methodType); 
  
  // copy r = b
  o_r.copyFrom(r);

  // subtract r = b - A*x
  ellipticScaledAdd(mesh, -1.f, o_Ax, 1.f, o_r);

  // gather-scatter 
  if(!strcmp(methodType, "H0"))
    ellipticParallelGatherScatter2D(mesh, ogs, o_r, o_r, dfloatString);
  
  if(!strcmp(iterationType,"PCG")){
    // Precon^{-1} (b-A*x)
    ellipticPreconditioner2D(mesh, precon, ogs, sendBuffer, recvBuffer,
			     o_r, o_zP, o_z, dfloatString, preconType, invDegree, invDegreeP); // r => rP => zP => z

    // p = z
    o_p.copyFrom(o_z); // PCG
  }
  else{
    // p = r
    o_p.copyFrom(o_r); // CG
  }

  //  MPI_Finalize();
  //  exit(0);
  
  // dot(r,r)
  dfloat rdotr0 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_r, o_tmp, tmp, methodType);
  dfloat rdotz0 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_z, o_tmp, tmp, methodType); 
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  iint Niter = 0;
  dfloat alpha, beta;
    
  do{
    
    // A*p
    ellipticOperator2D(mesh, ogs, lambda, o_p, o_gradp, o_Ap, methodType); // eventually add reduction in scatterKernel

    // dot(p,A*p)
    dfloat pAp = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_p, o_Ap, o_tmp, tmp, methodType);

    if(!strcmp(iterationType,"PCG"))
      // alpha = dot(r,z)/dot(p,A*p)
      alpha = rdotz0/pAp;
    else
      // alpha = dot(r,r)/dot(p,A*p)
      alpha = rdotr0/pAp;

    // x <= x + alpha*p
    ellipticScaledAdd(mesh,  alpha, o_p,  1.f, o_x);

    // r <= r - alpha*A*p
    ellipticScaledAdd(mesh, -alpha, o_Ap, 1.f, o_r);

    // dot(r,r)
    rdotr1 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_r, o_tmp, tmp, methodType);
    
    if(rdotr1 < tol*tol) break;

    if(!strcmp(iterationType,"PCG")){

      // z = Prvon^{-1} r
      ellipticPreconditioner2D(mesh, precon, ogs, sendBuffer, recvBuffer, o_r, o_zP, o_z,
			       dfloatString, preconType, invDegree, invDegreeP); 
      
      // dot(r,z)
      rdotz1 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_z, o_tmp, tmp, methodType);

      beta = rdotz1/rdotz0;

      // p = z + beta*p
      ellipticScaledAdd(mesh, 1.f, o_z, beta, o_p);

      // switch rdotz0 <= rdotz1
      rdotz0 = rdotz1;
    }
    else{
      beta = rdotr1/rdotr0;

      // p = r + beta*p
      ellipticScaledAdd(mesh, 1.f, o_r, beta, o_p);
    }

    // switch rdotr0 <= rdotr1
    rdotr0 = rdotr1;
    
    if(rank==0)
      printf("iter=%05d pAp = %g norm(r) = %g\n", Niter, pAp, sqrt(rdotr0));

    ++Niter;
    
  }while(rdotr0>(tol*tol));

  occa::printTimer();
  
  // copy solution from DEVICE to HOST
  o_x.copyTo(mesh->q);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn);
      dfloat error = fabs(exact-mesh->q[id]);
      
      maxError = mymax(maxError, error);

      mesh->q[id] -= exact;
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);
  
  meshPlotVTU2D(mesh, "foo", 0);
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
