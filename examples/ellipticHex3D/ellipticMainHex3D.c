#include "ellipticHex3D.h"

const int B = 256; // block size for reduction (hard coded)

void ellipticStartHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer){

  // count size of halo for this process
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  iint haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    
    // WARNING: uses dfloats
    mesh->haloExtractKernel(mesh->totalHaloPairs, mesh->Np, mesh->o_haloElementList,
			    o_q, mesh->o_haloBuffer);
    
    // copy extracted halo to HOST 
    mesh->o_haloBuffer.copyTo(sendBuffer);
    
    // start halo exchange HOST<>HOST
    meshHaloExchangeStart(mesh,
			  mesh->Np*sizeof(dfloat),
			  sendBuffer,
			  recvBuffer);
  }
}
    

void ellipticEndHaloExchange3D(mesh3D *mesh, occa::memory &o_q, dfloat *recvBuffer){
  
  // count size of halo for this process
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);
  iint haloOffset = mesh->Nelements*mesh->Np*sizeof(dfloat);
  
  // extract halo on DEVICE
  if(haloBytes){
    // finalize recv on HOST
    meshHaloExchangeFinish(mesh);
    
    // copy into halo zone of o_r  HOST>DEVICE
    o_q.copyFrom(recvBuffer, haloBytes, haloOffset);
  }
}

 

void ellipticParallelGatherScatter3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  mesh->device.finish();
  occa::tic("meshParallelGatherScatter3D");
  
  // use gather map for gather and scatter
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);

  mesh->device.finish();
  occa::toc("meshParallelGatherScatter3D");
  
}

void ellipticComputeDegreeVector(mesh3D *mesh, iint Ntotal, ogs_t *ogs, dfloat *deg){

  // build degree vector
  for(iint n=0;n<Ntotal;++n)
    deg[n] = 1;

  occa::memory o_deg = mesh->device.malloc(Ntotal*sizeof(dfloat), deg);
  
  o_deg.copyFrom(deg);
  
  ellipticParallelGatherScatter3D(mesh, ogs, o_deg, o_deg, dfloatString, "add");
  
  o_deg.copyTo(deg);

  mesh->device.finish();
  o_deg.free();
  
}

void ellipticOperator3D(mesh3D *mesh, dfloat *sendBuffer, dfloat *recvBuffer,
			ogs_t *ogs, dfloat lambda,
			occa::memory &o_q, occa::memory &o_gradq, occa::memory &o_Aq, const char *options){


  mesh->device.finish();
  occa::tic("AxKernel");
  
  // compute local element operations and store result in o_Aq
  if(strstr(options, "CONTINUOUS")){
    mesh->AxKernel(mesh->Nelements, mesh->o_ggeo, mesh->o_D, lambda, o_q, o_Aq);

    // parallel gather scatter
    ellipticParallelGatherScatter3D(mesh, ogs, o_Aq, o_Aq, dfloatString, "add");
  }
  else{

    ellipticStartHaloExchange3D(mesh, o_q, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_q, recvBuffer);

    iint allNelements = mesh->Nelements+mesh->totalHaloPairs;
    mesh->gradientKernel(allNelements, mesh->o_vgeo, mesh->o_D, o_q, o_gradq);

    dfloat tau = 2.f*mesh->Nq*mesh->Nq;
    mesh->ipdgKernel(mesh->Nelements,
		     mesh->o_vmapM,
		     mesh->o_vmapP,
		     lambda,
		     tau,
		     mesh->o_vgeo,
		     mesh->o_sgeo,
		     mesh->o_D,
		     o_gradq,
		     o_Aq);
        
  }

  mesh->device.finish();
  occa::toc("AxKernel");
}

dfloat ellipticScaledAdd(mesh3D *mesh, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b){

  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->device.finish();
  occa::tic("scaledAddKernel");
  
  // b[n] = alpha*a[n] + beta*b[n] n\in [0,Ntotal)
  mesh->scaledAddKernel(Ntotal, alpha, o_a, beta, o_b);

  mesh->device.finish();
  occa::toc("scaledAddKernel");
  
}

dfloat ellipticWeightedInnerProduct(mesh3D *mesh,
				    iint Nblock,
				    occa::memory &o_w,
				    occa::memory &o_a,
				    occa::memory &o_b,
				    occa::memory &o_tmp,
				    dfloat *tmp,
				    const char *options){

  mesh->device.finish();
  occa::tic("weighted inner product2");

  iint Ntotal = mesh->Nelements*mesh->Np;
  if(strstr(options,"CONTINUOUS"))
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

void ellipticProject(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_Pv){


  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->dotMultiplyKernel(Ntotal, mesh->o_projectL2, o_v, o_Pv);
  
  ellipticParallelGatherScatter3D(mesh, ogs, o_Pv, o_Pv, dfloatString, "add");
}

void ellipticPreconditioner3D(mesh3D *mesh,
			      precon_t *precon,
			      dfloat *sendBuffer,
			      dfloat *recvBuffer,
			      occa::memory &o_r,
			      occa::memory &o_zP,
			      occa::memory &o_z,
			      const char *type,
			      const char *options){

  if(strstr(options, "OAS")){

    ellipticStartHaloExchange3D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange3D(mesh, o_r, recvBuffer);
    
    mesh->device.finish();
    occa::tic("preconKernel");

    if(strstr(options, "CONTINUOUS")) {
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
      
      // gather-scatter precon blocks on DEVICE [ sum up all contributions ]
      // something goes wrong here:
      ellipticParallelGatherScatter3D(mesh, precon->ogsP, o_zP, o_zP, type, "add");
    }
    else{
      precon->preconKernel(mesh->Nelements,
			   mesh->o_vmapP,
			   precon->o_faceNodesP,
			   precon->o_oasForwardDg,
			   precon->o_oasDiagInvOpDg,
			   precon->o_oasBackDg,
			   o_r,
			   o_zP);

      iint NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;
      //      diagnostic(NtotalP, o_zP, "o_zP before GS"); // ok to here
      
      // OAS => additive Schwarz => sum up patch solutions
      ellipticParallelGatherScatter3D(mesh, precon->ogsDg, o_zP, o_zP, type, "add");

      //      diagnostic(NtotalP, o_zP, "o_zP after GS");

    }
    occa::toc("preconKernel");
	  
    // extract block interiors on DEVICE
    mesh->device.finish();
    occa::tic("restrictKernel");
    
    precon->restrictKernel(mesh->Nelements, o_zP, o_z);
    
    
    mesh->device.finish();
    occa::toc("restrictKernel");   
  }
  else if(strstr(options, "JACOBI")){

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

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be CONTINUOUS or IPDG
  char *options = strdup("solver=PCG preconditioner=OAS method=IPDG");
  //char *options = strdup("solver=PCG preconditioner=OAS method=CONTINUOUS"); 
  
  // set up mesh stuff
  mesh3D *meshSetupHex3D(char *, iint);
  mesh3D *mesh = meshSetupHex3D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // set up elliptic stuff

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 10;
  
  // set up
  ellipticSetupHex3D(mesh, &ogs, &precon, lambda);
  
  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NqP*mesh->NqP*mesh->NqP*mesh->Nelements;
  iint Nblock = (Ntotal+B-1)/B;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;
 
  dfloat *invDegree = (dfloat*) calloc(Ntotal, sizeof(dfloat));
  dfloat *degree = (dfloat*) calloc(Ntotal, sizeof(dfloat));

  occa::memory o_invDegree = mesh->device.malloc(Ntotal*sizeof(dfloat), invDegree);
  
  ellipticComputeDegreeVector(mesh, Ntotal, ogs, degree);

#if 0
  // superb it counts
  for(iint n=0;n<Ntotal;++n){
    if(degree[n]==0) printf("degree[%d]=%d\n", n, degree[n]);
    invDegree[n] = 1./degree[n];
  }
  o_invDegree.copyFrom(invDegree);
#else
  if(strstr(options, "CONTINUOUS")){
    for(iint n=0;n<Ntotal;++n){ // need to weight inner products{
      if(degree[n] == 0) printf("WARNING!!!!\n");
      invDegree[n] = 1./degree[n];
    }
  }
  else
    for(iint n=0;n<Ntotal;++n)
      invDegree[n] = 1.;
  
  o_invDegree.copyFrom(invDegree);
#endif
  
  dfloat *p   = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *r   = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *z   = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *zP  = (dfloat*) calloc(NallP, sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *Ap  = (dfloat*) calloc(Nall, sizeof(dfloat));
  dfloat *tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  dfloat *x4   = (dfloat*) calloc(4*(Ntotal+Nhalo), sizeof(dfloat));
  
  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  // load rhs into r
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){

      iint ggid = e*mesh->Np*mesh->Nggeo + n;
      dfloat wJ = mesh->ggeo[ggid+mesh->Np*GWJID];

      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat zn = mesh->z[id];

      dfloat f = -(3*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);

      r[id] = -wJ*f;

      x[id] = 0; // initial guess
    }
  }

  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  
  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  // need to rename o_r, o_x to avoid confusion
  occa::memory o_p   = mesh->device.malloc(Nall*sizeof(dfloat), p);
  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_z   = mesh->device.malloc(Nall*sizeof(dfloat), z);
  occa::memory o_zP  = mesh->device.malloc(NallP*sizeof(dfloat), zP); // CAUTION
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);
  occa::memory o_Ax  = mesh->device.malloc(Nall*sizeof(dfloat), x);
  occa::memory o_Ap  = mesh->device.malloc(Nall*sizeof(dfloat), Ap);
  occa::memory o_tmp = mesh->device.malloc(Nblock*sizeof(dfloat), tmp);

  occa::memory o_gradx  = mesh->device.malloc(Nall*4*sizeof(dfloat), x4);
  occa::memory o_gradp  = mesh->device.malloc(Nall*4*sizeof(dfloat), x4);

  // use this for OAS precon pairwise halo exchange
  iint haloBytes = mesh->totalHaloPairs*mesh->Np*sizeof(dfloat);

  dfloat *sendBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));
  dfloat *recvBuffer = (dfloat*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(dfloat));

  // copy initial guess for x to DEVICE
  o_x.copyFrom(x);

  // compute A*x
  ellipticOperator3D(mesh, sendBuffer, recvBuffer, ogs, lambda, o_x, o_gradx, o_Ax, options);
  
  // copy r = b
  o_r.copyFrom(r);

  // subtract r = b - A*x
  ellipticScaledAdd(mesh, -1.f, o_Ax, 1.f, o_r);

  // gather-scatter
  if(strstr(options,"CONTINUOUS"))
    ellipticParallelGatherScatter3D(mesh, ogs, o_r, o_r, dfloatString, "add");
  
  if(strstr(options,"PCG")){
    // Precon^{-1} (b-A*x)
    ellipticPreconditioner3D(mesh, precon, sendBuffer, recvBuffer,
			     o_r, o_zP, o_z, dfloatString, options); // r => rP => zP => z

    // p = z
    o_p.copyFrom(o_z); // PCG
  }
  else{
    // p = r
    o_p.copyFrom(o_r); // CG
  }

  // dot(r,r)
  dfloat rdotr0 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_r, o_tmp, tmp, options);
  dfloat rdotz0 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_z, o_tmp, tmp, options);
  
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  iint Niter = 0;
  dfloat alpha, beta;
    
  do{
    
    // A*p
    ellipticOperator3D(mesh, sendBuffer, recvBuffer, ogs, lambda, o_p, o_gradp, o_Ap, options); 

    // dot(p,A*p)
    dfloat pAp =
      ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_p, o_Ap, o_tmp, tmp, options);

    if(strstr(options,"PCG"))
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
    rdotr1 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_r, o_tmp, tmp, options);
    
    if(rdotr1 < tol*tol) break;

    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      ellipticPreconditioner3D(mesh, precon, sendBuffer, recvBuffer,
			       o_r, o_zP, o_z, dfloatString, options);
      
      // dot(r,z)
      rdotz1 = ellipticWeightedInnerProduct(mesh, Nblock, o_invDegree, o_r, o_z, o_tmp, tmp, options);
      
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
      dfloat zn = mesh->z[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
      dfloat error = fabs(exact-mesh->q[id]);
      
      maxError = mymax(maxError, error);

      mesh->q[id] -= exact;
    }
  }

  dfloat globalMaxError = 0;
  MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  if(rank==0)
    printf("globalMaxError = %g\n", globalMaxError);
  
  meshPlotVTU3D(mesh, "foo", 0);
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
