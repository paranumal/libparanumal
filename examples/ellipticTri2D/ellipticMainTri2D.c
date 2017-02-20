#include "ellipticTri2D.h"

const int B = 256; // block size for reduction (hard coded)

void ellipticStartHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *sendBuffer, dfloat *recvBuffer){

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
    

void ellipticEndHaloExchange2D(mesh2D *mesh, occa::memory &o_q, dfloat *recvBuffer){
  
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


void ellipticParallelGatherScatter2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  mesh->device.finish();
  occa::tic("meshParallelGatherScatter2D");
  
  // use gather map for gather and scatter
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);

  mesh->device.finish();
  occa::toc("meshParallelGatherScatter2D");
  
}

void ellipticOperator2D(mesh2D *mesh, dfloat *sendBuffer, dfloat *recvBuffer,
			ogs_t *ogs, dfloat lambda,
			occa::memory &o_q, occa::memory &o_gradq, occa::memory &o_Aq, const char *options){

  mesh->device.finish();
  occa::tic("AxKernel");

  // IPDG ONLY 
  ellipticStartHaloExchange2D(mesh, o_q, sendBuffer, recvBuffer);
  
  ellipticEndHaloExchange2D(mesh, o_q, recvBuffer);
  
  iint allNelements = mesh->Nelements+mesh->totalHaloPairs; 
  mesh->gradientKernel(allNelements,
		       mesh->o_vgeo,
		       mesh->o_DrT,
		       mesh->o_DsT,
		       o_q,
		       o_gradq);
  
  // TW NOTE WAS 2 !
  dfloat tau = 2.f*(mesh->N+1)*(mesh->N+1); // 1/h factor built into kernel 
  mesh->ipdgKernel(mesh->Nelements,
		   mesh->o_vmapM,
		   mesh->o_vmapP,
		   lambda,
		   tau,
		   mesh->o_vgeo,
		   mesh->o_sgeo,
		   mesh->o_DrT,
		   mesh->o_DsT,
		   mesh->o_LIFTT,
		   mesh->o_MM,
		   o_gradq,
		   o_Aq);

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

dfloat ellipticInnerProduct(mesh2D *mesh,
			    iint Nblock,
			    occa::memory &o_a,
			    occa::memory &o_b,
			    occa::memory &o_tmp,
			    dfloat *tmp,
			    const char *options){

  mesh->device.finish();
  occa::tic("weighted inner product2");

  iint Ntotal = mesh->Nelements*mesh->Np;

  mesh->innerProductKernel(Ntotal, o_a, o_b, o_tmp);
  
  mesh->device.finish();
  occa::toc("weighted inner product2");
  
  o_tmp.copyTo(tmp);

  dfloat ab = 0;
  for(iint n=0;n<Nblock;++n){
    ab += tmp[n];
  }
      
  dfloat globalab = 0;
  MPI_Allreduce(&ab, &globalab, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  return globalab;
}

void ellipticPreconditioner2D(mesh2D *mesh,
			      ogs_t *ogs,
			      precon_t *precon,
			      dfloat *sendBuffer,
			      dfloat *recvBuffer,
			      occa::memory &o_r,
			      occa::memory &o_zP,
			      occa::memory &o_z,
			      const char *type,
			      const char *options){

  if(strstr(options, "OAS")){
    
    ellipticStartHaloExchange2D(mesh, o_r, sendBuffer, recvBuffer);
    
    ellipticEndHaloExchange2D(mesh, o_r, recvBuffer);

    mesh->device.finish();
    occa::tic("preconKernel");

    precon->preconKernel(mesh->Nelements,
			 mesh->o_vmapP,
			 precon->o_oasForwardDgT,
			 precon->o_oasDiagInvOpDg,
			 precon->o_oasBackDgT,
			 o_r,
			 o_zP);
    
    iint NtotalP = mesh->NpP*mesh->Nelements;
    diagnostic(NtotalP, o_zP, "o_zP before GS"); // ok to here
    
    // OAS => additive Schwarz => sum up patch solutions
    ellipticParallelGatherScatter2D(mesh, precon->ogsDg, o_zP, o_zP, type, "add");
    
    diagnostic(NtotalP, o_zP, "o_zP after GS");

    mesh->device.finish();
    occa::toc("preconKernel");
    
    // extract block interiors on DEVICE
    mesh->device.finish();
    occa::tic("restrictKernel");

    precon->restrictKernel(mesh->Nelements, o_zP, o_z);

    if(strstr(options, "COARSEGRID")){ // should split into two parts

      // Z1*Z1'*PL1*(Z1*z1) = (Z1*rL)  HMMM
      precon->coarsenKernel(mesh->Nelements, precon->o_coarseInvDegree, precon->o_V1, o_r, precon->o_r1);

      // do we need to gather (or similar) here ?
      precon->o_r1.copyTo(precon->r1); 
      
      xxtSolve(precon->z1, precon->xxt,precon->r1);

      precon->o_z1.copyFrom(precon->z1);
      
      precon->prolongateKernel(mesh->Nelements, precon->o_V1, precon->o_z1, precon->o_ztmp);

      dfloat one = 1.;
      ellipticScaledAdd(mesh, one, precon->o_ztmp, one, o_z);
    }

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
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }
  
}



int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityTriH02.msh N\n");
    exit(-1);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // solver can be CG or PCG
  // preconditioner can be JACOBI, OAS, NONE
  // method can be IPDG
  //char *options = strdup("solver=PCG preconditioner=OAS method=IPDG");
  char *options = strdup("solver=PCG preconditioner=OAS method=IPDG coarse=COARSEGRID");
  //char *options = strdup("solver=PCG preconditioner=NONE method=IPDG");
  
  // set up mesh stuff
  mesh2D *meshSetupTri2D(char *, iint);
  mesh2D *mesh = meshSetupTri2D(argv[1], N);
  ogs_t *ogs;
  precon_t *precon;
  
  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda = 1;
  
  // set up
  ellipticSetupTri2D(mesh, &ogs, &precon, lambda);

  iint Ntotal = mesh->Np*mesh->Nelements;
  iint NtotalP = mesh->NpP*mesh->Nelements;
  iint Nblock = (Ntotal+B-1)/B;
  iint Nhalo = mesh->Np*mesh->totalHaloPairs;
  iint Nall   = Ntotal + Nhalo;
  iint NallP  = NtotalP;

  dfloat *p   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *z   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *zP  = (dfloat*) calloc(NallP,  sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *Ap  = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));

  dfloat *x4   = (dfloat*) calloc(4*(Ntotal+Nhalo), sizeof(dfloat));
  
  // convergence tolerance (currently absolute)
  const dfloat tol = 1e-6;

  // load rhs into r
  dfloat *cf = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){

#if 1
    for(iint n=0;n<mesh->cubNp;++n){
      dfloat cx = 0, cy = 0;
      for(iint m=0;m<mesh->Np;++m){
	cx += mesh->cubInterp[m+n*mesh->Np]*mesh->x[m+e*mesh->Np];
	cy += mesh->cubInterp[m+n*mesh->Np]*mesh->y[m+e*mesh->Np];
      }
      dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
      dfloat w = mesh->cubw[n];
      
      cf[n] = -J*w*(2*M_PI*M_PI+lambda)*cos(M_PI*cx)*cos(M_PI*cy);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->cubNp;++m){
	rhs += mesh->cubInterp[n+m*mesh->Np]*cf[m];
      }
      iint id = n+e*mesh->Np;
      r[id] = -rhs;
      x[id] = 0; // initial guess
    }
#else
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(iint n=0;n<mesh->Np;++n){
      dfloat xn = mesh->x[n+e*mesh->Np];
      dfloat yn = mesh->y[n+e*mesh->Np];
      nrhs[n] = -(2*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn);
    }
    for(iint n=0;n<mesh->Np;++n){
      dfloat rhs = 0;
      for(iint m=0;m<mesh->Np;++m){
	rhs += mesh->MM[n+m*mesh->Np]*nrhs[m];
      }
      iint id = n+e*mesh->Np;
      
      r[id] = rhs*J;
      x[id] = 0;
      mesh->q[id] = nrhs[n];
    }
#endif
  }
  free(nrhs);
  free(cf);

  meshPlotVTU2D(mesh, "foo1", 0);
  
  // placeholder conjugate gradient:
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method
  
  // placeholder preconditioned conjugate gradient
  // https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method

  // need to rename o_r, o_x to avoid confusion
  occa::memory o_p   = mesh->device.malloc(Nall*sizeof(dfloat), p);
  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_z   = mesh->device.malloc(Nall*sizeof(dfloat), z);
  occa::memory o_zP  = mesh->device.malloc(NallP*sizeof(dfloat),zP); // CAUTION
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
  ellipticOperator2D(mesh, sendBuffer, recvBuffer, ogs, lambda, o_x, o_gradx, o_Ax, options);
  
  // copy r = b
  o_r.copyFrom(r);

  // subtract r = b - A*x
  ellipticScaledAdd(mesh, -1.f, o_Ax, 1.f, o_r);

  if(strstr(options,"PCG")){

    // Precon^{-1} (b-A*x)
    ellipticPreconditioner2D(mesh, ogs, precon, sendBuffer, recvBuffer, 
			     o_r, o_zP, o_z, dfloatString, options); // r => rP => zP => z
    
    // p = z
    o_p.copyFrom(o_z); // PCG
  }
  else{
    // p = r
    o_p.copyFrom(o_r); // CG
  }


  // dot(r,r)
  dfloat rdotr0 = ellipticInnerProduct(mesh, Nblock, o_r, o_r, o_tmp, tmp, options);
  dfloat rdotz0 = ellipticInnerProduct(mesh, Nblock, o_r, o_z, o_tmp, tmp, options);
  dfloat rdotr1 = 0;
  dfloat rdotz1 = 0;
  iint Niter = 0;
  dfloat alpha, beta;

  if(rank==0)
    printf("rdotr0 = %g, rdotz0 = %g\n", rdotr0, rdotz0);
  
  do{

    diagnostic(Ntotal, o_p, "o_p");
    
    // A*p
    ellipticOperator2D(mesh, sendBuffer, recvBuffer, ogs, lambda, o_p, o_gradp, o_Ap, options); 

    diagnostic(Ntotal, o_p, "o_Ap");
    
    // dot(p,A*p)
    dfloat pAp =
      ellipticInnerProduct(mesh, Nblock, o_p, o_Ap, o_tmp, tmp, options);
    
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
    rdotr1 = ellipticInnerProduct(mesh, Nblock, o_r, o_r, o_tmp, tmp, options);
    
    if(rdotr1 < tol*tol) break;

    if(strstr(options,"PCG")){

      // z = Precon^{-1} r
      ellipticPreconditioner2D(mesh, ogs, precon, sendBuffer, recvBuffer, 
			       o_r, o_zP, o_z, dfloatString, options);
      
      // dot(r,z)
      rdotz1 = ellipticInnerProduct(mesh, Nblock, o_r, o_z, o_tmp, tmp, options);
      
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
