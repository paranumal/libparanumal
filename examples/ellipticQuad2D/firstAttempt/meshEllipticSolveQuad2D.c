#include "mpi.h"
#include "mesh2D.h"

void meshEllipticSolveQuad2D(mesh2D *mesh){

  // see p213: Mund, Deville, Fischer
  dfloat lambda = 1.;
  dfloat TOL = 1e-8;
  
  int NnodesL = mesh->Np*mesh->Nelements;

  // storage
  dfloat *W   = (dfloat*) calloc(NnodesL, sizeof(dfloat)); // reciprocal multiplicity of local nodes
  dfloat *ubL = (dfloat*) calloc(NnodesL, sizeof(dfloat)); // 
  dfloat *rL  = (dfloat*) calloc(NnodesL, sizeof(dfloat)); // 
  dfloat *fL  = (dfloat*) calloc(NnodesL, sizeof(dfloat));
  dfloat *zL  = (dfloat*) calloc(NnodesL, sizeof(dfloat));
  dfloat *xL  = (dfloat*) calloc(NnodesL, sizeof(dfloat));
  dfloat *qL  = (dfloat*) calloc(NnodesL, sizeof(dfloat));
  dfloat *wL  = (dfloat*) calloc(NnodesL, sizeof(dfloat));
  dfloat *ML  = (dfloat*) calloc(NnodesL, sizeof(dfloat));

  // mark gather nodes as boundary (or not) WILL break in parallel
  int *boundaryFlag = (int*) calloc(NnodesL, sizeof(int));
  for(int n=0;n<mesh->Nfp*mesh->Nfaces*mesh->Nelements;++n){
    if(mesh->vmapM[n] == mesh->vmapP[n]){
      int id = mesh->globalNumbering[mesh->vmapM[n]];
      boundaryFlag[id] = 1;
    }
  }

  // mark all boundary nodes and count WILL break in parallel
  int NnodesB = 0;
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    int id = mesh->globalNumbering[n];
    if(boundaryFlag[id]==1){
      boundaryFlag[n] = 1;
      ++NnodesB;
    }
  }
  
  // collect all boundary nodes
  int *nodesB = (int*) calloc(NnodesB, sizeof(int));
  int cnt = 0;
  for(int n=0;n<mesh->Np*mesh->Nelements;++n)
    if(boundaryFlag[n]==1)
      nodesB[cnt++] = n;

  // store reciprocal multiplicity of each node
  cnt = 0;
  for(int n=0;n<mesh->NgatherNodes;++n){  
    for(int m=0;m<mesh->gatherCounts[n];++m){
      int id = mesh->gatherIds[cnt];
      ML[id] = (1-boundaryFlag[id]);    // mask
      W[id] = 1./mesh->gatherCounts[n]; // reciprocal multiplicity 
      ++cnt;
    }
  }
  
  // rL = ML*GS*(MML*fL - HL*ubL)
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id = n + e*mesh->Np;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      int  gid = mesh->Np*mesh->Nvgeo*e + n + mesh->Np*JWID;
      dfloat Jw = mesh->vgeo[gid];
      // f = -(lambda+2*pi*pi)*cos(pi*x)*cos(pi*y)
      fL[id] = Jw*(lambda+2*M_PI*M_PI)*cos(M_PI*xn)*cos(M_PI*yn); // note this is massMatrixL
    }
  }

  for(int n=0;n<NnodesB;++n){
    int id = nodesB[n];
    dfloat xn = mesh->x[id];
    dfloat yn = mesh->y[id];
    ubL[id] = -cos(M_PI*xn)*cos(M_PI*yn); // -bc    
  }
  
  meshEllipticLocalOpQuad2D(mesh, ubL, lambda, rL);

  for(int n=0;n<NnodesL;++n)
    rL[n] += fL[n];
  
  meshGatherScatterQuad2D(mesh, rL);

  dfloat rho1 = 0;
  for(int n=0;n<NnodesL;++n){
    // mask boundary nodes
    rL[n] *= ML[n];
  
    // zL = PL\rL (precondition locally on elements + halo)
    zL[n] = rL[n];
  
    // wL = zL
    wL[n] = zL[n];
  
    // rho1 = rL.(W*zL)    
    rho1 += rL[n]*W[n]*zL[n];
  }
  
  // for j=0,1... until convergence do
  dfloat alpha, beta, rho0, wdotq;
  do{

    // qL = ML*(GS*(HL*wL)) (local op then gather-scatter then mask)
    {
      meshEllipticLocalOpQuad2D(mesh, wL, lambda, qL);    
      
      meshGatherScatterQuad2D(mesh, qL);
      
      for(int n=0;n<NnodesL;++n)
	qL[n] *= ML[n];
    }
    
    // alpha = rho1/(wL.*(W*qL)) (divide by multiplicity, compute inner product)
    wdotq = 0;
    for(int n=0;n<NnodesL;++n)
      wdotq += wL[n]*W[n]*qL[n];
    alpha = rho1/wdotq;

    dfloat rdotr = 0;
    for(int n=0;n<NnodesL;++n){
      // xL = xL + alpha*wL
      xL[n] += alpha*wL[n];
      
      // rL = rL - alpha*qL
      rL[n] -= alpha*qL[n];

      // 
      rdotr += rL[n]*rL[n]*W[n];
    }

    if(sqrt(rdotr)<TOL) break;
    
    // zL = PL\rL
    for(int n=0;n<NnodesL;++n)
      zL[n] = rL[n];
    
    // rho0 = rho1
    rho0 = rho1;

    // rho1 = rL.(W*zL) (divide by multiplicity, compute inner-product)
    rho1 = 0;
    for(int n=0;n<NnodesL;++n)
      rho1 += rL[n]*W[n]*zL[n];

    // beta = rho1/rho0
    beta = rho1/rho0;
    
    // wL = zL + beta*wL
    for(int n=0;n<NnodesL;++n)
      wL[n] = zL[n] + beta*wL[n];

    printf("alpha = %g, beta = %g, rho0 = %g, rho1 = %g, rdotr = %g\n",
	   alpha, beta, rho0, rho1, rdotr);
  }while(1);

  mesh->Nfields = 2;
  mesh->q = (dfloat*) calloc(mesh->Nfields*mesh->Np*mesh->Nelements, sizeof(dfloat));
  // solutionL = xL + ubL
  dfloat maxerr = 0;
  for(int n=0;n<NnodesL;++n){

    dfloat xn = mesh->x[n];
    dfloat yn = mesh->y[n];
    
    mesh->q[mesh->Nfields*n]   = xL[n] - ubL[n];
    
    dfloat err =  fabs(mesh->q[mesh->Nfields*n] - cos(M_PI*xn)*cos(M_PI*yn));
    mesh->q[mesh->Nfields*n+1] = err;

    maxerr = mymax(maxerr, err);
  }
  printf("maxerr=%g\n", maxerr);
  
  int fld = 0;
  meshPlotVTU2D(mesh, "foo", fld);
  
  free(W); free(ubL); free(rL); free(zL); free(xL); free(qL); free(wL); free(ML); 
  free(boundaryFlag);
  free(nodesB);
}
