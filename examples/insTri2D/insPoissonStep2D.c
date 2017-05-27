#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStep2D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer, 
				        char   * options, char * prSolverOptions){

mesh2D *mesh = ins->mesh; 

solver_t *solver = ins->prsolver;

dfloat t = tstep*ins->dt + ins->dt;


//Exctract Halo On Device

	if(mesh->totalHaloPairs>0){
	 
    ins->poissonHaloExtractKernel(mesh->Nelements,
                           mesh->totalHaloPairs,
                           mesh->o_haloElementList,
                           ins->o_U,
                           ins->o_velHaloBuffer);

    // copy extracted halo to HOST 
    ins->o_velHaloBuffer.copyTo(sendBuffer);            
    // start halo exchange
    meshHaloExchangeStart(mesh,
                          mesh->Np*ins->NVfields*sizeof(dfloat), 
                          sendBuffer,
                          recvBuffer);
  	}



   // computes div u^(n+1) volume term
   ins->poissonRhsVolumeKernel(mesh->Nelements,
                                 mesh->o_vgeo,
                                 mesh->o_DrT,
                                 mesh->o_DsT,
                                 mesh->o_MM,
                                 ins->o_U,  
                                 ins->o_rhsPr);


    // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
  // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);

    mesh->o_haloBuffer.copyFrom(recvBuffer); 

    ins->poissonHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  mesh->o_haloElementList,
                                  ins->o_U,
                                  ins->o_velHaloBuffer);
  }


   //computes div u^(n+1) surface term
  ins->poissonRhsSurfaceKernel(mesh->Nelements,
  	                          ins->dt,	
                              ins->g0,
                              mesh->o_sgeo,
                              mesh->o_FMMT,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              mesh->o_EToB,
                              t,
                              mesh->o_x,
                              mesh->o_y,
                              ins->o_U,
                              ins->o_rhsPr);



  // ins->poissonRhsIpdgBCKernel(mesh->Nelements,
  //                               mesh->o_sgeo,
  //                               mesh->o_vgeo,
  //                               mesh->o_DrT,
  //                               mesh->o_DsT,
  //                               mesh->o_FMMT,
  //                               mesh->o_vmapM,
  //                               mesh->o_vmapP,
  //                               mesh->o_EToB,
  //                               t,
  //                               ins->dt,
  //                               ins->tau,
  //                               mesh->o_x,
  //                               mesh->o_y,
  //                               ins->o_Pr,
  //                               ins->o_rhsPr
  //                               );



// SOLVE HELMHOLTZ EQUATION for ins->o_U
  #if 1
  dfloat lambda =1.0;
  iint Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  dfloat *r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  dfloat *x   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  
  // load rhs into r
  dfloat *cf = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat *nrhs = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  for(iint e=0;e<mesh->Nelements;++e){
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
      
      r[id] = -rhs*J;
      x[id] = 0;
      ins->Pr[id] = nrhs[n];
    }
  }
  free(nrhs);
  free(cf);

  occa::memory o_r   = mesh->device.malloc(Nall*sizeof(dfloat), r);
  occa::memory o_x   = mesh->device.malloc(Nall*sizeof(dfloat), x);

  ellipticSolveTri2D(solver,lambda , o_r, o_x, prSolverOptions);

  // copy solution from DEVICE to HOST
  o_x.copyTo(ins->Pr);

  dfloat maxError = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint   id = e*mesh->Np+n;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat exact = cos(M_PI*xn)*cos(M_PI*yn);
      dfloat error = fabs(exact-ins->Pr[id]);
      
      maxError = mymax(maxError, error);
    }
  }
  
  dfloat globalMaxError = maxError;
  //MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  //if(rank==0)
  printf("globalMaxError = %g\n", globalMaxError);
  
  meshPlotVTU2D(mesh, "foo", 0);

#else

ellipticSolveTri2D(solver, 0.0, ins->o_rhsPr, ins->o_PrI,  prSolverOptions);

#endif
   
}
