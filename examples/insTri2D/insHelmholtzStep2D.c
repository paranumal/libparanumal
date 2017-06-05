#include "ins2D.h"

// complete a time step using LSERK4
void insHelmholtzStep2D(ins_t *ins, iint tstep,  iint haloBytes,
				               dfloat * sendBuffer, dfloat * recvBuffer, 
                       char   * options){

	mesh2D *mesh = ins->mesh; 

  solver_t *solver = ins->velsolver; 

	dfloat t = tstep*ins->dt;


  // Exctract Halo On Device

	// if(mesh->totalHaloPairs>0){
	 
 //    ins->helmholtzHaloExtractKernel(mesh->Nelements,
 //                                    mesh->totalHaloPairs,
 //                                    mesh->o_haloElementList,
 //                                    ins->o_Ux,
 //                                    ins->o_Uy,
 //                                    ins->o_Pr,
 //                                    ins->o_totHaloBuffer);

 //    // copy extracted halo to HOST 
 //    ins->o_totHaloBuffer.copyTo(sendBuffer);            
 //    // start halo exchange
 //    meshHaloExchangeStart(mesh,
 //                          mesh->Np*(ins->NTfields)*sizeof(dfloat), // pressure also 
 //                          sendBuffer,
 //                          recvBuffer);
 //  	}

  	// Compute Volume Contribution
   ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            ins->o_Pr,  
                            ins->o_rhsUx,
                            ins->o_rhsUy);


  //   // COMPLETE HALO EXCHANGE
  // if(mesh->totalHaloPairs>0){
  // // wait for halo data to arrive
  //   meshHaloExchangeFinish(mesh);

  //   mesh->o_haloBuffer.copyFrom(recvBuffer); 

  //   ins->helmholtzHaloScatterKernel(mesh->Nelements,
  //                                   mesh->totalHaloPairs,
  //                                   mesh->o_haloElementList,
  //                                   ins->o_Ux,
  //                                   ins->o_Uy,
  //                                   ins->o_Pr,
  //                                   ins->o_totHaloBuffer);
  // }

 // Compute Surface Conribution
  ins->gradientSurfaceKernel(mesh->Nelements,
                              mesh->o_sgeo,
                              mesh->o_LIFTT,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              mesh->o_EToB,
                              mesh->o_x,
                              mesh->o_y,
                              t,
                              ins->dt, 
                              0,
                              ins->o_Pr,
                              ins->o_rhsUx,
                              ins->o_rhsUy);










 // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,	
                              ins->a0,
                              ins->a1,
                              ins->a2,
                              ins->b0,
                              ins->b1,
                              ins->b2,
                              ins->o_Ux,
                              ins->o_Uy,
                              ins->o_UO,
                              ins->o_NUx,
                              ins->o_NUy,
                              ins->o_NO,
                              ins->o_rhsUx,
                              ins->o_rhsUy);


  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->lamda,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsUx,
                                ins->o_rhsUy);

  
// USE STREAMING LATER!!!!!!!
  // SOLVE HELMHOLTZ EQUATION for ins->o_U
#if 0
 ins->o_Ux.copyTo(ins->Ux);
 ins->o_rhsUx.copyTo(ins->rhsUx);

 dfloat maxU = 0, minU = 1e9;
 dfloat maxR = 0, minR = 1e9;

 for(iint e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        iint id = n+e*mesh->Np;
        maxU  = mymax(maxU, fabs(ins->Ux[id]));
        minU  = mymin(minU, fabs(ins->Ux[id]));
        //
        maxR = mymax(maxR, fabs(ins->rhsUx[id]));
        minR = mymin(minR, fabs(ins->rhsUx[id]));
      }
    }

  printf("minU: %g maxU: %g minR: %g maxR: %g \n", minU, maxU, minR, maxR);   
#endif  

printf("Solving for Ux \n");
 ellipticSolveTri2D( solver, ins->lamda, ins->o_rhsUx, ins->o_Ux, ins->velsolverOptions);

#if 0
dfloat maxU = 0, minU = 1e9;
dfloat maxR = 0, minR = 1e9;

ins->o_Uy.copyTo(ins->Uy);
ins->o_rhsUy.copyTo(ins->rhsUy);

 for(iint e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        iint id = n+e*mesh->Np;
        maxU  = mymax(maxU, fabs(ins->Uy[id]));
        minU  = mymin(minU, fabs(ins->Uy[id]));
        //
        maxR = mymax(maxR, fabs(ins->rhsUy[id]));
        minR = mymin(minR, fabs(ins->rhsUy[id]));
      }
    }

  printf("minU: %g maxU: %g minR: %g maxR: %g \n", minU, maxU, minR, maxR);   

#endif



printf("Solving for Uy \n");
  ellipticSolveTri2D(solver, ins->lamda, ins->o_rhsUy, ins->o_Uy, ins->velsolverOptions);

}
