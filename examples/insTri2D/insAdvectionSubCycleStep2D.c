#include "ins2D.h"

#define SUBSTEP_METHOD 1 // Modified OIFS method of Maday
//#define SUBSTEP_METHOD 2  // Lagrangian Subscyling of Maday, Patera, Ronquist 

// complete a time step using LSERK4
void insAdvectionSubCycleStep2D(ins_t *ins, iint tstep, 
				dfloat * tsendBuffer, dfloat * trecvBuffer, 
				dfloat * vsendBuffer, dfloat *  vrecvBuffer, 
				char   * options){

#if SUBSTEP_METHOD==1
 
 //printf("SUBSTEP METHOD : SEMI-LAGRAGIAN OIFS METHOD\n");
 // Solve Stokes Problem if Nonlinear solver is deactivated
 dfloat activate_advection = 0.f; 
 if(ins->a0){activate_advection  = 1.f;} 	

 mesh2D *mesh = ins->mesh;

const iint NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);
const iint Ntotal         =  NtotalElements*mesh->Np;  

const iint voffset = 0; 
// field offset at this step
  iint offset0 = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  
  //Exctract Halo On Device, Assumes History is already done!
  if(mesh->totalHaloPairs>0){
		ins->totalHaloExtractKernel(mesh->Nelements,
				mesh->totalHaloPairs,
				mesh->o_haloElementList,
				offset0,
				ins->o_U,
				ins->o_V,
				ins->o_P,
				ins->o_tHaloBuffer);

		// copy extracted halo to HOST
		ins->o_tHaloBuffer.copyTo(tsendBuffer);

		// start halo exchange    
		meshHaloExchangeStart(mesh,
			  mesh->Np*(ins->NTfields)*sizeof(dfloat),
			  tsendBuffer,
			  trecvBuffer);
  }

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(trecvBuffer);
    ins->totalHaloScatterKernel(mesh->Nelements,
				mesh->totalHaloPairs,
				mesh->o_haloElementList,
				offset0,
				ins->o_U,
				ins->o_V,
				ins->o_P,
				ins->o_tHaloBuffer);
  }



if(activate_advection){
	// 
		const dfloat tn0 = (tstep-0)*ins->dt;
		const dfloat tn1 = (tstep-1)*ins->dt;
		const dfloat tn2 = (tstep-2)*ins->dt;
		// construct interpolating lagrange polynomial
		dfloat c0 = 0.f, c1 = 0.f, c2 = 0.f;
  
		// Solve for Each SubProblem
		for (iint torder=0; torder<ins->ExplicitOrder; torder++){
		  // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
			const iint sindex = (ins->index + 3 - torder)%3; 
			ins->o_Ud.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,sindex*Ntotal*sizeof(dfloat));
			ins->o_Vd.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,sindex*Ntotal*sizeof(dfloat));
			
			// SubProblem  starts from here from t^(n-torder)
			const dfloat tsub = tstep*ins->dt - torder*ins->dt;
			// Advance SubProblem to t^(n+1) 
			for(iint ststep = 0; ststep<ins->Nsubsteps*(torder+1);++ststep){

				const dfloat tstage = tsub + ststep*ins->sdt;			
				// LSERK4 stages
				for(iint rk=0;rk<mesh->Nrk;++rk){
				// Extrapolate velocity to subProblem stage time
					dfloat t = tstage +  ins->sdt*mesh->rkc[rk]; 

					switch(ins->ExplicitOrder){
						case 1:
							c0 = 1.f; c1 = 0.f; c2 = 0.f;
						break;
						case 2:
							c0 = (t-tn1)/(tn0-tn1);
							c1 = (t-tn0)/(tn1-tn0);
							c2 = 0.f; 
						break;
						case 3:
							c0 = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
							c1 = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
							c2 = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
						break;
					   }

					//printf("c0: %.5f, c1: %.5f c3: %.5f \n", c0,c1,c2);
					ins->subCycleExtKernel(NtotalElements,
						                     ins->index,
						                     NtotalElements,
					                       c0,
					                       c1,
					                       c2,
					                       ins->o_U,
					                       ins->o_V,
					                       ins->o_Ue,
					                       ins->o_Ve);

					//
					if(mesh->totalHaloPairs>0){

						ins->velocityHaloExtractKernel(mesh->Nelements,
						                         mesh->totalHaloPairs,
						                         mesh->o_haloElementList,
						                         voffset, // 0 offset
						                         ins->o_Ud,
						                         ins->o_Vd,
						                         ins->o_vHaloBuffer);

						// copy extracted halo to HOST 
						ins->o_vHaloBuffer.copyTo(vsendBuffer);            

						// start halo exchange
						meshHaloExchangeStart(mesh,
						                    mesh->Np*(ins->NVfields)*sizeof(dfloat), 
						                    vsendBuffer,
						                    vrecvBuffer);
						}
            

            occaTimerTic(mesh->device,"AdvectionVolume");
							// Compute Volume Contribution
						if(strstr(options, "CUBATURE")){

							ins->subCycleCubatureVolumeKernel(mesh->Nelements,
							           mesh->o_vgeo,
							           mesh->o_cubDrWT,
							           mesh->o_cubDsWT,
							           mesh->o_cubInterpT,
							           ins->o_Ue,
							           ins->o_Ve,
							           ins->o_Ud,
							           ins->o_Vd,
							           ins->o_rhsU,
							           ins->o_rhsV);
						}
						else{
							//Compute Volume
							ins->subCycleVolumeKernel(mesh->Nelements,
							                          mesh->o_vgeo,
							                          mesh->o_DrT,
							                          mesh->o_DsT,
							                          ins->o_Ue,
							                          ins->o_Ve,
							                          ins->o_Ud,
							                          ins->o_Vd,
							                          ins->o_rhsU,
							                          ins->o_rhsV);

						}

						occaTimerToc(mesh->device,"AdvectionVolume");




						if(mesh->totalHaloPairs>0){

							meshHaloExchangeFinish(mesh);

							ins->o_vHaloBuffer.copyFrom(vrecvBuffer); 

							ins->velocityHaloScatterKernel(mesh->Nelements,
							                          mesh->totalHaloPairs,
							                          mesh->o_haloElementList,
							                          voffset, //0 offset
							                          ins->o_Ud,
							                          ins->o_Vd,
							                          ins->o_vHaloBuffer);
						}

            occaTimerTic(mesh->device,"AdvectionSurface");
							  // Compute Volume Contribution
						if(strstr(options, "CUBATURE")){
							ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
							                                    mesh->o_sgeo,
							                                    mesh->o_intInterpT,
							                                    mesh->o_intLIFTT,
							                                    mesh->o_vmapM,
							                                    mesh->o_vmapP,
							                                    mesh->o_EToB,
							                                    t,
							                                    mesh->o_intx,
							                                    mesh->o_inty,
							                                    ins->o_Ue,
							                                    ins->o_Ve,
							                                    ins->o_Ud,
							                                    ins->o_Vd,
							                                    ins->o_rhsU,
							                                    ins->o_rhsV);
						}
						else{
						 	//Surface Kernel
							ins->subCycleSurfaceKernel(mesh->Nelements,
							                          mesh->o_sgeo,
							                          mesh->o_LIFTT,
							                          mesh->o_vmapM,
							                          mesh->o_vmapP,
							                          mesh->o_EToB,
							                          t,
							                          mesh->o_x,
							                          mesh->o_y,
							                          ins->o_Ue,
							                          ins->o_Ve,
							                          ins->o_Ud,
							                          ins->o_Vd,
							                          ins->o_rhsU,
							                          ins->o_rhsV);

						}

						occaTimerToc(mesh->device,"AdvectionSurface");
            



            occaTimerTic(mesh->device,"AdvectionUpdate");
						// Update Kernel
						ins->subCycleRKUpdateKernel(mesh->Nelements,
						                      activate_advection,
						                      ins->sdt,
						                      mesh->rka[rk],
						                      mesh->rkb[rk],
						                      ins->o_rhsU,
						                      ins->o_rhsV,
						                      ins->o_resU, 
						                      ins->o_resV,
						                      ins->o_Ud,
						                      ins->o_Vd);
					occaTimerToc(mesh->device,"AdvectionUpdate");
				}
  		}
  	// Finish SubProblem // Use NU to store Lagrange Velocities
		ins->o_Ud.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),sindex*Ntotal*sizeof(dfloat),0);
		ins->o_Vd.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),sindex*Ntotal*sizeof(dfloat),0);
  	}

}


dfloat tp1 = (tstep+1)*ins->dt;

occaTimerTic(mesh->device,"GradientVolume");
 
// Compute Volume Contribution for Pressure
ins->gradientVolumeKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_DrT,
                          mesh->o_DsT,
                          offset0,
                          ins->o_P,
                          ins->o_Px,
                          ins->o_Py);
occaTimerToc(mesh->device,"GradientVolume");

	
 
occaTimerTic(mesh->device,"GradientSurface");
const iint solverid = 0; // Pressure Solve
  // Compute Surface Conribution
  ins->gradientSurfaceKernel(mesh->Nelements,
												     mesh->o_sgeo,
												     mesh->o_LIFTT,
												     mesh->o_vmapM,
												     mesh->o_vmapP,
												     mesh->o_EToB,
												     mesh->o_x,
												     mesh->o_y,
												     tp1,
												     ins->dt,
												     ins->c0,
												     ins->c1,
												     ins->c2,
												     ins->index,
												     mesh->Nelements+mesh->totalHaloPairs,
												     solverid, // pressure BCs
													   ins->o_PI, //not used
												     ins->o_P,
												     ins->o_Px,
												     ins->o_Py);

occaTimerToc(mesh->device,"GradientSurface");
#endif




#if SUBSTEP_METHOD==2
 
 //printf("SUBSTEP METHOD : SEMI-LAGRAGIAN OIFS METHOD\n");
 // Solve Stokes Problem if Nonlinear solver is deactivated
 dfloat activate_advection = 0.f; 
 if(ins->a0){activate_advection  = 1.f;} 	

 mesh2D *mesh = ins->mesh;

const iint NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);
const iint Ntotal         =  NtotalElements*mesh->Np;  

const iint voffset = 0; 

   // field offset at this step
  iint offset0 = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  
  //Exctract Halo On Device, Assumes History is already done!
  if(mesh->totalHaloPairs>0){
		ins->totalHaloExtractKernel(mesh->Nelements,
				mesh->totalHaloPairs,
				mesh->o_haloElementList,
				offset0,
				ins->o_U,
				ins->o_V,
				ins->o_P,
				ins->o_tHaloBuffer);

		// copy extracted halo to HOST
		ins->o_tHaloBuffer.copyTo(tsendBuffer);

		// start halo exchange    
		meshHaloExchangeStart(mesh,
			  mesh->Np*(ins->NTfields)*sizeof(dfloat),
			  tsendBuffer,
			  trecvBuffer);
  }

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(trecvBuffer);
    ins->totalHaloScatterKernel(mesh->Nelements,
				mesh->totalHaloPairs,
				mesh->o_haloElementList,
				offset0,
				ins->o_U,
				ins->o_V,
				ins->o_P,
				ins->o_tHaloBuffer);
  }



if(activate_advection){
	// 
		const dfloat tn0 = (tstep-0)*ins->dt;
		const dfloat tn1 = (tstep-1)*ins->dt;
		const dfloat tn2 = (tstep-2)*ins->dt;
		// construct interpolating lagrange polynomial
		dfloat c0 = 0.f, c1 = 0.f, c2 = 0.f;
  
		// Solve for Each SubProblem
		for (iint torder=0; torder<ins->ExplicitOrder; torder++){
		  // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
			const iint sindex = (ins->index + 3 - torder)%3; 
			ins->o_Ud.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,sindex*Ntotal*sizeof(dfloat));
			ins->o_Vd.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,sindex*Ntotal*sizeof(dfloat));
			
			// SubProblem  starts from here from t^(n-torder)
			const dfloat tsub = tstep*ins->dt - torder*ins->dt;
			// Advance SubProblem to t^(n+1) 
			for(iint ststep = 0; ststep<ins->Nsubsteps*(torder+1);++ststep){

				const dfloat tstage = tsub + ststep*ins->sdt;			
				// LSERK4 stages
				for(iint rk=0;rk<mesh->Nrk;++rk){
				// Extrapolate velocity to subProblem stage time
					dfloat t = tstage +  ins->sdt*mesh->rkc[rk]; 
					//
					if(mesh->totalHaloPairs>0){

						ins->velocityHaloExtractKernel(mesh->Nelements,
						                         mesh->totalHaloPairs,
						                         mesh->o_haloElementList,
						                         voffset, // 0 offset
						                         ins->o_Ud,
						                         ins->o_Vd,
						                         ins->o_vHaloBuffer);

						// copy extracted halo to HOST 
						ins->o_vHaloBuffer.copyTo(vsendBuffer);            

						// start halo exchange
						meshHaloExchangeStart(mesh,
						                    mesh->Np*(ins->NVfields)*sizeof(dfloat), 
						                    vsendBuffer,
						                    vrecvBuffer);
						}

											// Compute Volume Contribution
						  if(strstr(options, "CUBATURE")){
						    ins->advectionCubatureVolumeKernel(mesh->Nelements,
										       mesh->o_vgeo,
										       mesh->o_cubDrWT,
										       mesh->o_cubDsWT,
										       mesh->o_cubInterpT,
										       voffset,
										       ins->o_Ud,
										       ins->o_Vd,
										       ins->o_rhsU,
										       ins->o_rhsV);
						  } else {
						    ins->advectionVolumeKernel(mesh->Nelements,
									       mesh->o_vgeo,
									       mesh->o_DrT,
									       mesh->o_DsT,
									       voffset,
									       ins->o_Ud,
									       ins->o_Vd,
									       ins->o_rhsU,
									       ins->o_rhsV);
						  }



						if(mesh->totalHaloPairs>0){

							meshHaloExchangeFinish(mesh);

							ins->o_vHaloBuffer.copyFrom(vrecvBuffer); 

							ins->velocityHaloScatterKernel(mesh->Nelements,
							                          mesh->totalHaloPairs,
							                          mesh->o_haloElementList,
							                          voffset, //0 offset
							                          ins->o_Ud,
							                          ins->o_Vd,
							                          ins->o_vHaloBuffer);
						}


									 if(strstr(options, "CUBATURE")){
							ins->advectionCubatureSurfaceKernel(mesh->Nelements,
									mesh->o_sgeo,
									mesh->o_intInterpT,
									mesh->o_intLIFTT,
									mesh->o_vmapM,
									mesh->o_vmapP,
									mesh->o_EToB,
									t,
									mesh->o_intx,
									mesh->o_inty,
									voffset,
									ins->o_Ud,
									ins->o_Vd,
									ins->o_rhsU,
									ins->o_rhsV);
							} else {
							ins->advectionSurfaceKernel(mesh->Nelements,
								mesh->o_sgeo,
								mesh->o_LIFTT,
								mesh->o_vmapM,
								mesh->o_vmapP,
								mesh->o_EToB,
								t,
								mesh->o_x,
								mesh->o_y,
								voffset,
								ins->o_Ud,
								ins->o_Vd,
								ins->o_rhsU,
								ins->o_rhsV);
							}


						// Update Kernel
						ins->subCycleRKUpdateKernel(mesh->Nelements,
						                      activate_advection,
						                      ins->sdt,
						                      mesh->rka[rk],
						                      mesh->rkb[rk],
						                      ins->o_rhsU,
						                      ins->o_rhsV,
						                      ins->o_resU, 
						                      ins->o_resV,
						                      ins->o_Ud,
						                      ins->o_Vd);
				}
  		}
  	// Finish SubProblem // Use NU to store Lagrange Velocities
		ins->o_Ud.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),sindex*Ntotal*sizeof(dfloat),0);
		ins->o_Vd.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),sindex*Ntotal*sizeof(dfloat),0);
  	}

}


dfloat tp1 = (tstep+1)*ins->dt;
// Compute Volume Contribution for Pressure
ins->gradientVolumeKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_DrT,
                          mesh->o_DsT,
                          offset0,
                          ins->o_P,
                          ins->o_Px,
                          ins->o_Py);

	
 

const iint solverid = 0; // Pressure Solve
  // Compute Surface Conribution
  ins->gradientSurfaceKernel(mesh->Nelements,
												     mesh->o_sgeo,
												     mesh->o_LIFTT,
												     mesh->o_vmapM,
												     mesh->o_vmapP,
												     mesh->o_EToB,
												     mesh->o_x,
												     mesh->o_y,
												     tp1,
												     ins->dt,
												     ins->c0,
												     ins->c1,
												     ins->c2,
												     ins->index,
												     mesh->Nelements+mesh->totalHaloPairs,
												     solverid, // pressure BCs
													   ins->o_PI, //not used
												     ins->o_P,
												     ins->o_Px,
												     ins->o_Py);


#endif


// #if SUBSTEP_METHOD==2

//   // field offset at this step
//   iint offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);

//   //Exctract Halo On Device
//   if(mesh->totalHaloPairs>0){
//     ins->totalHaloExtractKernel(mesh->Nelements,
// 				mesh->totalHaloPairs,
// 				mesh->o_haloElementList,
// 				offset,
// 				ins->o_U,
// 				ins->o_V,
// 				ins->o_P,
// 				ins->o_tHaloBuffer);

//     // copy extracted halo to HOST
//     ins->o_tHaloBuffer.copyTo(sendBuffer);

//     // start halo exchange    meshHaloExchangeStart(mesh,
// 			  mesh->Np*(ins->NTfields)*sizeof(dfloat),
// 			  sendBuffer,
// 			  recvBuffer);
//   }

//   // COMPLETE HALO EXCHANGE
//   if(mesh->totalHaloPairs>0){

//     meshHaloExchangeFinish(mesh);

//     ins->o_tHaloBuffer.copyFrom(recvBuffer);
//     ins->totalHaloScatterKernel(mesh->Nelements,
// 				mesh->totalHaloPairs,
// 				mesh->o_haloElementList,
// 				offset,
// 				ins->o_U,
// 				ins->o_V,
// 				ins->o_P,
// 				ins->o_tHaloBuffer);
//   }

//    // Solve Stokes Problem if Nonlinear solver is deactivated
//   dfloat activate_advection = 0.f; 
//   if(ins->a0){activate_advection  = 1.f;} 

//   const iint voffset = 0;  
//   // New subcycling // Initialize U^e = U^n
//   iint Ntotal =  (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
//   ins->o_Ue.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
//   ins->o_Ve.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  

// // Substepping 
// 	if(activate_advection){
//   	for(iint ststep = 0; ststep<ins->Nsubsteps;++ststep){
//     	dfloat tbase = t + ststep*ins->sdt;    
//     	// LSERK4 stages
//     		for(iint rk=0;rk<mesh->Nrk;++rk){
// 				// intermediate stage time
// 					dfloat tstage = tbase +  ins->sdt*mesh->rkc[rk]; 

// 				//
// 					if(mesh->totalHaloPairs>0){
// 						ins->velocityHaloExtractKernel(mesh->Nelements,
// 																					mesh->totalHaloPairs,
// 																					mesh->o_haloElementList,
// 																					voffset, // 0 offset
// 																					ins->o_Ue,
// 																					ins->o_Ve,
// 																					ins->o_vHaloBuffer);
// 						// copy extracted halo to HOST 
// 						ins->o_vHaloBuffer.copyTo(sendBuffer);            
// 						// start halo exchange
// 						meshHaloExchangeStart(mesh,
// 																	mesh->Np*(ins->NVfields)*sizeof(dfloat), 
// 																	sendBuffer,
// 																	recvBuffer);
// 					}

// 					// Compute Volume Contribution
// 					if(strstr(options, "CUBATURE")){
// 						ins->advectionCubatureVolumeKernel(mesh->Nelements,
// 															mesh->o_vgeo,
// 															mesh->o_cubDrWT,
// 															mesh->o_cubDsWT,
// 															mesh->o_cubInterpT,
// 															voffset, // zero offset
// 															ins->o_Ue,
// 															ins->o_Ve,
// 															ins->o_rhsU,
// 															ins->o_rhsV);
// 					} 
// 					else {
// 						ins->advectionVolumeKernel(mesh->Nelements,
// 																       mesh->o_vgeo,
// 																       mesh->o_DrT,
// 																       mesh->o_DsT,
// 																       voffset, //zero offset
// 																       ins->o_Ue,
// 																       ins->o_Ve,
// 																       ins->o_rhsU,
// 																       ins->o_rhsV);
// 						}


// 					if(mesh->totalHaloPairs>0){

// 						meshHaloExchangeFinish(mesh);

// 						ins->o_vHaloBuffer.copyFrom(recvBuffer); 

// 						ins->velocityHaloScatterKernel(mesh->Nelements,
// 																					 mesh->totalHaloPairs,
// 																					 mesh->o_haloElementList,
// 																					 voffset, //0 offset
// 																					 ins->o_Ue,
// 																					 ins->o_Ve,
// 																					 ins->o_vHaloBuffer);
// 					}

    
// 				  if(strstr(options, "CUBATURE")){
// 				  ins->advectionCubatureSurfaceKernel(mesh->Nelements,
// 								      mesh->o_sgeo,
// 								      mesh->o_intInterpT,
// 								      mesh->o_intLIFTT,
// 								      mesh->o_vmapM,
// 								      mesh->o_vmapP,
// 								      mesh->o_EToB,
// 								      tstage,
// 								      mesh->o_intx,
// 								      mesh->o_inty,
// 								      voffset, // 0
// 								      ins->o_Ue,
// 								      ins->o_Ve,
// 								      ins->o_rhsU,
// 								      ins->o_rhsV);
// 				} 
// 				else {
// 				  ins->advectionSurfaceKernel(mesh->Nelements,
// 							      mesh->o_sgeo,
// 							      mesh->o_LIFTT,
// 							      mesh->o_vmapM,
// 							      mesh->o_vmapP,
// 							      mesh->o_EToB,
// 							      tstage,
// 							      mesh->o_x,
// 							      mesh->o_y,
// 							      voffset, // 0
// 							      ins->o_Ue,
// 							      ins->o_Ve,
// 							      ins->o_rhsU,
// 							      ins->o_rhsV);
// 				}

// 			  // Update Kernel
// 				ins->subCycleRKUpdateKernel(mesh->Nelements,
// 							    activate_advection,
// 							    ins->sdt,
// 							    mesh->rka[rk],
// 							    mesh->rkb[rk],
// 							    ins->o_rhsU,
// 							    ins->o_rhsV,
// 							    ins->o_resU, 
// 							    ins->o_resV,
// 							    ins->o_Ue,
// 							    ins->o_Ve);
// 			}
// 		}
// 	}



//   // Compute Volume Contribution for Pressure
//   ins->gradientVolumeKernel(mesh->Nelements,
//                             mesh->o_vgeo,
//                             mesh->o_DrT,
//                             mesh->o_DsT,
//                             offset,
//                             ins->o_P,
//                             ins->o_Px,
//                             ins->o_Py);

	
 

// const iint solverid = 0; // Pressure Solve
//   // Compute Surface Conribution
//   ins->gradientSurfaceKernel(mesh->Nelements,
// 												     mesh->o_sgeo,
// 												     mesh->o_LIFTT,
// 												     mesh->o_vmapM,
// 												     mesh->o_vmapP,
// 												     mesh->o_EToB,
// 												     mesh->o_x,
// 												     mesh->o_y,
// 												     t,
// 												     ins->dt,
// 												     ins->c0,
// 												     ins->c1,
// 												     ins->c2,
// 												     ins->index,
// 												     mesh->Nelements+mesh->totalHaloPairs,
// 												     solverid, // pressure BCs
// 													   ins->o_PI, //not used
// 												     ins->o_P,
// 												     ins->o_Px,
// 												     ins->o_Py);


// //iint index1 = ins->index;
// // Use NU to store Ue
// ins->o_Ue.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat),0);
// ins->o_Ve.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),ins->index*Ntotal*sizeof(dfloat),0);

// #endif




// #if SUBSTEP_METHOD==3
//  mesh2D *mesh = ins->mesh;
//  dfloat     t = tstep*ins->dt;
 
//    // field offset at this step
//   iint offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);

//   //Exctract Halo On Device
//   if(mesh->totalHaloPairs>0){
//     ins->totalHaloExtractKernel(mesh->Nelements,
// 				mesh->totalHaloPairs,
// 				mesh->o_haloElementList,
// 				offset,
// 				ins->o_U,
// 				ins->o_V,
// 				ins->o_P,
// 				ins->o_tHaloBuffer);

//     // copy extracted halo to HOST
//     ins->o_tHaloBuffer.copyTo(sendBuffer);

//     // start halo exchange
//     meshHaloExchangeStart(mesh,
// 			  mesh->Np*(ins->NTfields)*sizeof(dfloat),
// 			  sendBuffer,
// 			  recvBuffer);
//   }

//   // COMPLETE HALO EXCHANGE
//   if(mesh->totalHaloPairs>0){

//     meshHaloExchangeFinish(mesh);

//     ins->o_tHaloBuffer.copyFrom(recvBuffer);
//     ins->totalHaloScatterKernel(mesh->Nelements,
// 				mesh->totalHaloPairs,
// 				mesh->o_haloElementList,
// 				offset,
// 				ins->o_U,
// 				ins->o_V,
// 				ins->o_P,
// 				ins->o_tHaloBuffer);
//   }






//   printf("SUBSTEP METHOD :2 \n");
//    // Solve Stokes Problem if Nonlinear solver is deactivated
//   dfloat activate_advection = 0.f; 
//   if(ins->a0){activate_advection  = 1.f;} 

//   const iint voffset = 0;  
//   //
//   iint Ntotal =  (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
// 	ins->o_Ue.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
// 	ins->o_Ve.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
// 	ins->o_Ud.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
// 	ins->o_Vd.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
  

// // Substepping 
// 	if(activate_advection){
//   	for(iint ststep = 0; ststep<ins->Nsubsteps;++ststep){
//     	dfloat tbase = t + ststep*ins->sdt;    
//     	// LSERK4 stages
//     		for(iint rk=0;rk<mesh->Nrk;++rk){
// 				// intermediate stage time
// 					dfloat tstage = tbase +  ins->sdt*mesh->rkc[rk]; 

// 				if(mesh->totalHaloPairs>0){
 
// 						ins->velocityHaloExtractKernel(mesh->Nelements,
// 						                         mesh->totalHaloPairs,
// 						                         mesh->o_haloElementList,
// 						                         voffset, // 0 offset
// 						                         ins->o_Ud,
// 						                         ins->o_Vd,
// 						                         ins->o_vHaloBuffer);

// 						// copy extracted halo to HOST 
// 						ins->o_vHaloBuffer.copyTo(sendBuffer);            

// 						// start halo exchange
// 						meshHaloExchangeStart(mesh,
// 						                    mesh->Np*(ins->NVfields)*sizeof(dfloat), 
// 						                    sendBuffer,
// 						                    recvBuffer);
// 					}

// 								// Compute Volume Contribution
// 					if(strstr(options, "CUBATURE")){

// 					  ins->subCycleCubatureVolumeKernel(mesh->Nelements,
// 					             mesh->o_vgeo,
// 					             mesh->o_cubDrWT,
// 					             mesh->o_cubDsWT,
// 					             mesh->o_cubInterpT,
// 					             ins->o_Ue,
// 					             ins->o_Ve,
// 					             ins->o_Ud,
// 					             ins->o_Vd,
// 					             ins->o_rhsU,
// 					             ins->o_rhsV);
// 					}
// 					else{
// 					  //Compute Volume
// 					  ins->subCycleVolumeKernel(mesh->Nelements,
// 					                            mesh->o_vgeo,
// 					                            mesh->o_DrT,
// 					                            mesh->o_DsT,
// 					                            ins->o_Ue,
// 					                            ins->o_Ve,
// 					                            ins->o_Ud,
// 					                            ins->o_Vd,
// 					                            ins->o_rhsU,
// 					                            ins->o_rhsV);

// 					}


// 					if(mesh->totalHaloPairs>0){
  
// 				    meshHaloExchangeFinish(mesh);

// 				    ins->o_vHaloBuffer.copyFrom(recvBuffer); 

// 				    ins->velocityHaloScatterKernel(mesh->Nelements,
// 				                                mesh->totalHaloPairs,
// 				                                mesh->o_haloElementList,
// 				                                voffset, //0 offset
// 				                                ins->o_Ud,
// 				                                ins->o_Vd,
// 				                                ins->o_vHaloBuffer);
// 					}


// 								  // Compute Volume Contribution
// 					if(strstr(options, "CUBATURE")){
// 					  ins->subCycleCubatureSurfaceKernel(mesh->Nelements,
// 					                                      mesh->o_sgeo,
// 					                                      mesh->o_intInterpT,
// 					                                      mesh->o_intLIFTT,
// 					                                      mesh->o_vmapM,
// 					                                      mesh->o_vmapP,
// 					                                      mesh->o_EToB,
// 					                                      t,
// 					                                      mesh->o_intx,
// 					                                      mesh->o_inty,
// 					                                      ins->o_Ue,
// 					                                      ins->o_Ve,
// 					                                      ins->o_Ud,
// 					                                      ins->o_Vd,
// 					                                      ins->o_rhsU,
// 					                                      ins->o_rhsV);
// 					  }
// 					else{
// 					   //Surface Kernel
// 					  ins->subCycleSurfaceKernel(mesh->Nelements,
// 					                            mesh->o_sgeo,
// 					                            mesh->o_LIFTT,
// 					                            mesh->o_vmapM,
// 					                            mesh->o_vmapP,
// 					                            mesh->o_EToB,
// 					                            t,
// 					                            mesh->o_x,
// 					                            mesh->o_y,
// 					                            ins->o_Ue,
// 					                            ins->o_Ve,
// 					                            ins->o_Ud,
// 					                            ins->o_Vd,
// 					                            ins->o_rhsU,
// 					                            ins->o_rhsV);

// 					}


// 			  // Update Kernel
//   ins->subCycleRKUpdateKernel(mesh->Nelements,
//                               activate_advection,
//                               ins->sdt,
//                               mesh->rka[rk],
//                               mesh->rkb[rk],
//                               ins->o_rhsU,
//                               ins->o_rhsV,
//                               ins->o_resU, 
//                               ins->o_resV,
//                               ins->o_Ud,
//                               ins->o_Vd);
  
// 			  //printf("Extrapolating Velocity to %d \n", ststep+1);
// 				// Extrapolate Velocity
// 				const dfloat t1 = tstep*ins->dt, t2 = (tstep-1)*ins->dt, t3 = (tstep-2)*ins->dt;
// 				// construct interpolating lagrange polynomial
// 				dfloat c0 = 0, c1 = 0, c2 = 0;
// 				#if 0
// 				c0 = 1;
// 				c1 = 0;
// 				c2 = 0;
// 				#endif
// 				#if 0
// 				c0 = (tstage-t2)/(t1-t2);
// 				c1 = (tstage-t1)/(t2-t1);
// 				c2 = 0;
// 				#endif
// 				#if 1

// 				c0 = (tstage-t2)*(tstage-t3)/((t1-t2)*(t1-t3)); 
// 				c1 = (tstage-t1)*(tstage-t3)/((t2-t1)*(t2-t3));
// 				c2 = (tstage-t1)*(tstage-t2)/((t3-t1)*(t3-t2));
// 				#endif

// 				iint offset0 = ((ins->index+0)%3)*Ntotal;
// 				iint offset1 = ((ins->index+2)%3)*Ntotal;
// 				iint offset2 = ((ins->index+1)%3)*Ntotal;

//         printf("c0: %.5f, c1: %.5f c3: %.5f \n", c0,c1,c2);
// 				ins->subCycleExtKernel(mesh->Nelements + mesh->totalHaloPairs,
// 					                     c0,
// 					                     c1,
// 					                     c2,
// 				                       offset0,
// 				                       offset1,
// 				                       offset2,
// 				                       ins->o_U,
// 				                       ins->o_V,
// 				                       ins->o_Ue,
// 				                       ins->o_Ve);

//       }
		


// 		  }
// 	}



//   // Compute Volume Contribution for Pressure
//   ins->gradientVolumeKernel(mesh->Nelements,
//                             mesh->o_vgeo,
//                             mesh->o_DrT,
//                             mesh->o_DsT,
//                             offset,
//                             ins->o_P,
//                             ins->o_Px,
//                             ins->o_Py);

	
 

// const iint solverid = 0; // Pressure Solve
//   // Compute Surface Conribution
//   ins->gradientSurfaceKernel(mesh->Nelements,
// 												     mesh->o_sgeo,
// 												     mesh->o_LIFTT,
// 												     mesh->o_vmapM,
// 												     mesh->o_vmapP,
// 												     mesh->o_EToB,
// 												     mesh->o_x,
// 												     mesh->o_y,
// 												     t,
// 												     ins->dt,
// 												     ins->c0,
// 												     ins->c1,
// 												     ins->c2,
// 												     ins->index,
// 												     mesh->Nelements+mesh->totalHaloPairs,
// 												     solverid, // pressure BCs
// 													   ins->o_PI, //not used
// 												     ins->o_P,
// 												     ins->o_Px,
// 												     ins->o_Py);


// iint index1 = ins->index;
// // Use NU to store Ue
// ins->o_Ud.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
// ins->o_Vd.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);

// #endif


// #if SUBSTEP_METHOD==3


//   mesh2D *mesh = ins->mesh; 
//   // field offset 
//   iint offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);  

//   // Subcycling with Burger Type, Modified

//   // May be modified : Solver may do smt while exchanging
//   if(mesh->totalHaloPairs>0){
//     ins->totalHaloExtractKernel(mesh->Nelements,
//                                 mesh->totalHaloPairs,
//                                 mesh->o_haloElementList,
//                                 offset,
//                                 ins->o_U,
//                                 ins->o_V,
//                                 ins->o_P,
//                                 ins->o_tHaloBuffer);
//     // copy extracted halo to HOST
//     ins->o_tHaloBuffer.copyTo(sendBuffer);
//     // start halo exchange
//     meshHaloExchangeStart(mesh,
//                           mesh->Np*(ins->NTfields)*sizeof(dfloat),
//                           sendBuffer,
//                           recvBuffer);
//   }


//   // Compute Volume Contribution of Pressure
//   ins->gradientVolumeKernel(mesh->Nelements,
// 			    mesh->o_vgeo,
// 			    mesh->o_DrT,
// 			    mesh->o_DsT,
// 			    offset,
// 			    ins->o_P,
// 			    ins->o_Px,
// 			    ins->o_Py);


//   if(mesh->totalHaloPairs>0){
//     meshHaloExchangeFinish(mesh);

//     ins->o_tHaloBuffer.copyFrom(recvBuffer);

//     ins->totalHaloScatterKernel(mesh->Nelements,
//                                 mesh->totalHaloPairs,
//                                 mesh->o_haloElementList,
//                                 offset,
//                                 ins->o_U,
//                                 ins->o_V,
//                                 ins->o_P,
//                                 ins->o_tHaloBuffer);
//   }


//   // Compute Pressure Surface Contribution
//   dfloat t = tstep*ins->dt;  
//   ins->gradientSurfaceKernel(mesh->Nelements,
// 			     mesh->o_sgeo,
// 			     mesh->o_LIFTT,
// 			     mesh->o_vmapM,
// 			     mesh->o_vmapP,
// 			     mesh->o_EToB,
// 			     mesh->o_x,
// 			     mesh->o_y,
// 			     t,
// 			     ins->dt,
// 			     ins->a0, // not used
// 			     ins->a1,
// 			     ins->a2,
// 			     ins->index,
// 			     mesh->Nelements+mesh->totalHaloPairs,
// 			     0, 
// 			     ins->o_PI, //not used
// 			     ins->o_P,
// 			     ins->o_Px,
// 			     ins->o_Py);  


//   // Solve Stokes Problem if Nonlinear solver is deactivated
//   dfloat activate_advection = 0.f; 
//   if(ins->a0){activate_advection  = 1.f;} // rhsU = -N(U) 

//   const iint voffset = 0;   // Use Nelements*Np size vector for Ue


//   //#if 1 // New subcycling
//   iint Ntotal =  (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

//   ins->o_Ue.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
//   ins->o_Ve.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));

//   if(activate_advection){
//     for(iint ststep = 0; ststep<ins->Nsubsteps;++ststep){
//       dfloat time = tstep*ins->dt + ststep*ins->sdt;    
//       // LSERK4 stages
//       for(iint rk=0;rk<mesh->Nrk;++rk){
// 	// intermediate stage time
// 	dfloat t = time +  ins->sdt*mesh->rkc[rk]; 

// 	if(mesh->totalHaloPairs>0){
// 	  ins->velocityHaloExtractKernel(mesh->Nelements,
// 					 mesh->totalHaloPairs,
// 					 mesh->o_haloElementList,
// 					 voffset, // 0 offset
// 					 ins->o_Ue,
// 					 ins->o_Ve,
// 					 ins->o_vHaloBuffer);
// 	  // copy extracted halo to HOST 
// 	  ins->o_vHaloBuffer.copyTo(sendBuffer);            
// 	  // start halo exchange
// 	  meshHaloExchangeStart(mesh,
// 				mesh->Np*(ins->NVfields)*sizeof(dfloat), 
// 				sendBuffer,
// 				recvBuffer);
// 	}


//         // Compute Volume Contribution
// 	if(strstr(options, "CUBATURE")){
// 	  ins->advectionCubatureVolumeKernel(mesh->Nelements,
// 					     mesh->o_vgeo,
// 					     mesh->o_cubDrWT,
// 					     mesh->o_cubDsWT,
// 					     mesh->o_cubInterpT,
// 					     voffset, // 0
// 					     ins->o_Ue,
// 					     ins->o_Ve,
// 					     ins->o_rhsU,
// 					     ins->o_rhsV);
// 	} else {
// 	  ins->advectionVolumeKernel(mesh->Nelements,
// 				     mesh->o_vgeo,
// 				     mesh->o_DrT,
// 				     mesh->o_DsT,
// 				     voffset, // 0
// 				     ins->o_Ue,
// 				     ins->o_Ve,
// 				     ins->o_rhsU,
// 				     ins->o_rhsV);
// 	}


// 	if(mesh->totalHaloPairs>0){

// 	  meshHaloExchangeFinish(mesh);

// 	  ins->o_vHaloBuffer.copyFrom(recvBuffer); 

// 	  ins->velocityHaloScatterKernel(mesh->Nelements,
// 					 mesh->totalHaloPairs,
// 					 mesh->o_haloElementList,
// 					 voffset, //0 offset
// 					 ins->o_Ue,
// 					 ins->o_Ve,
// 					 ins->o_vHaloBuffer);
// 	}


// 	if(strstr(options, "CUBATURE")){
// 	  ins->advectionCubatureSurfaceKernel(mesh->Nelements,
// 					      mesh->o_sgeo,
// 					      mesh->o_intInterpT,
// 					      mesh->o_intLIFTT,
// 					      mesh->o_vmapM,
// 					      mesh->o_vmapP,
// 					      mesh->o_EToB,
// 					      t,
// 					      mesh->o_intx,
// 					      mesh->o_inty,
// 					      voffset, // 0
// 					      ins->o_Ue,
// 					      ins->o_Ve,
// 					      ins->o_rhsU,
// 					      ins->o_rhsV);
// 	} else {
// 	  ins->advectionSurfaceKernel(mesh->Nelements,
// 				      mesh->o_sgeo,
// 				      mesh->o_LIFTT,
// 				      mesh->o_vmapM,
// 				      mesh->o_vmapP,
// 				      mesh->o_EToB,
// 				      t,
// 				      mesh->o_x,
// 				      mesh->o_y,
// 				      voffset, // 0
// 				      ins->o_Ue,
// 				      ins->o_Ve,
// 				      ins->o_rhsU,
// 				      ins->o_rhsV);
// 	}

// 	// t is current time
// 	const dfloat t1 = tstep*ins->dt, t2 = (tstep-1)*ins->dt, t3 = (tstep-2)*ins->dt;
	
// 	// construct interpolating lagrange polynomial
// 	dfloat c0 = 0, c1 = 0, c2 = 0;
// #if 0
// 	c0 = 1;
// 	c1 = 0;
// 	c2 = 0;
// #endif
// #if 0
// 	c0 = (t-t2)/(t1-t2);
// 	c1 = (t-t1)/(t2-t1);
// 	c2 = 0;
// #endif
// #if 0
// 	c0 = (t-t2)*(t-t3)/((t1-t2)*(t1-t3)); 
// 	c1 = (t-t1)*(t-t3)/((t2-t1)*(t2-t3));
// 	c2 = (t-t1)*(t-t2)/((t3-t1)*(t3-t2));
// #endif
	
// 	iint offset0 = ((ins->index+0)%3)*Ntotal;
// 	iint offset1 = ((ins->index+2)%3)*Ntotal;
// 	iint offset2 = ((ins->index+1)%3)*Ntotal;

// 	// Update Kernel
// 	ins->subCycleRKUpdateKernel(mesh->Nelements,
// 				    activate_advection,
// 				    ins->sdt,
// 				    mesh->rka[rk],
// 				    mesh->rkb[rk],
// 				    ins->o_rhsU,
// 				    ins->o_rhsV,
// 				    ins->o_resU, 
// 				    ins->o_resV,
// 				    offset0, c0,
// 				    offset1, c1,
// 				    offset2, c2,
// 				    ins->o_Px,
// 				    ins->o_Py,
// 				    ins->o_Ue,
// 				    ins->o_Ve);

      
//       }
//     }
//   }




//   // Now Compute N(Ue) Term in n+1

//    dfloat tn1 = (tstep+1)*ins->dt;

//   if(mesh->totalHaloPairs>0){

//     ins->velocityHaloExtractKernel(mesh->Nelements,
// 				   mesh->totalHaloPairs,
// 				   mesh->o_haloElementList,
// 				   voffset,  //voffset
// 				   ins->o_Ue,
// 				   ins->o_Ve,
// 				   ins->o_vHaloBuffer);

//     // copy extracted halo to HOST 
//     ins->o_vHaloBuffer.copyTo(sendBuffer);            

//     // start halo exchange
//     meshHaloExchangeStart(mesh,
// 			  mesh->Np*(ins->NVfields)*sizeof(dfloat), 
// 			  sendBuffer,
// 			  recvBuffer);
//   }


//   // Compute Volume Contribution
//   if(strstr(options, "CUBATURE")){
//     ins->advectionCubatureVolumeKernel(mesh->Nelements,
// 				       mesh->o_vgeo,
// 				       mesh->o_cubDrWT,
// 				       mesh->o_cubDsWT,
// 				       mesh->o_cubInterpT,
// 				       voffset, // 0
// 				       ins->o_Ue,
// 				       ins->o_Ve,
// 				       ins->o_rhsU,
// 				       ins->o_rhsV);
//   } else {
//     ins->advectionVolumeKernel(mesh->Nelements,
// 			       mesh->o_vgeo,
// 			       mesh->o_DrT,
// 			       mesh->o_DsT,
// 			       voffset, // 0
// 			       ins->o_Ue,
// 			       ins->o_Ve,
// 			       ins->o_rhsU,
// 			       ins->o_rhsV);
//   }


//   if(mesh->totalHaloPairs>0){

//     meshHaloExchangeFinish(mesh);

//     ins->o_vHaloBuffer.copyFrom(recvBuffer); 

//     ins->velocityHaloScatterKernel(mesh->Nelements,
// 				   mesh->totalHaloPairs,
// 				   mesh->o_haloElementList,
// 				   voffset, //0 offset
// 				   ins->o_Ue,
// 				   ins->o_Ve,
// 				   ins->o_vHaloBuffer);
//   }


//   if(strstr(options, "CUBATURE")){
//     ins->advectionCubatureSurfaceKernel(mesh->Nelements,
// 					mesh->o_sgeo,
// 					mesh->o_intInterpT,
// 					mesh->o_intLIFTT,
// 					mesh->o_vmapM,
// 					mesh->o_vmapP,
// 					mesh->o_EToB,
// 					tn1,
// 					mesh->o_intx,
// 					mesh->o_inty,
// 					voffset, // 0
// 					ins->o_Ue,
// 					ins->o_Ve,
// 					ins->o_rhsU,
// 					ins->o_rhsV);
//   } else {
//     ins->advectionSurfaceKernel(mesh->Nelements,
// 				mesh->o_sgeo,
// 				mesh->o_LIFTT,
// 				mesh->o_vmapM,
// 				mesh->o_vmapP,
// 				mesh->o_EToB,
// 				tn1,
// 				mesh->o_x,
// 				mesh->o_y,
// 				voffset, // 0
// 				ins->o_Ue,
// 				ins->o_Ve,
// 				ins->o_rhsU,
// 				ins->o_rhsV);
//   }

//   iint index1 = ins->index;
//   ins->o_rhsU.copyTo(ins->o_NU,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
//   ins->o_rhsV.copyTo(ins->o_NV,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);


//   //#endif
// #endif



}    
