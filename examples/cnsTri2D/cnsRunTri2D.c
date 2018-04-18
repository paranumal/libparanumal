#include "cnsTri2D.h"

void cnsRunTri2D(cns_t *cns, char *options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = cns->mesh;

  cnsReportTri2D(cns, 0, options);

  occa::timer timer;
  
  timer.initTimer(mesh->device);

  timer.tic("Run");
  
  if (strstr(options,"DOPRI5")) {
    // Dormand Prince -order (4) 5 with PID timestep control
    dfloat rkC[7] = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    dfloat rkA[7*7] ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                   0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                              3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                             44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                        19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                         9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0, 
                            35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
    dfloat rkE[7] = {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 }; 

    occa::memory o_rkA = mesh->device.malloc(7*7*sizeof(dfloat), rkA);
    occa::memory o_rkE = mesh->device.malloc(  7*sizeof(dfloat), rkE);
  
    dfloat dtMIN = 1E-7; //minumum allowed timestep
    dfloat ATOL = 1E-5;  //absolute error tolerance
    dfloat RTOL = 1E-3;  //relative error tolerance
    dfloat safe = 0.9;   //safety factor

    //error control parameters
    dfloat beta = 0.05;
    dfloat factor1 = 0.2;
    dfloat factor2 = 10.0;


    dfloat exp1 = 0.2 - 0.75*beta;
    dfloat invfactor1 = 1.0/factor1;
    dfloat invfactor2 = 1.0/factor2;
    dfloat facold = 1E-4;

    // hard code this for the moment
    dfloat outputInterval = .5;
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
    //initial time
    dfloat time = 0.0;
    int tstep=0, allStep = 0;

    int done =0;
    while (!done) {

      int advSwitch = 1;
      
      if (mesh->dt<dtMIN){
        printf("ERROR: Time step became too small at time step=%d\n", tstep);
        exit (-1);
      }
      if (isnan(mesh->dt)) {
        printf("ERROR: Solution became unstable at time step=%d\n", tstep);
        exit (-1);
      }

      // check for next output
      int isOutput = 0;
      if((time+mesh->dt > nextOutputTime) &&
	 (time<=nextOutputTime)){
	isOutput = 1;
	mesh->dt = nextOutputTime-time;
	
      }

      //check for final timestep
      if (time+mesh->dt > mesh->finalTime){
	mesh->dt = mesh->finalTime-time;
	done = 1;
	isOutput = 0;
      }

      
      //RK step
      for(int rk=0;rk<7;++rk){

        // t_rk = t + C_rk*dt
        dfloat currentTime = time + rkC[rk]*mesh->dt;
        
        //compute RK stage 
        // rkq = q + dt sum_{i=0}^{rk-1} a_{rk,i}*rhsq_i
        cns->rkStageKernel(mesh->Nelements,
                           rk,
                           mesh->dt,
                           o_rkA,
                           cns->o_q,
                           cns->o_rkrhsq,
                           cns->o_rkq);

        //compute RHS
        // rhsq = F(currentTIme, rkq)

        // extract q halo on DEVICE
        if(mesh->totalHaloPairs>0){
          int Nentries = mesh->Np*cns->Nfields;          
          mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_rkq, cns->o_haloBuffer);
          
          // copy extracted halo to HOST 
          cns->o_haloBuffer.copyTo(cns->sendBuffer);      
          
          // start halo exchange
          meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
        }

        // now compute viscous stresses
        cns->stressesVolumeKernel(mesh->Nelements, 
                                  mesh->o_vgeo, 
                                  mesh->o_DrT, 
                                  mesh->o_DsT, 
                                  cns->mu, 
                                  cns->o_rkq, 
                                  cns->o_viscousStresses);

        // wait for q halo data to arrive
        if(mesh->totalHaloPairs>0){
          meshHaloExchangeFinish(mesh);
          
          // copy halo data to DEVICE
          size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
          cns->o_rkq.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
        }
        
        cns->stressesSurfaceKernel(mesh->Nelements, 
                                   mesh->o_sgeo, 
                                   mesh->o_LIFTT,
                                   mesh->o_vmapM, 
                                   mesh->o_vmapP, 
                                   mesh->o_EToB, 
                                   currentTime,
                                   mesh->o_x, 
                                   mesh->o_y, 
                                   cns->mu, 
                                   cns->o_rkq, 
                                   cns->o_viscousStresses);
        
        // extract stresses halo on DEVICE
        if(mesh->totalHaloPairs>0){
          int Nentries = mesh->Np*cns->Nstresses;
          
          mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
          
          // copy extracted halo to HOST 
          cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
          
          // start halo exchange
          meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
        }
        
        // compute volume contribution to DG cns RHS
        if (strstr(options,"CUBATURE")) {
          cns->cubatureVolumeKernel(mesh->Nelements, 
                                    advSwitch, 
                                    mesh->o_vgeo, 
                                    mesh->o_cubDrWT,
                                    mesh->o_cubDsWT,
                                    mesh->o_cubInterpT,
                                    cns->o_viscousStresses, 
                                    cns->o_rkq, 
                                    cns->o_rhsq);
        } else {
          cns->volumeKernel(mesh->Nelements, 
                            advSwitch, 
                            mesh->o_vgeo, 
                            mesh->o_DrT,
                            mesh->o_DsT,
                            cns->o_viscousStresses, 
                            cns->o_rkq, 
                            cns->o_rhsq);
        }

        // wait for halo stresses data to arrive
        if(mesh->totalHaloPairs>0){
          meshHaloExchangeFinish(mesh);
          
          // copy halo data to DEVICE
          size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
          cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
        }
        
        // compute surface contribution to DG cns RHS (LIFTT ?)
        if (strstr(options,"CUBATURE")) {
          cns->cubatureSurfaceKernel(mesh->Nelements, 
                                     advSwitch, 
                                     mesh->o_sgeo, 
                                     mesh->o_intInterpT,
                                     mesh->o_intLIFTT, 
                                     mesh->o_vmapM, 
                                     mesh->o_vmapP, 
                                     mesh->o_EToB,
                                     currentTime, 
                                     mesh->o_intx, 
                                     mesh->o_inty, 
                                     cns->mu, 
                                     cns->o_rkq, 
                                     cns->o_viscousStresses, 
                                     cns->o_rhsq);
        } else {
          cns->surfaceKernel(mesh->Nelements, 
                             advSwitch, 
                             mesh->o_sgeo, 
                             mesh->o_LIFTT, 
                             mesh->o_vmapM, 
                             mesh->o_vmapP, 
                             mesh->o_EToB,
                             currentTime, 
                             mesh->o_x, 
                             mesh->o_y, 
                             cns->mu, 
                             cns->o_rkq, 
                             cns->o_viscousStresses, 
                             cns->o_rhsq);
        }
        
        // update solution using Runge-Kutta
        // rkrhsq_rk = rhsq
        // if rk==6 
        //   q = q + dt*sum_{i=0}^{rk} rkA_{rk,i}*rkrhs_i
        //   rkerr = dt*sum_{i=0}^{rk} rkE_{rk,i}*rkrhs_i
        cns->rkUpdateKernel(mesh->Nelements, 
                          rk,
                          mesh->dt, 
                          o_rkA, 
                          o_rkE, 
                          cns->o_q,
                          cns->o_rhsq, 
                          cns->o_rkrhsq, 
                          cns->o_rkq,
                          cns->o_rkerr);
      }

      //Error estimation 
      //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
      //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
      dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
      cns->rkErrorEstimateKernel(Ntotal, 
                                ATOL,
                                RTOL,
                                cns->o_q,
                                cns->o_rkq,
                                cns->o_rkerr,
                                cns->o_errtmp);
      cns->o_errtmp.copyTo(cns->errtmp);
      dfloat localerr = 0;
      dfloat err = 0;
      for(dlong n=0;n<cns->Nblock;++n){
        localerr += cns->errtmp[n];
      }
      MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      err = sqrt(err/cns->totalElements);

      dfloat fac1 = pow(err,exp1);
      dfloat fac = fac1/pow(facold,beta);

      fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted 
        time += mesh->dt;

        facold = mymax(err,1E-4);

        cns->o_q.copyFrom(cns->o_rkq);

        if(isOutput==1){

	  nextOutputTime += outputInterval;
	  printf("\r");
          cnsReportTri2D(cns, time, options);

        }
	//	printf("\r dt = %g accepted                      ", mesh->dt);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(invfactor1,fac1/safe));
	//	printf("\r dt = %g rejected, trying %g", mesh->dt, dtnew);
	done = 0;
      }
      mesh->dt = dtnew;
      allStep++;

    }

    o_rkA.free();
    o_rkE.free();

    mesh->device.finish();
    
    double elapsed  = timer.toc("Run");

    printf("run took %lg seconds for %d accepted steps and %d total steps\n", elapsed, tstep, allStep);
    
  } else if (strstr(options,"LSERK")) {
    // Low storage explicit Runge Kutta (5 stages, 4th order)
    for(int tstep=0;tstep<mesh->NtimeSteps;++tstep){

      int advSwitch = 1;//(tstep>100);
      
      for(int rk=0;rk<mesh->Nrk;++rk){

        dfloat currentTime = tstep*mesh->dt + mesh->rkc[rk]*mesh->dt;
        
        // extract q halo on DEVICE
        if(mesh->totalHaloPairs>0){
          int Nentries = mesh->Np*cns->Nfields;
          
          mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_q, cns->o_haloBuffer);
          
          // copy extracted halo to HOST 
          cns->o_haloBuffer.copyTo(cns->sendBuffer);      
          
          // start halo exchange
          meshHaloExchangeStart(mesh, mesh->Np*cns->Nfields*sizeof(dfloat), cns->sendBuffer, cns->recvBuffer);
        }

        // now compute viscous stresses
        cns->stressesVolumeKernel(mesh->Nelements, 
                                  mesh->o_vgeo, 
                                  mesh->o_DrT, 
                                  mesh->o_DsT, 
                                  cns->mu, 
                                  cns->o_q, 
                                  cns->o_viscousStresses);

        // wait for q halo data to arrive
        if(mesh->totalHaloPairs>0){
          meshHaloExchangeFinish(mesh);
          
          // copy halo data to DEVICE
          size_t offset = mesh->Np*cns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
          cns->o_q.copyFrom(cns->recvBuffer, cns->haloBytes, offset);
        }
        
        cns->stressesSurfaceKernel(mesh->Nelements, 
                                   mesh->o_sgeo, 
                                   mesh->o_LIFTT,
                                   mesh->o_vmapM, 
                                   mesh->o_vmapP, 
                                   mesh->o_EToB, 
                                   currentTime,
                                   mesh->o_x, 
                                   mesh->o_y, 
                                   cns->mu, 
                                   cns->o_q, 
                                   cns->o_viscousStresses);
        
        // extract stresses halo on DEVICE
        if(mesh->totalHaloPairs>0){
          int Nentries = mesh->Np*cns->Nstresses;
          
          mesh->haloExtractKernel(mesh->totalHaloPairs, Nentries, mesh->o_haloElementList, cns->o_viscousStresses, cns->o_haloStressesBuffer);
          
          // copy extracted halo to HOST 
          cns->o_haloStressesBuffer.copyTo(cns->sendStressesBuffer);      
          
          // start halo exchange
          meshHaloExchangeStart(mesh, mesh->Np*cns->Nstresses*sizeof(dfloat), cns->sendStressesBuffer, cns->recvStressesBuffer);
        }
        
        // compute volume contribution to DG cns RHS
        if (strstr(options,"CUBATURE")) {
          cns->cubatureVolumeKernel(mesh->Nelements, 
                                    advSwitch, 
                                    mesh->o_vgeo, 
                                    mesh->o_cubDrWT,
                                    mesh->o_cubDsWT,
                                    mesh->o_cubInterpT,
                                    cns->o_viscousStresses, 
                                    cns->o_q, 
                                    cns->o_rhsq);
        } else {
          cns->volumeKernel(mesh->Nelements, 
                            advSwitch, 
                            mesh->o_vgeo, 
                            mesh->o_DrT,
                            mesh->o_DsT,
                            cns->o_viscousStresses, 
                            cns->o_q, 
                            cns->o_rhsq);
        }

        // wait for halo stresses data to arrive
        if(mesh->totalHaloPairs>0){
          meshHaloExchangeFinish(mesh);
          
          // copy halo data to DEVICE
          size_t offset = mesh->Np*cns->Nstresses*mesh->Nelements*sizeof(dfloat); // offset for halo data
          cns->o_viscousStresses.copyFrom(cns->recvStressesBuffer, cns->haloStressesBytes, offset);
        }
        
        // compute surface contribution to DG cns RHS (LIFTT ?)
        if (strstr(options,"CUBATURE")) {
          cns->cubatureSurfaceKernel(mesh->Nelements, 
                                     advSwitch, 
                                     mesh->o_sgeo, 
                                     mesh->o_intInterpT,
                                     mesh->o_intLIFTT, 
                                     mesh->o_vmapM, 
                                     mesh->o_vmapP, 
                                     mesh->o_EToB,
                                     currentTime, 
                                     mesh->o_intx, 
                                     mesh->o_inty, 
                                     cns->mu, 
                                     cns->o_q, 
                                     cns->o_viscousStresses, 
                                     cns->o_rhsq);
        } else {
          cns->surfaceKernel(mesh->Nelements, 
                             advSwitch, 
                             mesh->o_sgeo, 
                             mesh->o_LIFTT, 
                             mesh->o_vmapM, 
                             mesh->o_vmapP, 
                             mesh->o_EToB,
                             currentTime, 
                             mesh->o_x, 
                             mesh->o_y, 
                             cns->mu, 
                             cns->o_q, 
                             cns->o_viscousStresses, 
                             cns->o_rhsq);
        }
        
        // update solution using Runge-Kutta
        cns->updateKernel(mesh->Nelements, 
                          mesh->dt, 
                          mesh->rka[rk], 
                          mesh->rkb[rk], 
                          cns->o_rhsq, 
                          cns->o_resq, 
                          cns->o_q);
      }
      
      if(((tstep+1)%mesh->errorStep)==0){
        dfloat time = (tstep+1)*mesh->dt;
        cnsReportTri2D(cns, time, options);
      }
    }
  }

}
