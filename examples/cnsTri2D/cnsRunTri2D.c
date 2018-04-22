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

    // hard code this for the moment
    dfloat outputInterval = .5;
    dfloat nextOutputTime = outputInterval;
    dfloat outputNumber = 0;
    
    //initial time
    dfloat time = 0.0;
    int tstep=0, allStep = 0;

    int done =0;
    while (!done) {

      cns->advSwitch = 1;
      
      if (mesh->dt<cns->dtMIN){
        printf("ERROR: Time step became too small at time step=%d\n", tstep);
        exit (-1);
      }
      if (isnan(mesh->dt)) {
        printf("ERROR: Solution became unstable at time step=%d\n", tstep);
        exit (-1);
      }

      //check for final timestep
      if (time+mesh->dt > mesh->finalTime){
	mesh->dt = mesh->finalTime-time;
	done = 1;
      }

      // try a step with the current time step
      cnsStepTri2D(cns, options, time);
      
      //Error estimation 
      //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
      //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
      dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
      cns->rkErrorEstimateKernel(Ntotal, 
				 cns->ATOL,
				 cns->RTOL,
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

      dfloat fac1 = pow(err,cns->exp1);
      dfloat fac = fac1/pow(cns->facold,cns->beta);

      fac = mymax(cns->invfactor2, mymin(cns->invfactor1,fac/cns->safe));
      dfloat dtnew = mesh->dt/fac;

      if (err<1.0) { //dt is accepted 
        time += mesh->dt;

        cns->facold = mymax(err,1E-4); // hard coded factor ?

        cns->o_q.copyFrom(cns->o_rkq);

	if(time-mesh->dt<nextOutputTime && time>nextOutputTime){
	  cnsReportTri2D(cns, time, options);
	  nextOutputTime += outputInterval;
	}
	
	printf("\r time = %g, dt = %g accepted                      ", time, mesh->dt);
        tstep++;
      } else {
        dtnew = mesh->dt/(mymax(cns->invfactor1,fac1/cns->safe));
	printf("\r time = %g, dt = %g rejected, trying %g", time, mesh->dt, dtnew);

	done = 0;
      }
      mesh->dt = dtnew;
      allStep++;

    }

#if 0
    cns->o_rkA.free();
    cns->o_rkE.free();
#endif
    
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
