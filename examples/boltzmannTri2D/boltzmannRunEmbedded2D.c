#include "boltzmann2D.h"

#define PID 1

void boltzmannRunEmbedded2D(bns_t *bns, int haloBytes, dfloat * sendBuffer,
	                                     dfloat *recvBuffer, setupAide &options){

	mesh2D *mesh = bns->mesh;

	
#if PID==1 // PID control
	printf("PID setup\n");
	dfloat safe  = 0.8;   //safety factor

	//error control parameters
	dfloat beta    = 0.05; 
	dfloat factor1 = 0.2;
	dfloat factor2 = 5.0;
	dfloat exp1       = 0.25 - 0.75*beta;
	dfloat invfactor1 = 1.0/factor1;
	dfloat invfactor2 = 1.0/factor2;
	dfloat facold     = 1E-4;
	
#else // PI control
	// dfloat dtMIN = 1E-8; //minumum allowed timestep
	// dfloat ATOL  = 1E-5;  //absolute error tolerance
	// dfloat RTOL  = 1E-4;  //relative error tolerance
	dfloat safe  = 0.9; //pow(0.25,0.25);   //safety factor

	//error control parameters
	dfloat facmin = 0.5;
	dfloat facmax = 2.0;
	dfloat exp1   = 0.2; 
#endif


	dfloat hmin = 1e9;
	for(dlong e=0;e<mesh->Nelements;++e){  

		for(int f=0;f<mesh->Nfaces;++f){
			dlong sid = mesh->Nsgeo*(mesh->Nfaces*e + f);
			dfloat sJ   = mesh->sgeo[sid + SJID];
			dfloat invJ = mesh->sgeo[sid + IJID];
			dfloat hest = .5/(sJ*invJ);
			hmin = mymin(hmin, hest);
		}
	}


	// hard code this for the moment
	dfloat outputInterval;
	options.getArgs("OUTPUT INTERVAL", outputInterval);

	dfloat nextOutputTime = outputInterval;
	dfloat outputNumber = 0;

	//initial time
	bns->time = 0.0;
	bns->tstep = 0;
	bns->atstep = 0; 
	bns->rtstep = 0;  

    if(bns->reportFlag)
	boltzmannReportAddaptive2D(bns, bns->time, options);

	int done =0;
	while (!done) {

		int advSwitch = 1;

		if (bns->dt<bns->dtMIN){
			printf("ERROR: Time step became too small at time step=%d\n", bns->tstep);
			exit (-1);
		}
		if (isnan(bns->dt)) {
			printf("ERROR: Solution became unstable at time step=%d\n", bns->tstep);
			exit (-1);
		}

		
		//check for final timestep
		if (bns->time+bns->dt > bns->finalTime){
			bns->dt = bns->finalTime-bns->time;
			done = 1;
		}

    
		if(options.compareArgs("TIME INTEGRATOR","IMEXRK")) // Kennedy-Carpanter Additve RK
		boltzmannIMEXRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

		if(options.compareArgs("TIME INTEGRATOR","SAADRK"))  // SA Adaptive RK 
		boltzmannSAADRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

		if(options.compareArgs("TIME INTEGRATOR","DOPRI5"))  // DORMAND-PRINCE explicit
		boltzmannDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

		if(options.compareArgs("TIME INTEGRATOR","XDOPRI")) // Integrating Factor DOPRI
		boltzmannXDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

		// // Error Control Nor completed yet !!!!!
        // boltzmannErrorControl2D(bns, options);
        // printf("\r Time: %.4e dt = %.4e Aiter= %d Riter = %d iter %d ", bns->time, bns->dt, bns->atstep, bns->rtstep, bns->tstep); fflush(stdout);


		//Error estimation 
		//E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
		//      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
		int Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;
		bns->errorEstimateKernel(Ntotal, 
		                    bns->ATOL,
		                    bns->RTOL,
		                    bns->o_q,
		                    bns->o_rkq,
		                    bns->o_rkerr,
		                    bns->o_errtmp);

		bns->o_errtmp.copyTo(bns->errtmp);
		dfloat localerr = 0;
		dfloat err = 0;
		for(int n=0;n<bns->Nblock;++n){
			localerr += bns->errtmp[n];
		}
		MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

		err = sqrt(err/(bns->totalElements*mesh->Np));


#if PID==1
		dfloat fac1 = pow(err,exp1);
		dfloat fac  = fac1/pow(facold,beta);

		fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
		dfloat dtnew = bns->dt/fac;

		if(err<1.0){
    
			if(bns->reportFlag){
			// check for output during this step and do a mini-step
				if(bns->time<nextOutputTime && bns->time+bns->dt>nextOutputTime){
					dfloat savedt = bns->dt;			  
					// save rkq
					bns->o_saveq.copyFrom(bns->o_rkq);

					// change dt to match output
					bns->dt = nextOutputTime-bns->time;

					// print
					printf("Taking output mini step: %g\n", bns->dt);

					if(options.compareArgs("TIME INTEGRATOR","IMEXRK")) // Kennedy-Carpanter Additve RK
					boltzmannIMEXRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","SAADRK"))  // SA Adaptive RK 
					boltzmannSAADRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","DOPRI5"))  // DORMAND-PRINCE explicit
					boltzmannDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","XDOPRI")) // Integrating Factor DOPRI
					boltzmannXDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

					// shift for output
					bns->o_rkq.copyTo(bns->o_q);

					// output  (print from rkq)
					boltzmannReportAddaptive2D(bns, nextOutputTime, options);

					// restore time step
					bns->dt = savedt;

					// increment next output time
					nextOutputTime += outputInterval;

					// accept saved rkq
					bns->o_q.copyFrom(bns->o_saveq);
				}
			}
		  

			bns->o_q.copyFrom(bns->o_rkq);
			bns->o_pmlqx.copyFrom(bns->o_rkqx);
			bns->o_pmlqy.copyFrom(bns->o_rkqy);

			facold = mymax(err,1E-4);
			bns->time += bns->dt;

			printf("\r time = %g (%d), dt = %g accepted (ratio dt/hmin = %g)               ", 
				         bns->time, bns->atstep, bns->dt, bns->dt/hmin);
			bns->tstep++;
		}else{
        bns->rtstep++; 
				dtnew = bns->dt/(mymax(invfactor1,fac1/safe));
				printf("\r time = %g (%d), dt = %g rejected (ratio dt/min = %g), trying %g", bns->time,bns->atstep, bns->dt, bns->dt/hmin, dtnew);
				done =0;
		}

		bns->dt = dtnew;
		bns->atstep++;


#else

		dfloat fac1 = pow(1.0/err,exp1);
		dfloat fac = mymin(facmax, mymax(facmin,safe*fac1));
		dfloat dtnew = bns->dt*fac;

		if(err<1.0){
    
			if(bns->reportFlag){
			// check for output during this step and do a mini-step
				if(bns->time<nextOutputTime && bns->time+bns->dt>nextOutputTime){
					dfloat savedt = bns->dt;			  
					// save rkq
					bns->o_saveq.copyFrom(bns->o_rkq);

					// change dt to match output
					bns->dt = nextOutputTime-bns->time;

					// print
					printf("Taking output mini step: %g\n", bns->dt);

					if(options.compareArgs("TIME INTEGRATOR","IMEXRK")) // Kennedy-Carpanter Additve RK
					boltzmannIMEXRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","SAADRK"))  // SA Adaptive RK 
					boltzmannSAADRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","DOPRI5"))  // DORMAND-PRINCE explicit
					boltzmannDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

					if(options.compareArgs("TIME INTEGRATOR","XDOPRI")) // Integrating Factor DOPRI
					boltzmannXDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

					// shift for output
					bns->o_rkq.copyTo(bns->o_q);

					// output  (print from rkq)
					boltzmannReportAddaptive2D(bns, nextOutputTime, options);

					// restore time step
					bns->dt = savedt;

					// increment next output time
					nextOutputTime += outputInterval;

					// accept saved rkq
					bns->o_q.copyFrom(bns->o_saveq);
				}
			}
		  

			bns->o_q.copyFrom(bns->o_rkq);
			bns->o_pmlqx.copyFrom(bns->o_rkqx);
			bns->o_pmlqy.copyFrom(bns->o_rkqy);

			bns->time += bns->dt;

			printf("\r time = %g (%d), dt = %g accepted (ratio dt/hmin = %g)               ", 
				         bns->time, bns->atstep, bns->dt, bns->dt/hmin);
			bns->tstep++;
		}else{

				dtnew = bns->dt*mymin(1.0, mymax(facmin,safe*fac1));
				printf("\r time = %g (%d), dt = %g rejected (ratio dt/min = %g), trying %g", 
					         bns->time,bns->atstep, bns->dt, bns->dt/hmin, dtnew);
				done =0;
		}

		bns->dt = dtnew;
		bns->atstep++;    
#endif
		
    if(options.compareArgs("TIME INTEGRATOR","SAADRK"))
		boltzmannSAADRKCoefficients(bns, options);

    char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt2D.dat");
    FILE *fp; fp = fopen(fname, "a");
    fprintf(fp, "%.5e %.5e\n", bns->time, bns->dt); 
    fclose(fp);


	}



}