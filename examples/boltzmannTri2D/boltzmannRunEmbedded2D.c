#include "boltzmann2D.h"

void boltzmannRunEmbedded2D(bns_t *bns, int haloBytes, dfloat * sendBuffer, dfloat *recvBuffer, char * options){

	mesh2D *mesh = bns->mesh;

	
	// #if PID==1 // PID control
	// printf("PID setup\n");
	// dfloat dtMIN = 1E-10; //minumum allowed timestep
	// dfloat ATOL  = 1E-5;  //absolute error tolerance
	// dfloat RTOL  = 1E-4;  //relative error tolerance
	// dfloat safe  = 0.8;   //safety factor

	// //error control parameters
	// dfloat beta    = 0.05; 
	// dfloat factor1 = 0.2;
	// dfloat factor2 = 5.0;


	// dfloat exp1       = 0.25 - 0.75*beta;
	// dfloat invfactor1 = 1.0/factor1;
	// dfloat invfactor2 = 1.0/factor2;
	// dfloat facold     = 1E-4;
	
	// #else // PI control
	// dfloat dtMIN = 1E-8; //minumum allowed timestep
	// dfloat ATOL  = 1E-5;  //absolute error tolerance
	// dfloat RTOL  = 1E-4;  //relative error tolerance
	// dfloat safe  = 0.9; //pow(0.25,0.25);   //safety factor

	// //error control parameters
	// dfloat facmin = 0.5;
	// dfloat facmax = 2.0;


	// dfloat exp1 = 0.2; 
	// dfloat fmin = 0.5; 
	// dfloat fmax = 2.0;
	// #endif


	// hard code this for the moment
	dfloat outputInterval = 0.5;
	dfloat nextOutputTime = outputInterval;
	dfloat outputNumber = 0;

	//initial time
	bns->time = 0.0;
	bns->tstep = 0;
	bns->atstep = 0; 
	bns->rtstep = 0;  

	// int tstep=0, allStep = 0;

	boltzmannReportAddaptive2D(bns, bns->time, options);

	int done =0;
	while (!done) {

		int advSwitch = 1;

		if (bns->dt<1E-10){
			printf("ERROR: Time step became too small at time step=%d\n", bns->tstep);
			exit (-1);
		}
		if (isnan(bns->dt)) {
			printf("ERROR: Solution became unstable at time step=%d\n", bns->tstep);
			exit (-1);
		}

		// // check for next output
		int isOutput = 0;
		// if((time+bns->dt > nextOutputTime) && (time<=nextOutputTime)){
		// 	isOutput = 1;
		// 	bns->dt = nextOutputTime-time;
		// }

		//check for final timestep
		if (bns->time+bns->dt > bns->finalTime){
			bns->dt = bns->finalTime-bns->time;
			done = 1;
			isOutput = 0;
		}

    
        if(strstr(options,"IMEXRK")) // Kennedy-Carpanter Additve RK
			boltzmannIMEXRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

	    if(strstr(options,"SAADRK"))  // SA Adaptive RK 
	    	boltzmannSAADRKStep2D(bns, bns->time, haloBytes, sendBuffer, recvBuffer, options);

	    if(strstr(options,"DOPRI5"))  // DORMAND-PRINCE explicit
	    	boltzmannDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);

	    if(strstr(options,"XDOPRI")) // Integrating Factor DOPRI
			boltzmannXDOPRIStep2D(bns,bns->time, haloBytes,sendBuffer, recvBuffer, options);





		// Error Control
        boltzmannErrorControl2D(bns, options);

        
        printf("\r Time: %.4e dt = %.4e Aiter= %d Riter = %d iter %d ", bns->time, bns->dt, bns->atstep, bns->rtstep, bns->tstep); fflush(stdout);


		// //Error estimation 
		// //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
		// //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
		// int Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;
		// bns->errorEstimateKernel(Ntotal, 
		//                     ATOL,
		//                     RTOL,
		//                     bns->o_q,
		//                     bns->o_rkq,
		//                     bns->o_rkerr,
		//                     bns->o_errtmp);

		// bns->o_errtmp.copyTo(bns->errtmp);
		// dfloat localerr = 0;
		// dfloat err = 0;
		// for(int n=0;n<bns->Nblock;++n){
		// 	localerr += bns->errtmp[n];
		// }
		// MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

		//  err = sqrt(err/(bns->totalElements*mesh->Np));


// #if PID==1
// 		dfloat fac1 = pow(err,exp1);
// 		dfloat fac  = fac1/pow(facold,beta);

// 		fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
// 		dfloat dtnew = bns->dt/fac;

// 		printf("\r err = %g ", err);

// 		if (err<1.0) { //dt is accepted 
// 			time += bns->dt;

// 			facold = mymax(err,1E-4);

// 			bns->o_q.copyFrom(bns->o_rkq);
// 			bns->o_pmlqx.copyFrom(bns->o_rkqx);
// 			bns->o_pmlqy.copyFrom(bns->o_rkqy);
// 		if(isOutput==1){

// 			nextOutputTime += outputInterval;
// 			// printf("\r");
// 			boltzmannReportAddaptive2D(bns, time, options);

// 		}
// 			printf("\r dt = %g  time: %g  Nstep:%d Nrejected: %d", bns->dt, time, tstep, allStep-tstep);
// 			tstep++;
// 		} else {
// 			dtnew = bns->dt/(mymax(invfactor1,fac1/safe));

// 		printf(" \r dt = %g rejected, trying %g at step %d  with rejected: %d", bns->dt, dtnew, tstep, allStep-tstep );
// 			done = 0;
// 		}
// 		bns->dt = dtnew;
// 		allStep++;
// #else
// 		// printf("\r err = %g ", err);
		
// 		dfloat fac1 = pow(1.0/err,exp1);

// 		dfloat fac = mymin(facmax, mymax(facmin,safe*fac1));

// 		dfloat dtnew = bns->dt*fac;


// 		printf("\r err = %g ", err);

// 		if (err<1.0) { //dt is accepted 
// 			time += bns->dt;
// 			bns->o_q.copyFrom(bns->o_rkq);
// 			bns->o_pmlqx.copyFrom(bns->o_rkqx);
// 			bns->o_pmlqy.copyFrom(bns->o_rkqy);
// 			if(isOutput==1){

// 				nextOutputTime += outputInterval;
// 				// printf("\r");
// 				boltzmannReportAddaptive2D(bns, time, options);

// 			}
// 		printf("\n dt = %g  time: %g  Nstep:%d Nrejected: %d", bns->dt, time, tstep, allStep-tstep);
// 		tstep++;
// 		} else {

// 		dtnew = bns->dt*mymin(1.0, mymax(facmin,safe*fac1));
		
// 		printf(" \r dt = %g rejected, trying %g at step %d  with rejected: %d", bns->dt, dtnew, tstep, allStep-tstep );
// 		done = 0;
// 		}

// 		bns->dt = dtnew;
// 		allStep++;
    
// #endif
		
    if(strstr(options,"SAADRK"))
		boltzmannSAADRKCoefficients(bns, options);

    char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt2D.dat");
    FILE *fp; fp = fopen(fname, "a");
    fprintf(fp, "%.5e %.5e\n", bns->time, bns->dt); 
    fclose(fp);


	}



}