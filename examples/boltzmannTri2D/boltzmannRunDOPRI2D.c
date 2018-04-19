#include "boltzmann2D.h"

void boltzmannRunDOPRI2D(bns_t *bns, int haloBytes, dfloat * sendBuffer, dfloat *recvBuffer, char * options){

	mesh2D *mesh = bns->mesh;

	dfloat dtMIN = 1E-7; //minumum allowed timestep
	dfloat ATOL  = 1E-6;  //absolute error tolerance
	dfloat RTOL  = 1E-4;  //relative error tolerance
	dfloat safe  = 0.9;   //safety factor

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

	boltzmannReportAddaptive2D(bns, time, options);

	int done =0;
	while (!done) {

		int advSwitch = 1;

		if (bns->dt<dtMIN){
			printf("ERROR: Time step became too small at time step=%d\n", tstep);
			exit (-1);
		}
		if (isnan(bns->dt)) {
			printf("ERROR: Solution became unstable at time step=%d\n", tstep);
			exit (-1);
		}

		// check for next output
		int isOutput = 0;
		if((time+bns->dt > nextOutputTime) && (time<=nextOutputTime)){
			isOutput = 1;
			bns->dt = nextOutputTime-time;
		}

		//check for final timestep
		if (time+bns->dt > bns->finalTime){
			bns->dt = bns->finalTime-time;
			done = 1;
			isOutput = 0;
		}

		// Perform an Full RK Step
		boltzmannDOPRIStep2D(bns,time, haloBytes,sendBuffer, recvBuffer, options);


		//Error estimation 
		//E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
		//      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
		int Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;
		bns->errorEstimateKernel(Ntotal, 
		                    ATOL,
		                    RTOL,
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

		// printf("\r err = %g ", err);
		 err = sqrt(err/(bns->totalElements*mesh->Np));

		dfloat fac1 = pow(err,exp1);
		dfloat fac = fac1/pow(facold,beta);

		fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
		dfloat dtnew = bns->dt/fac;

		printf("\r err = %g ", err);

		if (err<1.0) { //dt is accepted 
			time += bns->dt;

			facold = mymax(err,1E-4);

			bns->o_q.copyFrom(bns->o_rkq);
			bns->o_pmlqx.copyFrom(bns->o_rkqx);
			bns->o_pmlqy.copyFrom(bns->o_rkqy);

			if(isOutput==1){

				nextOutputTime += outputInterval;
				// printf("\r");
				boltzmannReportAddaptive2D(bns, time, options);

			}
		printf("\r dt = %g  %g  %d                       ", bns->dt, time, tstep+1);
		tstep++;
		} else {
		dtnew = bns->dt/(mymax(invfactor1,fac1/safe));
		//	printf("\r dt = %g rejected, trying %g", mesh->dt, dtnew);
		done = 0;
		}
		bns->dt = dtnew;
		allStep++;

    char fname[BUFSIZ]; sprintf(fname, "boltzmannAddaptiveDt2D.dat");
    FILE *fp; fp = fopen(fname, "a");
    fprintf(fp, "%.5e %.5e\n", time, bns->dt); 
    fclose(fp);


	}



}