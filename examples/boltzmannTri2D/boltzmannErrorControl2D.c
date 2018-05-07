#include "boltzmann2D.h"

void boltzmannErrorControl2D(bns_t *bns, setupAide &options){

	// mesh2D *mesh = bns->mesh;

	// //E. HAIRER, S.P. NORSETT AND G. WANNER, SOLVING ORDINARY
	// //      DIFFERENTIAL EQUATIONS I. NONSTIFF PROBLEMS. 2ND EDITION.
	// int Ntotal = mesh->Nelements*mesh->Np*bns->Nfields;

	// bns->errorEstimateKernel(Ntotal, 
	// 					        bns->ATOL,
	// 					        bns->RTOL,
	// 					        bns->o_q,
	// 					        bns->o_rkq,
	// 					        bns->o_rkerr,
	// 					        bns->o_errtmp);

	// bns->o_errtmp.copyTo(bns->errtmp);
	// dfloat localerr = 0;
	// dfloat err = 0;
	// for(int n=0;n<bns->Nblock;++n){
	// localerr += bns->errtmp[n];
	// }
	// MPI_Allreduce(&localerr, &err, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
    
 //    // root mean square of step error. i.e. (Q-Qem)/(ATOL + RTOL*Qmax);
	// err = sqrt(err/(bns->totalElements*mesh->Np));

 //  bns->tstepAccepted = (err<1.0) ? 1 : 0; 
  
 //      bns->ehist[0] = err;
 //  if(tstepAccepted){
 //    // update error history
 //    bns->ehist[2] = bns->ehist[1];
 //    bns->ehist[1] = bns->ehist[0];
 //  }


    
 //    // if(err<1.0){
 //    // 	bns->tstepAccepted = 1; 
 //    // 	bns->tstep++;
 //    // 	bns->time += bns->dt; 

 //    // }else{
 //    // 	bns->eflag = 0;
 //    //     bns->rtstep++;
 //    // }

 //    // bns->atstep++;
 


 //    if(bns->emethod ==0){ // Default PI controller with DORPI
	// 		dfloat safe       = 0.8;   //safety factor

	// 		//error control parameters
	// 		dfloat beta       = 0.05; 
	// 		dfloat factor1    = 0.2;
	// 		dfloat factor2    = 5.0;
	// 		dfloat exp1       = 0.25 - 0.75*beta;
	// 		dfloat invfactor1 = 1.0/factor1;
	// 		dfloat invfactor2 = 1.0/factor2;

	// 		dfloat fac1 = pow(bns->ehist[0],exp1);

	// 		dfloat facold  = mymax(bns->ehist[1],1E-4);
	// 		dfloat fac     = fac1/pow(facold,beta);

	// 		fac = mymax(invfactor2, mymin(invfactor1,fac/safe));
	// 		dfloat dtnew = bns->dt/fac;
			
 //      if(bns->tstepAccepted){

 //      	if(bns->reportFlag){


 //      	}








 //      }
			








	  
	//  }




	// 	// Use error control
	// 	if(bns->emethod ==1){ // PID control

 //    // coefficients are taken from ArkOde
	// 	dfloat safety    = 0.993;
	// 	dfloat facmin    = 0.2;
	// 	dfloat facmax    = 5.0;
	// 	dfloat k1        = -0.58/bns->rkp;
	// 	dfloat k2        = 0.21/bns->rkp;
	// 	dfloat k3        = -0.10/bns->rkp;

	// 	if(bns->eflag){
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1)*pow(bns->ehist[1],k2)*pow(bns->ehist[2],k3);
			
	// 		// check scaling limits
	// 		fac       = mymin(facmax, mymax(facmin,fac));

	// 		// printf("fac= %.4e err=%.4e  dt: %.4e\n", fac, err, fac*bns->dt);

	// 		dfloat dt = fac*bns->dt;

	// 		// will be changed
	// 		bns->dt  *=fac;  

	// 	}else{
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1)*pow(bns->ehist[1],k2)*pow(bns->ehist[2],k3);
	// 		fac           = mymin(1.0, mymax(facmin,fac));
	// 		// will be changed
	// 		bns->dt  *=fac;  
	// 	} 
	// }


	// if(bns->emethod ==2){ // PI control

 //    // coefficients are taken from ArkOde
	// 	dfloat safety    = 0.93;
	// 	dfloat facmin    = 0.2;
	// 	dfloat facmax    = 5.0;
	// 	dfloat k1        = -0.8/bns->rkp;
	// 	dfloat k2        = 0.35/bns->rkp;
		
	// 	if(bns->eflag){
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1)*pow(bns->ehist[1],k2);
			
	// 		// check scaling limits
	// 		fac       = mymin(facmax, mymax(facmin,fac));

	// 		printf("fac= %.4e err=%.4e  dt: %.4e\n", fac, err, fac*bns->dt);

	// 		dfloat dt = fac*bns->dt;

	// 		// will be changed
	// 		bns->dt  *=fac;  

	// 	}else{
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1)*pow(bns->ehist[1],k2);
	// 		fac           = mymin(1.0, mymax(facmin,fac));
	// 		// will be changed
	// 		bns->dt  *=fac;  
	// 	} 
	// }


	// if(bns->emethod ==3){ // I control

 //    // coefficients are taken from ArkOde
	// 	dfloat safety    = 0.93;
	// 	dfloat facmin    = 0.2;
	// 	dfloat facmax    = 5.0;
	// 	dfloat k1        = -1.0/bns->rkp;
				
	// 	if(bns->eflag){
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1);
			
	// 		// check scaling limits
	// 		fac       = mymin(facmax, mymax(facmin,fac));

	// 		// printf("fac= %.4e err=%.4e  dt: %.4e\n", fac, err, fac*bns->dt);

	// 		dfloat dt = fac*bns->dt;

	// 		// will be changed
	// 		bns->dt  *=fac;  

	// 	}else{
	// 		dfloat fac   = safety*pow(bns->ehist[0],k1);
	// 		fac           = mymin(1.0, mymax(facmin,fac));
	// 		// will be changed
	// 		bns->dt  *=fac;  
	// 	} 
	// }








    



 //    if(bns->eflag){
	// 		bns->o_q.copyFrom(bns->o_rkq);
	// 		if(options.compareArgs("ABSORBING LAYER","PML")){
	// 		bns->o_pmlqx.copyFrom(bns->o_rkqx);
	// 		bns->o_pmlqy.copyFrom(bns->o_rkqy);
	// 	}

 //    }



}