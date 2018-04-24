#include "boltzmann2D.h"

void boltzmannTimeStepperCoefficients(bns_t *bns, char *options){


mesh2D *mesh = bns->mesh; 

if(strstr(options,"MRSAAB") || strstr(options,"MRAB") || 
   strstr(options,"SRAB")   || strstr(options,"SAAB") ){

	int Nlevels = 0;

	if(mesh->MRABNlevels==0)
		Nlevels =1;
	else
		Nlevels = mesh->MRABNlevels;

	// Create circle on complex plane
	const int Nr = 32; 
	dfloat complex R[Nr]; 
	for(int ind =1; ind <= Nr; ++ind){
		const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr; 
		R[ind-1] = cexp(I*M_PI* theta);
	}


	bns->MRSAAB_A = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
	bns->MRSAAB_B = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
	bns->MRSAAB_C = (dfloat *) calloc(    Nlevels,sizeof(dfloat));
	bns->MRAB_A   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
	bns->MRAB_B   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
	bns->MRAB_C   = (dfloat *) calloc(    Nlevels,sizeof(dfloat));

	int MRABorder = 3; 

  for(int l = 0; l<Nlevels; ++l){
		// MRSAAB coefficients
		dfloat alpha = -bns->tauInv*bns->dt*pow(2,l);
		//dfloat alpha=0.0;
		dfloat h  = bns->dt * pow(2,l); 
    //
		for (int order=0; order<3; ++order){
      // computation of coefficients based on magnitude
      const int id = order*Nlevels*3 + l*3;
      if(order==0){

				double complex a1 = 0. + 0.* I; 
				double complex b1 = 0. + 0.* I; 

				for(int i = 0; i<Nr; ++i ){
					double complex lr = alpha  + R[i];
					a1 +=  h*(cexp(lr) - 1.)/lr;
					b1 +=  h*(cexp(lr/2.) - 1.)/lr;
				}
				// Full dt coeeficients
				bns->MRSAAB_A[id + 0] = creal(a1)/Nr;
				bns->MRSAAB_A[id + 1] = 0.f;
				bns->MRSAAB_A[id + 2] = 0.f;
				// Half coefficients
				bns->MRSAAB_B[id + 0] = creal(b1)/Nr;
				bns->MRSAAB_B[id + 1] = 0.f;
				bns->MRSAAB_B[id + 2] = 0.f;

				// MRAB coefficients
				bns->MRAB_A[id + 0]   =  h ;
				bns->MRAB_A[id + 1]   =  0.f ;
				bns->MRAB_A[id + 2]   =  0.f ;

				bns->MRAB_B[id+0]     =  h/2. ;
				bns->MRAB_B[id+1]     =  0.f ;
				bns->MRAB_B[id+2]     =  0.f ;
			}

			else if(order==1){

				double complex a1 = 0. + 0.* I; 
				double complex b1 = 0. + 0.* I; 
				double complex a2 = 0. + 0.* I; 
				double complex b2 = 0. + 0.* I; 

				for(int i = 0; i<Nr; ++i ){
					double complex lr = alpha  + R[i];
					a1 +=  h*(-2.*lr + (1.+lr)*cexp(lr) - 1.)/cpow(lr,2);
					a2 +=  h*(lr - cexp(lr) + 1.)/cpow(lr,2);
					b1 +=  h*(-1.5*lr + (1.+lr)*cexp(lr/2.) - 1.)/cpow(lr,2);
					b2 +=  h*(0.5*lr - cexp(lr/2.) + 1.)/cpow(lr,2);
				}
				// Full dt coeeficients
				bns->MRSAAB_A[id + 0] = creal(a1)/Nr;
				bns->MRSAAB_A[id + 1] = creal(a2)/Nr;
				bns->MRSAAB_A[id + 2] = 0.f;
				// Half coefficients
				bns->MRSAAB_B[id + 0] = creal(b1)/Nr;
				bns->MRSAAB_B[id + 1] = creal(b2)/Nr;
				bns->MRSAAB_B[id + 2] = 0.f;


				// MRAB coefficients
				bns->MRAB_A[id + 0]   =  3.*h/2. ;
				bns->MRAB_A[id + 1]   = -1.*h/2. ;
				bns->MRAB_A[id + 2]   =  0.f ;

				bns->MRAB_B[id + 0]   =  5.*h/8. ;
				bns->MRAB_B[id + 1]   = -1.*h/8. ;
				bns->MRAB_B[id + 2]   =   0.f ;
      }

			else{
				double complex a1 = 0. + 0.* I; 
				double complex b1 = 0. + 0.* I; 
				double complex a2 = 0. + 0.* I; 
				double complex b2 = 0. + 0.* I; 
				double complex a3 = 0. + 0.* I; 
				double complex b3 = 0. + 0.* I; 

				for(int i = 0; i<Nr; ++i ){
					double complex lr = alpha  + R[i];
					a1 += h*(-2.5*lr - 3.*cpow(lr,2) + (1.+cpow(lr,2)+1.5*lr)*cexp(lr) - 1.)/cpow(lr,3);
					a2 += h*(4.*lr + 3.*cpow(lr,2)- (2.*lr + 2.0)*cexp(lr) + 2.)/cpow(lr,3);
					a3 +=-h*(1.5*lr + cpow(lr,2)- (0.5*lr + 1.)*cexp(lr) + 1.)/cpow(lr,3);
					b1 += h*(cexp(lr/2.)- 2.*lr - (15.*cpow(lr,2))/8. + (cpow(lr,2) + 1.5*lr)*cexp(lr/2.) - 1.)/cpow(lr,3);
					b2 += h*(3.*lr - 2.*cexp(lr/2.0) + 1.25*cpow(lr,2) - 2.*lr*cexp(lr/2.) + 2.)/cpow(lr,3);
					b3 +=-h*(lr - cexp(lr/2.) + 0.375*cpow(lr,2) - 0.5*lr*cexp(lr/2.) + 1.)/cpow(lr,3);
				}

				// Full dt coeeficients
				bns->MRSAAB_A[id+0] = creal(a1)/Nr;
				bns->MRSAAB_A[id+1] = creal(a2)/Nr;
				bns->MRSAAB_A[id+2] = creal(a3)/Nr;
				// Half coefficients
				bns->MRSAAB_B[id+0] = creal(b1)/Nr;
				bns->MRSAAB_B[id+1] = creal(b2)/Nr;
				bns->MRSAAB_B[id+2] = creal(b3)/Nr;

				// MRAB coefficients
				bns->MRAB_A[id+0]   =  23.*h/12. ;
				bns->MRAB_A[id+1]   = -16.*h/12. ;
				bns->MRAB_A[id+2]   =  5. *h/12. ;

				bns->MRAB_B[id+0]   =  17.*h/24. ;
				bns->MRAB_B[id+1]   = - 7.*h/24. ;
				bns->MRAB_B[id+2]   =   2.*h/24. ;  
		}
  }

    // Exponential part
    bns->MRSAAB_C[l]    = exp(alpha);
    //printf("Coefficient: %.5f \n", bns->MRSAAB_C[l]);
    bns->MRAB_C[l]      =   h ;
 }

}


if(strstr(options, "LSERK")){
// 
printf("Relying the fact that  LSERK set on mesh structure\n");


}

if(strstr(options,"DOPRI5") || strstr(options,"XDOPRI")){

  dfloat rkC[bns->NrkStages]          = {0.0, 0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
  dfloat rkA[bns->NrkStages*bns->NrkStages]   ={             0.0,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                                   0.2,             0.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                              3.0/40.0,        9.0/40.0,            0.0,          0.0,             0.0,       0.0, 0.0,
                             44.0/45.0,      -56.0/15.0,       32.0/9.0,          0.0,             0.0,       0.0, 0.0,
                        19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0,             0.0,       0.0, 0.0,
                         9017.0/3168.0,     -355.0/33.0, 46732.0/5247.0,   49.0/176.0, -5103.0/18656.0,       0.0, 0.0, 
                            35.0/384.0,             0.0,   500.0/1113.0,  125.0/192.0,  -2187.0/6784.0, 11.0/84.0, 0.0 };
  dfloat rkE[bns->NrkStages]= {71.0/57600.0,  0.0, -71.0/16695.0, 71.0/1920.0, -17253.0/339200.0, 22.0/525.0, -1.0/40.0 }; 


  memcpy(bns->rkC, rkC, bns->NrkStages*sizeof(dfloat)); 
  memcpy(bns->rkA, rkA, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
  memcpy(bns->rkE, rkE, bns->NrkStages*sizeof(dfloat));
}


if(strstr(options,"SAADRK")){
	if(bns->NrkStages==5){ // SAARK43

		// first set non-semianalytic part of the integrator 
		dfloat rkC[bns->NrkStages]  = {0.0, 0.5, 0.5, 1.0, 1.0};
		dfloat rkA[bns->NrkStages*bns->NrkStages]
		                ={             0.0,             0.0,            0.0,          0.0,             0.0,
		                               0.5,             0.0,            0.0,          0.0,             0.0,
		                               0.0,             0.5,            0.0,          0.0,             0.0,
		                               0.0,             0.0,            1.0,          0.0,             0.0,
		                             1.0/6.0,         1.0/3.0,        1.0/3.0,       1.0/6.0,          0.0}; 
		dfloat rkE[bns->NrkStages]= {  0.0,             0.0,            0.0,        -1.0/6.0,        1.0/6.0}; 


		memcpy(bns->rkC, rkC, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkA, rkA, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
		memcpy(bns->rkE, rkE, bns->NrkStages*sizeof(dfloat));
    // Compute semi-analytic part of the integrator

    boltzmannSAADRKCoefficients(bns, options);


	}





}

if(strstr(options, "SARK")){
	//
	for(int i=0; i<5; i++){
		for(int j=0; j<5; j++){
			bns->SARK_A[i][j] = 0.0; bns->RK_A[i][j]   = 0.0; 
		}
		bns->SARK_B[i] = 0.0; bns->SARK_C[i] = 0.0;
		bns->RK_B[i]   = 0.0; bns->RK_C[i]   = 0.0;
	}


	dfloat a21 = 1.f/2.f;   dfloat a31 = -1.f ;    dfloat a32 = 2.f;
	dfloat b1 = 1.f/6.f;    dfloat b2 = 2./3.;       dfloat b3 = 1./6.; 
	dfloat c1 = 0.f;       dfloat c2 = 1./2.;        dfloat c3 = 1.; 

	// Base Method
	bns->RK_A[0][0] = 0.;  bns->RK_A[1][0] = a21;  bns->RK_A[2][0] = a31;   bns->RK_A[2][1] = a32; 
	bns->RK_B[0] = b1;     bns->RK_B[1] = b2;      bns->RK_B[2] = b3; 
	bns->RK_C[0] = c1;     bns->RK_C[1] = c2;      bns->RK_C[2] = c3; 

	dfloat alpha = -bns->tauInv*bns->dt, coef = -bns->tauInv, h   = bns->dt; 

	printf("alpha = %.16e\t h= %.16e\t coef=%.16e \n", alpha,bns->dt, -bns->tauInv);

	const int Nr = 32;   dfloat complex R[Nr]; 

	for(int ind =1; ind <= Nr; ++ind){
		const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr; 
		R[ind-1] = cexp(I*M_PI* theta);
	}

	double complex ca21 = 0. + 0.* I; 
	double complex cb1  = 0. + 0.* I; 
	double complex ca31 = 0. + 0.* I; 
	double complex cb2  = 0. + 0.* I; 
	double complex ca32 = 0. + 0.* I; 
	double complex cb3  = 0. + 0.* I; 

	for(int i = 0; i<Nr; ++i ){
		complex lr = alpha  + R[i];
		ca21 +=  (cexp(lr/2.) - 1.)/lr; 
		ca31 += -1.0*(cexp(lr)-1.0)/lr; 
		ca32 +=  2.0*(cexp(lr)-1.0)/lr;
		//
		cb1  +=  (-4. -lr + cexp(lr)*(4. -3.*lr + cpow(lr,2.)))/ cpow(lr,3.);
		cb2  +=  4.*(2. + lr + cexp(lr)*(-2. + lr)) / cpow(lr,3.) ;
		cb3  +=  (-4. -3.*lr - cpow(lr,2.)+ cexp(lr)*(4. - lr))/ cpow(lr,3.) ;
	}

	//  Exponential Coefficients
	bns->SARK_A[1][0] = creal(ca21)/ (double) Nr; 
	bns->SARK_A[2][0] = creal(ca31)/ (double) Nr; 
	bns->SARK_A[2][1] = creal(ca32)/ (double) Nr; 

	// If 1/tau*h is too small say <<1, need to write in terms of Taylor coefficients
	bns->SARK_B[0] = creal(cb1) / (double) Nr;
	bns->SARK_B[1] = creal(cb2) / (double) Nr;
	bns->SARK_B[2] = creal(cb3) / (double) Nr;
	//
	//
	bns->SARK_C[0] = exp(alpha*c2); 
	bns->SARK_C[1] = exp(alpha*c3); 
	bns->SARK_C[2] = exp(alpha*1.0);

	//printf("%.14e\t%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n",bns->SARK_A[1][0],bns->SARK_A[2][0],bns->SARK_A[2][1],
    //                                                       bns->SARK_B[0],  bns->SARK_B[1],  bns->SARK_B[2]);  
}
 



else if(strstr(options, "LSIMEX")){ 
	int Nimex = 4;
	dfloat ImB[4]   ={0.0, 673488652607.0 /2334033219546.0, 493801219040.0/853653026979.0, 184814777513.0/1389668723319.0 };
	dfloat ImC[4]   = { 0.0, 3375509829940.0/4525919076317.0, 272778623835.0/1039454778728.0, 1.0};
	dfloat ImAd[4]  = {0.0, 3375509829940.0/4525919076317.0, 566138307881.0/912153721139.0, 184814777513.0/1389668723319.0};

	dfloat ImAmBim[4] = {0.0, 0.0,
	        -11712383888607531889907.0/32694570495602105556248.0 - 673488652607.0 /2334033219546.0,
	         0.0};

	dfloat ImAmBex[4] = {0.0,
	          3375509829940.0/4525919076317.0,
	          272778623835.0/1039454778728.0 - 673488652607.0 /2334033219546.0,
	          1660544566939.0/2334033219546.0-493801219040.0/853653026979.0 };


	bns->Nimex = Nimex;
	memcpy(bns->LSIMEX_B, ImB, Nimex*sizeof(dfloat));
	memcpy(bns->LSIMEX_C, ImC, Nimex*sizeof(dfloat));
	memcpy(bns->LSIMEX_Ad, ImAd, Nimex*sizeof(dfloat));
	memcpy(bns->LSIMEX_ABi, ImAmBim, Nimex*sizeof(dfloat));
	memcpy(bns->LSIMEX_ABe, ImAmBex, Nimex*sizeof(dfloat));

	
}  







}
