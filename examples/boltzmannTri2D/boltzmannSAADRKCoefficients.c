#include "boltzmann2D.h"

void boltzmannSAADRKCoefficients(bns_t *bns, char *options){

mesh2D * mesh = bns->mesh; 

dfloat alpha = -bns->tauInv*bns->dt, coef = -bns->tauInv, h   = bns->dt; 

// printf("alpha = %.16e\t h= %.16e\t coef=%.16e \n", alpha,bns->dt, -bns->tauInv);

const int Nr = 32;   dfloat complex R[Nr]; 

for(int ind =1; ind <= Nr; ++ind){
	const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr; 
	R[ind-1] = cexp(I*M_PI* theta);
}



if(bns->NrkStages==5){
    // Initialize complex variables for contour integral
	double complex ca21 = 0. + 0.* I; 
	double complex ca31 = 0. + 0.* I; 
	double complex ca32 = 0. + 0.* I; 
	double complex ca41 = 0. + 0.* I; 
	// double complex ca42 = 0. + 0.* I; 
	double complex ca43 = 0. + 0.* I; 
	double complex ca51 = 0. + 0.* I; 
	double complex ca52 = 0. + 0.* I; 
	double complex ca53 = 0. + 0.* I; 
	double complex ca54 = 0. + 0.* I; 

	for(int i = 0; i<Nr; ++i ){
		complex lr = alpha  + R[i];

		ca21 +=  (cexp(lr/2.) - 1.)/lr; 

		ca31 +=  (4. + cexp(lr/2.)*(-4. + lr) + lr)/cpow(lr,2); 
		ca32 +=  (4.*cexp(lr/2.) -2.*lr -4.)/cpow(lr,2); 

		ca41 +=  ((cexp(lr) + 1.)*lr + 2. - 2.*cexp(lr))/cpow(lr,2);  
		ca43 +=  (2.*cexp(lr) - 2.*lr - 2.)/cpow(lr,2); 

		ca51 +=  (cexp(lr)*cpow(lr,2) + (- 3.*cexp(lr) - 1.)*lr + 4.*cexp(lr) - 4)/cpow(lr,3);  
		ca52 +=  ((2.*cexp(lr) + 2.)*lr + 4. - 4.*cexp(lr))/cpow(lr,3); 
		ca53 +=  ((2.*cexp(lr) + 2.)*lr + 4. - 4*cexp(lr))/cpow(lr,3);  
		ca54 +=  ((-cexp(lr) - 3.)*lr - cpow(lr,2) + 4.*cexp(lr) - 4.)/cpow(lr,3);  
	}
	

	dfloat a21=creal(ca21)/ (double) Nr; 

	dfloat a31=creal(ca31)/ (double) Nr; 
	dfloat a32=creal(ca32)/ (double) Nr; 

	dfloat a41=creal(ca41)/ (double) Nr; 
	dfloat a43=creal(ca43)/ (double) Nr; 

	dfloat a51=creal(ca51)/ (double) Nr; 
	dfloat a52=creal(ca52)/ (double) Nr; 
	dfloat a53=creal(ca53)/ (double) Nr; 
	dfloat a54=creal(ca54)/ (double) Nr; 
  //

// first set non-semianalytic part of the integrator 
	dfloat sarkC[bns->NrkStages]  = {0.0, exp(0.5*alpha), exp(0.5*alpha), exp(alpha), exp(alpha)  };
	dfloat sarkA[bns->NrkStages*bns->NrkStages]
	                ={             0.0,             0.0,            0.0,          0.0,             0.0,
	                               a21,             0.0,            0.0,          0.0,             0.0,
	                               a31,             a32,            0.0,          0.0,             0.0,
	                               a41,             0.0,            a43,          0.0,             0.0,
	                               a51,             a52,            a53,          a54,             0.0}; 
	dfloat sarkE[bns->NrkStages]= {  0.0,             0.0,            0.0,       -a54,             a54}; 


  // move data to device
	bns->o_sarkC = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), sarkC);
	bns->o_sarkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), sarkA);
  bns->o_sarkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), sarkE);
}







}