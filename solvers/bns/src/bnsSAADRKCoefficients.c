/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "bns.h"

void bnsSAADRKCoefficients(bns_t *bns, setupAide &options){

	mesh_t * mesh = bns->mesh; 

	dfloat alpha = -bns->tauInv*bns->dt, coef = -bns->tauInv, h   = bns->dt; 

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
	dfloat rkC[bns->NrkStages]  = {1.0, exp(0.5*alpha), exp(0.5*alpha), exp(alpha), exp(alpha)  };
	dfloat rkA[bns->NrkStages*bns->NrkStages]
	                ={             0.0,             0.0,            0.0,          0.0,             0.0,
	                               a21,             0.0,            0.0,          0.0,             0.0,
	                               a31,             a32,            0.0,          0.0,             0.0,
	                               a41,             0.0,            a43,          0.0,             0.0,
	                               a51,             a52,            a53,          a54,             0.0}; 
	dfloat rkE[bns->NrkStages]= {  0.0,             0.0,            0.0,       -a54,             a54}; 

	// move data to device
	bns->o_sarkC = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), rkC);
	bns->o_sarkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), rkA);
    bns->o_sarkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), rkE); 
}

else if(bns->NrkStages==7){

	   // Initialize complex variables for contour integral
	double complex ca21 = 0. + 0.* I; 

	double complex ca31 = 0. + 0.* I; 
	double complex ca32 = 0. + 0.* I; 

	double complex ca41 = 0. + 0.* I; 
	double complex ca43 = 0. + 0.* I; 


	double complex ca51 = 0. + 0.* I; 
	double complex ca52 = 0. + 0.* I; 
	double complex ca53 = 0. + 0.* I; 
	double complex ca54 = 0. + 0.* I; 

	double complex ca61 = 0. + 0.* I; 
	double complex ca62 = 0. + 0.* I; 	 
	double complex ca63 = 0. + 0.* I; 
	double complex ca64 = 0. + 0.* I; 
	double complex ca65 = 0. + 0.* I;

	double complex ca71 = 0. + 0.* I; 
	// double complex ca72 = 0. + 0.* I; 	 
	double complex ca73 = 0. + 0.* I; 
	double complex ca74 = 0. + 0.* I; 
	double complex ca75 = 0. + 0.* I;
	double complex ca76 = 0. + 0.* I;


    double complex cb1 = 0. + 0.* I; 
	// double complex cb2 = 0. + 0.* I; 
	double complex cb3 = 0. + 0.* I; 
	double complex cb4 = 0. + 0.* I; 
	// double complex cb5 = 0. + 0.* I;
	double complex cb6 = 0. + 0.* I;

	for(int i = 0; i<Nr; ++i ){
		complex lr = alpha  + R[i];

		ca21 +=  (cexp(lr/4.) - 1.)/lr; 

		ca31 +=  (cexp(lr/4.)*lr + 4. - 4.*cexp(lr/4))/cpow(lr,2); 
		ca32 +=  (4.*cexp(lr/4.) - lr - 4.)/cpow(lr,2);  

		ca41 +=  ((cexp(lr/2) + 1)*lr + 4 - 4*cexp(lr/2))/cpow(lr,2);
		ca43 +=  (4*cexp(lr/2) - 2*lr - 4)/cpow(lr,2);

		ca51 +=  ((2*cexp((3*lr)/4) + 1)*lr + 4 - 4*cexp((3*lr)/4))/(2*cpow(lr,2));   
		ca52 +=   -(cexp((3*lr)/4) - 1)/(2*lr);
		ca53 +=  (cexp((3*lr)/4) - 1)/(2*lr);  
		ca54 +=  (4*cexp((3*lr)/4) - 3*lr - 4)/(2*cpow(lr,2));  


		ca61 +=  ((- 77*cexp(lr) - 41)*lr + 118*cexp(lr) - 118)/(42*cpow(lr,2));   
		ca62 +=  (8*cexp(lr) - 8)/(7*lr);  
		ca63 +=  ((111*cexp(lr) + 63)*lr + 174 - 174*cexp(lr))/(28*cpow(lr,2));  
		ca64 += -(12*cexp(lr) - 12)/(7*lr);
		ca65 +=  ((- 47*cexp(lr) - 239)*lr + 286*cexp(lr) - 286)/(84*cpow(lr,2));


		ca71 +=  ((1799*cexp(lr) - 511)*cpow(lr,2) + (- 6958*cexp(lr) - 4382)*lr + 11340*cexp(lr) - 11340)/(2700*cpow(lr,3));   
		// ca72 +=  (8*cexp(lr) - 8)/(7*lr);  
		ca73 +=  ((1097*cexp(lr) + 287)*cpow(lr,2) + (1834 - 934*cexp(lr))*lr + 900 - 900*cexp(lr))/(1350*cpow(lr,3));  
		ca74 +=  ((112 - 98*cexp(lr))*cpow(lr,2) + (796*cexp(lr) + 824)*lr + 1620 - 1620*cexp(lr))/(225*cpow(lr,3));
		ca75 +=  ((- 313*cexp(lr) - 1183)*cpow(lr,2) + (1766*cexp(lr) - 1226)*lr + 540 - 540*cexp(lr))/(1350*cpow(lr,3));
		ca76 +=  ((509*cexp(lr) - 1741)*cpow(lr,2) + (- 4258*cexp(lr) - 6722)*lr + 10980*cexp(lr) - 10980)/(2700*cpow(lr,3));

    
    cb1 +=  ((313*cexp(lr) + 1183)*cpow(lr,2) + (1226 - 1766*cexp(lr))*lr + 540*cexp(lr) - 540)/(5400*cpow(lr,3));   
		// cb2 +=  (8*cexp(lr) - 8)/(7*lr);  
		cb3 +=  ((- 313*cexp(lr) - 1183)*cpow(lr,2) + (1766*cexp(lr) - 1226)*lr + 540 - 540*cexp(lr))/(1350*cpow(lr,3));  

		cb4 +=  ((313*cexp(lr) + 1183)*cpow(lr,2) + (1226 - 1766*cexp(lr))*lr + 540*cexp(lr) - 540)/(900*cpow(lr,3));
		// cb5 +=  ((- 313*cexp(lr) - 1183)*cpow(lr,2) + (1766*cexp(lr) - 1226)*lr + 540 - 540*cexp(lr))/(1350*cpow(lr,3));
		cb6 += ((313*cexp(lr) + 1183)*cpow(lr,2) + (1226 - 1766*cexp(lr))*lr + 540*cexp(lr) - 540)/(5400*cpow(lr,3));
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

	dfloat a61=creal(ca61)/ (double) Nr; 
	dfloat a62=creal(ca62)/ (double) Nr; 
	dfloat a63=creal(ca63)/ (double) Nr; 
	dfloat a64=creal(ca64)/ (double) Nr; 
	dfloat a65=creal(ca65)/ (double) Nr; 

	dfloat a71=creal(ca71)/ (double) Nr; 
	dfloat a73=creal(ca73)/ (double) Nr; 
	dfloat a74=creal(ca74)/ (double) Nr; 
	dfloat a75=creal(ca75)/ (double) Nr; 
	dfloat a76=creal(ca76)/ (double) Nr; 


	dfloat b1=creal(cb1)/ (double) Nr; 
	dfloat b3=creal(cb3)/ (double) Nr; 
	dfloat b4=creal(cb4)/ (double) Nr; 
	dfloat b6=creal(cb6)/ (double) Nr; 



  dfloat rkC[bns->NrkStages]  = {1.0, exp(0.25*alpha), exp(0.25*alpha),exp(0.5*alpha), exp(0.75*alpha), exp(alpha), exp(alpha)};
	dfloat rkA[bns->NrkStages*bns->NrkStages]
		                ={             0,    0,     0,     0,     0,    0,   0,
		                               a21,  0,     0,     0,     0,    0,   0,
		                               a31,  a32,   0,     0,     0,    0,   0,
		                               a41,  0,    a43,    0,     0,    0,   0,
		                               a51, a52,   a53,   a54,     0,   0,   0,
		                               a61, a62,   a63,   a64,   a65,   0,   0,
		                               a71,  0.,   a73,   a74,   a75,   a76, 0}; 

	dfloat rkE[bns->NrkStages]= {b1, 0, b3, b4, a75, b6, 0}; 


	// move data to device
	bns->o_sarkC = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), rkC);
	bns->o_sarkA = mesh->device.malloc(bns->NrkStages*bns->NrkStages*sizeof(dfloat), rkA);
	bns->o_sarkE = mesh->device.malloc(bns->NrkStages*sizeof(dfloat), rkE);


}

}