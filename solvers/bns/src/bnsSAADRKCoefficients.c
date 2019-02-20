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

	const int Nr = 32;
	std::complex<dfloat> R[Nr];

	for(int ind =1; ind <= Nr; ++ind){
		const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
		std::complex<dfloat> z = 0. + M_PI* theta*1i;
		R[ind-1] = exp(z);
	}



if(bns->NrkStages==5){
    // Initialize complex variables for contour integral
	std::complex<double> ca21 = 0. + 0.* 1i;
	std::complex<double> ca31 = 0. + 0.* 1i;
	std::complex<double> ca32 = 0. + 0.* 1i;
	std::complex<double> ca41 = 0. + 0.* 1i;
	// std::complex<double> ca42 = 0. + 0.* 1i;
	std::complex<double> ca43 = 0. + 0.* 1i;
	std::complex<double> ca51 = 0. + 0.* 1i;
	std::complex<double> ca52 = 0. + 0.* 1i;
	std::complex<double> ca53 = 0. + 0.* 1i;
	std::complex<double> ca54 = 0. + 0.* 1i;

	for(int i = 0; i<Nr; ++i ){
		std::complex<double> lr = alpha  + R[i];

		ca21 +=  (exp(lr/2.) - 1.)/lr;

		ca31 +=  (4. + exp(lr/2.)*(-4. + lr) + lr)/pow(lr,2);
		ca32 +=  (4.*exp(lr/2.) -2.*lr -4.)/pow(lr,2);

		ca41 +=  ((exp(lr) + 1.)*lr + 2. - 2.*exp(lr))/pow(lr,2);
		ca43 +=  (2.*exp(lr) - 2.*lr - 2.)/pow(lr,2);

		ca51 +=  (exp(lr)*pow(lr,2) + (- 3.*exp(lr) - 1.)*lr + 4.*exp(lr) - 4.)/pow(lr,3);
		ca52 +=  ((2.*exp(lr) + 2.)*lr + 4. - 4.*exp(lr))/pow(lr,3);
		ca53 +=  ((2.*exp(lr) + 2.)*lr + 4. - 4.*exp(lr))/pow(lr,3);
		ca54 +=  ((-exp(lr) - 3.)*lr - pow(lr,2) + 4.*exp(lr) - 4.)/pow(lr,3);
	}


	dfloat a21=real(ca21)/ (double) Nr;

	dfloat a31=real(ca31)/ (double) Nr;
	dfloat a32=real(ca32)/ (double) Nr;

	dfloat a41=real(ca41)/ (double) Nr;
	dfloat a43=real(ca43)/ (double) Nr;

	dfloat a51=real(ca51)/ (double) Nr;
	dfloat a52=real(ca52)/ (double) Nr;
	dfloat a53=real(ca53)/ (double) Nr;
	dfloat a54=real(ca54)/ (double) Nr;
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
	std::complex<double> ca21 = 0. + 0.* 1i;

	std::complex<double> ca31 = 0. + 0.* 1i;
	std::complex<double> ca32 = 0. + 0.* 1i;

	std::complex<double> ca41 = 0. + 0.* 1i;
	std::complex<double> ca43 = 0. + 0.* 1i;


	std::complex<double> ca51 = 0. + 0.* 1i;
	std::complex<double> ca52 = 0. + 0.* 1i;
	std::complex<double> ca53 = 0. + 0.* 1i;
	std::complex<double> ca54 = 0. + 0.* 1i;

	std::complex<double> ca61 = 0. + 0.* 1i;
	std::complex<double> ca62 = 0. + 0.* 1i;
	std::complex<double> ca63 = 0. + 0.* 1i;
	std::complex<double> ca64 = 0. + 0.* 1i;
	std::complex<double> ca65 = 0. + 0.* 1i;

	std::complex<double> ca71 = 0. + 0.* 1i;
	// std::complex<double> ca72 = 0. + 0.* 1i;
	std::complex<double> ca73 = 0. + 0.* 1i;
	std::complex<double> ca74 = 0. + 0.* 1i;
	std::complex<double> ca75 = 0. + 0.* 1i;
	std::complex<double> ca76 = 0. + 0.* 1i;


    std::complex<double> cb1 = 0. + 0.* 1i;
	// std::complex<double> cb2 = 0. + 0.* 1i;
	std::complex<double> cb3 = 0. + 0.* 1i;
	std::complex<double> cb4 = 0. + 0.* 1i;
	// std::complex<double> cb5 = 0. + 0.* 1i;
	std::complex<double> cb6 = 0. + 0.* 1i;

	for(int i = 0; i<Nr; ++i ){
		std::complex<double> lr = alpha  + R[i];

		ca21 +=  (exp(lr/4.) - 1.)/lr;

		ca31 +=  (exp(lr/4.)*lr + 4. - 4.*exp(lr/4.))/pow(lr,2);
		ca32 +=  (4.*exp(lr/4.) - lr - 4.)/pow(lr,2);

		ca41 +=  ((exp(lr/2.) + 1.)*lr + 4. - 4.*exp(lr/2.))/pow(lr,2);
		ca43 +=  (4.*exp(lr/2.) - 2.*lr - 4.)/pow(lr,2);

		ca51 +=  ((2.*exp((3.*lr)/4.) + 1.)*lr + 4. - 4.*exp((3.*lr)/4.))/(2.*pow(lr,2));
		ca52 +=   -(exp((3.*lr)/4.) - 1.)/(2.*lr);
		ca53 +=  (exp((3.*lr)/4.) - 1.)/(2.*lr);
		ca54 +=  (4.*exp((3.*lr)/4.) - 3.*lr - 4.)/(2.*pow(lr,2));


		ca61 +=  ((- 77.*exp(lr) - 41.)*lr + 118.*exp(lr) - 118.)/(42.*pow(lr,2));
		ca62 +=  (8.*exp(lr) - 8.)/(7.*lr);
		ca63 +=  ((111.*exp(lr) + 63.)*lr + 174. - 174.*exp(lr))/(28.*pow(lr,2));
		ca64 += -(12.*exp(lr) - 12.)/(7.*lr);
		ca65 +=  ((- 47.*exp(lr) - 239.)*lr + 286.*exp(lr) - 286.)/(84.*pow(lr,2));


		ca71 +=  ((1799.*exp(lr) - 511.)*pow(lr,2) + (- 6958.*exp(lr) - 4382.)*lr + 11340.*exp(lr) - 11340.)/(2700.*pow(lr,3));
		// ca72 +=  (8.*exp(lr) - 8.)/(7.*lr);
		ca73 +=  ((1097.*exp(lr) + 287.)*pow(lr,2) + (1834. - 934.*exp(lr))*lr + 900. - 900.*exp(lr))/(1350.*pow(lr,3));
		ca74 +=  ((112. - 98.*exp(lr))*pow(lr,2) + (796.*exp(lr) + 824.)*lr + 1620. - 1620.*exp(lr))/(225.*pow(lr,3));
		ca75 +=  ((- 313.*exp(lr) - 1183.)*pow(lr,2) + (1766.*exp(lr) - 1226.)*lr + 540. - 540.*exp(lr))/(1350.*pow(lr,3));
		ca76 +=  ((509.*exp(lr) - 1741.)*pow(lr,2) + (- 4258.*exp(lr) - 6722.)*lr + 10980.*exp(lr) - 10980.)/(2700.*pow(lr,3));


    cb1 +=  ((313.*exp(lr) + 1183.)*pow(lr,2) + (1226. - 1766.*exp(lr))*lr + 540.*exp(lr) - 540.)/(5400.*pow(lr,3));
		// cb2 +=  (8.*exp(lr) - 8.)/(7.*lr);
		cb3 +=  ((- 313.*exp(lr) - 1183.)*pow(lr,2) + (1766.*exp(lr) - 1226.)*lr + 540. - 540.*exp(lr))/(1350.*pow(lr,3));

		cb4 +=  ((313.*exp(lr) + 1183.)*pow(lr,2) + (1226. - 1766.*exp(lr))*lr + 540.*exp(lr) - 540.)/(900.*pow(lr,3));
		// cb5 +=  ((- 313.*exp(lr) - 1183.)*pow(lr,2) + (1766.*exp(lr) - 1226.)*lr + 540. - 540.*exp(lr))/(1350.*pow(lr,3));
		cb6 += ((313.*exp(lr) + 1183.)*pow(lr,2) + (1226. - 1766.*exp(lr))*lr + 540.*exp(lr) - 540.)/(5400.*pow(lr,3));
	}

  dfloat a21=real(ca21)/ (double) Nr;

	dfloat a31=real(ca31)/ (double) Nr;
	dfloat a32=real(ca32)/ (double) Nr;

	dfloat a41=real(ca41)/ (double) Nr;
	dfloat a43=real(ca43)/ (double) Nr;

	dfloat a51=real(ca51)/ (double) Nr;
	dfloat a52=real(ca52)/ (double) Nr;
	dfloat a53=real(ca53)/ (double) Nr;
	dfloat a54=real(ca54)/ (double) Nr;

	dfloat a61=real(ca61)/ (double) Nr;
	dfloat a62=real(ca62)/ (double) Nr;
	dfloat a63=real(ca63)/ (double) Nr;
	dfloat a64=real(ca64)/ (double) Nr;
	dfloat a65=real(ca65)/ (double) Nr;

	dfloat a71=real(ca71)/ (double) Nr;
	dfloat a73=real(ca73)/ (double) Nr;
	dfloat a74=real(ca74)/ (double) Nr;
	dfloat a75=real(ca75)/ (double) Nr;
	dfloat a76=real(ca76)/ (double) Nr;


	dfloat b1=real(cb1)/ (double) Nr;
	dfloat b3=real(cb3)/ (double) Nr;
	dfloat b4=real(cb4)/ (double) Nr;
	dfloat b6=real(cb6)/ (double) Nr;



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
