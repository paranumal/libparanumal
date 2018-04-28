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


if(strstr(options,"IMEXRK")){

	if(bns->NrkStages==6){ // Kennedy-Carpanter RK34

		dfloat rkCex[bns->NrkStages] ={0.0, 1.0/2.0, 83.0/250.0, 31.0/50.0, 17.0/20.0, 1.0};
		dfloat rkCim[bns->NrkStages] ={0.0, 1.0/2.0, 83.0/250.0, 31.0/50.0, 17.0/20.0, 1.0};
       
        // First fill the explicit part; 
        dfloat a11 = 0.0;

        dfloat a21 = 1.0/2.0; 
        dfloat a22 = 0.0; 

        dfloat a31 = 13861./62500.; 
        dfloat a32 = 6889./62500.; 
        dfloat a33 = 0.0; 
        
        dfloat a41 = -116923316275./2393684061468.; 
        dfloat a42 = -2731218467317./15368042101831.; 
        dfloat a43 = 9408046702089./11113171139209.; 
        dfloat a44 = 0.0; 

        dfloat a51 = -451086348788./2902428689909.; 
        dfloat a52 = -2682348792572./7519795681897.; 
        dfloat a53 = 12662868775082./11960479115383.; 
        dfloat a54 = 3355817975965./11060851509271.; 
        dfloat a55 = 0.0; 

        dfloat a61 = 647845179188./3216320057751.; 
        dfloat a62 = 73281519250./8382639484533.; 
        dfloat a63 = 552539513391./3454668386233.; 
        dfloat a64 = 3354512671639./8306763924573.;
        dfloat a65 = 4040./17871.; 
        dfloat a66 = 0.0; 


        dfloat b1  = 82889/524892;
        dfloat b2  = 0.0;
        dfloat b3  = 15625./83664.;
        dfloat b4  = 69875./102672.;
        dfloat b5  = -2260./8211.;
        dfloat b6  = 1.0/4.0;

        dfloat be1 = 4586570599./29645900160.;
        dfloat be2 = 0.0;
        dfloat be3 = 178811875.0/945068544.0;
        dfloat be4 = 814220225./1159782912.;
        dfloat be5 = -3700637./11593932.;
        dfloat be6 = 61727./225920.;



				dfloat rkAex[bns->NrkStages*bns->NrkStages] = {a11 ,  0., 0.,  0.,  0.,  0., 
		                                               a21,  a22, 0.,  0.,  0.,  0., 
		                                               a31,  a32, a33, 0.,  0.,  0.,
		                                               a41,  a42, a43, a44, 0.,  0., 
		                                               a51,  a52, a53, a54, a55, 0., 
		                                               a61,  a62, a63, a64, a65, a66 };

				dfloat rkBex[bns->NrkStages] = {b1, b2, b3, b4, b5, b6};

				dfloat rkEex[bns->NrkStages] = {b1-be1, b2-be2, b3-be3, b4-be4, b5-be5, b6-be6};



		// First fill the explicit part; 
        a11 = 0.0;

        a21 = 1.0/4.0; 
        a22 = 1.0/4.0; 
        
        a31 = 8611./62500.; 
        a32 = -1743./31250.; 
        a33 = 1.0/4.0; 
        
        a41 = 5012029./34652500.; 
        a42 = -654441./2922500.; 
        a43 = 174375./388108.; 
        a44 = 1.0/4.0;
        
        a51 = 15267082809./155376265600.; 
        a52 = -71443401./120774400.; 
        a53 = 730878875./902184768.; 
        a54 = 2285395./8070912.; 
        a55 = 1.0/4.0; 
        
        a61 = 82889./524892.; 
        a62 = 0.0; 
        a63 = 15625./83664.; 
        a64 = 69875./102672.;
        a65 = -2260./8211; 
        a66 = 1.0/4.0; 
        
        b1  = a61;
        b2  = a62;
        b3  = a63;
        b4  = a64;
        b5  = a65;
        b6  = a66;
        
        be1 = 4586570599./29645900160.0;
        be2 = 0.0;
        be3 = 178811875./945068544.;
        be4 = 814220225./1159782912.;
        be5 = -3700637./11593932.;
        be6 = 61727./225920.0;



        dfloat rkAim[bns->NrkStages*bns->NrkStages] = {a11 ,  0., 0.,  0.,  0.,  0., 
		                                               a21,  a22, 0.,  0.,  0.,  0., 
		                                               a31,  a32, a33, 0.,  0.,  0.,
		                                               a41,  a42, a43, a44, 0.,  0., 
		                                               a51,  a52, a53, a54, a55, 0., 
		                                               a61,  a62, a63, a64, a65, a66 };

		dfloat rkBim[bns->NrkStages] = {b1, b2, b3, b4, b5, b6};

		dfloat rkEim[bns->NrkStages] = {b1-be1, b2-be2, b3-be3, b4-be4, b5-be5, b6-be6};


		bns->rkCex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkBex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkAex    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
		bns->rkEex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));

		bns->rkCim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkBim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkAim    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
		bns->rkEim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));



		memcpy(bns->rkCex, rkCex, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkBex, rkBex, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkAex, rkAex, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
		memcpy(bns->rkEex, rkEex, bns->NrkStages*sizeof(dfloat));


		memcpy(bns->rkCim, rkCim, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkBim, rkBim, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkAim, rkAim, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
		memcpy(bns->rkEim, rkEim, bns->NrkStages*sizeof(dfloat));

	}




 	else if(bns->NrkStages==8){

 		printf("setting order 8 \n");
		dfloat c0 = 0.0; 	
		dfloat c1 = 41.0/100.0;
		dfloat c2 = 2935347310677.0/11292855782101.0;
		dfloat c3 = 1426016391358.0/7196633302097.0;
		dfloat c4 = 92.0/100.0;
		dfloat c5 = 24.0/100.0;
		dfloat c6 = 3.0/5.0;
		dfloat c7 = 1.0; 

		dfloat rkCex[bns->NrkStages] ={c0, c1, c2, c3, c4, c5, c6, c7};
		dfloat rkCim[bns->NrkStages] ={c0, c1, c2, c3, c4, c5, c6, c7};
       
		// First fill the explicit part; 
		dfloat a00 = 0.0; 

		dfloat a10 = 41.0/100.0;
		dfloat a11 = 0.0;

		dfloat a20 = 367902744464.0  /2072280473677.0;
		dfloat a21 = 677623207551.0  /8224143866563.0;
		dfloat a22 = 0.0;

		dfloat a30 = 1268023523408.0 /10340822734521.0;
		dfloat a31 = 0.0;
		dfloat a32 = 1029933939417.0 /13636558850479.0;
		dfloat a33 = 0.0;

		dfloat a40 = 14463281900351.0/6315353703477.0;
		dfloat a41 = 0.0;
		dfloat a42 = 66114435211212.0/5879490589093.0;
		dfloat a43 = -54053170152839.0/4284798021562.0;
		dfloat a44 = 0.0;

		dfloat a50 = 14090043504691.0/34967701212078.0;
		dfloat a51 = 0.0;
		dfloat a52 = 15191511035443.0/11219624916014.0;
		dfloat a53 = -18461159152457.0/12425892160975.0;
		dfloat a54 = -281667163811.0/9011619295870.0;
		dfloat a55 = 0.0;

		dfloat a60 = 19230459214898.0/13134317526959.0;
		dfloat a61 = 0.0;
		dfloat a62 = 21275331358303.0/2942455364971.0;
		dfloat a63 = -38145345988419.0/4862620318723.0;
		dfloat a64 = -1.0/8.0;
		dfloat a65 = -1.0/8.0;
		dfloat a66 = 0.0;

		dfloat a70 = -19977161125411.0/11928030595625.0;
		dfloat a71 = 0.0;
		dfloat a72 = -40795976796054.0/6384907823539.0;
		dfloat a73 = 177454434618887.0/12078138498510.0;
		dfloat a74 = 782672205425.0/8267701900261.0;
		dfloat a75 = -69563011059811.0/9646580694205.0;
		dfloat a76 = 7356628210526.0/4942186776405.0;
		dfloat a77 = 0.0;

		dfloat b0 = -872700587467.0/9133579230613.0;
		dfloat b1 = 0.0;
		dfloat b2 = 0.0;
		dfloat b3 = 22348218063261.0/9555858737531.0;
		dfloat b4 = -1143369518992.0/8141816002931.0;
		dfloat b5 = -39379526789629.0/19018526304540.0;
		dfloat b6 = 32727382324388.0/42900044865799.0;
		dfloat b7 = 41.0/200.0;

		dfloat be0 = -975461918565.0/9796059967033.0;
		dfloat be1 = 0.0;
		dfloat be2 = 0.0;
		dfloat be3 = 78070527104295.0/32432590147079.0;
		dfloat be4 = -548382580838.0/3424219808633.0;
		dfloat be5 = -33438840321285.0/15594753105479.0;
		dfloat be6 = 3629800801594.0/4656183773603.0;
		dfloat be7 = 4035322873751.0/18575991585200.0;


		dfloat rkAex[bns->NrkStages*bns->NrkStages] = {a00,  0.,  0.,  0.,  0.,  0.,  0., 0., 
		                                               a10,  a11, 0.,  0.,  0.,  0.,  0., 0., 
		                                               a20,  a21, a22, 0.,  0.,  0.,  0., 0., 
		                                               a30,  a31, a32, a33, 0.,  0.,  0., 0., 
		                                               a40,  a41, a42, a43, a44, 0.,  0., 0., 
		                                               a50,  a51, a52, a53, a54, a55, 0., 0., 
		                                               a60,  a61, a62, a63, a64, a65, a66, 0., 
		                                               a70,  a71, a72, a73, a74, a75, a76, a77};

		dfloat rkBex[bns->NrkStages] = {b0, b1, b2, b3, b4, b5, b6, b7};

		dfloat rkEex[bns->NrkStages] = {b0-be0, b1-be1, b2-be2, b3-be3, b4-be4, b5-be5, b6-be6, b7-be7};


		a00 = 0.0;

		a10 = 41.0/200.0;
		a11 = 41.0/200.0;

		a20 = 41.0/400.0;
		a21 = -567603406766.0/11931857230679.0;
		a22 = 41.0/200.0;

		a30 = 683785636431.0/9252920307686.0;
		a31 = 0.0; 
		a32 = -110385047103.0/1367015193373.0;
		a33 = 41.0/200.0;

		a40 = 3016520224154.0/10081342136671.0;
		a41 = 0.0; 
		a42 = 30586259806659.0/12414158314087.0;
		a43 = -22760509404356.0/11113319521817.0;
		a44 = 41.0/200.0;

		a50 = 218866479029.0/1489978393911.0;
		a51 = 0.0; 
		a52 = 638256894668.0/5436446318841.0;
		a53 = -1179710474555.0/5321154724896.0;
		a54 = -60928119172.0/8023461067671.0;
		a55 = 41.0/200.0;

		a60 = 1020004230633.0/5715676835656.0;
		a61 = 0.0; 
		a62 = 25762820946817.0/25263940353407.0;
		a63 = -2161375909145.0/9755907335909.0;
		a64 = -211217309593.0/5846859502534.0;
		a65 = -4269925059573.0/7827059040749.0;
		a66 = 41.0/200.0;

		a70 = -872700587467.0/9133579230613.0;
		a71 = 0.0;
		a72 = 0.0;  
		a73 = 22348218063261.0/9555858737531.0;
		a74 = -1143369518992.0/8141816002931.0;
		a75 = -39379526789629.0/19018526304540.0;
		a76 = 32727382324388.0/42900044865799.0;
		a77 = 41.0/200.0;


		b0 = -872700587467.0/9133579230613.0;
		b1 = 0.0;
		b2 = 0.0;  
		b3 = 22348218063261.0/9555858737531.0;
		b4 = -1143369518992.0/8141816002931.0;
		b5 = -39379526789629.0/19018526304540.0;
		b6 = 32727382324388.0/42900044865799.0;
		b7 = 41.0/200.0;

		be0 = -975461918565.0/9796059967033.0;
		be1 = 0.0; 
		be2 = 0.0; 
		be3 = 78070527104295.0/32432590147079.0;
		be4 = -548382580838.0/3424219808633.0;
		be5 = -33438840321285.0/15594753105479.0;
		be6 = 3629800801594.0/4656183773603.0;
		be7 = 4035322873751.0/18575991585200.0;


		dfloat rkAim[bns->NrkStages*bns->NrkStages] = {a00 ,  0., 0.,  0.,  0.,  0., 0., 0., 
		                                           a10,  a11, 0.,  0.,  0.,  0., 0., 0., 
		                                           a20,  a21, a22, 0.,  0.,  0., 0., 0., 
		                                           a30,  a31, a32, a33, 0.,  0., 0., 0., 
		                                           a40,  a41, a42, a43, a44, 0., 0., 0., 
		                                           a50,  a51, a52, a53, a54, a55, 0., 0., 
		                                           a60,  a61, a62, a63, a64, a65, a66, 0., 
		                                           a70,  a71, a72, a73, a74, a75, a76, a77};

		dfloat rkBim[bns->NrkStages] = {b0, b1, b2, b3, b4, b5, b6, b7};

		dfloat rkEim[bns->NrkStages] = {b0-be0, b1-be1, b2-be2, b3-be3, b4-be4, b5-be5, b6-be6, b7-be7};


		bns->rkCex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkBex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkAex    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
		bns->rkEex    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));

		bns->rkCim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkBim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
		bns->rkAim    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
		bns->rkEim    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));



		memcpy(bns->rkCex, rkCex, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkBex, rkBex, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkAex, rkAex, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
		memcpy(bns->rkEex, rkEex, bns->NrkStages*sizeof(dfloat));


		memcpy(bns->rkCim, rkCim, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkBim, rkBim, bns->NrkStages*sizeof(dfloat)); 
		memcpy(bns->rkAim, rkAim, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
		memcpy(bns->rkEim, rkEim, bns->NrkStages*sizeof(dfloat));
 	}




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

  
  bns->rkC    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
  bns->rkA    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
  bns->rkE    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));


  memcpy(bns->rkC, rkC, bns->NrkStages*sizeof(dfloat)); 
  memcpy(bns->rkA, rkA, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
  memcpy(bns->rkE, rkE, bns->NrkStages*sizeof(dfloat));
}


if(strstr(options,"SAADRK")){

	// dfloat rkC[bns->NrkStages], rkA[bns->NrkStages*bns->NrkStages], rkE[bns->NrkStages]; 
	if(bns->NrkStages==5){ // SAARK43

		// first set non-semianalytic part of the integrator 
		dfloat rkC[bns->NrkStages]  = {0.0, 0.5, 0.5, 1.0, 1.0};
		dfloat rkA[bns->NrkStages*bns->NrkStages]  
		   = {             0.0,             0.0,            0.0,          0.0,             0.0,
                           0.5,             0.0,            0.0,          0.0,             0.0,
                           0.0,             0.5,            0.0,          0.0,             0.0,
                           0.0,             0.0,            1.0,          0.0,             0.0,
                         1.0/6.0,         1.0/3.0,        1.0/3.0,       1.0/6.0,          0.0}; 
		dfloat rkE[bns->NrkStages]= {  0.0,             0.0,            0.0,        -1.0/6.0,        1.0/6.0}; 

	bns->rkC    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
	bns->rkA    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
	bns->rkE    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
	

	memcpy(bns->rkC, rkC, bns->NrkStages*sizeof(dfloat)); 
	memcpy(bns->rkA, rkA, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
	memcpy(bns->rkE, rkE, bns->NrkStages*sizeof(dfloat));
	// Compute semi-analytic part of the integrator

	}
	else if(bns->NrkStages==7){	

      printf("Numbe of stages in SAADRK is 7\n");

		dfloat rkC[bns->NrkStages]   = {0.0, 0.25, 0.25, 0.5, 0.75, 1.0, 1.0};
		dfloat rkA[bns->NrkStages*bns->NrkStages]  
		   =  {  0,    0,        0,     0,     0,    0,   0,
		         1/4,   0,        0,     0,     0,    0,   0,
                 1/8,  1/8,       0,     0,     0,    0,   0,
                  0,    0,        1/2,     0,     0,    0,   0,
               3./16., -3./8.,   3./8.,  9./16.,     0,    0, 0,
               -3./7.,  8./7.,   6./7., -12./7.,   8./7.,    0, 0,
               7./90.,    0.,    16./45.,  2./15., 16./45., 7./90., 0}; 

		dfloat rkE[bns->NrkStages]=  {-4./45., 0, 16./45., -8./15., 16./45., -4./45., 0.}; 

    bns->rkC    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
	bns->rkA    = (dfloat*) calloc(bns->NrkStages*bns->NrkStages, sizeof(dfloat));
	bns->rkE    = (dfloat*) calloc(bns->NrkStages, sizeof(dfloat));
	

	memcpy(bns->rkC, rkC, bns->NrkStages*sizeof(dfloat)); 
	memcpy(bns->rkA, rkA, bns->NrkStages*bns->NrkStages*sizeof(dfloat));
	memcpy(bns->rkE, rkE, bns->NrkStages*sizeof(dfloat));
		
	}

	

	boltzmannSAADRKCoefficients(bns, options);
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
