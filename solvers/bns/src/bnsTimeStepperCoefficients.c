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

void bnsTimeStepperCoefficients(bns_t *bns, setupAide &options){

  mesh_t *mesh = bns->mesh; 

  if(options.compareArgs("TIME INTEGRATOR","MRSAAB")){

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

        }else if(order==1){

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
        }else{
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


  if(options.compareArgs("TIME INTEGRATOR", "LSERK")){
  // printf("Relying the fact that  LSERK set on mesh structure\n");
  }

  if(options.compareArgs("TIME INTEGRATOR","SARK")){

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

    }else if(bns->NrkStages==7){     

    // printf("Numbe of stages in SAADRK is 7\n");

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
  bnsSAADRKCoefficients(bns, options);
  }

  }
