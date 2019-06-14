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

#include "ins.h"

// dfloat insL1Norm(ins_t* ins, dfloat*U); 
// dfloat insL2Norm(ins_t* ins, dfloat*U); 
dfloat insLInfNorm(ins_t* ins, dfloat*U); 

void insError(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;

  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  dfloat maxU = 0, minU = 1e9;
  dfloat maxV = 0, minV = 1e9;
  dfloat maxW = 0, minW = 1e9;
  dfloat maxP = 0, minP = 1e9; 
  
  if (ins->options.compareArgs("EXACT", "NONE")) {// just compute maximums
 
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        int id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        dfloat z = mesh->z[id];

        maxU = mymax(maxU, fabs(ins->U[id+0*offset]));
        minU = mymin(minU, fabs(ins->U[id+0*offset]));
        
        maxV = mymax(maxV, fabs(ins->U[id+1*offset]));
        minV = mymin(minV, fabs(ins->U[id+1*offset]));
        
        if (ins->dim==3) {
          maxW = mymax(maxW, fabs(ins->U[id+2*offset]));
          minW = mymin(minW, fabs(ins->U[id+2*offset]));  
        }

        maxP = mymax(maxP, fabs(ins->P[id]));
        minP = mymin(minP, fabs(ins->P[id]));
      }
    }

    // compute maximum over all processes
    dfloat gMaxU, gMinU;
    MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minU, &gMinU, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxV, gMinV;
    MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minV, &gMinV, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxW, gMinW;
    MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minW, &gMinW, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
    
    dfloat gMaxP, gMinP;
    MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minP, &gMinP, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    if(mesh->rank==0)
      if (ins->dim==3) {
        printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minW: %g maxW: %g minP: %g maxP: %g\n", 
           (int)((time-ins->startTime)/ins->dt)+1, time, gMinU, gMaxU, gMinV, gMaxV, gMinW, gMaxW, gMinP, gMaxP );
      } else {
        printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minP: %g maxP: %g\n", 
           (int)((time-ins->startTime)/ins->dt)+1, time, gMinU, gMaxU, gMinV, gMaxV, gMinP, gMaxP );
      }

    if( isnan(gMinU) || isnan(gMaxU) || 
        isnan(gMinV) || isnan(gMaxV) || 
        isnan(gMinW) || isnan(gMaxW) || 
        isnan(gMinP) || isnan(gMaxP) )
      exit(EXIT_FAILURE);
  } 
 else { //compare to an exact solution
    // 
    dfloat *uErr, *vErr, *wErr, *pErr; 
    uErr = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
    vErr = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));
    if(ins->dim==3)
       wErr = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));

    pErr = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));

    if (ins->options.compareArgs("EXACT","VORTEX")) { //2D Taylor vortex

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          int id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];

          dfloat uExact = -sin(2.f*M_PI*y)*exp(-ins->nu*4.f*M_PI*M_PI*time);
          dfloat vExact =  sin(2.f*M_PI*x)*exp(-ins->nu*4.f*M_PI*M_PI*time);
          dfloat pExact = -cos(2.f*M_PI*x)*cos(2.f*M_PI*y)*exp(-ins->nu*8.f*M_PI*M_PI*time);

           // Store for L2 Norm of Error
          uErr[id] = fabs(ins->U[id+0*offset]-uExact);
          vErr[id] = fabs(ins->U[id+1*offset]-vExact);
          pErr[id] = fabs(ins->P[id+0*offset]-pExact);
          #if 1
            ins->U[id+0*offset] -= uExact;
            ins->U[id+1*offset] -= vExact;
            ins->P[id] -= pExact;
          #endif
        }
      }
    } else if (ins->options.compareArgs("EXACT","KOVASZNAY")) { //2D Kovasznay flow

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          int id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];

          dfloat lambda = 0.5f*1.0/ins->nu - sqrt(1.f/(4.0f*ins->nu*ins->nu) + 4.f*M_PI*M_PI);

          dfloat uExact = 1.0 - exp(lambda*x)*cos(2.f*M_PI*y);
          dfloat vExact = 0.5f*lambda/M_PI*exp(lambda*x)*sin(2.f*M_PI*y);
          dfloat pExact = 0.5f*(1.f - exp(2.f*lambda*x));

           // Store for L2 Norm of Error
          uErr[id] = fabs(ins->U[id+0*offset]-uExact);
          vErr[id] = fabs(ins->U[id+1*offset]-vExact);
          pErr[id] = fabs(ins->P[id+0*offset]-pExact);
          #if 1
            ins->U[id+0*offset] = uErr[id];
            ins->U[id+1*offset] = vErr[id];
            ins->P[id]          = pErr[id];
          #endif
        }
      }
    } else if (ins->options.compareArgs("EXACT","BELTRAMI")) { //3D Beltrami flow
      
      
      dfloat lpam = 0.0; // approximated pr 
      dfloat lpem = 0.0; // exact pressure 
      dfloat lpvl = 0.0; // volume of the domain

      // Compute Exact Pressure and project
      dfloat *pExact  = (dfloat *)calloc(mesh->Nelements*mesh->Np, sizeof(dfloat));

      for(dlong e=0;e<mesh->Nelements;++e){

        for(int n=0;n<mesh->Np;++n){
          dlong id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

          dfloat a = M_PI/4.f;
          dfloat d = M_PI/2.f;

          pExact[id] = -0.5f*a*a*exp(-2.f*ins->nu*d*d*time)*(
                            2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+ 
                            2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+ 
                            2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+
                            exp(2.f*a*z) + exp(2.f*a*y) +exp(2.f*a*x) );
        }
      }
      occa::memory o_pExact = mesh->device.malloc(mesh->Nelements*mesh->Np*sizeof(dfloat), pExact);

      ellipticZeroMean(ins->pSolver, o_pExact); 
      ellipticZeroMean(ins->pSolver, ins->o_P); 

      o_pExact.copyTo(pExact); 
      ins->o_P.copyTo(ins->P);

      dfloat *pEx  = (dfloat *)calloc(mesh->Np, sizeof(dfloat));
      dfloat *pPr  = (dfloat *)calloc(mesh->Np, sizeof(dfloat));
      dfloat *pVl  = (dfloat *)calloc(mesh->Np, sizeof(dfloat));
     
      for(dlong e=0;e<mesh->Nelements;++e){

        for(int n=0;n<mesh->Np;++n){
          dlong id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];
          //
          dfloat J  = mesh->vgeo[e*mesh->Np*mesh->Nvgeo+mesh->Np*JID + n]; 
#if 0
          dfloat a = M_PI/4.f;
          dfloat d = M_PI/2.f;
          dfloat pExact = -0.5f*a*a*exp(-2.f*ins->nu*d*d*time)*(
                            2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+ 
                            2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+ 
                            2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+
                            exp(2.f*a*z) + exp(2.f*a*y) +exp(2.f*a*x) );

          pEx[n]  = J*pExact; 
          pPr[n]  = J*ins->P[id]; 
          pVl[n]  = J; 
#else
          pEx[n]  = J*pExact[id]; 
          pPr[n]  = J*ins->P[id]; 
          pVl[n]  = J; 
#endif
        }

        for(int n=0;n<mesh->Np;++n){    
          dfloat sumex  = 0.0; 
          dfloat sumpr  = 0.0; 
          dfloat sumvol = 0.0; 
          for(int i=0; i<mesh->Np; i++){
            sumex  += mesh->MM[n*mesh->Np +i]*pEx[i]; 
            sumpr  += mesh->MM[n*mesh->Np +i]*pPr[i]; 
            sumvol += mesh->MM[n*mesh->Np +i]*pVl[i]; 
          }
          //
          lpem  += sumex; 
          lpam  += sumpr; 
          lpvl  += sumvol; 
        }
      }

      free(pEx); free(pPr); free(pVl); 

        dfloat gpam  = 0.0; 
        dfloat gpem  = 0.0;
        dfloat gpvl  = 0.0;

        MPI_Allreduce(&lpam, &gpam, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
        MPI_Allreduce(&lpem, &gpem, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
        MPI_Allreduce(&lpvl, &gpvl, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

        gpam /= (gpvl);   
        gpem /= (gpvl); 

       printf("gpmean: %.4e and gpexmean: %.4e gVol = %.4e\n ", gpam, gpem, gpvl);

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          dlong id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

          dfloat a = M_PI/4.f;
          dfloat d = M_PI/2.f;

          dfloat uExact = -a*exp(-ins->nu*d*d*time)*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y));
          dfloat vExact = -a*exp(-ins->nu*d*d*time)*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z));
          dfloat wExact = -a*exp(-ins->nu*d*d*time)*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x));

#if 0

           dfloat pExact = -0.5f*a*a*exp(-2.f*ins->nu*d*d*time)*(
                            2.f*exp(a*(z+y))*cos(a*z+d*x)*sin(a*x+d*y)+ 
                            2.f*exp(a*(z+x))*cos(a*x+d*y)*sin(a*y+d*z)+ 
                            2.f*exp(a*(y+x))*cos(a*y+d*z)*sin(a*z+d*x)+
                            exp(2.f*a*z) + exp(2.f*a*y) +exp(2.f*a*x) );


          // move exact solution to the same mean of approximated pressure
          pExact = (pExact-gpem) + gpam; 
          // Store for L2 Norm of Error
          uErr[id] = fabs(ins->U[id+0*offset]-uExact);
          vErr[id] = fabs(ins->U[id+1*offset]-vExact);
          wErr[id] = fabs(ins->U[id+2*offset]-wExact);
          pErr[id] = fabs(ins->P[id+0*offset]-pExact);
#else
           // Store for L2 Norm of Error
          uErr[id] = fabs(ins->U[id+0*offset]-uExact);
          vErr[id] = fabs(ins->U[id+1*offset]-vExact);
          wErr[id] = fabs(ins->U[id+2*offset]-wExact);
          pErr[id] = fabs(ins->P[id+0*offset]-((pExact[id] -gpem) + gpam));
#endif
        
          #if 0
            ins->U[id+0*offset] -= uExact;
            ins->U[id+1*offset] -= vExact;
            ins->U[id+2*offset] -= wExact;
            ins->P[id]          -= pExact;
          #endif
        }
      }
    }


      // compute maximum over all processes
      dfloat uLiNorm, vLiNorm, wLiNorm, pLiNorm; 
      dfloat uL2norm, vL2norm, wL2norm, pL2norm; 
      
      uLiNorm = insLInfNorm(ins, uErr); 
      vLiNorm = insLInfNorm(ins, vErr); 
      if(ins->dim==3) wLiNorm = insLInfNorm(ins, wErr); 
      pLiNorm = insLInfNorm(ins, pErr); 

      // uL2norm  = insL2Norm(ins, uErr); 
      // vL2norm  = insL2Norm(ins, vErr); 
      // if(ins->dim==3) wL2norm  = insL2Norm(ins, wErr); 

      // pL2norm  = insL2Norm(ins, pErr); 


      if(ins->dim==2){

      if( isnan(uLiNorm) || isnan(vLiNorm) || isnan(pLiNorm) )
        exit(EXIT_FAILURE);
       if(mesh->rank==0){
       printf("Linf error-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorP: %.4e\n", 
           (int)(time/ins->dt), time, uLiNorm, vLiNorm, pLiNorm);

        // printf("L2 error-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorP: %.4e\n", 
           // (int)(time/ins->dt), time, uL2norm, vL2norm, pL2norm);
      }

      }else if(ins->dim=3){
         if( isnan(uLiNorm) || isnan(vLiNorm) || isnan(wLiNorm) || isnan(pLiNorm) )
           exit(EXIT_FAILURE);

         if(mesh->rank==0){
        printf("Linf error\t-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorW: %.4e ErrorP: %.4e\n", 
           (int)(time/ins->dt), time, uLiNorm, vLiNorm, wLiNorm, pLiNorm);

        // printf("L2 error\t-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorW: %.4e ErrorP: %.4e\n", 
           // (int)(time/ins->dt), time, uL2norm, vL2norm, wL2norm, pL2norm);
      }


      }
      free(uErr);  
      free(vErr);  
      if(ins->dim==3) 
        free(wErr);  
      free(pErr); 
  }
}


// // l2 norm
// dfloat insL2Norm(ins_t *ins, dfloat *U){

// mesh_t *mesh = ins->mesh; 
// dfloat l2norm = 0.0; 

//  for(dlong e=0;e<mesh->Nelements;++e){
//     for(int n=0;n<mesh->Np;++n){  
//         dfloat sum = 0.0;  
//           for(int i=0;i<mesh->Np;++i){  
//             dfloat J  = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + mesh->Np*JID + i]; 
//             dfloat ui = U[e*mesh->Np+i]; 
//             sum += mesh->MM[n*mesh->Np +i]*J*ui*ui; 
//           }
        
//         l2norm += sum; 
//       }
//     }

// dfloat gl2norm; 
// MPI_Allreduce(&l2norm, &gl2norm, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

// return sqrt(gl2norm);

// }


// // l2 norm
// dfloat insL1Norm(ins_t *ins, dfloat *U){

// mesh_t *mesh = ins->mesh; 
// dfloat l1norm = 0.0; 

//  for(dlong e=0;e<mesh->Nelements;++e){
//     for(int n=0;n<mesh->Np;++n){  
//         dfloat sum = 0.0;  
//           for(int i=0;i<mesh->Np;++i){  
//             dfloat J  = mesh->vgeo[e*mesh->Np*mesh->Nvgeo + mesh->Np*JID + i]; 
//             dfloat ui = U[e*mesh->Np+i]; 
//             sum += mesh->MM[n*mesh->Np +i]*J*ui; 
//           }
        
//         l1norm += sum; 
//       }
//     }

// dfloat gl1norm; 
// MPI_Allreduce(&l1norm, &gl1norm, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

// return gl1norm;

// }


// l2 norm
dfloat insLInfNorm(ins_t *ins, dfloat *U){

mesh_t *mesh = ins->mesh; 
dfloat linfnorm = 0.0; 

 for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){  
      linfnorm = mymax(U[e*mesh->Np+n],linfnorm);                   
      }
    }

dfloat glinfnorm; 
MPI_Allreduce(&linfnorm, &glinfnorm, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

return glinfnorm;

}
