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
dfloat insL2Norm(ins_t* ins, dlong offset, dfloat*U); 
dfloat insLInfNorm(ins_t* ins, dlong offset, dfloat*U); 
void insError(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;
  
  if (ins->options.compareArgs("EXACT", "NONE")) {// just compute maximums
    
    dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
    dfloat maxU = 0, minU = 1e9;
    dfloat maxV = 0, minV = 1e9;
    dfloat maxW = 0, minW = 1e9;
    dfloat maxP = 0, minP = 1e9; 
 
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

      dfloat uLi, vLi, wLi, pLi; 
      dfloat uL2, vL2, wL2, pL2; 
      dfloat uexL2, vexL2, wexL2, pexL2; 

      dfloat mexct = 0.f; 
      dfloat msoln = 0.f; 

      if(ins->pSolver->allNeumann){

       occa::memory o_P = mesh->device.malloc(ins->Ntotal*sizeof(dfloat), ins->rkP);

       ins->setFlowFieldKernel(mesh->Nelements, 
                            time, 
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            ins->fieldOffset, 
                            ins->o_rkU, // dummy 
                            o_P);
      mexct = insMean(ins, o_P);
      msoln = insMean(ins, ins->o_P);
      printf("exact_mean : %.8e solution_mean: %.8e\n", mexct, msoln);
      }

      dfloat shift = msoln - mexct; 

    // Compute exact solution in cubature nodes
    const dlong offset = mesh->Nelements*mesh->cubNp; 
      ins->setFlowFieldCubKernel(mesh->Nelements, 
                              time, 
                              mesh->o_cubInterpT,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              offset, 
                              ins->o_Uex, 
                              ins->o_Pex);

   
    ins->o_Uex.copyTo(ins->Uer); 
    ins->o_Pex.copyTo(ins->Per);
  

    // Compute the L2 norm of exact soln
    uexL2  = insL2Norm(ins, 0*offset, ins->Uer); 
    vexL2  = insL2Norm(ins, 1*offset, ins->Uer); 
      if(ins->dim==3) wexL2  = insL2Norm(ins, 2*offset, ins->Uer); 
    pexL2  = insL2Norm(ins, 0*offset, ins->Per);  
    
    // Compute the error on cubature nodes i.e. err = abs(Interp*U - Uexact)
    ins->errorKernel(mesh->Nelements,
                     mesh->o_cubInterpT,
                     ins->fieldOffset,
                     offset,
                     shift, 
                     ins->o_U,
                     ins->o_P,
                     ins->o_Uex,
                     ins->o_Pex);

    ins->o_Uex.copyTo(ins->Uer); 
    ins->o_Pex.copyTo(ins->Per);

   
      // Compute Inf Norm of error  
      uLi = insLInfNorm(ins, 0*offset, ins->Uer); 
      vLi = insLInfNorm(ins, 1*offset, ins->Uer); 
      if(ins->dim==3) wLi = insLInfNorm(ins, 2*offset, ins->Uer); 
      pLi = insLInfNorm(ins, 0*offset, ins->Per); 
      // Compute L2 Norm of error
      uL2  = insL2Norm(ins, 0*offset, ins->Uer); 
      vL2  = insL2Norm(ins, 1*offset, ins->Uer); 
      if(ins->dim==3) wL2  = insL2Norm(ins, 2*offset, ins->Uer); 

      pL2  = insL2Norm(ins, 0*offset, ins->Per);  


    if(ins->dim==2){

    if( isnan(uLi) || isnan(vLi) || isnan(pLi) )
      exit(EXIT_FAILURE);
     if(mesh->rank==0){
     printf("Linf error-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorP: %.4e\n", 
         (int)(time/ins->dt), time, uLi, vLi, pLi);

      printf("L2   error-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorP: %.4e\n", 
         (int)(time/ins->dt), time, uL2/uexL2, vL2/vexL2, pL2/pexL2);
    }

    }else if(ins->dim=3){
       if( isnan(uLi) || isnan(vLi) || isnan(wLi) || isnan(pLi) )
         exit(EXIT_FAILURE);

       if(mesh->rank==0){
      printf("Linf error\t-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorW: %.4e ErrorP: %.4e\n", 
         (int)(time/ins->dt), time, uLi, vLi, wLi, pLi);

      printf("L2   error\t-->\tStep: %d Time: %g ErrorU: %.4e ErrorV: %.4e ErrorW: %.4e ErrorP: %.4e\n", 
         (int)(time/ins->dt), time, uL2/uexL2, vL2/vexL2, wL2/wexL2, pL2/pexL2);
        }


    }
  }
}

// l2 norm
dfloat insL2Norm(ins_t *ins, dlong offset, dfloat *U){

mesh_t *mesh = ins->mesh; 
dfloat l2norm = 0.0; 

if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA){

 for(dlong e=0;e<mesh->Nelements;++e){
      dfloat sum = 0.0;  
    for(int j=0;j<mesh->cubNq;++j){
      for(int i=0;i<mesh->cubNq;++i){
          dlong vbase = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq;
          dlong nbase = mesh->cubNp*e + i + j*mesh->cubNq;
         dfloat JW  = mesh->cubvgeo[vbase + mesh->cubNp*JWID]; 
         dfloat ui = U[nbase + offset]; 
         sum +=ui*JW*ui;
       }
     }

      l2norm += sum; 
  }
}else if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA){
for(dlong e=0;e<mesh->Nelements;++e){
      dfloat sum = 0.0;  
  for(int k=0;k<mesh->cubNq;++k){
    for(int j=0;j<mesh->cubNq;++j){
      for(int i=0;i<mesh->cubNq;++i){
          dlong vbase = mesh->Nvgeo*mesh->cubNp*e + i + j*mesh->cubNq + k*mesh->cubNq*mesh->cubNq;
          dlong nbase = mesh->cubNp*e + i + j*mesh->cubNq + k*mesh->cubNq*mesh->cubNq;
         dfloat JW  = mesh->cubvgeo[vbase + mesh->cubNp*JWID]; 
         dfloat ui = U[nbase + offset]; 
         sum +=ui*JW*ui;
       }
     }
   }
  l2norm += sum; 
  }
}

dfloat gl2norm; 
MPI_Allreduce(&l2norm, &gl2norm, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

return sqrt(gl2norm);
}

// linf norm
dfloat insLInfNorm(ins_t *ins, dlong offset, dfloat *U){

  mesh_t *mesh = ins->mesh; 
  dfloat linfnorm = 0.0; 

  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->cubNp;++n){  
      linfnorm = mymax(U[e*mesh->cubNp+ n + offset],linfnorm);                   
    }
  }

  dfloat glinfnorm; 
  MPI_Allreduce(&linfnorm, &glinfnorm, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

  return glinfnorm;

}
