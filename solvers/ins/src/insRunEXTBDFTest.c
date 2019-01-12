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

void extbdfCoefficents(ins_t *ins, int order);

void insRunEXTBDFTest(ins_t *ins, int maxiter){
   printf("RUNNING.........\n");   
  mesh_t *mesh = ins->mesh;
  
  for(int tstep=0; tstep<maxiter;  ++tstep){

    if(tstep<1) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      extbdfCoefficents(ins,tstep+1);

    dfloat time = ins->startTime + tstep*ins->dt;

    hlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

    printf("ADVECTION STEP.........");   
    if(ins->Nsubsteps) {
      insSubCycle(ins, time, ins->Nstages, ins->o_U, ins->o_NU);
    } else {
      insAdvection(ins, time, ins->o_U, ins->o_NU);
    }

    printf("done\n");   

    insGradient (ins, time, ins->o_P, ins->o_GP);
    
    
    insVelocityRhs  (ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW);
    insVelocitySolve(ins, time+ins->dt, ins->Nstages, ins->o_rhsU, ins->o_rhsV, ins->o_rhsW, ins->o_rkU);


    insPressureRhs  (ins, time+ins->dt, ins->Nstages);
    // insPressureSolve(ins, time+ins->dt, ins->Nstages); 

    insPressureUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkP);
    insGradient(ins, time+ins->dt, ins->o_rkP, ins->o_rkGP);

    
    //cycle history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_U.copyFrom(ins->o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_P.copyFrom(ins->o_P, ins->Ntotal*sizeof(dfloat), 
                                  (s-1)*ins->Ntotal*sizeof(dfloat), 
                                  (s-2)*ins->Ntotal*sizeof(dfloat));
    }

    //copy updated pressure
    ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

    //update velocity
    insVelocityUpdate(ins, time+ins->dt, ins->Nstages, ins->o_rkGP, ins->o_rkU);

#if 0
    if(mesh->dim==3){
      if(ins->elementType == QUADRILATERALS ||
	 ins->elementType == TRIANGLES){
	ins->constrainKernel(mesh->Nelements,
			     offset,
			     mesh->o_x,
			     mesh->o_y,
			     mesh->o_z,
			     ins->o_rkU);
      }
    }
#endif
      
    

    
    //copy updated pressure
    ins->o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

    //cycle rhs history
    for (int s=ins->Nstages;s>1;s--) {
      ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
      ins->o_GP.copyFrom(ins->o_GP, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
			 (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    }



#if 0
    if(((tstep+1)%10)==0){
      char fname[BUFSIZ];
      sprintf(fname, "Iterations_Quad3D_%d", ins->SNrk);
      FILE *fp; 
      fp = fopen(fname, "a");
      fprintf(fp, "%d %d %d %d %d\n", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP);
      fclose(fp);
    }
#endif

    // if (ins->dim==2 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterP); fflush(stdout);
    // if (ins->dim==3 && mesh->rank==0) printf("\rtstep = %d, solver iterations: U - %3d, V - %3d, W - %3d, P - %3d", tstep+1, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); fflush(stdout);
    
 }


//   dfloat finalTime = ins->NtimeSteps*ins->dt;
//   printf("\n");
  
  // if(mesh->rank==0) occa::printTimer();
}

