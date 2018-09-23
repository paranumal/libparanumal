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
          maxW = mymax(maxW, fabs(ins->U[id+3*offset]));
          minW = mymin(minW, fabs(ins->U[id+3*offset]));  
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
  } else { //compare to an exact solution

    if (ins->options.compareArgs("EXACT","VORTEX")) { //2D Taylor vortex

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          int id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];

          dfloat uExact = -sin(2.f*M_PI*y)*exp(-ins->nu*4.f*M_PI*M_PI*time);
          dfloat vExact =  sin(2.f*M_PI*x)*exp(-ins->nu*4.f*M_PI*M_PI*time);
          dfloat pExact = -cos(2.f*M_PI*x)*cos(2.f*M_PI*y)*exp(-ins->nu*8.f*M_PI*M_PI*time);

          maxU = mymax(maxU, fabs(ins->U[id+0*offset]-uExact));
          maxV = mymax(maxV, fabs(ins->U[id+1*offset]-vExact));
          maxP = mymax(maxP, fabs(ins->P[id]-pExact));

          #if 0
            ins->U[id+0*offset] -= uExact;
            ins->U[id+1*offset] -= vExact;
            ins->P[id] -= pExact;
          #endif
        }
      }

      // compute maximum over all processes
      dfloat gMaxU;
      MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      dfloat gMaxV;
      MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
      
      dfloat gMaxP;
      MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      if(mesh->rank==0)
        printf("Step: %d Time: %g ErrorU: %g ErrorV: %g ErrorP: %g \n", 
           (int)(time/ins->dt), time, gMaxU, gMaxV, gMaxP);

      if( isnan(gMaxU) || 
          isnan(gMaxV) || 
          isnan(gMaxP) )
        exit(EXIT_FAILURE);
    } else if (ins->options.compareArgs("EXACT","BELTRAMI")) { //3D Beltrami flow

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          dlong id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

          dfloat a = M_PI/4.f;
          dfloat d = M_PI/2.f;

          dfloat uExact = -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-d*d*time);
          dfloat vExact = -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-d*d*time);
          dfloat wExact = -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-d*d*time);
          dfloat pExact = -a*a*exp(-2.f*d*d*time)*(exp(2.f*a*x)+exp(2.f*a*y)+exp(2.f*a*z))
                                                  *( sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))
                                                    +sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))
                                                    +sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))); 

          maxU = mymax(maxU, fabs(ins->U[id+0*offset]-uExact));
          maxV = mymax(maxV, fabs(ins->U[id+1*offset]-vExact));
          maxW = mymax(maxW, fabs(ins->U[id+2*offset]-wExact));
          maxP = mymax(maxP, fabs(ins->P[id]-pExact));

          #if 0
            ins->U[id+0*offset] -= uExact;
            ins->U[id+1*offset] -= vExact;
            ins->U[id+2*offset] -= wExact;
            ins->P[id] -= pExact;
          #endif
        }
      }

      // compute maximum over all processes
      dfloat gMaxU;
      MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      dfloat gMaxV;
      MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      dfloat gMaxW;
      MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
      
      dfloat gMaxP;
      MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      if(mesh->rank==0)
        printf("Step: %d Time: %g ErrorU: %g ErrorV: %g ErrorW: %g ErrorP: %g \n", 
           (int)(time/ins->dt), time, gMaxU, gMaxV, gMaxW, gMaxP);

      if( isnan(gMaxU) || 
          isnan(gMaxV) || 
          isnan(gMaxW) || 
          isnan(gMaxP) )
        exit(EXIT_FAILURE);
    }
  }
}
