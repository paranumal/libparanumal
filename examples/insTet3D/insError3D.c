#include "ins3D.h"

// currently maximum
void insError3D(ins_t *ins, dfloat time, char *options){

    mesh3D *mesh = ins->mesh;

    dfloat maxU = 0, minU = 1e9;
    dfloat maxV = 0, minV = 1e9;
    dfloat maxW = 0, minW = 1e9;
    dfloat maxP = 0, minP = 1e9; 
    for(iint e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        iint id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        dfloat z = mesh->z[id];

        id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
        maxU = mymax(maxU, fabs(ins->U[id]));
        minU = mymin(minU, fabs(ins->U[id]));
        
        maxV = mymax(maxV, fabs(ins->V[id]));
        minV = mymin(minV, fabs(ins->V[id]));
        
        maxW = mymax(maxW, fabs(ins->W[id]));
        minW = mymin(minW, fabs(ins->W[id]));  
        //
        maxP = mymax(maxP, fabs(ins->P[id]));
        minP = mymin(minP, fabs(ins->P[id]));
      }
    }

    // compute maximum over all processes
    dfloat gMaxU, gMinU;
    MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minU, &gMinU, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    dfloat gMaxV, gMinV;
    MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minV, &gMinV, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    dfloat gMaxW, gMinW;
    MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minW, &gMinW, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    
    dfloat gMaxP, gMinP;
    MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minP, &gMinP, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0)
      printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minW: %g maxW: %g minP: %g maxP: %g\n", 
	         (int)(time/ins->dt), time, gMinU, gMaxU, gMinV, gMaxV, gMinW, gMaxW, gMinP, gMaxP );

    if( isnan(gMinU) || isnan(gMaxU) || 
        isnan(gMinV) || isnan(gMaxV) || 
        isnan(gMinW) || isnan(gMaxW) || 
        isnan(gMinP) || isnan(gMaxP) )
      exit(EXIT_FAILURE);

 
  
}
