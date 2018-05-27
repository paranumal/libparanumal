#include "mns.h"


void mnsError(mns_t *mns, dfloat time){

  mesh_t *mesh = mns->mesh;

  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  dfloat maxU   = 0, minU   = 1e9;
  dfloat maxV   = 0, minV   = 1e9;
  dfloat maxW   = 0, minW   = 1e9;
  dfloat maxP   = 0, minP   = 1e9; 
  dfloat maxPhi = 0, minPhi = 1e9; 


  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];

      maxU  = mymax(maxU, fabs(mns->U[id+0*offset]));
      minU  = mymin(minU, fabs(mns->U[id+0*offset]));
      
      maxV  = mymax(maxV, fabs(mns->U[id+1*offset]));
      minV  = mymin(minV, fabs(mns->U[id+1*offset]));
      

      if (mns->dim==3) {
        maxW = mymax(maxW, fabs(mns->U[id+3*offset]));
        minW = mymin(minW, fabs(mns->U[id+3*offset]));  
      }

      maxP = mymax(maxP, fabs(mns->P[id]));
      minP = mymin(minP, fabs(mns->P[id]));

      maxPhi = mymax(maxPhi, fabs(mns->Phi[id]));
      minPhi = mymin(minPhi, fabs(mns->Phi[id]));
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

  dfloat gMaxPhi, gMinPhi;
  MPI_Allreduce(&maxPhi, &gMaxPhi, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&minPhi, &gMinPhi, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
    if (mns->dim==3) {
      printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minW: %g maxW: %g minP: %g maxP: %g minPhi: %g maxPhi: %g\n", 
         (int)(time/mns->dt), time, gMinU, gMaxU, gMinV, gMaxV, gMinW, gMaxW, gMinP, gMaxP, gMinPhi, gMaxPhi );
    } else {
      printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minP: %g maxP: %g minPhi: %g maxPhi: %g\n", 
         (int)(time/mns->dt), time, gMinU, gMaxU, gMinV, gMaxV, gMinP, gMaxP, gMinPhi, gMaxPhi);
    }

  if( isnan(gMinU)   || isnan(gMaxU)  || 
      isnan(gMinV)   || isnan(gMaxV)  || 
      isnan(gMinW)   || isnan(gMaxW)  || 
      isnan(gMinP)   || isnan(gMaxP)  ||
      isnan(gMinPhi) || isnan(gMaxPhi) )
    exit(EXIT_FAILURE);

}
