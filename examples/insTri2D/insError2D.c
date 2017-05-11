#include "ins2D.h"

// currently maximum
void insError2D(solver_t *ins, dfloat time,char *options){

    mesh2D *mesh = ins->mesh;

    dfloat maxUx = 0, minUx = 1e9;
    dfloat maxUy = 0, minUy = 1e9;
    dfloat maxPr = 0, minPr = 1e9; 
    for(iint e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        iint id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];

        maxUx = mymax(maxUx, fabs(ins->U[id + UXID*mesh->Nelements*mesh->Np]));
        minUx = mymin(minUx, fabs(ins->U[id + UXID*mesh->Nelements*mesh->Np]));
        //
        maxUy = mymax(maxUy, fabs(ins->U[id + UYID*mesh->Nelements*mesh->Np]));
        minUy = mymin(minUy, fabs(ins->U[id + UYID*mesh->Nelements*mesh->Np]));
        //
        maxPr = mymax(maxPr, fabs(ins->Pr[id]));
        minPr = mymin(minPr, fabs(ins->Pr[id]));
      }
    }

    // compute maximum over all processes
    dfloat gMaxUx, gMinUx;
    MPI_Allreduce(&maxUx, &gMaxUx, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minUx, &gMinUx, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    //
    dfloat gMaxUy, gMinUy;
    MPI_Allreduce(&maxUy, &gMaxUy, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minUy, &gMinUy, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    //
    dfloat gMaxPr, gMinPr;
    MPI_Allreduce(&maxPr, &gMaxPr, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minPr, &gMinPr, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0)
      printf("Time: %g minUx: %g maxUx: %g minUy: %g maxUy: %g minPr: %g maxPr: %g\n", 
              time, gMinUx, gMaxUx, gMinUy, gMaxUy,gMinPr, gMaxPr );

    if( isnan(gMinUx) || isnan(gMaxUx) || 
        isnan(gMinUy) || isnan(gMaxUy) || 
        isnan(gMinPr) || isnan(gMaxPr) )
      exit(EXIT_FAILURE);

 
  
}
