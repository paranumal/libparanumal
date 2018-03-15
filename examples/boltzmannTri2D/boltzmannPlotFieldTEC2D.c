#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "mesh2D.h"

// interpolate data to plot nodes and save to file (one per process
void boltzmannPlotFieldTEC2D(mesh2D *mesh, char *fileName, dfloat time, iint fld){

    
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *fp;
  
  fp = fopen(fileName, "a");

  fprintf(fp,"#ELEMENT NUMBER        =%d\n",mesh->Nelements);
  fprintf(fp,"#SOLUTION ORDER        =%d\n",mesh->N);
  fprintf(fp,"#POSTPROCESS NODES     =%d\n",mesh->plotNp);
  fprintf(fp,"#POSTPROCESS Elements  =%d\n",mesh->plotNelements);
  fprintf(fp,"#SPEED of SOUND        =%f\n",mesh->sqrtRT);

  fprintf(fp,"VARIABLES=x,y,field\n");

  iint TotalPoints = mesh->Nelements*mesh->plotNp;
  iint TotalCells  = mesh->Nelements*mesh->plotNelements;

  fprintf(fp, "ZONE  N = %d  E = %d  F=FEPOINT , ET=TRIANGLE , STRANDID=1, SOLUTIONTIME = %.8f \n", TotalPoints,TotalCells,time);

  // Write data
  // compute plot node coordinates on the fly
  for(iint e=0;e<mesh->Nelements;++e){
    // First compute the curl and div
    
    for(iint n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotfld=0;
      for(iint m=0;m<mesh->Np;++m){
        plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];

        iint base    = mesh->Nfields*(m + e*mesh->Np);
        dfloat field = mesh->q[base + fld];
        plotfld     += mesh->plotInterp[n*mesh->Np+m]*field;
     }

      fprintf(fp,"%.10e\t%.10e\t%.10e\n",plotxn, plotyn, plotfld);
    }
}


  // Write Connectivity
   for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->plotNelements;++n){
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%9d\t ", 1 + e*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
}
