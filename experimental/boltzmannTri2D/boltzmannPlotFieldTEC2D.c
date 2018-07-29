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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "mesh2D.h"

// interpolate data to plot nodes and save to file (one per process
void boltzmannPlotFieldTEC2D(mesh2D *mesh, char *fileName, dfloat time, int fld){

    
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

  int TotalPoints = mesh->Nelements*mesh->plotNp;
  int TotalCells  = mesh->Nelements*mesh->plotNelements;

  fprintf(fp, "ZONE  N = %d  E = %d  F=FEPOINT , ET=TRIANGLE , STRANDID=1, SOLUTIONTIME = %.8f \n", TotalPoints,TotalCells,time);

  // Write data
  // compute plot node coordinates on the fly
  for(int e=0;e<mesh->Nelements;++e){
    // First compute the curl and div
    
    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotfld=0;
      for(int m=0;m<mesh->Np;++m){
        plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];

        int base    = mesh->Nfields*(m + e*mesh->Np);
        dfloat field = mesh->q[base + fld];
        plotfld     += mesh->plotInterp[n*mesh->Np+m]*field;
     }

      fprintf(fp,"%.10e\t%.10e\t%.10e\n",plotxn, plotyn, plotfld);
    }
}


  // Write Connectivity
   for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->plotNelements;++n){
      for(int m=0;m<mesh->plotNverts;++m){
        fprintf(fp, "%9d\t ", 1 + e*mesh->plotNp + mesh->plotEToV[n*mesh->plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
}
