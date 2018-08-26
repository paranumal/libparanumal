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

#include "boltzmann2D.h"

// interpolate data to plot nodes and save to file (one per process
void boltzmannPlotTEC2D(bns_t *bns, char *fileName, dfloat time){

  mesh2D *mesh = bns->mesh; 

  FILE *fp;
  
  fp = fopen(fileName, "a");

  fprintf(fp,"#ELEMENT NUMBER        =%d\n",mesh->Nelements);
  fprintf(fp,"#SOLUTION ORDER        =%d\n",mesh->N);
  fprintf(fp,"#POSTPROCESS NODES     =%d\n",mesh->plotNp);
  fprintf(fp,"#POSTPROCESS Elements  =%d\n",mesh->plotNelements);
  fprintf(fp,"#SPEED of SOUND        =%f\n",bns->sqrtRT);

  fprintf(fp,"VARIABLES=x,y,u,v,p,w,div\n");

  int TotalPoints = mesh->Nelements*mesh->plotNp;
  int TotalCells  = mesh->Nelements*mesh->plotNelements;

  fprintf(fp, "ZONE  N = %d  E = %d  F=FEPOINT , ET=TRIANGLE , STRANDID=1, SOLUTIONTIME = %.8f \n", TotalPoints,TotalCells,time);

  


  dfloat *curlU   = (dfloat *) calloc(mesh->Np, sizeof(dfloat));
  dfloat *divU    = (dfloat *) calloc(mesh->Np, sizeof(dfloat));
  // Write data
  // compute plot node coordinates on the fly
  for(int e=0;e<mesh->Nelements;++e){
    // First compute the curl and div
    for(int n=0;n<mesh->Np;++n){
      dfloat dUdr = 0, dUds = 0, dVdr = 0, dVds = 0;
      for(int m=0;m<mesh->Np;++m){
        int base = bns->Nfields*(m + e*mesh->Np);
        dfloat rho = bns->q[base + 0];
        dfloat u = bns->q[1 + base]*bns->sqrtRT/rho;
        dfloat v = bns->q[2 + base]*bns->sqrtRT/rho;
        dUdr += mesh->Dr[n*mesh->Np+m]*u;
        dUds += mesh->Ds[n*mesh->Np+m]*u;
        dVdr += mesh->Dr[n*mesh->Np+m]*v;
        dVds += mesh->Ds[n*mesh->Np+m]*v;
      }

      dfloat rx = mesh->vgeo[e*mesh->Nvgeo+RXID];
      dfloat ry = mesh->vgeo[e*mesh->Nvgeo+RYID];
      dfloat sx = mesh->vgeo[e*mesh->Nvgeo+SXID];
      dfloat sy = mesh->vgeo[e*mesh->Nvgeo+SYID];

      dfloat dUdx = rx*dUdr + sx*dUds;
      dfloat dUdy = ry*dUdr + sy*dUds;
      dfloat dVdx = rx*dVdr + sx*dVds;
      dfloat dVdy = ry*dVdr + sy*dVds;

      curlU[n] = dVdx-dUdy;
      divU[n]  = dUdx+dVdy;
    }   

    for(int n=0;n<mesh->plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotun=0;
      dfloat plotvn = 0, plotpn =0,  plotwn=0, plotdn=0;

      for(int m=0;m<mesh->Np;++m){

        plotxn += mesh->plotInterp[n*mesh->Np+m]*mesh->x[m+e*mesh->Np];
        plotyn += mesh->plotInterp[n*mesh->Np+m]*mesh->y[m+e*mesh->Np];
        //
        int base = bns->Nfields*(m + e*mesh->Np);
        dfloat rho = bns->q[base];
        dfloat pm = bns->sqrtRT*bns->sqrtRT*rho; 
        dfloat um = bns->q[1 + base]*bns->sqrtRT/rho;
        dfloat vm = bns->q[2 + base]*bns->sqrtRT/rho;
        //
        dfloat wz = curlU[m];
        dfloat du = divU[m];
        //
        plotpn  += mesh->plotInterp[n*mesh->Np+m]*pm;
        plotun  += mesh->plotInterp[n*mesh->Np+m]*um;
        plotvn  += mesh->plotInterp[n*mesh->Np+m]*vm;
        plotwn  += mesh->plotInterp[n*mesh->Np+m]*wz;
        plotdn  += mesh->plotInterp[n*mesh->Np+m]*du;
     }

      fprintf(fp,"%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",plotxn,plotyn,plotun,plotvn,plotpn,plotwn,plotdn);
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


  free(curlU);
  free(divU);

}
