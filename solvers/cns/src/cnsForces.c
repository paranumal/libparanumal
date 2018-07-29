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

#include "cns.h"


void cnsForces(cns_t *cns, dfloat time){


  mesh_t *mesh = cns->mesh; 
  // Hard coded weights i.e. sum(MM1D,1), MM1D = inv(V1)'*V1
  dfloat W[mesh->Nfp];
  if(mesh->N==1){
    W[0] = 1.000000000000000; W[1] = 1.000000000000000; 
  }

 if(mesh->N==2){
    W[0] = 0.333333333333333; W[1] = 1.333333333333333;
    W[2] = 0.333333333333333; 
  }

  if(mesh->N==3){
    W[0] = 0.166666666666666; W[1] = 0.833333333333334;
    W[2] = 0.833333333333334; W[3] = 0.166666666666666;
  }
  if(mesh->N==4){
    W[0] = 0.100000000000000; W[1] = 0.544444444444444;
    W[2] = 0.711111111111111; W[3] = 0.544444444444444;
    W[4] = 0.100000000000000; 
  }
  if(mesh->N==5){
    W[0] = 0.066666666666666; W[1] = 0.378474956297847;
    W[2] = 0.554858377035487; W[3] = 0.554858377035486;
    W[4] = 0.378474956297847; W[5] = 0.066666666666667;
  }

  if(mesh->N==6){
    W[0] = 0.047619047619047; W[1] = 0.276826047361566;
    W[2] = 0.431745381209863; W[3] = 0.487619047619048;
    W[4] = 0.431745381209863; W[5] = 0.276826047361566;
    W[6] = 0.047619047619047;
  }

   

  dfloat Fx =0. , Fy = 0. ;
  dfloat *dUdx = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dUdy = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dVdx = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dVdy = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *Pr   = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  //
   for(int e=0;e<mesh->Nelements;++e){
     int flag = 0; 
    for(int f=0;f<mesh->Nfaces;++f){
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if(bc == 1){ flag = 1; }
    }

    if(flag){

      dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
      dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
      dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];

      for(int n=0;n<mesh->Np;++n){
        dfloat dudr = 0, duds = 0, dvdr = 0, dvds = 0;     
        for(int i=0;i<mesh->Np;++i){
          // load data at node i of element e (note Nfields==4)
          int id = e*mesh->Np*cns->Nfields + i;
          dfloat r = mesh->q[id+0*mesh->Np];
          dfloat u = mesh->q[id+1*mesh->Np]/r;
          dfloat v = mesh->q[id+2*mesh->Np]/r;
          //  
          dfloat Drni = mesh->Dr[n*mesh->Np+i];
          dfloat Dsni = mesh->Ds[n*mesh->Np+i];
          // differentiate (u,v,p) with respect to 'r' and 's'
          dudr += Drni*u;
          duds += Dsni*u;
          dvdr += Drni*v;
          dvds += Dsni*v;   
        }
        dUdx[n] = drdx*dudr + dsdx*duds;
        dUdy[n] = drdy*dudr + dsdy*duds;
        dVdx[n] = drdx*dvdr + dsdx*dvds;
        dVdy[n] = drdy*dvdr + dsdy*dvds;
        Pr[n]   = mesh->q[e*mesh->Np*cns->Nfields + n]*cns->RT;
      }

   
      for(int f=0;f<mesh->Nfaces;++f){
        int bc = mesh->EToB[e*mesh->Nfaces+f];
        if(bc==1){
          for(int n=0;n<mesh->Nfp; n++){
            int sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
            dfloat nx = mesh->sgeo[sid+0];
            dfloat ny = mesh->sgeo[sid+1];
            dfloat sJ = mesh->sgeo[sid+2];

            int vid  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
            int idM  = mesh->vmapM[vid];
            //
            dfloat dudx = dUdx[idM-e*mesh->Np]; 
            dfloat dudy = dUdy[idM-e*mesh->Np]; 
            dfloat dvdx = dVdx[idM-e*mesh->Np]; 
            dfloat dvdy = dVdy[idM-e*mesh->Np];
            dfloat p    = Pr[idM-e*mesh->Np];  // Pressure
                       
            Fx += W[n]*(sJ*(-p*nx + cns->mu*(nx*( 2.0*dudx - 2.0/3.0 * (dudx + dvdy)) + ny*(dvdx + dudy))));
            Fy += W[n]*(sJ*(-p*ny + cns->mu*(nx*(dvdx + dudy) + ny*(2.0*dvdy- 2.0/3.0 * (dudx + dvdy)) )));
          }
        }
      }
    }
  }


  dfloat gFx = 0., gFy = 0.; 
  // Add all processors force, 
  MPI_Allreduce(&Fx, &gFx, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&Fy, &gFy, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  //
  if(mesh->rank==0){
    char fname[BUFSIZ];
    sprintf(fname, "CNSForceData_N%d.dat", mesh->N);

    FILE *fp; 
    fp = fopen(fname, "a");  

    fprintf(fp, "%.4e %.8e %.8e \n", time, gFx, gFy); 

    fclose(fp);
  }


}
