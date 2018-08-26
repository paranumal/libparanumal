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

#include "bns.h"


void bnsForces(bns_t *bns, dfloat time, setupAide &options){


  mesh_t *mesh = bns->mesh; 

  
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
  //
  for(int e=0;e<mesh->Nelements;++e){
    int flag = 0; 
    for(int f=0;f<mesh->Nfaces;++f){
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if(bc == 1){ flag = 1; }
    }
    if(flag){

      for(int f=0;f<mesh->Nfaces;++f){
	int bc = mesh->EToB[e*mesh->Nfaces+f];
	if(bc==1){
	  for(int n=0;n<mesh->Nfp; n++){
	    int sid  = mesh->Nsgeo*(e*mesh->Nfaces+f);
	    dfloat nx = mesh->sgeo[sid+NXID];
	    dfloat ny = mesh->sgeo[sid+NYID];
	    dfloat sJ = mesh->sgeo[sid+SJID];
	    //
	    int id   = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
            int idM  = mesh->vmapM[id];

            const int vidM = idM%mesh->Np;
            const int qidM = e*mesh->Np*bns->Nfields + vidM;
            		    	 
	    dfloat q1  = bns->q[qidM + 0*mesh->Np];
	    dfloat q2  = bns->q[qidM + 1*mesh->Np];
	    dfloat q3  = bns->q[qidM + 2*mesh->Np];
	    dfloat q4  = bns->q[qidM + 3*mesh->Np];
	    dfloat q5  = bns->q[qidM + 4*mesh->Np];
	    dfloat q6  = bns->q[qidM + 5*mesh->Np];

	    // Compute Stress Tensor
	    dfloat s11 = -bns->RT*(sqrt(2.0)*q5 - q2*q2/q1);
	    dfloat s12 = -bns->RT*(          q4 - q2*q3/q1);
	    dfloat s22 = -bns->RT*(sqrt(2.0)*q6 - q3*q3/q1);

	    dfloat P   = q1*bns->RT; // rho*RT 
	    
            Fx  += W[n]*sJ*(-P*nx  + (s11*nx + s12*ny) );
            Fy  += W[n]*sJ*(-P*ny  + (s12*nx + s22*ny) );
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
    sprintf(fname, "BNSForceData_N%d.dat", mesh->N);

    FILE *fp; 
    fp = fopen(fname, "a");  

    fprintf(fp, "%.4e %.8e %.8e \n", time, gFx, gFy); 

    fclose(fp);
  }
}
