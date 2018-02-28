#include "acoustics3D.h"

void acousticsError3D(mesh3D *mesh, dfloat time){

  iint NMax = mesh->NMax;

  dfloat maxErrorP = 0;
  for(iint e=0;e<mesh->Nelements;++e){    
    iint id = e*mesh->Nverts;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];
    dfloat xe4 = mesh->EX[id+3];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    dfloat ye4 = mesh->EY[id+3];
    
    dfloat ze1 = mesh->EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = mesh->EZ[id+1];
    dfloat ze3 = mesh->EZ[id+2];
    dfloat ze4 = mesh->EZ[id+3];

    for(iint n=0;n<mesh->NpMax;++n){
      /* (r,s,t) coordinates of plot nodes*/
      dfloat rn = mesh->r[NMax][n]; 
      dfloat sn = mesh->s[NMax][n];
      dfloat tn = mesh->t[NMax][n];

      /* physical coordinate of interpolation node */
      dfloat x = -0.5*(1+rn+sn+tn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3 + 0.5*(1+tn)*xe4;
      dfloat y = -0.5*(1+rn+sn+tn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3 + 0.5*(1+tn)*ye4;
      dfloat z = -0.5*(1+rn+sn+tn)*ze1 + 0.5*(1+rn)*ze2 + 0.5*(1+sn)*ze3 + 0.5*(1+tn)*ze4;
      
      dfloat u,v,w,p;
      acousticsCavitySolution3D(x, y, z, time, &u, &v, &w, &p);

      dfloat qm = mesh->q[(e*mesh->NpMax+n)*mesh->Nfields+3];

      maxErrorP = mymax(maxErrorP, fabs(p-qm));
    }
  }

  // compute maximum over all processes
  dfloat globalMaxErrorP;
  MPI_Allreduce(&maxErrorP, &globalMaxErrorP, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0)
    printf("%g, %g (time,maxError(pressure)\n", time, globalMaxErrorP);
  
}
