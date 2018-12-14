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

#include "cds.h"


void cdsError(cds_t *cds, dfloat time){

  mesh_t *mesh = cds->mesh;

  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  dfloat maxS = 0, minS = 1e9;
 
  if (cds->options.compareArgs("EXACT", "NONE")) {// just compute maximums
 
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        int id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        dfloat z = mesh->z[id];

        maxS = mymax(maxS, fabs(cds->S[id+0*cds->sOffset]));
        minS = mymin(minS, fabs(cds->S[id+0*cds->sOffset]));
      }
    }

    // compute maximum over all processes
    dfloat gMaxS, gMinS;
    MPI_Allreduce(&maxS, &gMaxS, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minS, &gMinS, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    if(mesh->rank==0){
      int cstep = (int)((time-cds->startTime)/cds->dt)+1; 
      printf("Step: %d Time: %g minS: %g maxS: %g\n",cstep, time, gMinS, gMaxS);
    }

    if( isnan(gMinS) || isnan(gMaxS) )  exit(EXIT_FAILURE);
  } else { //compare to an exact solution

    if (cds->options.compareArgs("EXACT","2D EXACT")) { 
      // NADA 
    } else if (cds->options.compareArgs("EXACT","UCD")) { //3D Unsteady Convection Diffusion

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          dlong id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];
          dfloat z = mesh->z[id];

	  dfloat ax = 0.01;
	  dfloat ay = 0.01;
	  dfloat az = 0.01; 

	  dfloat cx = 0.8;
	  dfloat cy = 0.8;
	  dfloat cz = 0.8; 

          dfloat alf   = 1.0/pow((4.0*time +1.0),3.0/2.0);
          dfloat xterm = pow((x - cx*time - 0.5),2.0) / ( ax*(4.0*time+1.0) );
	  dfloat yterm = pow((y - cy*time - 0.5),2.0) / ( ay*(4.0*time+1.0) );
	  dfloat zterm = pow((z - cz*time - 0.5),2.0) / ( az*(4.0*time+1.0) );


          dfloat sExact = alf*exp(-xterm -yterm - zterm);
          maxS = mymax(maxS, fabs(cds->S[id+0*cds->sOffset]-sExact));
#if 1
	  cds->S[id+0*cds->sOffset] -= sExact;
#endif
        }
      }

      // compute maximum over all processes
      dfloat gMaxS;
      MPI_Allreduce(&maxS, &gMaxS, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      if(mesh->rank==0)
        printf("Step: %d Time: %g ErrorS: %g\n", (int)(time/cds->dt), time, gMaxS);
      if( isnan(gMaxS)) exit(EXIT_FAILURE);
    }
  }
}
