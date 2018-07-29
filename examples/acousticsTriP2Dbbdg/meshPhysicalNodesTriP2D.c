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

#include <stdio.h>
#include <stdlib.h>
#include "mesh2D.h"

void meshPhysicalNodesTriP2D(mesh2D *mesh){
  
  mesh->x = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  mesh->y = (dfloat*) calloc(mesh->Nelements*mesh->NpMax,sizeof(dfloat));
  
  for(int e=0;e<mesh->Nelements;++e){ /* for each element */
    int N = mesh->N[e];

    int id = e*mesh->Nverts+0;

    dfloat xe1 = mesh->EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = mesh->EX[id+1];
    dfloat xe3 = mesh->EX[id+2];

    dfloat ye1 = mesh->EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = mesh->EY[id+1];
    dfloat ye3 = mesh->EY[id+2];
    
    for(int n=0;n<mesh->Np[N];++n){ /* for each node */
      
      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = mesh->r[N][n]; 
      dfloat sn = mesh->s[N][n];

      /* physical coordinate of interpolation node */
      mesh->x[e*mesh->NpMax+n] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      mesh->y[e*mesh->NpMax+n] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
    }
  }

  // create halo extension for x,y arrays
  int totalHaloNodes = mesh->totalHaloPairs*mesh->NpMax;
  int localNodes     = mesh->Nelements*mesh->NpMax;

  // temporary send buffer
  dfloat *sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  mesh->x = (dfloat*) realloc(mesh->x, (localNodes+totalHaloNodes)*sizeof(dfloat));
  mesh->y = (dfloat*) realloc(mesh->y, (localNodes+totalHaloNodes)*sizeof(dfloat));
  
  // send halo data and recv into extended part of arrays
  meshHaloExchange(mesh, mesh->NpMax*sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
  meshHaloExchange(mesh, mesh->NpMax*sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);

  // grab EX,EY,EZ from halo
  mesh->EX = (dfloat*) realloc(mesh->EX, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts*sizeof(dfloat));
  mesh->EY = (dfloat*) realloc(mesh->EY, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts*sizeof(dfloat));

  // send halo data and recv into extended part of arrays
  meshHaloExchange(mesh, mesh->Nverts*sizeof(dfloat), mesh->EX, sendBuffer, mesh->EX + mesh->Nverts*mesh->Nelements);
  meshHaloExchange(mesh, mesh->Nverts*sizeof(dfloat), mesh->EY, sendBuffer, mesh->EY + mesh->Nverts*mesh->Nelements);

  if (totalHaloNodes) free(sendBuffer);
}
