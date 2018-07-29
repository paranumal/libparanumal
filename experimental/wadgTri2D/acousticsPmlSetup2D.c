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

#include "acoustics2D.h"

void acousticsPmlSetup2D(mesh2D *mesh,
			 dfloat xmin, dfloat xmax, // bounding box for non-pml sub-domain
			 dfloat ymin, dfloat ymax,
			 dfloat xsigma, dfloat ysigma){

  // absorbtion coefficient value per element
  mesh->pmlNfields = 4;
  mesh->pmlSigmaX      = (dfloat*) calloc(mesh->Nelements, sizeof(dfloat));
  mesh->pmlSigmaY      = (dfloat*) calloc(mesh->Nelements, sizeof(dfloat));
  mesh->pmlElementList = (int*) calloc(mesh->Nelements, sizeof(int));
  
  // find elements with center inside PML zone
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(int n=0;n<mesh->Nverts;++n){
      cx += mesh->EX[e*mesh->Nverts+n];
      cy += mesh->EY[e*mesh->Nverts+n];
    }
    cx /= mesh->Nverts;
    cy /= mesh->Nverts;

    // add element outside [xmin,xmax]x[ymin,ymax] to pml
    if(cx>xmax || cx<xmin || cy>ymax || cy<ymin){
      mesh->pmlElementList[cnt] = e;
      if(cx<xmin || cx>xmax)
	mesh->pmlSigmaX[cnt] = xsigma;
      if(cy<xmin || cy>xmax)
	mesh->pmlSigmaY[cnt] = ysigma;
#if 0
      for(int n=0;n<mesh->Np;++n){
	dfloat x = mesh->x[n + e*mesh->Np];
	dfloat y = mesh->y[n + e*mesh->Np];
	if(cx>xmax)
	  mesh->pmlSigmaX[mesh->Np*cnt + n] = xsigma*pow(x-xmax,2);
	if(cx<xmin)
	  mesh->pmlSigmaX[mesh->Np*cnt + n] = xsigma*pow(x-xmin,2);
	if(cy>ymax)
	  mesh->pmlSigmaY[mesh->Np*cnt + n] = ysigma*pow(y-ymax,2);
	if(cy<ymin)
	  mesh->pmlSigmaY[mesh->Np*cnt + n] = ysigma*pow(y-ymin,2);
      }
#endif
      
      ++cnt;
    }
  }
  mesh->pmlNelements = cnt;
  mesh->pmlElementList = (int*)   realloc(mesh->pmlElementList, cnt*sizeof(int));
  mesh->pmlSigmaX      = (dfloat*) realloc(mesh->pmlSigmaX,      cnt*sizeof(dfloat));
  mesh->pmlSigmaY      = (dfloat*) realloc(mesh->pmlSigmaY,      cnt*sizeof(dfloat));

  // assume quiescent pml
  mesh->pmlq    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
  mesh->pmlrhsq = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
  mesh->pmlresq = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
  
  printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
	 mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

  // set up PML on DEVICE
  mesh->o_pmlq = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlq);
  mesh->o_pmlrhsq = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsq);
  mesh->o_pmlresq = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlresq);
  mesh->o_pmlSigmaX = mesh->device.malloc(mesh->pmlNelements*sizeof(dfloat), mesh->pmlSigmaX);
  mesh->o_pmlSigmaY = mesh->device.malloc(mesh->pmlNelements*sizeof(dfloat), mesh->pmlSigmaY);
  mesh->o_pmlElementList = mesh->device.malloc(mesh->pmlNelements*sizeof(int), mesh->pmlElementList);
}

