#include <math.h>
#include <stdlib.h>

#include "mesh2D.h"

void meshAcousticsPmlSetup2D(mesh2D *mesh,
			     dfloat xmin, dfloat xmax, // bounding box for non-pml sub-domain
			     dfloat ymin, dfloat ymax,
			     dfloat xsigma, dfloat ysigma){

  // absorbtion coefficient value per element
  mesh->pmlNfields = 4;
  mesh->pmlSigmaX      = (dfloat*) calloc(mesh->Nelements, sizeof(dfloat));
  mesh->pmlSigmaY      = (dfloat*) calloc(mesh->Nelements, sizeof(dfloat));
  mesh->pmlElementList = (iint*) calloc(mesh->Nelements, sizeof(iint));
  
  // find elements with center inside PML zone
  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat cx = 0, cy = 0;
    for(iint n=0;n<mesh->Nverts;++n){
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
      for(iint n=0;n<mesh->Np;++n){
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
  mesh->pmlElementList = (iint*)   realloc(mesh->pmlElementList, cnt*sizeof(iint));
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
  mesh->o_pmlElementList = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlElementList);
}

