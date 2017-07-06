#include "acoustics2D.h"

void acousticsPmlSetup2D(mesh2D *mesh){

  //constant pml absorption coefficient
  dfloat xsigma = 40, ysigma = 40;

  // absorbtion coefficient value per element
  mesh->pmlSigmaX      = (dfloat*) calloc(mesh->pmlNelements, sizeof(dfloat));
  mesh->pmlSigmaY      = (dfloat*) calloc(mesh->pmlNelements, sizeof(dfloat));
  mesh->pmlElementList = (iint*) calloc(mesh->pmlNelements, sizeof(iint));

  // assume quiescent pml
  mesh->pmlNfields = 4;
  mesh->pmlq    = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));
  mesh->pmlrhsq = (dfloat*) calloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields, sizeof(dfloat));

  iint cnt =0;
  for (iint e=0;e<mesh->Nelements;e++) {
    int type = mesh->elementInfo[e];
    if (type==100) { //X Pml
      mesh->pmlSigmaX[cnt] = xsigma;
      mesh->pmlElementList[cnt] = e;
      cnt++;
    } else if (type==200) { //Y Pml
      mesh->pmlSigmaY[cnt] = ysigma;
      mesh->pmlElementList[cnt] = e;
      cnt++;
    } else if (type==300) { //XY Pml
      mesh->pmlSigmaX[cnt] = xsigma;
      mesh->pmlSigmaY[cnt] = ysigma;
      mesh->pmlElementList[cnt] = e;
      cnt++;
    }
  }

  printf("PML: found %d elements inside absorbing layers and %d elements outside\n",
  mesh->pmlNelements, mesh->Nelements-mesh->pmlNelements);

  // set up PML on DEVICE
  mesh->o_pmlq = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlq);
  mesh->o_pmlrhsq = mesh->device.malloc(mesh->pmlNelements*mesh->Np*mesh->pmlNfields*sizeof(dfloat), mesh->pmlrhsq);

  mesh->o_pmlSigmaX = mesh->device.malloc(mesh->pmlNelements*sizeof(dfloat), mesh->pmlSigmaX);
  mesh->o_pmlSigmaY = mesh->device.malloc(mesh->pmlNelements*sizeof(dfloat), mesh->pmlSigmaY);
  mesh->o_pmlElementList = mesh->device.malloc(mesh->pmlNelements*sizeof(iint), mesh->pmlElementList);
}