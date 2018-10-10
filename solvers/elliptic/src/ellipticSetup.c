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

#include "elliptic.h"
#include "omp.h"
#include <unistd.h>

void reportMemoryUsage(occa::device &device, const char *mess);

elliptic_t *ellipticSetup(mesh_t *mesh, dfloat lambda, occa::properties &kernelInfo, setupAide options){

  elliptic_t *elliptic = (elliptic_t*) calloc(1, sizeof(elliptic_t));

  options.getArgs("MESH DIMENSION", elliptic->dim);
  options.getArgs("ELEMENT TYPE", elliptic->elementType);

  elliptic->mesh = mesh;
  elliptic->options = options;

  mesh->Nfields = 1;

  // compute samples of q at interpolation nodes
  mesh->q = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));

  if(elliptic->dim==3){
    if(elliptic->elementType != QUADRILATERALS)
      meshOccaSetup3D(mesh, options, kernelInfo);
    else
      meshOccaSetupQuad3D(mesh, options, kernelInfo); 
  } 
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  if (mesh->rank==0)
    reportMemoryUsage(mesh->device, "after occa setup");

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};
  elliptic->BCType = (int*) calloc(3,sizeof(int));
  memcpy(elliptic->BCType,BCType,3*sizeof(int));

  // build trilinear geometric factors for hexes (do before solve setup)
  if(elliptic->elementType==HEXAHEDRA){
    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      if(options.compareArgs("ELEMENT MAP", "TRILINEAR")){
	printf("mesh->dim = %d, mesh->Nverts = %d\n", mesh->dim, mesh->Nverts);

	// pack gllz, gllw, and elementwise EXYZ
	hlong Nxyz = mesh->Nelements*mesh->dim*mesh->Nverts;
	dfloat *EXYZ = (dfloat*) calloc(Nxyz, sizeof(dfloat));
	dfloat *gllzw = (dfloat*) calloc(2*mesh->Nq, sizeof(dfloat));

	int sk = 0;
	for(int n=0;n<mesh->Nq;++n)
	  gllzw[sk++] = mesh->gllz[n];
	for(int n=0;n<mesh->Nq;++n)
	  gllzw[sk++] = mesh->gllw[n];

	sk = 0;
	for(hlong e=0;e<mesh->Nelements;++e){
	  for(int v=0;v<mesh->Nverts;++v)
	    EXYZ[sk++] = mesh->EX[e*mesh->Nverts+v];
	  for(int v=0;v<mesh->Nverts;++v)
	    EXYZ[sk++] = mesh->EY[e*mesh->Nverts+v];
	  for(int v=0;v<mesh->Nverts;++v)
	    EXYZ[sk++] = mesh->EZ[e*mesh->Nverts+v];
	}

	// nodewise ggeo with element coordinates and gauss node info
	elliptic->o_EXYZ = mesh->device.malloc(Nxyz*sizeof(dfloat), EXYZ);
	elliptic->o_gllzw = mesh->device.malloc(2*mesh->Nq*sizeof(dfloat), gllzw);

	free(EXYZ);
	free(gllzw);
      }
    }
  }

  //
  ellipticSolveSetup(elliptic, lambda, kernelInfo);


  dlong Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  elliptic->r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load forcing into r
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){

      dfloat J;
      if (elliptic->elementType==TRIANGLES || elliptic->elementType==TETRAHEDRA) {
        J = mesh->vgeo[e*mesh->Nvgeo+JID];
      } else {
        J = mesh->vgeo[mesh->Np*(e*mesh->Nvgeo + JID) + n];
      }
      dlong id = n+e*mesh->Np;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      dfloat zn = mesh->z[id];
      
      dfloat forcing; 
      if(elliptic->dim==2)
        elliptic->r[id] = J*(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
      else{
        if(elliptic->elementType==QUADRILATERALS){
          // dfloat alpha = 1.0;
          // dfloat beta  = 3.0;

          // // dfloat rad = sqrt(xn*xn + yn*yn + zn*zn); // has to be one !!!

          // dfloat theta = acos(zn);
          // dfloat phi   = atan2(xn,yn);
          // forcing      = J*sin(alpha*phi)*sin(beta*theta)*(alpha*alpha/(sin(theta)*sin(theta))  
          //                                          +beta*beta-beta*cos(theta)*cos(beta*theta)/(sin(theta)*sin(beta*theta)));

          // if(isnan(forcing))
          //    forcing =J*sin(alpha*phi)*sin(beta*theta); 

          forcing = -(2*(pow(xn,3) - pow(xn,2)*yn - 2*pow(xn,2)*zn - 2*xn*pow(yn,2) + xn*pow(zn,2) - pow(yn,3) + pow(yn,2)*zn + 2*yn*pow(zn,2) + pow(zn,3)))/(pow(xn,2) + pow(yn,2) + pow(zn,2));
          // Added for making well-posed lamda
          forcing += lambda*(xn*pow(yn,2) - yn*pow(zn,2) + pow(xn,2)*zn);
          elliptic->r[id] = J*forcing; 

        }
        else
        elliptic->r[id] = J*(3*M_PI*M_PI+lambda)*cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);

      }
      elliptic->x[id] = 0;
    }
  }


#if 0
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "AAA%s_%04d.vtu",(char*)outName.c_str(), mesh->rank);
    if(elliptic->dim==3)
      meshPlotVTU3D(mesh, fname, 0);
    else
      meshPlotVTU2D(mesh, fname, 0);
#endif

  //Apply some element matrix ops to r depending on our solver
  if (options.compareArgs("BASIS","BERN"))   meshApplyElementMatrix(mesh,mesh->invVB,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","BERN"))   meshApplyElementMatrix(mesh,mesh->BBMM,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","NODAL"))  meshApplyElementMatrix(mesh,mesh->MM,elliptic->r,elliptic->r);

  //copy to occa buffers
  elliptic->o_r   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->r);
  elliptic->o_x   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->x);


  string boundaryHeaderFileName;
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo["includes"] += (char*)boundaryHeaderFileName.c_str();

  // set kernel name suffix
  char *suffix;

  if(elliptic->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(elliptic->elementType==QUADRILATERALS){
    if(elliptic->dim==2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D"); 
  }
  if(elliptic->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  //add boundary condition contribution to rhs
  if (options.compareArgs("DISCRETIZATION","IPDG") && 
      !(elliptic->dim==3 && elliptic->elementType==QUADRILATERALS) ) {
    for(int r=0;r<mesh->size;++r){
      if(r==mesh->rank){
	sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBCIpdg%s.okl", suffix);
	sprintf(kernelName, "ellipticRhsBCIpdg%s", suffix);

	elliptic->rhsBCIpdgKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);
      }
      MPI_Barrier(mesh->comm);
    }
    dfloat zero = 0.f;
    elliptic->rhsBCIpdgKernel(mesh->Nelements,
			      mesh->o_vmapM,
			      elliptic->tau,
			      zero,
			      mesh->o_x,
			      mesh->o_y,
			      mesh->o_z,
			      mesh->o_vgeo,
			      mesh->o_sgeo,
			      elliptic->o_EToB,
			      mesh->o_Dmatrices,
			      mesh->o_LIFTT,
			      mesh->o_MM,
			      elliptic->o_r);
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS") &&
       !(elliptic->dim==3 && elliptic->elementType==QUADRILATERALS) ) {
    for(int r=0;r<mesh->size;++r){
      if(r==mesh->rank){
      sprintf(fileName, DELLIPTIC "/okl/ellipticRhsBC%s.okl", suffix);
        sprintf(kernelName, "ellipticRhsBC%s", suffix);

        elliptic->rhsBCKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticAddBC%s.okl", suffix);
        sprintf(kernelName, "ellipticAddBC%s", suffix);

        elliptic->addBCKernel = mesh->device.buildKernel(fileName,kernelName, kernelInfo);
      }
      MPI_Barrier(mesh->comm);
    }

    dfloat zero = 0.f;
    elliptic->rhsBCKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        mesh->o_sgeo,
                        mesh->o_Dmatrices,
                        mesh->o_Smatrices,
                        mesh->o_MM,
                        mesh->o_vmapM,
                        mesh->o_sMT,
                        lambda,
                        zero,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_z,
                        elliptic->o_mapB,
                        elliptic->o_r);
  }

  // gather-scatter
 if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ogsGatherScatter(elliptic->o_r, ogsDfloat, ogsAdd, mesh->ogs);
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_r);
  }

  return elliptic;
}
