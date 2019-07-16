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
    if(elliptic->elementType == TRIANGLES)
      meshOccaSetupTri3D(mesh, options, kernelInfo);
    else if(elliptic->elementType == QUADRILATERALS)
      meshOccaSetupQuad3D(mesh, options, kernelInfo);
    else
      meshOccaSetup3D(mesh, options, kernelInfo);
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

    // Set Timer
  elliptic->profiler = new timer; 
  elliptic->profiler->setTimer(elliptic->options);
  elliptic->profiler->initTimer(mesh->device);

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

#if 0
	  dfloat exact = pow(xn,2);
	  dfloat forcing = -2*(- 2*pow(xn,2) + pow(yn,2) + pow(zn,2));
#endif
	  dfloat a = 1, b = 2, c = 3;
	  dfloat pi = M_PI;

#if 0
	  dfloat exact = sin(pi*xn)*sin(pi*yn)*sin(pi*zn);
	  dfloat forcing =
	    - 2*pi*pi*sin(pi*xn)*sin(pi*yn)*sin(pi*zn)
	    - 2*xn*pi*cos(pi*xn)*sin(pi*yn)*sin(pi*zn)
	    - 2*yn*pi*cos(pi*yn)*sin(pi*xn)*sin(pi*zn)
	    - 2*zn*pi*cos(pi*zn)*sin(pi*xn)*sin(pi*yn)
	    - 2*xn*yn*pi*pi*cos(pi*xn)*cos(pi*yn)*sin(pi*zn)
	    - 2*xn*zn*pi*pi*cos(pi*xn)*cos(pi*zn)*sin(pi*yn)
	    - 2*yn*zn*pi*pi*cos(pi*yn)*cos(pi*zn)*sin(pi*xn);
#endif

	  dfloat exact = sin(a*xn)*sin(b*yn)*sin(c*zn);
	  
	  dfloat forcing = 
	      b*b*yn*yn*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    - c*c*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    - a*a*yn*yn*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    - a*a*zn*zn*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    - b*b*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    + c*c*zn*zn*sin(a*xn)*sin(b*yn)*sin(c*zn)
	    - 2*a*xn*cos(a*xn)*sin(b*yn)*sin(c*zn)
	    - 2*b*yn*cos(b*yn)*sin(a*xn)*sin(c*zn)
	    - 2*c*zn*cos(c*zn)*sin(a*xn)*sin(b*yn)
	    - 2*a*c*xn*zn*cos(a*xn)*cos(c*zn)*sin(b*yn)
	    - 2*b*c*yn*zn*cos(b*yn)*cos(c*zn)*sin(a*xn)
	    - 2*a*b*xn*yn*cos(a*xn)*cos(b*yn)*sin(c*zn);
	  
	  
	  forcing = -forcing + lambda*exact;
	  
          elliptic->r[id] = J*forcing; 

        }
        else{
	  int mode = 1;
	  elliptic->r[id] =
	    J*(3*mode*mode*M_PI*M_PI+lambda)*cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);
	  //	  elliptic->r[id] += 0.1*2*(drand48()-0.5);
	}

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
  if (options.compareArgs("BASIS","NODAL")){
    if(options.compareArgs("ELLIPTIC INTEGRATION", "NODAL")){
      printf("MASS APPLY NODAL\n");
      meshApplyElementMatrix(mesh,mesh->MM,elliptic->r,elliptic->r);
    }
    else{
      printf("MASS APPLY CUBATURE\n");
      dfloat *cubx = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
      dfloat *cuby = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
      dfloat *cubz = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
      dfloat *cubrhs = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
      
      dfloat *cubInterpT = (dfloat*) calloc(mesh->cubNq*mesh->Nq, sizeof(dfloat));
      for(int n=0;n<mesh->Nq;++n){
	for(int m=0;m<mesh->cubNq;++m){        
	  cubInterpT[m+n*mesh->cubNq] = mesh->cubInterp[m*mesh->Nq+n];
	  printf("%g ", cubInterpT[m+n*mesh->cubNq]);
	}
	printf("\n");
      }
      
      for(hlong e=0;e<mesh->Nelements;++e){
	
	interpolateHex3D(mesh->cubInterp, mesh->x+mesh->Np*e, mesh->Nq, cubx, mesh->cubNq);
	interpolateHex3D(mesh->cubInterp, mesh->y+mesh->Np*e, mesh->Nq, cuby, mesh->cubNq);
	interpolateHex3D(mesh->cubInterp, mesh->z+mesh->Np*e, mesh->Nq, cubz, mesh->cubNq);
	
	for(int n=0;n<mesh->cubNp;++n){
	  dfloat JW = mesh->cubggeo[e*mesh->cubNp*mesh->Nggeo + n + GWJID*mesh->cubNp];
	  //	  cubrhs[n] = JW*(3*M_PI*M_PI+lambda)*cos(M_PI*cubx[n])*cos(M_PI*cuby[n])*cos(M_PI*cubz[n]);
	  int mode = 1;
	  cubrhs[n] = JW*(3*mode*mode*M_PI*M_PI+lambda)*cos(mode*M_PI*cubx[n])*cos(mode*M_PI*cuby[n])*cos(mode*M_PI*cubz[n]);
	  //	  cubrhs[n] += 0.1*2*(drand48()-0.5);
	}
	
	interpolateHex3D(cubInterpT, cubrhs, mesh->cubNq, elliptic->r+e*mesh->Np, mesh->Nq);

	//	for(int n=0;n<mesh->Np;++n){
	//	  printf("elliptic->r[%d]=%g\n", e*mesh->Np+n, elliptic->r[e*mesh->Np+n]);
	//	}	  
      }
    }	
  }	

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

    dfloat zero = 0.f, mone = -1.0f, one = 1.0f;
    if(options.compareArgs("ELLIPTIC INTEGRATION", "NODAL")){
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
    }else{

      // first set Dirichlet bc in tmp field
      hlong Etotal = (mesh->Nelements+mesh->totalHaloPairs);
      dfloat *qbc = (dfloat*) calloc(Etotal*mesh->Np, sizeof(dfloat));
      // note the zeroing
      occa::memory o_qbc = mesh->device.malloc(Etotal*mesh->Np*sizeof(dfloat), qbc);
      occa::memory o_rbc = mesh->device.malloc(Etotal*mesh->Np*sizeof(dfloat), qbc);
      
      elliptic->addBCKernel(mesh->Nelements,
			    zero, // time
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    elliptic->o_mapB,
			    o_qbc);

      // A*qbc
      elliptic->partialCubatureAxKernel(mesh->NlocalGatherElements,
					mesh->o_localGatherElementList,
					mesh->o_cubggeo,
					mesh->o_cubD,
					mesh->o_cubInterpT,
					lambda,
					o_qbc,
					o_rbc);
      // r -= A*q_bc
      elliptic->scaledAddKernel(Etotal*mesh->Np, mone, o_rbc, one, elliptic->o_r);

      // add Neumann fluxes later
    }
    
    
  }

  // gather-scatter
 if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ogsGatherScatter(elliptic->o_r, ogsDfloat, ogsAdd, mesh->ogs);
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_r);
  }

  return elliptic;
}
