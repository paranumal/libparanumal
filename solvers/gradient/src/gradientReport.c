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

#include "gradient.h"

void gradientReport(gradient_t *gradient, dfloat time, setupAide &options){

  mesh3D *mesh = gradient->mesh;

  int isoField = 3;
  int isoNlevels = 4;
  int isoMaxNtris = 1E8;
  
  dfloat *isoLevels = (dfloat*) calloc(isoNlevels, sizeof(dfloat));
  int *isoNtris = (int*) calloc(1, sizeof(int));

  size_t isoMax = (mesh->dim + gradient->Nfields)*3*isoMaxNtris;
  dfloat *isoq = (dfloat*) calloc(isoMax, sizeof(dfloat));

  dfloat isoMinVal = -.8;
  dfloat isoMaxVal = .8;
  
  for(int l=0;l<isoNlevels;++l)
    isoLevels[l] = isoMinVal + (isoMaxVal-isoMinVal)*l/(dfloat)(isoNlevels-1);

  occa::memory o_isoLevels = mesh->device.malloc(isoNlevels*sizeof(dfloat), isoLevels);
  occa::memory o_isoq      = mesh->device.malloc(isoMax*sizeof(dfloat), isoq);

  // note that this should be zero before calling isoSurface kernel
  occa::memory o_isoNtris  = mesh->device.malloc(1*sizeof(int), isoNtris);
    
						 
  gradient->isoSurfaceKernel(mesh->Nelements,    // number of elements
			     isoField,     // which field to use for isosurfacing
			     isoNlevels,   // number of isosurface levels
			     o_isoLevels,  // array of isosurface levels
			     isoMaxNtris,  // maximum number of generated triangles
			     mesh->o_x,
			     mesh->o_y,
			     mesh->o_z,
			     gradient->o_q,
			     gradient->o_plotInterp,
			     gradient->o_plotEToV,
			     o_isoNtris,  // output: number of generated triangles
			     o_isoq       // output: (p_dim+p_Nfields)*3*isoNtris[0] values (x,y,z,q0,q1..)
			     );
  

  // find number of generated triangles
  o_isoNtris.copyTo(isoNtris);
  isoNtris[0] = mymin(isoNtris[0], isoMaxNtris);
  
  printf("generated %d triangles\n", isoNtris[0]);
  
  //
  int offset = 0;
  o_isoq.copyTo(isoq, isoNtris[0]*(mesh->dim+gradient->Nfields)*3*sizeof(dfloat), offset);

#if 0
  for(int n=0;n<isoNtris[0];++n){
    for(int v=0;v<3;++v){
      printf("[%d,%d]: [", n, v);
      for(int fld=0;fld<(mesh->dim+gradient->Nfields);++fld){
	int id = (n*3+v)*(mesh->dim+gradient->Nfields) + fld;
	printf(" %f", isoq[id]);
      }
      printf("]\n");
    }
  }
#endif			  
  
#if 1
  // output field files
  char fname[BUFSIZ];
  string outName;
  options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, gradient->frame++);

  gradientPlotVTU(gradient, isoNtris[0], isoq, fname);
#endif

}
