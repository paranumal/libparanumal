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

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./ellipticMain setupfile\n");

    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  options.getArgs("MESH FILE", fileName);
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  // set up mesh
   // set up mesh
  mesh_t *mesh;
  switch(elementType){
    case TRIANGLES:
      mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
    case QUADRILATERALS:{
      if(dim==2){
        mesh = meshSetupQuad2D((char*)fileName.c_str(), N);
      }
      else{
        dfloat radius = 1;
        options.getArgs("SPHERE RADIUS", radius);
        mesh = meshSetupQuad3D((char*)fileName.c_str(), N, radius);
      }
    break;
  }
  case TETRAHEDRA:
  mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
  mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  if(mesh->Nelements<10)
  meshPrint3D(mesh);
#if 0
  char fname[BUFSIZ];
  sprintf(fname,"meshQuad3D.vtu");
  meshVTU3D(mesh, fname);
#endif
// parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  // set up
  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  elliptic_t *elliptic = ellipticSetup(mesh, lambda, kernelInfo, options);

  if(options.compareArgs("BENCHMARK", "BK5") ||
     options.compareArgs("BENCHMARK", "BP5")){

    // test Ax throughput
    occa::streamTag startAx = mesh->device.tagStream();

    int NAx = 1;

    for(int it=0;it<NAx;++it){
      // include gather-scatter
      if(options.compareArgs("BENCHMARK", "BP5"))
        ellipticOperator(elliptic, lambda, elliptic->o_x, elliptic->o_Ax, dfloatString); // standard precision

      if(options.compareArgs("BENCHMARK", "BK5")){
        if(!options.compareArgs("ELEMENT MAP", "TRILINEAR")){
          elliptic->partialAxKernel(mesh->NlocalGatherElements,
                                    mesh->o_localGatherElementList,
                                    mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM,
                                    lambda, elliptic->o_x, elliptic->o_Ax);
        }
        else{
          elliptic->partialAxKernel(mesh->NlocalGatherElements,
                                    mesh->o_localGatherElementList,
                                    elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM,
                                    lambda, elliptic->o_x, elliptic->o_Ax);
        }
      }
    }

    occa::streamTag stopAx = mesh->device.tagStream();

    mesh->device.finish();

    double elapsedAx = mesh->device.timeBetween(startAx, stopAx);
    elapsedAx /= NAx;


    if (mesh->rank==0)
      printf("%d, %d, %g, %d, %g, %g; \%\%elemental: N, dofs, elapsed, dummy, time per node, nodes/time %s\n",
           mesh->N,
           mesh->NlocalGatherElements*mesh->Np,
           elapsedAx,
           0,
           elapsedAx/(mesh->Np*mesh->Nelements),
           mesh->Nelements*mesh->Np/elapsedAx,
           (char*) options.getArgs("DISCRETIZATION").c_str());

  }
  else{

    // convergence tolerance
    dfloat tol = 1e-8;

    occa::streamTag startTag = mesh->device.tagStream();

    int it = ellipticSolve(elliptic, lambda, tol, elliptic->o_r, elliptic->o_x);
    occa::streamTag stopTag = mesh->device.tagStream();
    mesh->device.finish();

    double elapsed = mesh->device.timeBetween(startTag, stopTag);

    if (mesh->rank==0)
      printf("%d, %d, %g, %d, %g, %g; \%\%global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
           mesh->N,
           mesh->Nelements*mesh->Np,
           elapsed,
           it,
           elapsed/(mesh->Np*mesh->Nelements),
           mesh->Nelements*(it*mesh->Np/elapsed),
           (char*) options.getArgs("PRECONDITIONER").c_str());

    if(options.compareArgs("DISCRETIZATION","CONTINUOUS") && 
       !(elliptic->dim==3 && elliptic->elementType==QUADRILATERALS)){
      dfloat zero = 0.;
      elliptic->addBCKernel(mesh->Nelements,
                            zero,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            elliptic->o_mapB,
                            elliptic->o_x);
    }

    // copy solution from DEVICE to HOST
    elliptic->o_x.copyTo(mesh->q);

    if (options.compareArgs("BASIS","BERN"))
      meshApplyElementMatrix(mesh,mesh->VB,mesh->q,mesh->q);

    dfloat maxError = 0;
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        dlong   id = e*mesh->Np+n;
        dfloat xn = mesh->x[id];
        dfloat yn = mesh->y[id];
        dfloat zn = mesh->z[id];

        dfloat exact;
        if (elliptic->dim==2)
          exact = sin(M_PI*xn)*sin(M_PI*yn);
        else{
          if(elliptic->elementType==QUADRILATERALS){
#if 0
	    exact = xn*xn;
#endif

#if 0
	    exact = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
#endif
	    dfloat a = 1, b = 2, c = 3;
	    exact = sin(a*xn)*sin(b*yn)*sin(c*zn);
	  }
	  else
	    exact = cos(M_PI*xn)*cos(M_PI*yn)*cos(M_PI*zn);
        }

        dfloat error = fabs(exact-mesh->q[id]);

	// store error
	mesh->q[id] = fabs(mesh->q[id] - exact);
        maxError = mymax(maxError, error);
      }
    }

    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    if(mesh->rank==0)
      printf("globalMaxError = %g\n", globalMaxError);

#if 1
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d.vtu",(char*)outName.c_str(), mesh->rank);
    if(elliptic->dim==3)
      meshPlotVTU3D(mesh, fname, 0);
    else
      meshPlotVTU2D(mesh, fname, 0);
#endif
  }

#if 0
  {
    dfloat maxError = 0;
    for(int n=0;n<mesh->Np*mesh->Nelements;++n){
      dfloat sc = M_PI;
      dfloat exact = cos(sc*mesh->x[n])*cos(sc*mesh->y[n])*cos(sc*mesh->z[n]);
      dfloat error =  exact - mesh->q[n];
      maxError = mymax(maxError, fabs(error));
    }
    printf("maxError = %g\n", maxError);

    printf("PRINTING VTU\n");
    ellipticPlotVTUHex3D(mesh, "bah", 0);
  }
#endif

  // close down MPI
  MPI_Finalize();

  return 0;
}
