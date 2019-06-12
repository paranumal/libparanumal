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

#include "elliptic.hpp"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(argc!=2)
    LIBP_ABORT(string("Usage: ./ellipticMain setupfile"));

  ellipticSettings_t settings(comm); //sets default settings
  settings.readSettingsFromFile(argv[1]);
  if (!rank) settings.report();

  // set up occa device
  occa::device device;
  occa::properties props;
  occaDeviceConfig(device, comm, settings, props);

  // set up mesh
  mesh_t& mesh = mesh_t::Setup(device, comm, settings, props);

  // set up elliptic solver
  elliptic_t& elliptic = elliptic_t::Setup(mesh);

  // run
  elliptic.Run();

  // close down MPI
  MPI_Finalize();
  return LIBP_SUCCESS;

  {
    occa::memory o_r = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), elliptic->o_r);
    occa::memory o_x = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), elliptic->o_x);

    // convergence tolerance
    dfloat tol = 1e-8;

    // warm up
    //    int it = ellipticSolve(elliptic, lambda, tol, elliptic->o_r, elliptic->o_x);
    int it;

    MPI_Barrier(mesh->comm);

    occa::streamTag startTag = mesh->device.tagStream();
    int Ntests = 1;
    it = 0;
    for(int test=0;test<Ntests;++test){
      o_r.copyTo(elliptic->o_r);
      o_x.copyTo(elliptic->o_x);
      it += ellipticSolve(elliptic, lambda, tol, elliptic->o_r, elliptic->o_x);
    }

    MPI_Barrier(mesh->comm);

    occa::streamTag stopTag = mesh->device.tagStream();
    mesh->device.finish();

    double elapsed = mesh->device.timeBetween(startTag, stopTag);

    double globalElapsed;
    hlong globalNelements;

    MPI_Reduce(&elapsed, &globalElapsed, 1, MPI_DOUBLE, MPI_MAX, 0, mesh->comm);
    MPI_Reduce(&(mesh->Nelements), &globalNelements, 1, MPI_HLONG, MPI_SUM, 0, mesh->comm);

    if (mesh->rank==0)
      printf("%d, " hlongFormat ", %g, %d, %g, %g; \%\%global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
	     mesh->N,
	     globalNelements*mesh->Np,
	     globalElapsed,
	     it,
	     globalElapsed/(mesh->Np*globalNelements),
	     globalNelements*((dfloat)it*mesh->Np/globalElapsed),
	     (char*) options.getArgs("PRECONDITIONER").c_str());

    if (options.compareArgs("VERBOSE", "TRUE")){
      fflush(stdout);
      MPI_Barrier(mesh->comm);
      printf("rank %d has %d internal elements and %d non-internal elements\n",
	     mesh->rank,
	     mesh->NinternalElements,
	     mesh->NnotInternalElements);
      MPI_Barrier(mesh->comm);
    }

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
          exact = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
//           if(elliptic->elementType==QUADRILATERALS){
// #if 0
// 	    exact = xn*xn;
// #endif

// #if 0
// 	    exact = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
// #endif
// 	    dfloat a = 1, b = 2, c = 3;
// 	    exact = sin(a*xn)*sin(b*yn)*sin(c*zn);
// 	  }
// 	  else{
// 	    double mode = 1.0;
// 	    exact = cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);
// 	  }
        }

        dfloat error = fabs(exact-mesh->q[id]);

	//	mesh->q[id] -= exact;

        // store error
        // mesh->q[id] = fabs(mesh->q[id] - exact);
        maxError = mymax(maxError, error);
      }
    }

    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    if(mesh->rank==0)
      printf("globalMaxError = %g\n", globalMaxError);

#if 0
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    // original
    elliptic->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d",(char*)outName.c_str(), mesh->rank);
    ellipticPlotVTUHex3D(mesh, fname, 0);
#endif

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

  //  cout << kernelInfo;

  // build one-ring ( to rule them all )
  //  ellipticBuildOneRing(elliptic, kernelInfo);

  // close down MPI
  MPI_Finalize();

  return 0;
}
