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

#include "adaptive.h"

/** Initialize the libraries that we are using
 *
 */
static int init_libs(MPI_Comm comm, int verbosity)
{
  int rank;
  ASD_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  int logpriority = ASD_MAX(SC_LP_STATISTICS - verbosity, SC_LP_ALWAYS);
  sc_init(comm, 0, 0, NULL, logpriority);
  p4est_init(NULL, logpriority);

  int loglevel = ASD_MAX(ASD_LL_INFO - verbosity, ASD_LL_ALWAYS);
  asd_log_init(rank, stdout, loglevel);

  // add signal handler to get backtrace on abort
  asd_signal_handler_set();

  return loglevel;
}

static void print_precision()
{
  const char *comp = (sizeof(double) == sizeof(dfloat_t)) ? "double" : "single";
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("----- Precision ------------------------------------------");
  ASD_ROOT_INFO("compute precision = %s", comp);
  ASD_ROOT_INFO("----------------------------------------------------------");
}

int main(int argc, char **argv){
  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./adaptiveMain setupfile\n");

    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  int verbosity = 0;
  options.getArgs("VERBOSITY", verbosity);
  init_libs(comm, verbosity);
  print_precision();

  app_t *app = app_new(options, comm);

  app_free(app);

#if 0
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);
  if(dim!=3 || elementType!=HEXAHEDRA){ printf("EXITING: HEX3D ONLY\n"); MPI_Finalize(); exit(0); }

  mesh_t *mesh;

  // set up mesh
  if(options.getArgs("MESH FILE", fileName)){
    mesh = meshSetup((char*) fileName.c_str(), N, options);
  }
  else if(options.compareArgs("BOX DOMAIN", "TRUE")){
    mesh = adaptiveSetupBoxHex3D(N, options);
  }


  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  // set up
  occa::properties kernelInfo;
  kernelInfo["defines"].asObject();
  kernelInfo["includes"].asArray();
  kernelInfo["header"].asArray();
  kernelInfo["flags"].asObject();

  // configure mesh on OCCA Device
  meshOccaSetup3D(mesh, options, kernelInfo);
  
  adaptive_t *adaptive = adaptiveSetup(mesh, lambda, kernelInfo, options);
  
  {    
    occa::memory o_r = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), adaptive->o_r);
    occa::memory o_x = mesh->device.malloc(mesh->Np*mesh->Nelements*sizeof(dfloat), adaptive->o_x);    

    // convergence tolerance
    dfloat tol = 1e-8;

    // warm up
    //    int it = adaptiveSolve(adaptive, lambda, tol, adaptive->o_r, adaptive->o_x);
    int it;

    MPI_Barrier(mesh->comm);
    
    occa::streamTag startTag = mesh->device.tagStream();
    int Ntests = 1;
    it = 0;
    for(int test=0;test<Ntests;++test){
      o_r.copyTo(adaptive->o_r);
      o_x.copyTo(adaptive->o_x);
      it += adaptiveSolve(adaptive, lambda, tol, adaptive->o_r, adaptive->o_x);
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
      printf("%d, %d, %g, %d, %g, %g; \%\%global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
	     mesh->N,
	     globalNelements*mesh->Np,
	     globalElapsed,
	     it,
	     globalElapsed/(mesh->Np*globalNelements),
	     globalNelements*(it*mesh->Np/globalElapsed),
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
    
    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      dfloat zero = 0.;
      adaptive->addBCKernel(mesh->Nelements,
                            zero,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            adaptive->o_mapB,
                            adaptive->o_x);
    }

    // copy solution from DEVICE to HOST
    adaptive->o_x.copyTo(mesh->q);

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
	double mode = 1.0;
	exact = cos(mode*M_PI*xn)*cos(mode*M_PI*yn)*cos(mode*M_PI*zn);

        dfloat error = fabs(exact-mesh->q[id]);

        maxError = mymax(maxError, error);
      }
    }

    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    if(mesh->rank==0)
      printf("globalMaxError = %g\n", globalMaxError);

    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    // original
    adaptive->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d",(char*)outName.c_str(), mesh->rank);
    adaptivePlotVTUHex3D(mesh, fname, 0);
  }
#endif

  // close down MPI
  MPI_Finalize();

  return 0;
}
