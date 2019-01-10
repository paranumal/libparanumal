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
#include "cuda.h"
#include "cudaProfiler.h"
#include "ins.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);

  // set up mesh stuff

  int dim, elementType;

  string fileName = argv[2];
  int    N        = atoi(argv[3]);
  int    maxiter  = atoi(argv[4]);

  int NTEST      = 1; // number of test for 
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);
 
  // set up mesh
  mesh_t *mesh;
  switch(elementType){
  case TRIANGLES:
    mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
  case QUADRILATERALS:{
    if(dim==2)
      mesh = meshSetupQuad2D((char*)fileName.c_str(), N);
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

  ins_t *ins = insSetup(mesh,options);

  if (!ins->options.compareArgs("TIME INTEGRATOR", "EXTBDF")){ 
    printf("currently only EXTBDF.....\n");
    exit(EXIT_FAILURE);
  }
    
  // if (ins->options.compareArgs("TIME INTEGRATOR", "ARK"))  insRunARK(ins);


  double start = 0.0, end =0.0;
  mesh->device.finish();
  MPI_Barrier(mesh->comm);
  
  start = MPI_Wtime();
  
  cuProfilerStart();
  insRunEXTBDFTest(ins, maxiter);
  cuProfilerStop();


  mesh->device.finish();
  end = MPI_Wtime();
  double localElapsed = end-start;


  int size = mesh->size;

  hlong   localDofs     = (hlong) mesh->Np*mesh->Nelements;
  hlong   localElements = (hlong) mesh->Nelements;
  double  globalElapsed;
  hlong   globalDofs;
  hlong   globalElements;

  MPI_Reduce(&localElapsed, &globalElapsed, 1, MPI_DOUBLE,  MPI_MAX, 0, mesh->comm );
  MPI_Reduce(&localDofs,    &globalDofs,    1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );
  MPI_Reduce(&localElements,&globalElements,1, MPI_HLONG,   MPI_SUM, 0, mesh->comm );

  globalElapsed /= maxiter;
  int ppn = 1;
  if (mesh->rank==0){
#if 0
    char fname[BUFSIZ];
    sprintf(fname,"BP5_Scaling_M%02d.dat", mesh->size);
    FILE *fp;
    fp = fopen(fname, "a");

    fprintf(fp, "%02d %02d "hlongFormat" "hlongFormat" %d %17.15lg %g %g \n",
      mesh->size, mesh->N, globalElements, globalDofs, maxIter, globalElapsed, globalElapsed/(globalDofs),1.0/(globalElapsed/(globalDofs)));

    fprintf(fp, "lipP,CUDA,BP5,%02d,%02d,%d,"hlongFormat",%17.15lg\n",mesh->N, mesh->size, ppn, globalElements,globalElapsed);
    fclose(fp);
#else
    printf("%02d %02d "hlongFormat" "hlongFormat" %d %17.15lg %g %g\t [ RANKS N NELEMENTS DOFS ITERATIONS ELAPSEDTIME TIME_PER_DOF, DOF_PER_TIME] \n",
     mesh->size, mesh->N, globalElements, globalDofs, maxiter, globalElapsed, globalElapsed/(globalDofs),
     1.0/(globalElapsed/(globalDofs)));
#endif
  }







  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
