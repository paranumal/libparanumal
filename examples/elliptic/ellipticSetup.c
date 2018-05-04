#include "elliptic.h"
#include "omp.h"
#include <unistd.h>

elliptic_t *ellipticSetup(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, setupAide options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  int totalDevices = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
  }
  for (int r=0;r<size;r++) {
    if (hostIds[r]==hostId) totalDevices++;
  }

  // read thread model/device/platform from options
  if(options.compareArgs("THREAD MODEL", "CUDA")){
    int dev;
    options.getArgs("DEVICE NUMBER" ,dev);
    sprintf(deviceConfig, "mode = CUDA, deviceID = %d",dev);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenCL")){
    int dev, plat;
    options.getArgs("DEVICE NUMBER", dev);
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode = OpenCL, deviceID = %d, platformID = %d", dev, plat);
  }
  else if(options.compareArgs("THREAD MODEL", "OpenMP")){
    sprintf(deviceConfig, "mode = OpenMP");
  }
  else{
    sprintf(deviceConfig, "mode = Serial");
  }
   
  //set number of omp threads to use
  int Ncores = sysconf(_SC_NPROCESSORS_ONLN);
  int Nthreads = Ncores/totalDevices;
  omp_set_num_threads(Nthreads);
  
  if (rank==0 && options.compareArgs("VERBOSE","TRUE")) 
    printf("Rank %d: Ncores = %d, Nthreads = %d\n", rank, Ncores, Nthreads);

  elliptic_t *elliptic = (elliptic_t*) calloc(1, sizeof(elliptic_t));

  options.getArgs("MESH DIMENSION", elliptic->dim);
  options.getArgs("ELEMENT TYPE", elliptic->elementType);

  elliptic->mesh = mesh;
  elliptic->options = options;

  mesh->Nfields = 1;

  // compute samples of q at interpolation nodes
  mesh->q = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));

  if(elliptic->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  // allocate unified derivative matrices
  if(elliptic->elementType == HEXAHEDRA || elliptic->elementType == QUADRILATERALS){
    elliptic->o_Dmatrices = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D);
    elliptic->o_Smatrices = mesh->device.malloc(mesh->Nq*mesh->Nq*sizeof(dfloat), mesh->D); //dummy
  } else if(elliptic->elementType == TRIANGLES){
    // build Dr, Ds transposes
    dfloat *DrsT = (dfloat*) calloc(2*mesh->Np*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        DrsT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
        DrsT[n+m*mesh->Np+mesh->Np*mesh->Np] = mesh->Ds[n*mesh->Np+m];
      }
    }

    dfloat *ST = (dfloat*) calloc(3*mesh->Np*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        ST[n+m*mesh->Np+0*mesh->Np*mesh->Np] = mesh->Srr[n*mesh->Np+m];
        ST[n+m*mesh->Np+1*mesh->Np*mesh->Np] = mesh->Srs[n*mesh->Np+m]+mesh->Ssr[n*mesh->Np+m];
        ST[n+m*mesh->Np+2*mesh->Np*mesh->Np] = mesh->Sss[n*mesh->Np+m];
      }
    }
   
    elliptic->o_Dmatrices = mesh->device.malloc(2*mesh->Np*mesh->Np*sizeof(dfloat), DrsT);
    elliptic->o_Smatrices = mesh->device.malloc(3*mesh->Np*mesh->Np*sizeof(dfloat), ST);
    free(DrsT); free(ST);
  } else {
    // build Dr, Ds transposes
    dfloat *DrstT = (dfloat*) calloc(3*mesh->Np*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        DrstT[n+m*mesh->Np] = mesh->Dr[n*mesh->Np+m];
        DrstT[n+m*mesh->Np+mesh->Np*mesh->Np] = mesh->Ds[n*mesh->Np+m];
        DrstT[n+m*mesh->Np+2*mesh->Np*mesh->Np] = mesh->Dt[n*mesh->Np+m];
      }
    }

    dfloat *ST = (dfloat*) calloc(6*mesh->Np*mesh->Np, sizeof(dfloat));
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        ST[n+m*mesh->Np+0*mesh->Np*mesh->Np] = mesh->Srr[n*mesh->Np+m];
        ST[n+m*mesh->Np+1*mesh->Np*mesh->Np] = mesh->Srs[n*mesh->Np+m]+mesh->Ssr[n*mesh->Np+m];
        ST[n+m*mesh->Np+2*mesh->Np*mesh->Np] = mesh->Srt[n*mesh->Np+m]+mesh->Str[n*mesh->Np+m];
        ST[n+m*mesh->Np+3*mesh->Np*mesh->Np] = mesh->Sss[n*mesh->Np+m];
        ST[n+m*mesh->Np+4*mesh->Np*mesh->Np] = mesh->Sst[n*mesh->Np+m]+mesh->Sts[n*mesh->Np+m];
        ST[n+m*mesh->Np+5*mesh->Np*mesh->Np] = mesh->Stt[n*mesh->Np+m];
      }
    }

    elliptic->o_Dmatrices = mesh->device.malloc(3*mesh->Np*mesh->Np*sizeof(dfloat), DrstT);
    elliptic->o_Smatrices = mesh->device.malloc(6*mesh->Np*mesh->Np*sizeof(dfloat), ST);
    free(DrstT); free(ST);
  }

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};
  elliptic->BCType = (int*) calloc(3,sizeof(int));
  memcpy(elliptic->BCType,BCType,3*sizeof(int));

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

      if(elliptic->dim==2)
        elliptic->r[id] = J*(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
      else 
        elliptic->r[id] = J*(3*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
      elliptic->x[id] = 0;
    }
  }

  //Apply some element matrix ops to r depending on our solver
  if (options.compareArgs("BASIS","BERN"))   meshApplyElementMatrix(mesh,mesh->invVB,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","BERN"))   meshApplyElementMatrix(mesh,mesh->BBMM,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","NODAL"))  meshApplyElementMatrix(mesh,mesh->MM,elliptic->r,elliptic->r);

  //copy to occa buffers
  elliptic->o_r   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->r);
  elliptic->o_x   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->x);


  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo.addInclude((char*)boundaryHeaderFileName.c_str());

  // set kernel name suffix
  char *suffix;
  
  if(elliptic->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(elliptic->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(elliptic->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  //add boundary condition contribution to rhs
  if (options.compareArgs("DISCRETIZATION","IPDG")) {

    sprintf(fileName, "okl/ellipticRhsBCIpdg%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBCIpdg%s", suffix);

    elliptic->rhsBCIpdgKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

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
                            elliptic->o_Dmatrices,
                            mesh->o_LIFTT,
                            mesh->o_MM,
                            elliptic->o_r);
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {

    sprintf(fileName, "okl/ellipticRhsBC%s.okl", suffix);
    sprintf(kernelName, "ellipticRhsBC%s", suffix);

    elliptic->rhsBCKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

    sprintf(fileName, "okl/ellipticAddBC%s.okl", suffix);
    sprintf(kernelName, "ellipticAddBC%s", suffix);

    elliptic->addBCKernel = mesh->device.buildKernelFromSource(fileName,kernelName, kernelInfo);

    dfloat zero = 0.f;
    elliptic->rhsBCKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        mesh->o_sgeo,
                        elliptic->o_Dmatrices,
                        elliptic->o_Smatrices,
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
    ellipticParallelGatherScatter(mesh, mesh->ogs, elliptic->o_r, dfloatString, "add");  
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_r);
  }

  return elliptic;
}
