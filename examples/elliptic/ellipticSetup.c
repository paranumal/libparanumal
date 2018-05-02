#include "elliptic.h"

elliptic_t *ellipticSetup(mesh_t *mesh, occa::kernelInfo &kernelInfo, setupAide options){

  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  long int hostId = gethostid();

  long int* hostIds = (long int*) calloc(size,sizeof(long int));
  MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,MPI_COMM_WORLD);

  int deviceID = 0;
  for (int r=0;r<rank;r++) {
    if (hostIds[r]==hostId) deviceID++;
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
        
  elliptic_t *elliptic = (elliptic_t*) calloc(1, sizeof(elliptic_t));

  options.getArgs("MESH DIMENSION", elliptic->dim);
  options.getArgs("ELEMENT TYPE", elliptic->elementType);

  elliptic->options = options;

  mesh->Nfields = 1;

  // compute samples of q at interpolation nodes
  mesh->q = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields, sizeof(dfloat));

  if(elliptic->dim==3)
    meshOccaSetup3D(mesh, deviceConfig, kernelInfo);
  else
    meshOccaSetup2D(mesh, deviceConfig, kernelInfo);

  // Boundary Type translation. Just default from the mesh file.
  int BCType[3] = {0,1,2};
  elliptic->BCType = (int*) calloc(3*sizeof(int));
  memcpy(elliptic->BCType,BCType,3*sizeof(int));

  //tau
  elliptic->tau = 2.0*(mesh->N+1)*(mesh->N+2)/2.0;
  
  dlong Nall = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);
  elliptic->r   = (dfloat*) calloc(Nall,   sizeof(dfloat));
  elliptic->x   = (dfloat*) calloc(Nall,   sizeof(dfloat));

  // load forcing into r
  for(dlong e=0;e<mesh->Nelements;++e){
    dfloat J = mesh->vgeo[e*mesh->Nvgeo+JID];
    for(int n=0;n<mesh->Np;++n){
      dlong id = n+e*mesh->Np;
      dfloat xn = mesh->x[id];
      dfloat yn = mesh->y[id];
      elliptic->r[id] = J*(2*M_PI*M_PI+lambda)*sin(M_PI*xn)*sin(M_PI*yn);
      elliptic->x[id] = 0;
    }
  }

  //Apply some element matrix ops to r depending on our solver
  if (options.compareArgs("BASIS","BERN"))   applyElementMatrix(mesh,mesh->invVB,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","BERN"))   applyElementMatrix(mesh,mesh->BBMM,elliptic->r,elliptic->r);
  if (options.compareArgs("BASIS","NODAL"))  applyElementMatrix(mesh,mesh->MM,elliptic->r,elliptic->r);

  //copy to occa buffers
  elliptic->o_r   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->r);
  elliptic->o_x   = mesh->device.malloc(Nall*sizeof(dfloat), elliptic->x);


  ellipticSolveSetupTri2D(elliptic, kernelInfo, options);

  string boundaryHeaderFileName; 
  options.getArgs("DATA FILE", boundaryHeaderFileName);
  kernelInfo.addInclude((char*)boundaryHeaderFileName.c_str());

  // set kernel name suffix
  char *suffix;
  
  if(cns->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(cns->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(cns->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(cns->elementType==HEXAHEDRA)
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
                            mesh->o_vmapP,
                            elliptic->tau,
                            zero,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vgeo,
                            mesh->o_sgeo,
                            elliptic->o_EToB,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            mesh->o_DtT,
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
                        mesh->o_SrrT,
                        mesh->o_SrsT,
                        mesh->o_SsrT,
                        mesh->o_SssT,
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
    ellipticParallelGatherScatterTri2D(mesh, mesh->ogs, o_r, dfloatString, "add");  
    if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, elliptic->o_r);
  }

  return elliptic;
}
