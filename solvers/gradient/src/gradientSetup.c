#include "gradient.h"

gradient_t *gradientSetup(mesh_t *mesh, setupAide &options){

  gradient_t *gradient = (gradient_t*) calloc(1, sizeof(gradient_t));

  gradient->mesh = mesh;
  
  options.getArgs("MESH DIMENSION", gradient->dim);
  options.getArgs("ELEMENT TYPE", gradient->elementType);
  
  mesh->Nfields = 1;
  gradient->Nfields = mesh->Nfields;
  
  dlong Ntotal = mesh->Nelements*mesh->Np*mesh->Nfields;
  gradient->Nblock = (Ntotal+blockSize-1)/blockSize;
  
  hlong localElements = (hlong) mesh->Nelements;
  MPI_Allreduce(&localElements, &(gradient->totalElements), 1, MPI_HLONG, MPI_SUM, MPI_COMM_WORLD);

  int check;

  // compute samples of q at interpolation nodes
  mesh->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
                                sizeof(dfloat));

  gradient->q    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*mesh->Nfields,
				    sizeof(dfloat));
  
  gradient->gradientq = (dfloat*) calloc(mesh->Nelements*mesh->Np*gradient->dim,
					 sizeof(dfloat));
  
  gradient->frame = 0;


  // OCCA stuff
  occa::kernelInfo kernelInfo;
  if(gradient->dim==3)
    meshOccaSetup3D(mesh, options, kernelInfo);
  else
    meshOccaSetup2D(mesh, options, kernelInfo);

  //add boundary data to kernel info  
  gradient->o_q =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->Nfields*sizeof(dfloat), mesh->q);

  gradient->o_gradientq =
    mesh->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*mesh->dim*sizeof(dfloat));

  dfloat *plotInterp = (dfloat*) calloc(mesh->plotNp*mesh->Np, sizeof(dfloat));
  for(int n=0;n<mesh->plotNp;++n){
    for(int m=0;m<mesh->Np;++m){
      plotInterp[n+m*mesh->plotNp] = mesh->plotInterp[n*mesh->Np+m];
    }
  }
  gradient->o_plotInterp = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(dfloat), plotInterp);

  int *plotEToV = (int*) calloc(mesh->plotNp*mesh->Np, sizeof(int));
  for(int n=0;n<mesh->plotNelements;++n){
    for(int m=0;m<mesh->plotNverts;++m){
      plotEToV[n+m*mesh->plotNelements] = mesh->plotEToV[n*mesh->plotNverts+m];
    }
  }
  gradient->o_plotEToV = mesh->device.malloc(mesh->plotNp*mesh->Np*sizeof(int), plotEToV);
  
  
  if(mesh->totalHaloPairs>0){
    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      mesh->device.malloc(mesh->totalHaloPairs*mesh->Np*mesh->Nfields*sizeof(dfloat));

    // MPI send buffer
    gradient->haloBytes = mesh->totalHaloPairs*mesh->Np*gradient->Nfields*sizeof(dfloat);
    occa::memory o_sendBuffer = mesh->device.mappedAlloc(gradient->haloBytes, NULL);
    occa::memory o_recvBuffer = mesh->device.mappedAlloc(gradient->haloBytes, NULL);
    gradient->o_haloBuffer = mesh->device.malloc(gradient->haloBytes);
    gradient->sendBuffer = (dfloat*) o_sendBuffer.getMappedPointer();
    gradient->recvBuffer = (dfloat*) o_recvBuffer.getMappedPointer();

  }

#if 0
  occa::setVerboseCompilation(false);
#endif
  
  kernelInfo.addDefine("p_Nfields", mesh->Nfields);
  kernelInfo.addDefine("p_dim", mesh->dim);
  kernelInfo.addDefine("p_plotNp", mesh->plotNp);
  kernelInfo.addDefine("p_plotNelements", mesh->plotNelements);

  int plotNthreads = mymax(mesh->Np, mymax(mesh->plotNp, mesh->plotNelements));
  kernelInfo.addDefine("p_plotNthreads", plotNthreads);
  
  const dfloat p_one = 1.0, p_two = 2.0, p_half = 1./2., p_third = 1./3., p_zero = 0;

  kernelInfo.addDefine("p_two", p_two);
  kernelInfo.addDefine("p_one", p_one);
  kernelInfo.addDefine("p_half", p_half);
  kernelInfo.addDefine("p_third", p_third);
  kernelInfo.addDefine("p_zero", p_zero);
  
  int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
  kernelInfo.addDefine("p_maxNodes", maxNodes);

  int NblockV = 512/mesh->Np; // works for CUDA
  kernelInfo.addDefine("p_NblockV", NblockV);

  int NblockS = 512/maxNodes; // works for CUDA
  kernelInfo.addDefine("p_NblockS", NblockS);

  int cubMaxNodes = mymax(mesh->Np, (mesh->intNfp*mesh->Nfaces));
  kernelInfo.addDefine("p_cubMaxNodes", cubMaxNodes);

  int cubMaxNodes1 = mymax(mesh->Np, (mesh->intNfp));
  kernelInfo.addDefine("p_cubMaxNodes1", cubMaxNodes1);


  kernelInfo.addDefine("p_blockSize", blockSize);

  kernelInfo.addParserFlag("automate-add-barriers", "disabled");

  // set kernel name suffix
  char *suffix;
  
  if(gradient->elementType==TRIANGLES)
    suffix = strdup("Tri2D");
  if(gradient->elementType==QUADRILATERALS)
    suffix = strdup("Quad2D");
  if(gradient->elementType==TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(gradient->elementType==HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  for (int r=0;r<mesh->size;r++) {

    MPI_Barrier(MPI_COMM_WORLD);
    if (r==mesh->rank) {

      // kernels from volume file
      if(mesh->dim==3){
	sprintf(fileName, DHOLMES "/okl/meshIsoSurface3D.okl");
	sprintf(kernelName, "meshIsoSurface3D");
	
	gradient->isoSurfaceKernel =
	  mesh->device.buildKernelFromSource(fileName, kernelName, kernelInfo);
      }
      
      // kernels from volume file
      sprintf(fileName, DGRADIENT "/okl/gradientVolume%s.okl",
	      suffix);
      sprintf(kernelName, "gradientVolume%s", suffix);

      gradient->gradientKernel =
	mesh->device.buildKernelFromSource(fileName,
					   kernelName,
					   kernelInfo);

#if 0
      // fix this later
      mesh->haloExtractKernel =
        mesh->device.buildKernelFromSource(DHOLMES "/okl/meshHaloExtract3D.okl",
                                           "meshHaloExtract3D",
					   kernelInfo);
#endif
    }
  }

  return gradient;
}
