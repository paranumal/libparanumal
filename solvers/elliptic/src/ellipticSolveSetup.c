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


void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::properties &kernelInfo){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;




  for (int r=0;r<2;r++){
    if ((r==0 && mesh->rank==0) || (r==1 && mesh->rank>0)) {

      mesh->maskKernel =
        mesh->device.buildKernel(LIBP_DIR "/okl/mask.okl",
                   "mask",
                   kernelInfo);


      kernelInfo["defines/" "p_blockSize"]= blockSize;

      // add custom defines
      kernelInfo["defines/" "p_NpP"]= (mesh->Np+mesh->Nfp*mesh->Nfaces);
      kernelInfo["defines/" "p_Nverts"]= mesh->Nverts;

      //sizes for the coarsen and prolongation kernels. degree N to degree 1
      kernelInfo["defines/" "p_NpFine"]= mesh->Np;
      kernelInfo["defines/" "p_NpCoarse"]= mesh->Nverts;


      if (elliptic->elementType==QUADRILATERALS || elliptic->elementType==HEXAHEDRA) {
        kernelInfo["defines/" "p_NqFine"]= mesh->N+1;
        kernelInfo["defines/" "p_NqCoarse"]= 2;
      }

      kernelInfo["defines/" "p_NpFEM"]= mesh->NpFEM;

      int Nmax = mymax(mesh->Np, mesh->Nfaces*mesh->Nfp);
      kernelInfo["defines/" "p_Nmax"]= Nmax;

      int maxNodes = mymax(mesh->Np, (mesh->Nfp*mesh->Nfaces));
      kernelInfo["defines/" "p_maxNodes"]= maxNodes;

      int NblockV = mymax(1,maxNthreads/mesh->Np); // works for CUDA
      int NnodesV = 1; //hard coded for now
      kernelInfo["defines/" "p_NblockV"]= NblockV;
      kernelInfo["defines/" "p_NnodesV"]= NnodesV;
      kernelInfo["defines/" "p_NblockVFine"]= NblockV;
      kernelInfo["defines/" "p_NblockVCoarse"]= NblockV;

      int NblockS = mymax(1,maxNthreads/maxNodes); // works for CUDA
      kernelInfo["defines/" "p_NblockS"]= NblockS;

      int NblockP = mymax(1,maxNthreads/(4*mesh->Np)); // get close to maxNthreads threads
      kernelInfo["defines/" "p_NblockP"]= NblockP;

      int NblockG;
      if(mesh->Np<=32) NblockG = ( 32/mesh->Np );
      else NblockG = maxNthreads/mesh->Np;
      kernelInfo["defines/" "p_NblockG"]= NblockG;

      kernelInfo["defines/" "p_halfC"]= (int)((mesh->cubNq+1)/2);
      kernelInfo["defines/" "p_halfN"]= (int)((mesh->Nq+1)/2);

      kernelInfo["defines/" "p_NthreadsUpdatePCG"] = (int) NthreadsUpdatePCG; // WARNING SHOULD BE MULTIPLE OF 32
      kernelInfo["defines/" "p_NwarpsUpdatePCG"] = (int) (NthreadsUpdatePCG/32); // WARNING: CUDA SPECIFIC

      //      cout << kernelInfo ;

      //add standard boundary functions
      char *boundaryHeaderFileName;
      if (elliptic->dim==2)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary2D.h");
      else if (elliptic->dim==3)
        boundaryHeaderFileName = strdup(DELLIPTIC "/data/ellipticBoundary3D.h");
      kernelInfo["includes"] += boundaryHeaderFileName;


      sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
      sprintf(kernelName, "ellipticAx%s", suffix);

      occa::properties dfloatKernelInfo = kernelInfo;
      occa::properties floatKernelInfo = kernelInfo;
      floatKernelInfo["defines/" "pfloat"]= "float";
      dfloatKernelInfo["defines/" "pfloat"]= dfloatString;

      elliptic->AxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);

      if(elliptic->elementType!=HEXAHEDRA){
        sprintf(kernelName, "ellipticPartialAx%s", suffix);
      }
      else{
        if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR")){
          sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
        }else{
          sprintf(kernelName, "ellipticPartialAx%s", suffix);
        }
      }

      elliptic->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      elliptic->partialAxKernel2 = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      elliptic->partialFloatAxKernel = mesh->device.buildKernel(fileName,kernelName,floatKernelInfo);

      // only for Hex3D - cubature Ax
      if(elliptic->elementType==HEXAHEDRA){
	//	printf("BUILDING partialCubatureAxKernel\n");
	sprintf(fileName,  DELLIPTIC "/okl/ellipticCubatureAx%s.okl", suffix);

	sprintf(kernelName, "ellipticCubaturePartialAx%s", suffix);
	elliptic->partialCubatureAxKernel = mesh->device.buildKernel(fileName,kernelName,dfloatKernelInfo);
      }




      // combined update for Non-blocking PCG
      elliptic->update1NBPCGKernel =
	mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
				 "ellipticUpdate1NBPCG", dfloatKernelInfo);

      elliptic->update2NBPCGKernel =
	mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
				 "ellipticUpdate2NBPCG", dfloatKernelInfo);


      // combined update for Non-blocking flexible PCG
      elliptic->update0NBFPCGKernel =
	mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
				 "ellipticUpdate0NBFPCG", dfloatKernelInfo);

      elliptic->update1NBFPCGKernel =
	mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
				 "ellipticUpdate1NBFPCG", dfloatKernelInfo);


        sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
        sprintf(kernelName, "ellipticGradient%s", suffix);

        elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialGradient%s", suffix);
        elliptic->partialGradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
        sprintf(kernelName, "ellipticAxIpdg%s", suffix);
        elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
        elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);


      // Use the same kernel with quads for the following kenels
      if(elliptic->dim==3){
	if(elliptic->elementType==QUADRILATERALS)
	  suffix = strdup("Quad2D");
	else if(elliptic->elementType==TRIANGLES)
	  suffix = strdup("Tri2D");
      }

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
      elliptic->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
      sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
      elliptic->precon->prolongateKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);



      sprintf(fileName, DELLIPTIC "/okl/ellipticBlockJacobiPrecon.okl");
      sprintf(kernelName, "ellipticBlockJacobiPrecon");
      elliptic->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
      elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      sprintf(fileName, DELLIPTIC "/okl/ellipticPatchSolver.okl");
      sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
      elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

      if (   elliptic->elementType == TRIANGLES
          || elliptic->elementType == TETRAHEDRA) {
        elliptic->precon->SEMFEMInterpKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                     "ellipticSEMFEMInterp",
                     kernelInfo);

        elliptic->precon->SEMFEMAnterpKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                     "ellipticSEMFEMAnterp",
                     kernelInfo);
      }
    }

    MPI_Barrier(mesh->comm);
  }

  // TW: WARNING C0 appropriate only
  mesh->sumKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, elliptic->o_tmp);
  elliptic->o_tmp.copyTo(elliptic->tmp);

  dfloat nullProjectWeightLocal = 0;
  dfloat nullProjectWeightGlobal = 0;
  for(dlong n=0;n<elliptic->Nblock;++n)
    nullProjectWeightLocal += elliptic->tmp[n];

  MPI_Allreduce(&nullProjectWeightLocal, &nullProjectWeightGlobal, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

  elliptic->nullProjectWeightGlobal = 1./nullProjectWeightGlobal;



  long long int pre = mesh->device.memoryAllocated();

  ellipticPreconditionerSetup(elliptic, elliptic->ogs, lambda, kernelInfo);

  long long int usedBytes = mesh->device.memoryAllocated()-pre;

  elliptic->precon->preconBytes = usedBytes;

}
