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

#include "lss.hpp"

void lss_t::SetupStabilizer(){

  if(advection){
    if(mesh.rank==0)
      printf("Level-Set Advection doesn't need stabilization mechanism: Skipping setup \n");
  }else{
    if(subcellStabilization){
      subcell = NULL;
      subcell = &(subcell_t::Setup(*this));

     
    // Storage required for subcell stabilization
    const dlong sNlocal = mesh.Nelements*subcell->Nsubcells; 
    // this->sq = (dfloat*) calloc(sNlocal, sizeof(dfloat));
    // this->o_sq = mesh.device.malloc(sNlocal*sizeof(dfloat), this->sq);
    
    this->ssgnq = (dfloat*) calloc(sNlocal, sizeof(dfloat));
    this->o_ssgnq = mesh.device.malloc(sNlocal*sizeof(dfloat), this->ssgnq);
    
    this->sface = (dfloat *) calloc(sNlocal*mesh.Nfaces, sizeof(dfloat)); 
    this->o_sface = mesh.device.malloc(sNlocal*mesh.Nfaces*sizeof(dfloat), this->sface); 

      char *suffix;
      if(mesh.elementType==TRIANGLES)
      suffix = strdup("Tri2D");
      if(mesh.elementType==QUADRILATERALS)
      suffix = strdup("Quad2D");
      if(mesh.elementType==TETRAHEDRA)
      suffix = strdup("Tet3D");
      if(mesh.elementType==HEXAHEDRA)
      suffix = strdup("Hex3D");

      char fileName[BUFSIZ], kernelName[BUFSIZ];
  
    occa::properties kernelInfo = props; // copy base props


    sprintf(fileName, DLSS "/okl/lssStabilize%s.okl", suffix);

    sprintf(kernelName, "skyline%s", suffix);
    skylineKernel = buildKernel(device, fileName, kernelName, kernelInfo, comm);


    sprintf(kernelName, "lssProject%s", suffix);
    projectKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 

    sprintf(kernelName, "lssReconstruct%s", suffix);
    reconstructKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 

    sprintf(kernelName, "lssPartialRedistanceVolume%s", suffix);
    partialRedistanceVolumeKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 

    sprintf(kernelName, "lssReconstructInternalFace%s", suffix);
    reconstructInternalFaceKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 


    sprintf(kernelName, "lssProjectDG%s", suffix);
    projectDGKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 


    sprintf(kernelName, "lssReconstructExternalFace%s", suffix);
    reconstructExternalFaceKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm);

    sprintf(kernelName, "lssPartialRedistanceSurface%s", suffix);
    partialRedistanceSurfaceKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 

    // sprintf(kernelName, "lssSubcellReconstructFace%s", suffix);
    // subcellReconstructFaceKernel = buildKernel(mesh.device, fileName, kernelName, 
    //                                     kernelInfo, mesh.comm); 

    sprintf(kernelName, "lssSubcellCompute%s", suffix);
    subcellComputeKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 

    
    sprintf(fileName, DLSS "/okl/lssSubcellRegularizedSign2D.okl");
    sprintf(kernelName, "lssSubcellSign2D");
    subcellSignKernel = buildKernel(mesh.device, fileName, kernelName, 
                                     kernelInfo, mesh.comm); 


  }else{
    if(mesh.rank==0 ){
      stringstream ss;
      ss << "Currently fv-subcell stabilization is supported: exitting...";
      LIBP_ABORT(ss.str());
    }
  }
}

}
