/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "maxwell.hpp"

void maxwell_t::Setup(platform_t& _platform, mesh_t& _mesh,
                        maxwellSettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  // type of material
  materialType = (settings.compareSetting("MATERIAL TYPE", "ISOTROPIC")) ? ISOTROPIC:HETEROGENEOUS;

  if(materialType==HETEROGENEOUS){
    mesh.CubatureSetup();
    //    mesh.CubaturePhysicalNodes();
  }
  
  Nfields = (mesh.dim==3) ? 6:3;

  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd", "sum", "norm2"});

  /*setup trace halo exchange */
  traceHalo = mesh.HaloTraceSetup(Nfields);

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  // set penalty parameter
  dfloat Lambda2 = 0.5;

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(Nlocal+Nhalo);

  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 1024;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

  if (materialType==HETEROGENEOUS){
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;
    int cubMaxNp = std::max(mesh.Np, (mesh.cubNp));
    kernelInfo["defines/" "p_cubMaxNp"]= cubMaxNp;
    int NblockC = std::max(1, blockMax/cubMaxNp);
    kernelInfo["defines/" "p_NblockC"]= NblockC;
  }
  
  // set kernel name suffix
  std::string suffix = mesh.elementSuffix();
  std::string oklFilePrefix = DMAXWELL "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "maxwellVolume" + suffix + oklFileSuffix;
  kernelName = "maxwellVolume" + suffix;

  volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  // kernels from surface file
  fileName   = oklFilePrefix + "maxwellSurface" + suffix + oklFileSuffix;
  kernelName = "maxwellSurface" + suffix;

  surfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);

  // kernels from error file
  fileName   = oklFilePrefix + "maxwellError" + suffix + oklFileSuffix;
  kernelName = "maxwellError" + suffix;

  errorKernel = platform.buildKernel(fileName, kernelName,
				     kernelInfo);

  if (materialType==HETEROGENEOUS){
    // kernels from heterogeneous surface file
    fileName   = oklFilePrefix + "maxwellHeterogeneousSurface" + suffix + oklFileSuffix;
    kernelName = "maxwellHeterogeneousSurface" + suffix;
    
    heterogeneousSurfaceKernel = platform.buildKernel(fileName, kernelName,
						      kernelInfo);

    // kernels from heterogeneous project file
    fileName   = oklFilePrefix + "maxwellHeterogeneousProject" + suffix + oklFileSuffix;
    kernelName = "maxwellHeterogeneousProject" + suffix;
    
    heterogeneousProjectKernel = platform.buildKernel(fileName, kernelName,
						      kernelInfo);
    
  }
  
  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "maxwellInitialCondition2D" + oklFileSuffix;
    kernelName = "maxwellInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "maxwellInitialCondition3D" + oklFileSuffix;
    kernelName = "maxwellInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);

  if (materialType==HETEROGENEOUS){
    
    materialNfields = 2;
    
    materialNlocal = mesh.Nelements*mesh.Np*materialNfields;
    materialNhalo  = mesh.totalHaloPairs*mesh.Np*materialNfields;
    materialNtotal = materialNhalo+materialNlocal;
    
    materialHalo = mesh.HaloTraceSetup(materialNfields);
    materialCoefficients.malloc(materialNlocal+materialNhalo);

    for(dlong e=0;e<mesh.Nelements;++e){
      for(dlong n=0;n<mesh.Np;++n){
	dlong id = e*mesh.Np*materialNfields + n;
	dfloat xn = mesh.x[e*mesh.Np+n];
	dfloat yn = mesh.y[e*mesh.Np+n];
	dfloat zn = (mesh.dim==3) ? mesh.z[e*mesh.Np+n] :0;

	materialCoefficients[id] = 1.0; // mu
	materialCoefficients[id+mesh.Np] = 1.0 + 0.03*sin(xn*M_PI)*sin(yn*M_PI); // epsilon
      }
    }

    // exchange halo  (so that we can access traces of the coefficients from elements on other ranks)
    materialHalo.Exchange(materialCoefficients, 1);

    // compute integration trace of sqrt coefficient ratios
    materialUpwindWeights.malloc(mesh.Nelements*mesh.intNfp*mesh.Nfaces*4); // dfloat4 not available here

    // for each element extract traces of the coefficients and compute upwind coefficients at each interpolation node
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int f=0;f<mesh.Nfaces;++f){
	for(dlong n=0;n<mesh.intNfp;++n){
	  dfloat fmuM = 0, fmuP = 0, fepsilonM = 0, fepsilonP = 0;
	  // interpolate traces of material coefficients to flux quadrature/cubature
	  for(dlong m=0;m<mesh.Nfp;++m){
	    dlong fid  = m+f*mesh.Nfp+e*mesh.Nfp*mesh.Nfaces;
	    dlong vidM = mesh.vmapM[fid];
	    dlong vidP = mesh.vmapP[fid];
	    dlong idM  = materialNfields*mesh.Np*(vidM/mesh.Np) + vidM%mesh.Np;
	    dlong idP  = materialNfields*mesh.Np*(vidP/mesh.Np) + vidP%mesh.Np;

	    dfloat muM      = materialCoefficients[idM];
	    dfloat epsilonM = materialCoefficients[idM+mesh.Np];
	    dfloat muP      = materialCoefficients[idP];
	    dfloat epsilonP = materialCoefficients[idP+mesh.Np];
	    
            dfloat Inm = mesh.intInterp[m+n*mesh.Nfp+f*mesh.intNfp*mesh.Nfp];
	    fmuM += Inm*muM;
	    fepsilonM += Inm*epsilonM;
	    fmuP += Inm*muP;
	    fepsilonP += Inm*epsilonP;
	  }
	  dfloat ZM = sqrt(fmuM/fepsilonM);
	  dfloat ZP = sqrt(fmuP/fepsilonP);
	  dfloat YM = sqrt(fepsilonM/fmuM);
	  dfloat YP = sqrt(fepsilonP/fmuP);

	  dlong cid = 4*(e*mesh.intNfp*mesh.Nfaces + f*mesh.intNfp+n);

	  materialUpwindWeights[cid+0] = YP/(YM+YP);
	  materialUpwindWeights[cid+1] = 1./(YM+YP);
	  materialUpwindWeights[cid+2] = ZP/(ZM+ZP);
	  materialUpwindWeights[cid+3] = 1./(ZM+ZP);
	}
      }
    }

    o_materialUpwindWeights = platform.malloc<dfloat>(materialUpwindWeights);
    
    // interpolate to cubature nodes
    materialInverseWeights.malloc(mesh.Nelements*mesh.cubNp*materialNfields);
    for(dlong e=0;e<mesh.Nelements;++e){
      for(dlong n=0;n<mesh.cubNp;++n){
	dfloat cubMu = 0, cubEpsilon = 0;
	for(dlong m=0;m<mesh.Np;++m){
	  dfloat cubInm = mesh.cubInterp[n*mesh.Np+m];
	  dlong id = e*mesh.Np*materialNfields + m;
	  cubMu += cubInm*materialCoefficients[id];
	  cubEpsilon += cubInm*materialCoefficients[id+mesh.Np];
	}
	dlong cubId = e*mesh.cubNp*materialNfields + n;
	// NOTE: cubProjectT on DEVICE has cubature weights built in
	materialInverseWeights[cubId] = 1./cubMu;
	materialInverseWeights[cubId+mesh.cubNp] = 1./cubEpsilon;
      }
    }

    o_materialInverseWeights = platform.malloc<dfloat>(materialInverseWeights);
  }  
}
