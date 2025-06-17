
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

#include "ins.hpp"
#include "stress.hpp"


//  Solves gamma*U - nu*Laplacian*U = rhs
void ins_t::StressSolve(deviceMemory<dfloat>& o_U, deviceMemory<dfloat>& o_RHS,
			const dfloat gamma, const dfloat T) {

  dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*NVfields;

  deviceMemory<dfloat> o_UH = platform.reserve<dfloat>(Ntotal);
  deviceMemory<dfloat> o_rhsU = platform.reserve<dfloat>(Ntotal);

  // compute RHS = MM*RHS/nu + BCdata
  // and split fields to separate arrays
  dfloat nuInv = 1./nu;
  stressRhsKernel(mesh.Nelements,
		  mesh.o_wJ,
		  mesh.o_vgeo,
		  mesh.o_sgeo,
		  mesh.o_ggeo,
		  mesh.o_S,
		  mesh.o_D,
		  mesh.o_LIFT,
		  mesh.o_MM,
		  mesh.o_sM,
		  mesh.o_vmapM,
		  mesh.o_EToB,
		  mesh.o_mapB,
		  vTau,
		  T,
		  mesh.o_x,
		  mesh.o_y,
		  mesh.o_z,
		  gamma/nu,
		  nuInv,
		  stressSolver.o_nut,
		  o_U,
		  o_RHS,
		  o_UH,
		  o_rhsU);
  
  int maxIter = 5000;
  int verbose = 0;

  stressSolver.lambda = gamma/nu;

  //  Solve lambda*U - Laplacian*U = rhs
  if (vDisc_c0){

    //    printf("SOLVING U *********************\n");
    // gather, solve, scatter
    deviceMemory<dfloat> o_GrhsU = platform.reserve<dfloat>(NVfields*(uSolver.Ndofs+uSolver.Nhalo));
    deviceMemory<dfloat> o_GUH   = platform.reserve<dfloat>(NVfields*(uSolver.Ndofs+uSolver.Nhalo));

    stressSolver.ogsMasked.Gather(o_GrhsU, o_rhsU, NVfields, ogs::Add, ogs::Trans);
    NiterU = stressSolver.Solve(stressLinearSolver, o_GUH, o_GrhsU, velTOL, maxIter, verbose);
    stressSolver.ogsMasked.Scatter(o_UH, o_GUH, NVfields, ogs::NoTrans);

    o_GUH.free(); o_GrhsU.free();

    
  } else {
    std::cout << " Not implemented " << std::endl;
    exit(-1);
  }

  // merge arrays back, and enter BCs if C0
  stressBCKernel(mesh.Nelements,
		 mesh.o_sgeo,
		 mesh.o_vmapM,
		 mesh.o_mapB,
		 T,
		 mesh.o_x,
		 mesh.o_y,
		 mesh.o_z,
		 nu,
		 vDisc_c0,
		 o_UH,
		 o_U);
}
