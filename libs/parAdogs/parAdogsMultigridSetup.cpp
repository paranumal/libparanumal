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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"

namespace libp {

namespace paradogs {

/****************************************/
/* Construct Multigrid Hierarchy        */
/****************************************/
void graph_t::MultigridSetup() {

  CreateLaplacian();

  /*Target size for coarsest graph*/
  const int coarseSize = 100;

  /*Coarsening tolerance. If a coarse graph isn't at least
    coarseTol times smaller than a fine graph, we consider the coarsening
    to be stalling*/
  const float coarseTol = 0.8;

  /*Stength threashold*/
  dfloat theta=0.08;

  do{
    /*Get coarsest level*/
    mgLevel_t& Lf = L[Nlevels-1];

    /*If the graph is small enough, we're done*/
    if (Lf.Nglobal <= coarseSize) {
      coarseSolver.Setup(Lf.A, Lf.null);
      break;
    }

    LIBP_ABORT("Paradogs: Max levels exceeded in coarse graph creation. Increase MAX_LEVELS.",
               Nlevels>=MAX_LEVELS);

    Lf.SetupSmoother();

    /*Construct next level via coarsening*/
    mgLevel_t& Lc = L[Nlevels];    
    Lc.CoarsenLevel(Lf, theta);
    Nlevels++;
    
    // Increase coarsening rate as we add levels.
    //See: Algebraic Multigrid On Unstructured Meshes, P Vanek, J. Mandel, M. Brezina.
    theta=theta/2;

    /*Check for stalls*/
    if (Lc.Nglobal > coarseTol*Lf.Nglobal) {
      LIBP_FORCE_WARNING("Paradogs: Graph coarsening stalling. Coarse graph has " << Lc.Nglobal << " nodes.");
      coarseSolver.Setup(Lc.A, Lc.null);
      break;
    }
  } while(true);

  for (int l=0;l<Nlevels;++l) {
    L[l].AllocateScratch(l);
  }
}

void mgLevel_t::SetupSmoother() {

  // estimate rho(invD * A)
  A.rho = A.rhoDinvA(null);

  /*Smoothing params*/
  lambda1 = A.rho;
  lambda0 = A.rho/10.;

}

void mgLevel_t::AllocateScratch(const int l) {

  /*Space for Fiedler*/
  Fiedler.malloc(Ncols);

  RES.malloc(Ncols);

  if (l>0) {
    /*Multigrid buffers*/
    RHS.malloc(Nrows);
    X.malloc(Ncols);
  }

  /*Scratch space*/
  scratch.malloc(2*Ncols);
}



/*Coarsen a graph using an aggregation*/
void mgLevel_t::CoarsenLevel(mgLevel_t& Lf, const dfloat theta) {

  /*Create a FineToCoarse mapping*/
  const dlong Nf = Lf.Nrows;

  /*Create a vertex matching*/
  dlong Nc=0;
  memory<hlong> FineToCoarse(Lf.Ncols);
  Lf.A.Aggregate(Nc, theta, FineToCoarse);

  /* Tentative prolongation operator*/
  parCSR T = TentativeProlongator(Nf, Nc,
                                  Lf.A.platform, Lf.A.comm,
                                  FineToCoarse,
                                  Lf.null, null);
  FineToCoarse.free();

  /* Smoothed prologontion */
  Lf.P = SmoothProlongator(Lf.A, T);
  T = parCSR(); //Free T

  /* R = P^T*/
  Lf.R = Transpose(Lf.P);
  Lf.Ncols = std::max(Lf.Ncols, Lf.R.Ncols);

  /*Galerkin product*/
  parCSR AP = SpMM(Lf.A, Lf.P);
  A = SpMM(Lf.R, AP);
  // A.GalerkinProduct(Lf.A, Lf.P);
  AP= parCSR(); //Free AP

  /*fill diagonal*/
  A.diagA.malloc(A.Ncols);
  A.diagInv.malloc(A.Nrows);

  #pragma omp parallel for
  for (dlong i=0;i<A.Nrows;i++) {
    const dlong start = A.diag.rowStarts[i];
    const dlong end   = A.diag.rowStarts[i+1];

    for (dlong j=start;j<end;j++) {
      //record the diagonal
      if (A.diag.cols[j]==i) {
        A.diagA[i] = A.diag.vals[j];
        A.diagInv[i] = 1.0/A.diagA[i];
        break;
      }
    }
  }

  //fill the halo region
  A.halo.Exchange(A.diagA, 1);

  Nrows = A.Nrows;
  Ncols = std::max(A.Ncols, Lf.P.Ncols);

  Nglobal = static_cast<hlong>(Nrows);
  A.comm.Allreduce(Nglobal);
}

/*Free coarse levels of hierarchy*/
void graph_t::MultigridDestroy() {
  colIds.free();
  coarseSolver = coarseSolver_t();
  for (int n=Nlevels-1;n>=0;--n) L[n] = mgLevel_t();
  Nlevels=0;
}

} //namespace paradogs

} //namespace libp
