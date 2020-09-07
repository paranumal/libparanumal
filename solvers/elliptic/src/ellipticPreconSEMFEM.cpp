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

#include "ellipticPrecon.hpp"

// Cast problem into spectrally-equivalent N=1 FEM space and precondition with AMG
void SEMFEMPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  if (mesh.elementType==TRIANGLES) {

    dlong Ntotal = mesh.Np*mesh.Nelements;

    // Mr = invDegree.*r
    elliptic.linAlg.amxpy(Ntotal, 1.0, elliptic.o_weight, o_r, 0.0, o_Mr);

    SEMFEMInterpKernel(mesh.Nelements,mesh.o_SEMFEMAnterp,o_Mr,o_rFEM);
    FEMogs->Gather(o_GrFEM, o_rFEM, ogs_dfloat, ogs_add, ogs_trans);

    parAlmond.Operator(o_GrFEM, o_GzFEM);

    FEMogs->Scatter(o_zFEM, o_GzFEM, ogs_dfloat, ogs_add, ogs_notrans);
    SEMFEMAnterpKernel(mesh.Nelements,mesh.o_SEMFEMAnterp,o_zFEM,o_Mr);

    // Mr = invDegree.*Mr
    elliptic.linAlg.amx(Ntotal, 1.0, elliptic.o_weight, o_Mr);

    elliptic.ogsMasked->GatherScatter(o_Mr, ogs_dfloat, ogs_add, ogs_sym);

  } else {

    FEMogs->Gather(o_rhsG, o_r, ogs_dfloat, ogs_add, ogs_notrans);
    parAlmond.Operator(o_rhsG, o_xG);
    FEMogs->Scatter(o_Mr, o_xG, ogs_dfloat, ogs_add, ogs_notrans);
  }

  if (elliptic.Nmasked)
      elliptic.maskKernel(elliptic.Nmasked, elliptic.o_maskIds, o_Mr);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}

SEMFEMPrecon::SEMFEMPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic), mesh(_elliptic.mesh), settings(_elliptic.settings),
  parAlmond(elliptic.platform, settings) {

  //sanity checking
  if (!settings.compareSetting("DISCRETIZATION", "CONTINUOUS") )
    LIBP_ABORT(string("SEMFEM is supported for CONTINUOUS only"));

  //make a low-order fem mesh from the sem mesh (also return globalIds of the enriched sem nodes, and faceNode mapping)
  int Nfp = 0;
  int *faceNodes = NULL;
  hlong *globalIds = NULL;
  femMesh = mesh.SetupSEMFEM(&globalIds, &Nfp, &faceNodes);

  //use the BCs to make a maskedGlobalIds array
  dlong Ntotal = mesh.NpFEM*mesh.Nelements;
  hlong* maskedGlobalIds = (hlong *) calloc(Ntotal,sizeof(hlong));
  memcpy(maskedGlobalIds, globalIds, Ntotal*sizeof(hlong));
  if (mesh.elementType==TRIANGLES) { //build a new mask for NpFEM>Np node sets
    // gather-scatter
    int verbose = 0;
    ogs_t *ogs = ogs_t::Setup(Ntotal, globalIds, mesh.comm, verbose, elliptic.platform);

    //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
    const int largeNumber = 1<<20;
    int *mapB = (int *) calloc(Ntotal,sizeof(int));
    for (dlong e=0;e<mesh.Nelements;e++) {
      for (int n=0;n<mesh.NpFEM;n++) mapB[n+e*mesh.NpFEM] = largeNumber;

      for (int f=0;f<mesh.Nfaces;f++) {
        int bc = elliptic.EToB[f+e*mesh.Nfaces];
        if (bc>0) {
          for (int n=0;n<Nfp;n++) {
            int fid = faceNodes[n+f*Nfp];
            mapB[fid+e*mesh.NpFEM] = mymin(bc,mapB[fid+e*mesh.Np]);
          }
        }
      }
    }
    ogs->GatherScatter(mapB, ogs_int, ogs_min, ogs_sym);

    //use the bc flags to find masked ids
    for (dlong n=0;n<mesh.Nelements*mesh.NpFEM;n++)
      if (mapB[n] == 1) //Dirichlet boundary
        maskedGlobalIds[n] = 0;

    free(mapB);
    ogs->Free();
  } else {
    //mask using the original mask
    for (dlong n=0;n<elliptic.Nmasked;n++)
      maskedGlobalIds[elliptic.maskIds[n]] = 0;
  }

  //build masked gs handle to gather from enriched sem nodes to assembled fem problem
  int verbose = 0;
  ogs_t::Unique(maskedGlobalIds, Ntotal, mesh.comm);     //flag a unique node in every gather node
  FEMogs = ogs_t::Setup(Ntotal, maskedGlobalIds, mesh.comm, verbose, elliptic.platform);

  //make a map from the fem mesh's nodes to the (enriched) sem nodes
  dlong *localIds = (dlong *) calloc(femMesh->Nelements*femMesh->Nverts,sizeof(dlong));
  for(dlong e=0;e<mesh.Nelements;++e){
    for (int n=0;n<mesh.NelFEM;n++) {
      dlong id[femMesh->Nverts];

      //local ids in the subelement fem grid
      for (int i=0;i<mesh.Nverts;i++)
        id[i] = e*mesh.NpFEM + mesh.FEMEToV[n*femMesh->Nverts+i];

      dlong femId = e*mesh.NelFEM*femMesh->Nverts+n*mesh.Nverts;
      switch(mesh.elementType){
      case TRIANGLES:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[2];
        break;
      case QUADRILATERALS:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[3];  //need to swap this as the Np nodes are ordered [0,1,3,2] in a degree 1 element
        localIds[femId+3] = id[2];
        break;
      case TETRAHEDRA:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[2];
        localIds[femId+3] = id[3];
        break;
      case HEXAHEDRA:
        localIds[femId+0] = id[0];
        localIds[femId+1] = id[1];
        localIds[femId+2] = id[3];  //need to swap this as the Np nodes are ordered [0,1,3,2,4,5,7,6] in a degree 1 element
        localIds[femId+3] = id[2];
        localIds[femId+4] = id[4];
        localIds[femId+5] = id[5];
        localIds[femId+6] = id[7];
        localIds[femId+7] = id[6];
        break;
      }
    }
  }

  //make a fem elliptic solver
  femElliptic = new elliptic_t(elliptic.platform, *femMesh,
                               elliptic.settings, elliptic.lambda);
  femElliptic->ogsMasked = FEMogs; //only for getting Ngather when building matrix

  // number of degrees of freedom on this rank (after gathering)
  hlong Ngather = FEMogs->Ngather;

  // create a global numbering system
  hlong *globalIds2 = (hlong *) calloc(Ngather,sizeof(hlong));
  int   *owner     = (int *) calloc(Ngather,sizeof(int));

  // every gathered degree of freedom has its own global id
  hlong *globalStarts = (hlong *) calloc(mesh.size+1,sizeof(hlong));
  MPI_Allgather(&Ngather, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh.comm);
  for(int r=0;r<mesh.size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  //use the offsets to set a consecutive global numbering
  for (dlong n =0;n<FEMogs->Ngather;n++) {
    globalIds2[n] = n + globalStarts[mesh.rank];
    owner[n] = mesh.rank;
  }
  free(globalStarts);

  //scatter this numbering to the original nodes
  Ntotal = mesh.NpFEM*mesh.Nelements;
  hlong* maskedGlobalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  int  * maskedGlobalOwners    = (int *) calloc(Ntotal,sizeof(int));
  for (dlong n=0;n<Ntotal;n++) maskedGlobalNumbering[n] = -1;
  FEMogs->Scatter(maskedGlobalNumbering, globalIds2, ogs_hlong, ogs_add, ogs_notrans);
  FEMogs->Scatter(maskedGlobalOwners, owner, ogs_int, ogs_add, ogs_notrans);
  free(globalIds2); free(owner);

  //transfer the consecutive global numbering to the fem mesh
  Ntotal = femMesh->Np*femMesh->Nelements;
  femElliptic->maskedGlobalNumbering = (hlong *) calloc(Ntotal,sizeof(hlong));
  femElliptic->maskedGlobalOwners    = (int *) calloc(Ntotal,sizeof(int));

  for (dlong e=0;e<femMesh->Nelements;e++) {
    for (int n=0;n<femMesh->Np;n++) {
      dlong id = e*femMesh->Np + n;
      dlong localId = localIds[id];
      femElliptic->maskedGlobalNumbering[id] = maskedGlobalNumbering[localId];
      femElliptic->maskedGlobalOwners[id] = maskedGlobalOwners[localId];
    }
  }
  free(localIds); free(maskedGlobalNumbering); free(maskedGlobalOwners);

  //finally, build the fem matrix and pass to parAlmond
  parAlmond::parCOO A(elliptic.platform);
  femElliptic->BuildOperatorMatrixContinuous(A);

  //populate null space unit vector
  int rank = elliptic.platform.rank;
  int size = elliptic.platform.size;
  hlong TotalRows = A.globalStarts[size];
  dlong numLocalRows = (dlong) (A.globalStarts[rank+1]-A.globalStarts[rank]);
  dfloat *null = (dfloat *) malloc(numLocalRows*sizeof(dfloat));
  for (dlong i=0;i<numLocalRows;i++) null[i] = 1.0/sqrt(TotalRows);

  parAlmond.AMGSetup(A, elliptic.allNeumann, null, elliptic.allNeumannPenalty);
  free(null);

  parAlmond.Report();

  if (mesh.elementType==TRIANGLES) {
    // build interp and anterp
    dfloat *SEMFEMAnterp = (dfloat*) calloc(mesh.NpFEM*mesh.Np, sizeof(dfloat));
    for(int n=0;n<mesh.NpFEM;++n){
      for(int m=0;m<mesh.Np;++m){
        SEMFEMAnterp[n+m*mesh.NpFEM] = mesh.SEMFEMInterp[n*mesh.Np+m];
      }
    }

    mesh.o_SEMFEMInterp = elliptic.platform.malloc(mesh.NpFEM*mesh.Np*sizeof(dfloat),mesh.SEMFEMInterp);
    mesh.o_SEMFEMAnterp = elliptic.platform.malloc(mesh.NpFEM*mesh.Np*sizeof(dfloat),SEMFEMAnterp);

    free(SEMFEMAnterp);

    dfloat *dummy = (dfloat*) calloc(mesh.Nelements*mesh.NpFEM,sizeof(dfloat)); //need this to avoid uninitialized memory warnings
    o_rFEM = elliptic.platform.malloc(mesh.Nelements*mesh.NpFEM*sizeof(dfloat), dummy);
    o_zFEM = elliptic.platform.malloc(mesh.Nelements*mesh.NpFEM*sizeof(dfloat), dummy);
    free(dummy);

    dlong Ncols = parAlmond.getNumCols(0);
    dummy = (dfloat*) calloc(Ncols,sizeof(dfloat));
    o_GrFEM = elliptic.platform.malloc(Ncols*sizeof(dfloat),dummy);
    o_GzFEM = elliptic.platform.malloc(Ncols*sizeof(dfloat),dummy);
    free(dummy);

    //build kernels
    occa::properties kernelInfo = mesh.props;

    kernelInfo["defines/" "p_Np"]= mesh.Np;
    kernelInfo["defines/" "p_NpFEM"]= mesh.NpFEM;

    int NblockV = 512/mesh.NpFEM;
    kernelInfo["defines/" "p_NblockV"]= NblockV;

    SEMFEMInterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                                     "ellipticSEMFEMInterp", kernelInfo);

    SEMFEMAnterpKernel = elliptic.platform.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                                     "ellipticSEMFEMAnterp", kernelInfo);
  } else {
    dlong Ncols = parAlmond.getNumCols(0);
    // rhsG = (dfloat*) calloc(Ncols,sizeof(dfloat));
    // xG   = (dfloat*) calloc(Ncols,sizeof(dfloat));
    dfloat *dummy = (dfloat*) calloc(Ncols,sizeof(dfloat));
    o_rhsG = elliptic.platform.malloc(Ncols*sizeof(dfloat),dummy);
    o_xG   = elliptic.platform.malloc(Ncols*sizeof(dfloat),dummy);
    free(dummy);
  }
}

SEMFEMPrecon::~SEMFEMPrecon() {
  femElliptic->ogsMasked->Free();

  femMesh->halo->Free();

  SEMFEMInterpKernel.free();
  SEMFEMAnterpKernel.free();
}