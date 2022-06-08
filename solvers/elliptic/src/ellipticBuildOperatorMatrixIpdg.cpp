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

#include "elliptic.hpp"

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

void elliptic_t::BuildOperatorMatrixIpdg(parAlmond::parCOO& A){

  switch(mesh.elementType){
  case Mesh::TRIANGLES:
  {
    if(mesh.dim==2)
      BuildOperatorMatrixIpdgTri2D(A);
    else
      BuildOperatorMatrixIpdgTri3D(A);
    break;
  }
  case Mesh::QUADRILATERALS:{
    if(mesh.dim==2)
      BuildOperatorMatrixIpdgQuad2D(A);
    else
      BuildOperatorMatrixIpdgQuad3D(A);
    break;
  }
  case Mesh::TETRAHEDRA:
    BuildOperatorMatrixIpdgTet3D(A); break;
  case Mesh::HEXAHEDRA:
    BuildOperatorMatrixIpdgHex3D(A); break;
  }
}

void elliptic_t::BuildOperatorMatrixIpdgTri2D(parAlmond::parCOO& A){

  int Np = mesh.Np;
  int Nfp = mesh.Nfp;
  int Nfaces = mesh.Nfaces;
  dlong Nelements = mesh.Nelements;

  // number of degrees of freedom on this rank
  hlong Nnum = Np*Nelements;

  // create a global numbering system
  memory<hlong> globalIds((Nelements+mesh.totalHaloPairs)*Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = Nelements;
  hlong globalElementOffset = Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np + n] = n + (e + globalElementOffset)*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, Np);

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(Nfaces*Nfp*Nfp);
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<Nfp;n++) {
      int fn = mesh.faceNodes[f*Nfp+n];

      for (int m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<Np;i++){
          MSnm += mesh.MM[fn+i*Np]*mesh.LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }


  // reset non-zero counter
  dlong nnz = 0;

  A.entries.malloc(nnzLocalBound);

  memory<dfloat> SM(Np*Np);
  memory<dfloat> SP(Np*Np);

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // loop over all elements
  for(dlong eM=0;eM<Nelements;++eM){

    dlong vbase = eM*mesh.Nvgeo;
    dfloat drdx = mesh.vgeo[vbase+mesh.RXID];
    dfloat drdy = mesh.vgeo[vbase+mesh.RYID];
    dfloat dsdx = mesh.vgeo[vbase+mesh.SXID];
    dfloat dsdy = mesh.vgeo[vbase+mesh.SYID];
    dfloat J = mesh.vgeo[vbase+mesh.JID];

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        SM[n*Np+m]  = J*lambda*mesh.MM[n*Np+m];
        SM[n*Np+m] += J*drdx*drdx*mesh.Srr[n*Np+m];
        SM[n*Np+m] += J*drdx*dsdx*mesh.Srs[n*Np+m];
        SM[n*Np+m] += J*dsdx*dsdx*mesh.Sss[n*Np+m];

        SM[n*Np+m] += J*drdy*drdy*mesh.Srr[n*Np+m];
        SM[n*Np+m] += J*drdy*dsdy*mesh.Srs[n*Np+m];
        SM[n*Np+m] += J*dsdy*dsdy*mesh.Sss[n*Np+m];
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {

      for (int n=0;n<Np*Np;n++) SP[n] =0;

      // load surface geofactors for this face
      dlong sid = mesh.Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh.sgeo[sid+mesh.NXID];
      dfloat ny = mesh.sgeo[sid+mesh.NYID];
      dfloat sJ = mesh.sgeo[sid+mesh.SJID];
      dfloat hinv = mesh.sgeo[sid+mesh.IHID];
      dfloat penalty = tau*hinv;

      dlong eP = mesh.EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      dlong vbaseP = eP*mesh.Nvgeo;
      dfloat drdxP = mesh.vgeo[vbaseP+mesh.RXID];
      dfloat drdyP = mesh.vgeo[vbaseP+mesh.RYID];
      dfloat dsdxP = mesh.vgeo[vbaseP+mesh.SXID];
      dfloat dsdyP = mesh.vgeo[vbaseP+mesh.SYID];

      int bcD = 0, bcN =0;
      int bc = mesh.EToB[fM+Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // reset eP
      eP = mesh.EToE[eM*Nfaces+fM];

      // mass matrix for this face
      memory<dfloat> MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh.faceNodes[fM*Nfp+n];
          int mM = mesh.faceNodes[fM*Nfp+m];
          int mP  = (int) (mesh.vmapP[idM]%Np);

          dfloat MSfnm = sJ*MSf[n*Nfp+m];

          SM[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          SP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh.faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh.faceNodes[fM*Nfp+i];
            int iP = (int) (mesh.vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = drdx*mesh.Dr[iM*Np+m] + dsdx*mesh.Ds[iM*Np+m];
            dfloat DyMim = drdy*mesh.Dr[iM*Np+m] + dsdy*mesh.Ds[iM*Np+m];

            dfloat DxPim = drdxP*mesh.Dr[iP*Np+m] + dsdxP*mesh.Ds[iP*Np+m];
            dfloat DyPim = drdyP*mesh.Dr[iP*Np+m] + dsdyP*mesh.Ds[iP*Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            SM[nM*Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

            SP[nM*Np+m] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            SP[nM*Np+m] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh.faceNodes[fM*Nfp+m];
          int mP = (int) (mesh.vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

          for(int i=0;i<Nfp;++i){
            int iM = mesh.faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh.Dr[iM*Np+n] + dsdx*mesh.Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh.Dr[iM*Np+n] + dsdy*mesh.Ds[iM*Np+n];

            SM[n*Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;

            SP[n*Np+mP] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            SP[n*Np+mP] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }

      // store non-zeros for off diagonal block
      for(int n=0;n<Np;++n){
        for(int m=0;m<Np;++m){
          dfloat val = SP[n*Np+m];
          if(std::abs(val)>tol){
            A.entries[nnz].row = globalIds[eM*Np + n];
            A.entries[nnz].col = globalIds[eP*Np + m];
            A.entries[nnz].val = val;
            ++nnz;
          }
        }
      }
    }
    // store non-zeros for diagonal block
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        dfloat val = SM[n*Np+m];
        if(std::abs(val)>tol){
          A.entries[nnz].row = globalIds[eM*Np + n];
          A.entries[nnz].col = globalIds[eM*Np + m];
          A.entries[nnz].val = val;
          ++nnz;
        }
      }
    }
  }

  //printf("nnz = %d\n", nnz);

  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");

#if 0
  dfloat* Ap = (dfloat *) calloc(Np*Np*Nelements*Nelements,sizeof(dfloat));
  for (int n=0;n<nnz;n++) {
    int row = A.entries[n].row;
    int col = A.entries[n].col;

    Ap[col+row*Np*Nelements] = A.entries[n].val;
  }

  for (int i=0;i<Np*Nelements;i++) {
    for (int j =0;j<Nelements*Np;j++) {
      printf("%4.2f \t", Ap[j+i*Np*Nelements]);
    }
    printf("\n");
  }
#endif
}

void elliptic_t::BuildOperatorMatrixIpdgTri3D(parAlmond::parCOO& A){

  int Np = mesh.Np;
  int Nfp = mesh.Nfp;
  int Nfaces = mesh.Nfaces;
  dlong Nelements = mesh.Nelements;

  // number of degrees of freedom on this rank
  hlong Nnum = Np*Nelements;

  // create a global numbering system
  memory<hlong> globalIds((Nelements+mesh.totalHaloPairs)*Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = Nelements;
  hlong globalElementOffset = Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np + n] = n + (e + globalElementOffset)*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, Np);

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(Nfaces*Nfp*Nfp);
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<Nfp;n++) {
      int fn = mesh.faceNodes[f*Nfp+n];

      for (int m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<Np;i++){
          MSnm += mesh.MM[fn+i*Np]*mesh.LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }


  // reset non-zero counter
  dlong nnz = 0;

  A.entries.malloc(nnzLocalBound);

  memory<dfloat> SM(Np*Np);
  memory<dfloat> SP(Np*Np);

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // loop over all elements
  for(dlong eM=0;eM<Nelements;++eM){

    dlong vbase = eM*mesh.Nvgeo;
    dfloat drdx = mesh.vgeo[vbase+mesh.RXID];
    dfloat drdy = mesh.vgeo[vbase+mesh.RYID];
    dfloat drdz = mesh.vgeo[vbase+mesh.RZID];
    dfloat dsdx = mesh.vgeo[vbase+mesh.SXID];
    dfloat dsdy = mesh.vgeo[vbase+mesh.SYID];
    dfloat dsdz = mesh.vgeo[vbase+mesh.SZID];
    dfloat J = mesh.vgeo[vbase+mesh.JID];

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        SM[n*Np+m]  = J*lambda*mesh.MM[n*Np+m];
        SM[n*Np+m] += J*drdx*drdx*mesh.Srr[n*Np+m];
        SM[n*Np+m] += J*drdx*dsdx*mesh.Srs[n*Np+m];
        SM[n*Np+m] += J*dsdx*dsdx*mesh.Sss[n*Np+m];

        SM[n*Np+m] += J*drdy*drdy*mesh.Srr[n*Np+m];
        SM[n*Np+m] += J*drdy*dsdy*mesh.Srs[n*Np+m];
        SM[n*Np+m] += J*dsdy*dsdy*mesh.Sss[n*Np+m];

        SM[n*Np+m] += J*drdz*drdz*mesh.Srr[n*Np+m];
        SM[n*Np+m] += J*drdz*dsdz*mesh.Srs[n*Np+m];
        SM[n*Np+m] += J*dsdz*dsdz*mesh.Sss[n*Np+m];

      }
    }

    for (int fM=0;fM<Nfaces;fM++) {

      for (int n=0;n<Np*Np;n++) SP[n] =0;

      // load surface geofactors for this face
      dlong sid = mesh.Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh.sgeo[sid+mesh.NXID];
      dfloat ny = mesh.sgeo[sid+mesh.NYID];
      dfloat nz = mesh.sgeo[sid+mesh.NZID];
      dfloat sJ = mesh.sgeo[sid+mesh.SJID];
      dfloat hinv = mesh.sgeo[sid+mesh.IHID];
      dfloat penalty = tau*hinv;

      dlong eP = mesh.EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      dlong vbaseP = eP*mesh.Nvgeo;
      dfloat drdxP = mesh.vgeo[vbaseP+mesh.RXID];
      dfloat drdyP = mesh.vgeo[vbaseP+mesh.RYID];
      dfloat drdzP = mesh.vgeo[vbaseP+mesh.RZID];
      dfloat dsdxP = mesh.vgeo[vbaseP+mesh.SXID];
      dfloat dsdyP = mesh.vgeo[vbaseP+mesh.SYID];
      dfloat dsdzP = mesh.vgeo[vbaseP+mesh.SZID];

      int bcD = 0, bcN =0;
      int bc = mesh.EToB[fM+Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // reset eP
      eP = mesh.EToE[eM*Nfaces+fM];

      // mass matrix for this face
      memory<dfloat> MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh.faceNodes[fM*Nfp+n];
          int mM = mesh.faceNodes[fM*Nfp+m];
          int mP  = (int) (mesh.vmapP[idM]%Np);

          dfloat MSfnm = sJ*MSf[n*Nfp+m];

          SM[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          SP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh.faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh.faceNodes[fM*Nfp+i];
            int iP = (int) (mesh.vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = drdx*mesh.Dr[iM*Np+m] + dsdx*mesh.Ds[iM*Np+m];
            dfloat DyMim = drdy*mesh.Dr[iM*Np+m] + dsdy*mesh.Ds[iM*Np+m];
            dfloat DzMim = drdz*mesh.Dr[iM*Np+m] + dsdz*mesh.Ds[iM*Np+m];

            dfloat DxPim = drdxP*mesh.Dr[iP*Np+m] + dsdxP*mesh.Ds[iP*Np+m];
            dfloat DyPim = drdyP*mesh.Dr[iP*Np+m] + dsdyP*mesh.Ds[iP*Np+m];
            dfloat DzPim = drdzP*mesh.Dr[iP*Np+m] + dsdzP*mesh.Ds[iP*Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            SM[nM*Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
            SM[nM*Np+m] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;

            SP[nM*Np+m] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            SP[nM*Np+m] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
            SP[nM*Np+m] += -0.5*nz*(1-bcD)*(1-bcN)*MSfni*DzPim;
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh.faceNodes[fM*Nfp+m];
          int mP = (int) (mesh.vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

          for(int i=0;i<Nfp;++i){
            int iM = mesh.faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh.Dr[iM*Np+n] + dsdx*mesh.Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh.Dr[iM*Np+n] + dsdy*mesh.Ds[iM*Np+n];
            dfloat DzMin = drdz*mesh.Dr[iM*Np+n] + dsdz*mesh.Ds[iM*Np+n];

            SM[n*Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
            SM[n*Np+mM] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;

            SP[n*Np+mP] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            SP[n*Np+mP] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
            SP[n*Np+mP] +=  +0.5*nz*(1-bcD)*(1-bcN)*DzMin*MSfim;
          }
        }
      }

      // store non-zeros for off diagonal block
      for(int n=0;n<Np;++n){
        for(int m=0;m<Np;++m){
          dfloat val = SP[n*Np+m];
          if(std::abs(val)>tol){
            A.entries[nnz].row = globalIds[eM*Np + n];
            A.entries[nnz].col = globalIds[eP*Np + m];
            A.entries[nnz].val = val;
            ++nnz;
          }
        }
      }
    }
    // store non-zeros for diagonal block
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        dfloat val = SM[n*Np+m];
        if(std::abs(val)>tol){
          A.entries[nnz].row = globalIds[eM*Np + n];
          A.entries[nnz].col = globalIds[eM*Np + m];
          A.entries[nnz].val = val;
          ++nnz;
        }
      }
    }
  }

  //printf("nnz = %d\n", nnz);

  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");
}



void elliptic_t::BuildOperatorMatrixIpdgQuad2D(parAlmond::parCOO& A){

  int Np = mesh.Np;
  int Nfaces = mesh.Nfaces;
  dlong Nelements = mesh.Nelements;

  hlong Nnum = mesh.Np*mesh.Nelements;

  // create a global numbering system
  memory<hlong> globalIds((Nelements+mesh.totalHaloPairs)*Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = Nelements;
  hlong globalElementOffset = Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np + n] = n + (e + globalElementOffset)*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, Np);

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  memory<dfloat> B (mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Br(mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Bs(mesh.Np*mesh.Np, 0.0);

  int mode = 0;
  for(int nj=0;nj<mesh.N+1;++nj){
    for(int ni=0;ni<mesh.N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh.N+1;++j){
        for(int i=0;i<mesh.N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh.Np+node] = 1;
          if(nj==j)
            Br[mode*mesh.Np+node] = mesh.D[ni+mesh.Nq*i];
          if(ni==i)
            Bs[mode*mesh.Np+node] = mesh.D[nj+mesh.Nq*j];

          ++node;
        }
      }

      ++mode;
    }
  }

  A.entries.malloc(nnzLocalBound);

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // reset non-zero counter
  dlong nnz = 0;

  // loop over all elements
  for(dlong eM=0;eM<mesh.Nelements;++eM){

    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh.Np;++n){
      for(int m=0;m<mesh.Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh.Np;++i){
          dlong base = eM*mesh.Np*mesh.Nvgeo + i;
          dfloat drdx = mesh.vgeo[base+mesh.Np*mesh.RXID];
          dfloat drdy = mesh.vgeo[base+mesh.Np*mesh.RYID];
          dfloat dsdx = mesh.vgeo[base+mesh.Np*mesh.SXID];
          dfloat dsdy = mesh.vgeo[base+mesh.Np*mesh.SYID];
          dfloat JW   = mesh.vgeo[base+mesh.Np*mesh.JWID];

          int idn = n*mesh.Np+i;
          int idm = m*mesh.Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];
          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh.Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh.Nfp;++i){
            int vidM = mesh.faceNodes[i+fM*mesh.Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh.Np*mesh.Nvgeo + vidM;
            dfloat drdxM = mesh.vgeo[baseM+mesh.Np*mesh.RXID];
            dfloat drdyM = mesh.vgeo[baseM+mesh.Np*mesh.RYID];
            dfloat dsdxM = mesh.vgeo[baseM+mesh.Np*mesh.SXID];
            dfloat dsdyM = mesh.vgeo[baseM+mesh.Np*mesh.SYID];

            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh.Nfp*mesh.Nfaces+fM*mesh.Nfp+i;
            int vidP      = (int) (mesh.vmapP[idM]%mesh.Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh.vmapP[idM]/mesh.Np;
            dlong baseP   = localEP*mesh.Np*mesh.Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh.vgeo[baseP+mesh.Np*mesh.RXID];
            dfloat drdyP = mesh.vgeo[baseP+mesh.Np*mesh.RYID];
            dfloat dsdxP = mesh.vgeo[baseP+mesh.Np*mesh.SXID];
            dfloat dsdyP = mesh.vgeo[baseP+mesh.Np*mesh.SYID];

            // grab surface geometric factors
            dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
            dfloat nx = mesh.sgeo[base+mesh.NXID];
            dfloat ny = mesh.sgeo[base+mesh.NYID];
            dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
            dfloat hinv = mesh.sgeo[base+mesh.IHID];

            // form negative trace terms in IPDG
            int idnM = n*mesh.Np+vidM;
            int idmM = m*mesh.Np+vidM;
            int idmP = m*mesh.Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
            dfloat lmM = B[idmM];

            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP;
            dfloat lmP = B[idmP];

            dfloat penalty = tau*hinv;

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            dlong eP    = mesh.EToE[eM*mesh.Nfaces+fM];
            if (eP < 0) {
              int qSgn, gradqSgn;
              int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag
              int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              Anm += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              Anm += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              Anm += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-)
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
          }
          if(std::abs(AnmP)>tol){
            // remote info
            dlong eP    = mesh.EToE[eM*mesh.Nfaces+fM];
            A.entries[nnz].row = globalIds[eM*mesh.Np + n];
            A.entries[nnz].col = globalIds[eP*mesh.Np + m];
            A.entries[nnz].val = AnmP;
            ++nnz;
          }
        }
        if(std::abs(Anm)>tol){
          // local block
          A.entries[nnz].row = globalIds[eM*mesh.Np+n];
          A.entries[nnz].col = globalIds[eM*mesh.Np+m];
          A.entries[nnz].val = Anm;
          ++nnz;
        }
      }
    }
  }

  // sort received non-zero entries by row block
  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");
}


void elliptic_t::BuildOperatorMatrixIpdgQuad3D(parAlmond::parCOO& A){

  int Np = mesh.Np;
  int Nfaces = mesh.Nfaces;
  dlong Nelements = mesh.Nelements;

  hlong Nnum = mesh.Np*mesh.Nelements;

  // create a global numbering system
  memory<hlong> globalIds((Nelements+mesh.totalHaloPairs)*Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = Nelements;
  hlong globalElementOffset = Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np + n] = n + (e + globalElementOffset)*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, Np);

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  memory<dfloat> B (mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Br(mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Bs(mesh.Np*mesh.Np, 0.0);

  int mode = 0;
  for(int nj=0;nj<mesh.N+1;++nj){
    for(int ni=0;ni<mesh.N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh.N+1;++j){
        for(int i=0;i<mesh.N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh.Np+node] = 1;
          if(nj==j)
            Br[mode*mesh.Np+node] = mesh.D[ni+mesh.Nq*i];
          if(ni==i)
            Bs[mode*mesh.Np+node] = mesh.D[nj+mesh.Nq*j];

          ++node;
        }
      }

      ++mode;
    }
  }

  A.entries.malloc(nnzLocalBound);

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // reset non-zero counter
  dlong nnz = 0;

  // loop over all elements
  for(dlong eM=0;eM<mesh.Nelements;++eM){

    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh.Np;++n){
      for(int m=0;m<mesh.Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh.Np;++i){
          dlong base = eM*mesh.Np*mesh.Nvgeo + i;
          dfloat drdx = mesh.vgeo[base+mesh.Np*mesh.RXID];
          dfloat drdy = mesh.vgeo[base+mesh.Np*mesh.RYID];
          dfloat drdz = mesh.vgeo[base+mesh.Np*mesh.RZID];
          dfloat dsdx = mesh.vgeo[base+mesh.Np*mesh.SXID];
          dfloat dsdy = mesh.vgeo[base+mesh.Np*mesh.SYID];
          dfloat dsdz = mesh.vgeo[base+mesh.Np*mesh.SZID];
          // dfloat dtdx = mesh.vgeo[base+mesh.Np*mesh.TXID];
          // dfloat dtdy = mesh.vgeo[base+mesh.Np*mesh.TYID];
          // dfloat dtdz = mesh.vgeo[base+mesh.Np*mesh.TZID];
          dfloat JW   = mesh.vgeo[base+mesh.Np*mesh.JWID];

          int idn = n*mesh.Np+i;
          int idm = m*mesh.Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];// + dtdx;
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];// + dtdy;
          dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn];// + dtdz;

          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];// + dtdx;
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];// + dtdy;
          dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm];// + dtdz;

          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy + dlndz*dlmdz);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh.Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh.Nfp;++i){
            int vidM = mesh.faceNodes[i+fM*mesh.Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh.Np*mesh.Nvgeo + vidM;
            dfloat drdxM = mesh.vgeo[baseM+mesh.Np*mesh.RXID];
            dfloat drdyM = mesh.vgeo[baseM+mesh.Np*mesh.RYID];
            dfloat drdzM = mesh.vgeo[baseM+mesh.Np*mesh.RZID];

            dfloat dsdxM = mesh.vgeo[baseM+mesh.Np*mesh.SXID];
            dfloat dsdyM = mesh.vgeo[baseM+mesh.Np*mesh.SYID];
            dfloat dsdzM = mesh.vgeo[baseM+mesh.Np*mesh.SZID];

            // dfloat dtdxM = mesh.vgeo[baseM+mesh.Np*mesh.TXID];
            // dfloat dtdyM = mesh.vgeo[baseM+mesh.Np*mesh.TYID];
            // dfloat dtdzM = mesh.vgeo[baseM+mesh.Np*mesh.TZID];


            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh.Nfp*mesh.Nfaces+fM*mesh.Nfp+i;
            int vidP      = (int) (mesh.vmapP[idM]%mesh.Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh.vmapP[idM]/mesh.Np;
            dlong baseP   = localEP*mesh.Np*mesh.Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh.vgeo[baseP+mesh.Np*mesh.RXID];
            dfloat drdyP = mesh.vgeo[baseP+mesh.Np*mesh.RYID];
            dfloat drdzP = mesh.vgeo[baseP+mesh.Np*mesh.RZID];

            dfloat dsdxP = mesh.vgeo[baseP+mesh.Np*mesh.SXID];
            dfloat dsdyP = mesh.vgeo[baseP+mesh.Np*mesh.SYID];
            dfloat dsdzP = mesh.vgeo[baseP+mesh.Np*mesh.SZID];

            // dfloat dtdxP = mesh.vgeo[baseP+mesh.Np*mesh.TXID];
            // dfloat dtdyP = mesh.vgeo[baseP+mesh.Np*mesh.TYID];
            // dfloat dtdzP = mesh.vgeo[baseP+mesh.Np*mesh.TZID];

            // grab surface geometric factors
            dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
            dfloat nx = mesh.sgeo[base+mesh.NXID];
            dfloat ny = mesh.sgeo[base+mesh.NYID];
            dfloat nz = mesh.sgeo[base+mesh.NZID];
            dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
            dfloat hinv = mesh.sgeo[base+mesh.IHID];

            // form negative trace terms in IPDG
            int idnM = n*mesh.Np+vidM;
            int idmM = m*mesh.Np+vidM;
            int idmP = m*mesh.Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];// + dtdxM;
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];// + dtdyM;
            dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM];// + dtdzM;

            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM + nz*dlndzM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];// + dtdxM;
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];// + dtdyM;
            dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM];// + dtdzM;
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM + nz*dlmdzM;
            dfloat lmM = B[idmM];

            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP];// + dtdxP;
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP];// + dtdyP;
            dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP];// + dtdzP;
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;
            dfloat lmP = B[idmP];

            dfloat penalty = tau*hinv;

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
            AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
            AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
          }
          if(std::abs(AnmP)>tol){
            // remote info
            dlong eP    = mesh.EToE[eM*mesh.Nfaces+fM];
            A.entries[nnz].row = globalIds[eM*mesh.Np + n];
            A.entries[nnz].col = globalIds[eP*mesh.Np + m];
            A.entries[nnz].val = AnmP;
            ++nnz;
          }
        }

        if(std::abs(Anm)>tol){
          // local block
          A.entries[nnz].row = globalIds[eM*mesh.Np+n];
          A.entries[nnz].col = globalIds[eM*mesh.Np+m];
          A.entries[nnz].val = Anm;
          ++nnz;
        }
      }
    }
  }

  // sort received non-zero entries by row block
  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");

#if 0
  {
    FILE *fp = fopen("DGS.dat", "w");
    for(int n=0;n<nnz;++n){
      fprintf(fp, "%d %d %17.15lf\n",
              A.entries[n].row+1,
              A.entries[n].col+1,
              A.entries[n].val);
    }
    fclose(fp);
  }
#endif
}






void elliptic_t::BuildOperatorMatrixIpdgTet3D(parAlmond::parCOO& A){

  // number of degrees of freedom on this rank
  hlong Nnum = mesh.Np*mesh.Nelements;

  // create a global numbering system
  memory<hlong> globalIds((mesh.Nelements+mesh.totalHaloPairs)*mesh.Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = mesh.Nelements;
  hlong globalElementOffset = mesh.Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - mesh.Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<mesh.Nelements;e++) {
    for (int n=0;n<mesh.Np;n++) {
      globalIds[e*mesh.Np + n] = n + (e + globalElementOffset)*mesh.Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, mesh.Np);

  dlong nnzLocalBound = mesh.Np*mesh.Np*(1+mesh.Nfaces)*mesh.Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(mesh.Nfaces*mesh.Np*mesh.Nfp);
  for (int f=0;f<mesh.Nfaces;f++) {
    for (int n=0;n<mesh.Np;n++) {
      for (int m=0;m<mesh.Nfp;m++) {
        dfloat MSnm = 0;
        for (int i=0;i<mesh.Np;i++)
          MSnm += mesh.MM[n+i*mesh.Np]*mesh.LIFT[i*mesh.Nfp*mesh.Nfaces+f*mesh.Nfp+m];

        MS[m+n*mesh.Nfp + f*mesh.Nfp*mesh.Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS, DtT*MS
  memory<dfloat> DrTMS(mesh.Nfaces*mesh.Np*mesh.Nfp);
  memory<dfloat> DsTMS(mesh.Nfaces*mesh.Np*mesh.Nfp);
  memory<dfloat> DtTMS(mesh.Nfaces*mesh.Np*mesh.Nfp);
  for (int f=0;f<mesh.Nfaces;f++) {
    for (int n=0;n<mesh.Np;n++) {
      for (int i=0;i<mesh.Nfp;i++) {
        DrTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np] = 0.;
        DsTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np] = 0.;
        DtTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np] = 0.;
        for (int m=0;m<mesh.Np;m++) {
          DrTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np]
            += mesh.Dr[n+m*mesh.Np]*MS[i+m*mesh.Nfp+f*mesh.Nfp*mesh.Np];
          DsTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np]
            += mesh.Ds[n+m*mesh.Np]*MS[i+m*mesh.Nfp+f*mesh.Nfp*mesh.Np];
          DtTMS[i+n*mesh.Nfp + f*mesh.Nfp*mesh.Np]
            += mesh.Dt[n+m*mesh.Np]*MS[i+m*mesh.Nfp+f*mesh.Nfp*mesh.Np];
        }
      }
    }
  }

  A.entries.malloc(nnzLocalBound);

  // reset non-zero counter
  dlong nnz = 0;

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // loop over all elements
  //#pragma omp parallel
{

  memory<dfloat> BM(mesh.Np*mesh.Np);

  memory<dfloat> qmP(mesh.Nfp);
  memory<dfloat> qmM(mesh.Nfp);
  memory<dfloat> ndotgradqmM(mesh.Nfp);
  memory<dfloat> ndotgradqmP(mesh.Nfp);

  //#pragma omp for
  for(dlong eM=0;eM<mesh.Nelements;++eM){

    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase+mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase+mesh.G01ID];
    dfloat Grt = mesh.ggeo[gbase+mesh.G02ID];
    dfloat Gss = mesh.ggeo[gbase+mesh.G11ID];
    dfloat Gst = mesh.ggeo[gbase+mesh.G12ID];
    dfloat Gtt = mesh.ggeo[gbase+mesh.G22ID];
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      for(int m=0;m<mesh.Np;++m){
        BM[m+n*mesh.Np]  = J*lambda*mesh.MM[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Grr*mesh.Srr[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Grs*mesh.Srs[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Grt*mesh.Srt[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Gss*mesh.Sss[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Gst*mesh.Sst[m+n*mesh.Np];
        BM[m+n*mesh.Np] += Gtt*mesh.Stt[m+n*mesh.Np];
      }
    }

    dlong vbase = eM*mesh.Nvgeo;
    dfloat drdx = mesh.vgeo[vbase+mesh.RXID];
    dfloat drdy = mesh.vgeo[vbase+mesh.RYID];
    dfloat drdz = mesh.vgeo[vbase+mesh.RZID];
    dfloat dsdx = mesh.vgeo[vbase+mesh.SXID];
    dfloat dsdy = mesh.vgeo[vbase+mesh.SYID];
    dfloat dsdz = mesh.vgeo[vbase+mesh.SZID];
    dfloat dtdx = mesh.vgeo[vbase+mesh.TXID];
    dfloat dtdy = mesh.vgeo[vbase+mesh.TYID];
    dfloat dtdz = mesh.vgeo[vbase+mesh.TZID];

    for (int m=0;m<mesh.Np;m++) {
      for (int fM=0;fM<mesh.Nfaces;fM++) {
        // load surface geofactors for this face
        dlong sid = mesh.Nsgeo*(eM*mesh.Nfaces+fM);
        dfloat nx = mesh.sgeo[sid+mesh.NXID];
        dfloat ny = mesh.sgeo[sid+mesh.NYID];
        dfloat nz = mesh.sgeo[sid+mesh.NZID];
        dfloat sJ = mesh.sgeo[sid+mesh.SJID];
        dfloat hinv = mesh.sgeo[sid+mesh.IHID];

        dlong eP = mesh.EToE[eM*mesh.Nfaces+fM];
        if (eP < 0) eP = eM;
        dlong vbaseP = eP*mesh.Nvgeo;
        dfloat drdxP = mesh.vgeo[vbaseP+mesh.RXID];
        dfloat drdyP = mesh.vgeo[vbaseP+mesh.RYID];
        dfloat drdzP = mesh.vgeo[vbaseP+mesh.RZID];
        dfloat dsdxP = mesh.vgeo[vbaseP+mesh.SXID];
        dfloat dsdyP = mesh.vgeo[vbaseP+mesh.SYID];
        dfloat dsdzP = mesh.vgeo[vbaseP+mesh.SZID];
        dfloat dtdxP = mesh.vgeo[vbaseP+mesh.TXID];
        dfloat dtdyP = mesh.vgeo[vbaseP+mesh.TYID];
        dfloat dtdzP = mesh.vgeo[vbaseP+mesh.TZID];

        // extract trace nodes
        for (int i=0;i<mesh.Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          int idM    = eM*mesh.Nfp*mesh.Nfaces+fM*mesh.Nfp+i;
          int vidM   = mesh.faceNodes[i+fM*mesh.Nfp];
          int vidP   = (int) (mesh.vmapP[idM]%mesh.Np); // only use this to identify location of positive trace vgeo

          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;

          ndotgradqmM[i] = (nx*drdx+ny*drdy+nz*drdz)*mesh.Dr[m+vidM*mesh.Np]
                          +(nx*dsdx+ny*dsdy+nz*dsdz)*mesh.Ds[m+vidM*mesh.Np]
                          +(nx*dtdx+ny*dtdy+nz*dtdz)*mesh.Dt[m+vidM*mesh.Np];
          ndotgradqmP[i] = (nx*drdxP+ny*drdyP+nz*drdzP)*mesh.Dr[m+vidP*mesh.Np]
                          +(nx*dsdxP+ny*dsdyP+nz*dsdzP)*mesh.Ds[m+vidP*mesh.Np]
                          +(nx*dtdxP+ny*dtdyP+nz*dtdzP)*mesh.Dt[m+vidP*mesh.Np];
        }

        dfloat penalty = tau*hinv;
        eP = mesh.EToE[eM*mesh.Nfaces+fM];
        for (int n=0;n<mesh.Np;n++) {
          for (int i=0;i<mesh.Nfp;i++) {
            BM[m+n*mesh.Np] += -0.5*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*ndotgradqmM[i];
            BM[m+n*mesh.Np] += -0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i]
                                -0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i]
                                -0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i];
            BM[m+n*mesh.Np] += +0.5*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*penalty*qmM[i];
          }

          dfloat AnmP = 0;
          if (eP < 0) {
            int qSgn, gradqSgn;
            int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag
            int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
            if(bcType==1){ // Dirichlet
              qSgn     = -1;
              gradqSgn =  1;
            } else if (bcType==2){ // Neumann
              qSgn     =  1;
              gradqSgn = -1;
            } else { // Neumann for now
              qSgn     =  1;
              gradqSgn = -1;
            }

            for (int i=0;i<mesh.Nfp;i++) {
              BM[m+n*mesh.Np] += -0.5*gradqSgn*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*ndotgradqmM[i];
              BM[m+n*mesh.Np] += +0.5*qSgn*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmM[i];
              BM[m+n*mesh.Np] += -0.5*qSgn*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*penalty*qmM[i];
            }
          } else {
            for (int i=0;i<mesh.Nfp;i++) {
              AnmP += -0.5*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*ndotgradqmP[i];
              AnmP += +0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmP[i]
                      +0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmP[i]
                      +0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*qmP[i];
              AnmP += -0.5*sJ*MS[i+n*mesh.Nfp+fM*mesh.Nfp*mesh.Np]*penalty*qmP[i];
            }
          }

          if(std::abs(AnmP)>tol){
            //#pragma omp critical
            {
              // remote info
              A.entries[nnz].row = globalIds[eM*mesh.Np+n];
              A.entries[nnz].col = globalIds[eP*mesh.Np+m];
              A.entries[nnz].val = AnmP;
              ++nnz;
            }
          }
        }
      }
    }

    for (int n=0;n<mesh.Np;n++) {
      for (int m=0;m<mesh.Np;m++) {
        dfloat Anm = BM[m+n*mesh.Np];

        if(std::abs(Anm)>tol){
          //#pragma omp critical
          {
            A.entries[nnz].row = globalIds[eM*mesh.Np+n];
            A.entries[nnz].col = globalIds[eM*mesh.Np+m];
            A.entries[nnz].val = Anm;
            ++nnz;
          }
        }
      }
    }
  }
}

  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });
  // free up unused storage
  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");
}

void elliptic_t::BuildOperatorMatrixIpdgHex3D(parAlmond::parCOO& A){

  int Np = mesh.Np;
  int Nfaces = mesh.Nfaces;
  dlong Nelements = mesh.Nelements;

  hlong Nnum = mesh.Np*mesh.Nelements;

  // create a global numbering system
  memory<hlong> globalIds((Nelements+mesh.totalHaloPairs)*Np);

  // every degree of freedom has its own global id
  A.globalRowStarts.malloc(mesh.size+1,0);
  A.globalColStarts.malloc(mesh.size+1,0);
  mesh.comm.Allgather(Nnum, A.globalRowStarts+1);
  for(int r=0;r<mesh.size;++r) {
    A.globalRowStarts[r+1] = A.globalRowStarts[r]+A.globalRowStarts[r+1];
    A.globalColStarts[r+1] = A.globalRowStarts[r+1];
  }

  /* so find number of elements on each rank */
  hlong gNelements = Nelements;
  hlong globalElementOffset = Nelements;
  mesh.comm.Scan(gNelements, globalElementOffset);
  globalElementOffset = globalElementOffset - Nelements;
  //use the offsets to set a global id
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np + n] = n + (e + globalElementOffset)*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  mesh.halo.Exchange(globalIds, Np);

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  memory<dfloat> B (mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Br(mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Bs(mesh.Np*mesh.Np, 0.0);
  memory<dfloat> Bt(mesh.Np*mesh.Np, 0.0);

  int mode = 0;
  for(int nk=0;nk<mesh.N+1;++nk){
    for(int nj=0;nj<mesh.N+1;++nj){
      for(int ni=0;ni<mesh.N+1;++ni){

        int node = 0;

        for(int k=0;k<mesh.N+1;++k){
          for(int j=0;j<mesh.N+1;++j){
            for(int i=0;i<mesh.N+1;++i){

              if(nk==k && nj==j && ni==i)
                B[mode*mesh.Np+node] = 1;
              if(nj==j && nk==k)
                Br[mode*mesh.Np+node] = mesh.D[ni+mesh.Nq*i];
              if(ni==i && nk==k)
                Bs[mode*mesh.Np+node] = mesh.D[nj+mesh.Nq*j];
              if(ni==i && nj==j)
                Bt[mode*mesh.Np+node] = mesh.D[nk+mesh.Nq*k];

              ++node;
            }
          }
        }

        ++mode;
      }
    }
  }

  A.entries.malloc(nnzLocalBound);

  if(Comm::World().rank()==0) {printf("Building full IPDG matrix...");fflush(stdout);}

  // reset non-zero counter
  dlong nnz = 0;

  // loop over all elements
  //#pragma omp parallel for
  for(dlong eM=0;eM<mesh.Nelements;++eM){

    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(int n=0;n<mesh.Np;++n){
      for(int m=0;m<mesh.Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh.Np;++i){
          dlong base = eM*mesh.Np*mesh.Nvgeo + i;
          dfloat drdx = mesh.vgeo[base+mesh.Np*mesh.RXID];
          dfloat drdy = mesh.vgeo[base+mesh.Np*mesh.RYID];
          dfloat drdz = mesh.vgeo[base+mesh.Np*mesh.RZID];
          dfloat dsdx = mesh.vgeo[base+mesh.Np*mesh.SXID];
          dfloat dsdy = mesh.vgeo[base+mesh.Np*mesh.SYID];
          dfloat dsdz = mesh.vgeo[base+mesh.Np*mesh.SZID];
          dfloat dtdx = mesh.vgeo[base+mesh.Np*mesh.TXID];
          dfloat dtdy = mesh.vgeo[base+mesh.Np*mesh.TYID];
          dfloat dtdz = mesh.vgeo[base+mesh.Np*mesh.TZID];
          dfloat JW   = mesh.vgeo[base+mesh.Np*mesh.JWID];

          int idn = n*mesh.Np+i;
          int idm = m*mesh.Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
          dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm] + dtdx*Bt[idm];
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm] + dtdy*Bt[idm];
          dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm] + dtdz*Bt[idm];
          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy+dlndz*dlmdz);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh.Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh.Nfp;++i){
            int vidM = mesh.faceNodes[i+fM*mesh.Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh.Np*mesh.Nvgeo + vidM;
            dfloat drdxM = mesh.vgeo[baseM+mesh.Np*mesh.RXID];
            dfloat drdyM = mesh.vgeo[baseM+mesh.Np*mesh.RYID];
            dfloat drdzM = mesh.vgeo[baseM+mesh.Np*mesh.RZID];
            dfloat dsdxM = mesh.vgeo[baseM+mesh.Np*mesh.SXID];
            dfloat dsdyM = mesh.vgeo[baseM+mesh.Np*mesh.SYID];
            dfloat dsdzM = mesh.vgeo[baseM+mesh.Np*mesh.SZID];
            dfloat dtdxM = mesh.vgeo[baseM+mesh.Np*mesh.TXID];
            dfloat dtdyM = mesh.vgeo[baseM+mesh.Np*mesh.TYID];
            dfloat dtdzM = mesh.vgeo[baseM+mesh.Np*mesh.TZID];

            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh.Nfp*mesh.Nfaces+fM*mesh.Nfp+i;
            int vidP    = (int) (mesh.vmapP[idM]%mesh.Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh.vmapP[idM]/mesh.Np;
            dlong baseP   = localEP*mesh.Np*mesh.Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh.vgeo[baseP+mesh.Np*mesh.RXID];
            dfloat drdyP = mesh.vgeo[baseP+mesh.Np*mesh.RYID];
            dfloat drdzP = mesh.vgeo[baseP+mesh.Np*mesh.RZID];
            dfloat dsdxP = mesh.vgeo[baseP+mesh.Np*mesh.SXID];
            dfloat dsdyP = mesh.vgeo[baseP+mesh.Np*mesh.SYID];
            dfloat dsdzP = mesh.vgeo[baseP+mesh.Np*mesh.SZID];
            dfloat dtdxP = mesh.vgeo[baseP+mesh.Np*mesh.TXID];
            dfloat dtdyP = mesh.vgeo[baseP+mesh.Np*mesh.TYID];
            dfloat dtdzP = mesh.vgeo[baseP+mesh.Np*mesh.TZID];

            // grab surface geometric factors
            dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
            dfloat nx = mesh.sgeo[base+mesh.NXID];
            dfloat ny = mesh.sgeo[base+mesh.NYID];
            dfloat nz = mesh.sgeo[base+mesh.NZID];
            dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
            dfloat hinv = mesh.sgeo[base+mesh.IHID];

            // form negative trace terms in IPDG
            int idnM = n*mesh.Np+vidM;
            int idmM = m*mesh.Np+vidM;
            int idmP = m*mesh.Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
            dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM] + dtdxM*Bt[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM] + dtdyM*Bt[idmM];
            dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM] + dtdzM*Bt[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
            dfloat lmM = B[idmM];

            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP] + dtdxP*Bt[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP] + dtdyP*Bt[idmP];
            dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP] + dtdzP*Bt[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;
            dfloat lmP = B[idmP];

            dfloat penalty = tau*hinv;

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            dlong eP    = mesh.EToE[eM*mesh.Nfaces+fM];
            if (eP<0) {
              int qSgn, gradqSgn;
              int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag
              int bcType = BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              Anm += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              Anm += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              Anm += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-)
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
          }
          if(std::abs(AnmP)>tol){
            //#pragma omp critical
            {
              // remote info
              dlong eP    = mesh.EToE[eM*mesh.Nfaces+fM];
              A.entries[nnz].row = globalIds[eM*mesh.Np + n];
              A.entries[nnz].col = globalIds[eP*mesh.Np + m];
              A.entries[nnz].val = AnmP;
              ++nnz;
            }
          }
        }
        if(std::abs(Anm)>tol){
          //#pragma omp critical
          {
            // local block
            A.entries[nnz].row = globalIds[eM*mesh.Np+n];
            A.entries[nnz].col = globalIds[eM*mesh.Np+m];
            A.entries[nnz].val = Anm;
            ++nnz;
          }
        }
      }
    }
  }

  // sort received non-zero entries by row block
  sort(A.entries.ptr(), A.entries.ptr()+nnz,
      [](const parAlmond::parCOO::nonZero_t& a,
         const parAlmond::parCOO::nonZero_t& b) {
        if (a.row < b.row) return true;
        if (a.row > b.row) return false;

        return a.col < b.col;
      });

  //*A = (parAlmond::parCOO::nonZero_t*) realloc(*A, nnz*sizeof(parAlmond::parCOO::nonZero_t));
  A.nnz = nnz;

  if(Comm::World().rank()==0) printf("done.\n");
}
