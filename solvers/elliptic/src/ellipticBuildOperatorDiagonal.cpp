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

void elliptic_t::BuildOperatorDiagonal(memory<dfloat>& diagA){

  if(Comm::World().rank()==0) {printf("Building diagonal...");fflush(stdout);}

  if (settings.compareSetting("DISCRETIZATION","IPDG")) {
    switch(mesh.elementType){
      case Mesh::TRIANGLES:
      {
        if(mesh.dim==2)
          BuildOperatorDiagonalIpdgTri2D(diagA);
        else
          BuildOperatorDiagonalIpdgTri3D(diagA);
        break;
      }
      case Mesh::QUADRILATERALS:
        BuildOperatorDiagonalIpdgQuad2D(diagA);
        break;
      case Mesh::TETRAHEDRA:
        BuildOperatorDiagonalIpdgTet3D(diagA);
        break;
      case Mesh::HEXAHEDRA:
        BuildOperatorDiagonalIpdgHex3D(diagA);
        break;
    }
  } else if (settings.compareSetting("DISCRETIZATION","CONTINUOUS")) {
    memory<dfloat> diagAL(mesh.Np*mesh.Nelements);

    switch(mesh.elementType){
      case Mesh::TRIANGLES:
        BuildOperatorDiagonalContinuousTri2D(diagAL);
        break;
      case Mesh::QUADRILATERALS:
      {
        if(mesh.dim==2)
          BuildOperatorDiagonalContinuousQuad2D(diagAL);
        else
          BuildOperatorDiagonalContinuousQuad3D(diagAL);
        break;
      }
      case Mesh::TETRAHEDRA:
        BuildOperatorDiagonalContinuousTet3D(diagAL);
        break;
      case Mesh::HEXAHEDRA:
        BuildOperatorDiagonalContinuousHex3D(diagAL);
        break;
    }

    //gather the diagonal to assemble it
    ogsMasked.Gather(diagA, diagAL, 1, ogs::Add, ogs::Trans);
  }
  if(Comm::World().rank()==0) printf("done.\n");
}

void elliptic_t::BuildOperatorDiagonalIpdgTri2D(memory<dfloat>& A) {

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(mesh.Nfaces*mesh.Nfp*mesh.Nfp);
  for (int f=0;f<mesh.Nfaces;f++) {
    for (int n=0;n<mesh.Nfp;n++) {
      int fn = mesh.faceNodes[f*mesh.Nfp+n];

      for (int m=0;m<mesh.Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh.Np;i++){
          MSnm += mesh.MM[fn+i*mesh.Np]*mesh.LIFT[i*mesh.Nfp*mesh.Nfaces+f*mesh.Nfp+m];
        }

        MS[m+n*mesh.Nfp + f*mesh.Nfp*mesh.Nfp]  = MSnm;
      }
    }
  }

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong vbase = eM*mesh.Nvgeo;
    dfloat drdx = mesh.vgeo[vbase+mesh.RXID];
    dfloat drdy = mesh.vgeo[vbase+mesh.RYID];
    dfloat dsdx = mesh.vgeo[vbase+mesh.SXID];
    dfloat dsdy = mesh.vgeo[vbase+mesh.SYID];
    dfloat J = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n]  = J*lambda*mesh.MM[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdx*drdx*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdx*dsdx*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*dsdx*dsdx*mesh.Sss[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdy*drdy*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdy*dsdy*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*dsdy*dsdy*mesh.Sss[n*mesh.Np+n];
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
      }
    }

    for (int fM=0;fM<mesh.Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh.Nsgeo*(eM*mesh.Nfaces+fM);
      dfloat nx = mesh.sgeo[sid+mesh.NXID];
      dfloat ny = mesh.sgeo[sid+mesh.NYID];
      dfloat sJ = mesh.sgeo[sid+mesh.SJID];
      dfloat hinv = mesh.sgeo[sid+mesh.IHID];

      int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag

      dfloat penalty = tau*hinv;

      int bcD = 0, bcN =0;
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

      // mass matrix for this face
      memory<dfloat> MSf = MS+fM*mesh.Nfp*mesh.Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<mesh.Nfp;++n){
        int nM = mesh.faceNodes[fM*mesh.Nfp+n];

        for(int m=0;m<mesh.Nfp;++m){
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];
          if (mM == nM) {
            // OP11 = OP11 + 0.5*( gtau*mmE )
            dfloat MSfnm = sJ*MSf[n*mesh.Nfp+m];
            A[eM*mesh.Np+nM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          }
        }
      }

      // now add differential surface terms
      for(int n=0;n<mesh.Nfp;++n){
        int nM = mesh.faceNodes[fM*mesh.Nfp+n];

        for(int i=0;i<mesh.Nfp;++i){
          int iM = mesh.faceNodes[fM*mesh.Nfp+i];

          dfloat MSfni = sJ*MSf[n*mesh.Nfp+i]; // surface Jacobian built in

          dfloat DxMim = drdx*mesh.Dr[iM*mesh.Np+nM] + dsdx*mesh.Ds[iM*mesh.Np+nM];
          dfloat DyMim = drdy*mesh.Dr[iM*mesh.Np+nM] + dsdy*mesh.Ds[iM*mesh.Np+nM];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[eM*mesh.Np+nM] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          A[eM*mesh.Np+nM] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
        }
      }

      for(int n=0;n<mesh.Np;++n){
        for(int m=0;m<mesh.Nfp;++m){
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];

          if (mM==n) {
            for(int i=0;i<mesh.Nfp;++i){
              int iM = mesh.faceNodes[fM*mesh.Nfp+i];

              dfloat MSfim = sJ*MSf[i*mesh.Nfp+m];

              dfloat DxMin = drdx*mesh.Dr[iM*mesh.Np+n] + dsdx*mesh.Ds[iM*mesh.Np+n];
              dfloat DyMin = drdy*mesh.Dr[iM*mesh.Np+n] + dsdy*mesh.Ds[iM*mesh.Np+n];

              // OP11 = OP11 + (- Dn1'*mmE );
              A[eM*mesh.Np+n] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
              A[eM*mesh.Np+n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
            }
          }
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalIpdgTri3D(memory<dfloat>& A) {

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(mesh.Nfaces*mesh.Nfp*mesh.Nfp);
  for (int f=0;f<mesh.Nfaces;f++) {
    for (int n=0;n<mesh.Nfp;n++) {
      int fn = mesh.faceNodes[f*mesh.Nfp+n];

      for (int m=0;m<mesh.Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh.Np;i++){
          MSnm += mesh.MM[fn+i*mesh.Np]*mesh.LIFT[i*mesh.Nfp*mesh.Nfaces+f*mesh.Nfp+m];
        }

        MS[m+n*mesh.Nfp + f*mesh.Nfp*mesh.Nfp]  = MSnm;
      }
    }
  }

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong vbase = eM*mesh.Nvgeo;
    dfloat drdx = mesh.vgeo[vbase+mesh.RXID];
    dfloat drdy = mesh.vgeo[vbase+mesh.RYID];
    dfloat drdz = mesh.vgeo[vbase+mesh.RZID];
    dfloat dsdx = mesh.vgeo[vbase+mesh.SXID];
    dfloat dsdy = mesh.vgeo[vbase+mesh.SYID];
    dfloat dsdz = mesh.vgeo[vbase+mesh.SZID];
    dfloat J = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n]  = J*lambda*mesh.MM[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdx*drdx*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdx*dsdx*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*dsdx*dsdx*mesh.Sss[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdy*drdy*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdy*dsdy*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*dsdy*dsdy*mesh.Sss[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdz*drdz*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*drdz*dsdz*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*dsdz*dsdz*mesh.Sss[n*mesh.Np+n];
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
      }
    }

    for (int fM=0;fM<mesh.Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh.Nsgeo*(eM*mesh.Nfaces+fM);
      dfloat nx = mesh.sgeo[sid+mesh.NXID];
      dfloat ny = mesh.sgeo[sid+mesh.NYID];
      dfloat nz = mesh.sgeo[sid+mesh.NZID];
      dfloat sJ = mesh.sgeo[sid+mesh.SJID];
      dfloat hinv = mesh.sgeo[sid+mesh.IHID];

      int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag

      dfloat penalty = tau*hinv;

      int bcD = 0, bcN =0;
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

      // mass matrix for this face
      memory<dfloat> MSf = MS+fM*mesh.Nfp*mesh.Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<mesh.Nfp;++n){
        int nM = mesh.faceNodes[fM*mesh.Nfp+n];

        for(int m=0;m<mesh.Nfp;++m){
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];
          if (mM == nM) {
            // OP11 = OP11 + 0.5*( gtau*mmE )
            dfloat MSfnm = sJ*MSf[n*mesh.Nfp+m];
            A[eM*mesh.Np+nM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          }
        }
      }

      // now add differential surface terms
      for(int n=0;n<mesh.Nfp;++n){
        int nM = mesh.faceNodes[fM*mesh.Nfp+n];

        for(int i=0;i<mesh.Nfp;++i){
          int iM = mesh.faceNodes[fM*mesh.Nfp+i];

          dfloat MSfni = sJ*MSf[n*mesh.Nfp+i]; // surface Jacobian built in

          dfloat DxMim = drdx*mesh.Dr[iM*mesh.Np+nM] + dsdx*mesh.Ds[iM*mesh.Np+nM];
          dfloat DyMim = drdy*mesh.Dr[iM*mesh.Np+nM] + dsdy*mesh.Ds[iM*mesh.Np+nM];
          dfloat DzMim = drdz*mesh.Dr[iM*mesh.Np+nM] + dsdz*mesh.Ds[iM*mesh.Np+nM];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[eM*mesh.Np+nM] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          A[eM*mesh.Np+nM] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
          A[eM*mesh.Np+nM] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;
        }
      }

      for(int n=0;n<mesh.Np;++n){
        for(int m=0;m<mesh.Nfp;++m){
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];

          if (mM==n) {
            for(int i=0;i<mesh.Nfp;++i){
              int iM = mesh.faceNodes[fM*mesh.Nfp+i];

              dfloat MSfim = sJ*MSf[i*mesh.Nfp+m];

              dfloat DxMin = drdx*mesh.Dr[iM*mesh.Np+n] + dsdx*mesh.Ds[iM*mesh.Np+n];
              dfloat DyMin = drdy*mesh.Dr[iM*mesh.Np+n] + dsdy*mesh.Ds[iM*mesh.Np+n];
              dfloat DzMin = drdz*mesh.Dr[iM*mesh.Np+n] + dsdz*mesh.Ds[iM*mesh.Np+n];

              // OP11 = OP11 + (- Dn1'*mmE );
              A[eM*mesh.Np+n] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
              A[eM*mesh.Np+n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
              A[eM*mesh.Np+n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DzMin*MSfim;
            }
          }
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalContinuousTri2D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n] = J*lambda*mesh.MM[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalIpdgQuad2D(memory<dfloat>& A) {

  // build some monolithic basis arrays (for quads and hexes)
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

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n] = 0;

      // (grad phi_n, grad phi_m)_{D^e}
      for(int i=0;i<mesh.Np;++i){
        dlong base = eM*mesh.Np*mesh.Nvgeo + i;
        dfloat drdx = mesh.vgeo[base+mesh.Np*mesh.RXID];
        dfloat drdy = mesh.vgeo[base+mesh.Np*mesh.RYID];
        dfloat dsdx = mesh.vgeo[base+mesh.Np*mesh.SXID];
        dfloat dsdy = mesh.vgeo[base+mesh.Np*mesh.SYID];
        dfloat JW   = mesh.wJ[eM*mesh.Np + i];

        int idn = n*mesh.Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
        A[eM*mesh.Np+n] += JW*(dlndx*dlndx+dlndy*dlndy);
        A[eM*mesh.Np+n] += lambda*JW*B[idn]*B[idn];
      }

      for (int fM=0;fM<mesh.Nfaces;fM++) {
        // accumulate flux terms for negative and positive traces
        for(int i=0;i<mesh.Nfp;++i){
          int vidM = mesh.faceNodes[i+fM*mesh.Nfp];

          // grab vol geofacs at surface nodes
          dlong baseM = eM*mesh.Np*mesh.Nvgeo + vidM;
          dfloat drdxM = mesh.vgeo[baseM+mesh.Np*mesh.RXID];
          dfloat drdyM = mesh.vgeo[baseM+mesh.Np*mesh.RYID];
          dfloat dsdxM = mesh.vgeo[baseM+mesh.Np*mesh.SXID];
          dfloat dsdyM = mesh.vgeo[baseM+mesh.Np*mesh.SYID];

          // grab surface geometric factors
          dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
          dfloat nx = mesh.sgeo[base+mesh.NXID];
          dfloat ny = mesh.sgeo[base+mesh.NYID];
          dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
          dfloat hinv = mesh.sgeo[base+mesh.IHID];

          // form negative trace terms in IPDG
          int idnM = n*mesh.Np+vidM;

          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
          dfloat lnM = B[idnM];

          dfloat penalty = tau*hinv;
          int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag

          int bcD = 0, bcN =0;
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

          A[eM*mesh.Np+n] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
          A[eM*mesh.Np+n] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
          A[eM*mesh.Np+n] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalContinuousQuad2D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    for (int ny=0;ny<mesh.Nq;ny++) {
      for (int nx=0;nx<mesh.Nq;nx++) {
        int iid = nx+ny*mesh.Nq;
        if (mapB[nx+ny*mesh.Nq+eM*mesh.Np]!=1) {
          A[eM*mesh.Np+iid] = 0;

          for (int k=0;k<mesh.Nq;k++) {
            int id = k+ny*mesh.Nq;
            dfloat Grr = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G00ID*mesh.Np];
            A[eM*mesh.Np+iid] += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[nx+k*mesh.Nq];
          }

          for (int k=0;k<mesh.Nq;k++) {
            int id = nx+k*mesh.Nq;
            dfloat Gss = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G11ID*mesh.Np];
            A[eM*mesh.Np+iid] += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[ny+k*mesh.Nq];
          }

          int id = nx+ny*mesh.Nq;
          dfloat Grs = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
          A[eM*mesh.Np+iid] += 2*Grs*mesh.D[nx+nx*mesh.Nq]*mesh.D[ny+ny*mesh.Nq];

          dfloat JW = mesh.wJ[eM*mesh.Np + iid];
          A[eM*mesh.Np+iid] += JW*lambda;

        } else {
          A[eM*mesh.Np+iid] = 1; //just put a 1 so A is invertable
        }
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}


void elliptic_t::BuildOperatorDiagonalIpdgQuad3D(memory<dfloat>& A) {

  // build some monolithic basis arrays (for quads and hexes)
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

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n] = 0;

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
        dfloat JW   = mesh.wJ[eM*mesh.Np + i];

        int idn = n*mesh.Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx;
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy;
        dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz;
        A[eM*mesh.Np+n] += JW*(dlndx*dlndx + dlndy*dlndy + dlndz*dlndz);
        A[eM*mesh.Np+n] += lambda*JW*B[idn]*B[idn];
      }

      for (int fM=0;fM<mesh.Nfaces;fM++) {
        // accumulate flux terms for negative and positive traces
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

          // grab surface geometric factors
          dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
          dfloat nx = mesh.sgeo[base+mesh.NXID];
          dfloat ny = mesh.sgeo[base+mesh.NYID];
          dfloat nz = mesh.sgeo[base+mesh.NZID];
          dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
          dfloat hinv = mesh.sgeo[base+mesh.IHID];

          // form negative trace terms in IPDG
          int idnM = n*mesh.Np+vidM;

          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM;
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM;
          dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM;
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM ;
          dfloat lnM = B[idnM];

          dfloat penalty = tau*hinv;

          A[eM*mesh.Np+n] += -0.5*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
          A[eM*mesh.Np+n] += -0.5*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
          A[eM*mesh.Np+n] += +0.5*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalContinuousQuad3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    for (int ny=0;ny<mesh.Nq;ny++) {
      for (int nx=0;nx<mesh.Nq;nx++) {
        int iid = nx+ny*mesh.Nq;
        A[eM*mesh.Np+iid] = 0;

        for (int k=0;k<mesh.Nq;k++) {
          int id = k+ny*mesh.Nq;
          dfloat Grr = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G00ID*mesh.Np];
          A[eM*mesh.Np+iid] += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[nx+k*mesh.Nq];
        }

        for (int k=0;k<mesh.Nq;k++) {
          int id = nx+k*mesh.Nq;
          dfloat Gss = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G11ID*mesh.Np];
          A[eM*mesh.Np+iid] += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[ny+k*mesh.Nq];
        }

        int id = nx+ny*mesh.Nq;
        dfloat Grs = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G01ID*mesh.Np];
        A[eM*mesh.Np+iid] += 2*Grs*mesh.D[nx+nx*mesh.Nq]*mesh.D[ny+ny*mesh.Nq];

        // id = nx+ny*mesh.Nq;
        // dfloat Grt = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G02ID*mesh.Np];
        // A[eM*mesh.Np+iid] += 2*Grt*mesh.D[nx+nx*mesh.Nq];

        // id = nx+ny*mesh.Nq;
        // dfloat Gst = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G12ID*mesh.Np];
        // A[eM*mesh.Np+iid] += 2*Gst*mesh.D[ny+ny*mesh.Nq];

        // dfloat Gtt = mesh.ggeo[eM*mesh.Np*mesh.Nggeo + id + mesh.G22ID*mesh.Np];
        // A[eM*mesh.Np+iid] += Gtt;



        dfloat JW  = mesh.wJ[eM*mesh.Np + iid];
        A[eM*mesh.Np+iid] += JW*lambda;
      }
    }
  }
}



void elliptic_t::BuildOperatorDiagonalIpdgTet3D(memory<dfloat>& A) {

  // surface mass matrices MS = MM*LIFT
  memory<dfloat> MS(mesh.Nfaces*mesh.Nfp*mesh.Nfp);
  for (int f=0;f<mesh.Nfaces;f++) {
    for (int n=0;n<mesh.Nfp;n++) {
      int fn = mesh.faceNodes[f*mesh.Nfp+n];

      for (int m=0;m<mesh.Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<mesh.Np;i++){
          MSnm += mesh.MM[fn+i*mesh.Np]*mesh.LIFT[i*mesh.Nfp*mesh.Nfaces+f*mesh.Nfp+m];
        }

        MS[m+n*mesh.Nfp + f*mesh.Nfp*mesh.Nfp]  = MSnm;
      }
    }
  }

  for(dlong eM=0;eM<mesh.Nelements;++eM){
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
    dfloat J = mesh.wJ[eM];

    dfloat G00 = drdx*drdx + drdy*drdy + drdz*drdz;
    dfloat G01 = drdx*dsdx + drdy*dsdy + drdz*dsdz;
    dfloat G02 = drdx*dtdx + drdy*dtdy + drdz*dtdz;

    dfloat G11 = dsdx*dsdx + dsdy*dsdy + dsdz*dsdz;
    dfloat G12 = dsdx*dtdx + dsdy*dtdy + dsdz*dtdz;
    dfloat G22 = dtdx*dtdx + dtdy*dtdy + dtdz*dtdz;


    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n]  = J*lambda*mesh.MM[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G00*mesh.Srr[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G01*mesh.Srs[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G02*mesh.Srt[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G11*mesh.Sss[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G12*mesh.Sst[n*mesh.Np+n];
      A[eM*mesh.Np+n] += J*G22*mesh.Stt[n*mesh.Np+n];
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
      }
    }

    for (int fM=0;fM<mesh.Nfaces;fM++) {
      // load surface geofactors for this face
      dlong sid = mesh.Nsgeo*(eM*mesh.Nfaces+fM);
      dfloat nx = mesh.sgeo[sid+mesh.NXID];
      dfloat ny = mesh.sgeo[sid+mesh.NYID];
      dfloat nz = mesh.sgeo[sid+mesh.NZID];
      dfloat sJ = mesh.sgeo[sid+mesh.SJID];
      dfloat hinv = mesh.sgeo[sid+mesh.IHID];

      int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag

      dfloat penalty = tau*hinv;

      int bcD = 0, bcN =0;
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

      // mass matrix for this face
      memory<dfloat> MSf = MS+fM*mesh.Nfp*mesh.Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<mesh.Nfp;++n){
        for(int m=0;m<mesh.Nfp;++m){
          int nM = mesh.faceNodes[fM*mesh.Nfp+n];
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];

          if (mM==nM) {
            // OP11 = OP11 + 0.5*( gtau*mmE )
            dfloat MSfnm = sJ*MSf[n*mesh.Nfp+m];
            A[eM*mesh.Np+nM] += 0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          }
        }
      }

      // now add differential surface terms
      for(int n=0;n<mesh.Nfp;++n){
        int nM = mesh.faceNodes[fM*mesh.Nfp+n];

        for(int i=0;i<mesh.Nfp;++i){
          int iM = mesh.faceNodes[fM*mesh.Nfp+i];

          dfloat MSfni = sJ*MSf[n*mesh.Nfp+i]; // surface Jacobian built in

          dfloat DxMim = drdx*mesh.Dr[iM*mesh.Np+nM] + dsdx*mesh.Ds[iM*mesh.Np+nM] + dtdx*mesh.Dt[iM*mesh.Np+nM];
          dfloat DyMim = drdy*mesh.Dr[iM*mesh.Np+nM] + dsdy*mesh.Ds[iM*mesh.Np+nM] + dtdy*mesh.Dt[iM*mesh.Np+nM];
          dfloat DzMim = drdz*mesh.Dr[iM*mesh.Np+nM] + dsdz*mesh.Ds[iM*mesh.Np+nM] + dtdz*mesh.Dt[iM*mesh.Np+nM];

          // OP11 = OP11 + 0.5*( - mmE*Dn1)
          A[eM*mesh.Np+nM] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
          A[eM*mesh.Np+nM] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
          A[eM*mesh.Np+nM] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;
        }
      }

      for(int n=0;n<mesh.Np;++n){
        for(int m=0;m<mesh.Nfp;++m){
          int mM = mesh.faceNodes[fM*mesh.Nfp+m];

          if (mM==n) {
            for(int i=0;i<mesh.Nfp;++i){
              int iM = mesh.faceNodes[fM*mesh.Nfp+i];

              dfloat MSfim = sJ*MSf[i*mesh.Nfp+m];

              dfloat DxMin = drdx*mesh.Dr[iM*mesh.Np+n] + dsdx*mesh.Ds[iM*mesh.Np+n] + dtdx*mesh.Dt[iM*mesh.Np+n];
              dfloat DyMin = drdy*mesh.Dr[iM*mesh.Np+n] + dsdy*mesh.Ds[iM*mesh.Np+n] + dtdy*mesh.Dt[iM*mesh.Np+n];
              dfloat DzMin = drdz*mesh.Dr[iM*mesh.Np+n] + dsdz*mesh.Ds[iM*mesh.Np+n] + dtdz*mesh.Dt[iM*mesh.Np+n];

              // OP11 = OP11 + (- Dn1'*mmE );
              A[eM*mesh.Np+n] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
              A[eM*mesh.Np+n] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
              A[eM*mesh.Np+n] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;
            }
          }
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalContinuousTet3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    dlong gbase = eM*mesh.Nggeo;
    dfloat Grr = mesh.ggeo[gbase + mesh.G00ID];
    dfloat Grs = mesh.ggeo[gbase + mesh.G01ID];
    dfloat Grt = mesh.ggeo[gbase + mesh.G02ID];
    dfloat Gss = mesh.ggeo[gbase + mesh.G11ID];
    dfloat Gst = mesh.ggeo[gbase + mesh.G12ID];
    dfloat Gtt = mesh.ggeo[gbase + mesh.G22ID];
    dfloat J   = mesh.wJ[eM];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
        A[eM*mesh.Np+n] = J*lambda*mesh.MM[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grr*mesh.Srr[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grs*mesh.Srs[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Grt*mesh.Srt[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gss*mesh.Sss[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gst*mesh.Sst[n+n*mesh.Np];
        A[eM*mesh.Np+n] += Gtt*mesh.Stt[n+n*mesh.Np];
      } else {
        A[eM*mesh.Np+n] = 1; //just put a 1 so A is invertable
      }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalIpdgHex3D(memory<dfloat>& A) {

  // build some monolithic basis arrays (for quads and hexes)
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

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    /* start with stiffness matrix  */
    for(int n=0;n<mesh.Np;++n){
      A[eM*mesh.Np+n] = 0;

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
        dfloat JW   = mesh.wJ[eM*mesh.Np + i];

        int idn = n*mesh.Np+i;
        dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
        dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
        dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];
        A[eM*mesh.Np+n] += JW*(dlndx*dlndx+dlndy*dlndy+dlndz*dlndz);
        A[eM*mesh.Np+n] += lambda*JW*B[idn]*B[idn];
      }

      for (int fM=0;fM<mesh.Nfaces;fM++) {
        // accumulate flux terms for negative and positive traces
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

          // grab surface geometric factors
          dlong base = mesh.Nsgeo*(eM*mesh.Nfp*mesh.Nfaces + fM*mesh.Nfp + i);
          dfloat nx = mesh.sgeo[base+mesh.NXID];
          dfloat ny = mesh.sgeo[base+mesh.NYID];
          dfloat nz = mesh.sgeo[base+mesh.NZID];
          dfloat wsJ = mesh.sgeo[base+mesh.WSJID];
          dfloat hinv = mesh.sgeo[base+mesh.IHID];

          // form negative trace terms in IPDG
          int idnM = n*mesh.Np+vidM;
          dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
          dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
          dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
          dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
          dfloat lnM = B[idnM];

          dfloat penalty = tau*hinv;
          int bc = mesh.EToB[fM+mesh.Nfaces*eM]; //raw boundary flag

          int bcD = 0, bcN =0;
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

          A[eM*mesh.Np+n] += -0.5*(1+bcD)*(1-bcN)*wsJ*lnM*ndotgradlnM;  // -(ln^-, N.grad lm^-)
          A[eM*mesh.Np+n] += -0.5*(1+bcD)*(1-bcN)*wsJ*ndotgradlnM*lnM;  // -(N.grad ln^-, lm^-)
          A[eM*mesh.Np+n] += +0.5*(1+bcD)*(1-bcN)*wsJ*penalty*lnM*lnM; // +((tau/h)*ln^-,lm^-)
        }
      }
    }
  }
}

void elliptic_t::BuildOperatorDiagonalContinuousHex3D(memory<dfloat>& A) {

  for(dlong eM=0;eM<mesh.Nelements;++eM){
    for (int nz=0;nz<mesh.Nq;nz++) {
    for (int ny=0;ny<mesh.Nq;ny++) {
    for (int nx=0;nx<mesh.Nq;nx++) {
      int idn = nx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
      if (mapB[idn+eM*mesh.Np]!=1) {
        A[eM*mesh.Np+idn] = 0;

        int id = nx+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
        dlong base = eM*mesh.Np*mesh.Nggeo;


        dfloat Grs = mesh.ggeo[base + id + mesh.G01ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Grs*mesh.D[nx+nx*mesh.Nq]*mesh.D[ny+ny*mesh.Nq];

        dfloat Grt = mesh.ggeo[base + id + mesh.G02ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Grt*mesh.D[nx+nx*mesh.Nq]*mesh.D[nz+nz*mesh.Nq];

        dfloat Gst = mesh.ggeo[base + id + mesh.G12ID*mesh.Np];
        A[eM*mesh.Np+idn] += 2*Gst*mesh.D[ny+ny*mesh.Nq]*mesh.D[nz+nz*mesh.Nq];

        for (int k=0;k<mesh.Nq;k++) {
          int iid = k+ny*mesh.Nq+nz*mesh.Nq*mesh.Nq;
          dfloat Grr = mesh.ggeo[base + iid + mesh.G00ID*mesh.Np];
          A[eM*mesh.Np+idn] += Grr*mesh.D[nx+k*mesh.Nq]*mesh.D[nx+k*mesh.Nq];
        }

        for (int k=0;k<mesh.Nq;k++) {
          int iid = nx+k*mesh.Nq+nz*mesh.Nq*mesh.Nq;
          dfloat Gss = mesh.ggeo[base + iid + mesh.G11ID*mesh.Np];
          A[eM*mesh.Np+idn] += Gss*mesh.D[ny+k*mesh.Nq]*mesh.D[ny+k*mesh.Nq];
        }

        for (int k=0;k<mesh.Nq;k++) {
          int iid = nx+ny*mesh.Nq+k*mesh.Nq*mesh.Nq;
          dfloat Gtt = mesh.ggeo[base + iid + mesh.G22ID*mesh.Np];
          A[eM*mesh.Np+idn] += Gtt*mesh.D[nz+k*mesh.Nq]*mesh.D[nz+k*mesh.Nq];
        }

        dfloat JW = mesh.wJ[eM*mesh.Np + idn];
        A[eM*mesh.Np+idn] += JW*lambda;
      } else {
        A[eM*mesh.Np+idn] = 1; //just put a 1 so A is invertable
      }
    }
    }
    }

    //add the rank boost for the allNeumann Poisson problem
    if (allNeumann) {
      for(int n=0;n<mesh.Np;++n){
        if (mapB[n+eM*mesh.Np]!=1) { //dont fill rows for masked nodes
          A[eM*mesh.Np+n] += allNeumannPenalty*allNeumannScale*allNeumannScale;
        }
      }
    }
  }
}
