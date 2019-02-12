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

#include "adaptive.h"

int parallelCompareRowColumn(const void *a, const void *b);

void adaptiveBuildIpdgHex3D(adaptive_t *adaptive, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);



void adaptiveBuildIpdg(adaptive_t *adaptive, int basisNp, dfloat *basis,
                       dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  adaptiveBuildIpdgHex3D(adaptive, basisNp, basis,lambda, A, nnzA, globalStarts);

}

void adaptiveBuildIpdgHex3D(adaptive_t *adaptive, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = adaptive->mesh;
  setupAide options = adaptive->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
    rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nk=0;nk<mesh->N+1;++nk){
    for(int nj=0;nj<mesh->N+1;++nj){
      for(int ni=0;ni<mesh->N+1;++ni){

        int node = 0;

        for(int k=0;k<mesh->N+1;++k){
          for(int j=0;j<mesh->N+1;++j){
            for(int i=0;i<mesh->N+1;++i){

              if(nk==k && nj==j && ni==i)
                B[mode*mesh->Np+node] = 1;
              if(nj==j && nk==k)
                Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
              if(ni==i && nk==k)
                Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
              if(ni==i && nj==j)
                Bt[mode*mesh->Np+node] = mesh->D[nk+mesh->Nq*k]; 
              
              ++node;
            }
          }
        }
        
        ++mode;
      }
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));
  
  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // reset non-zero counter
  dlong nnz = 0;
  
  // loop over all elements
  //#pragma omp parallel for
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          dlong base = eM*mesh->Np*mesh->Nvgeo + i;
          dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
          dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
          dfloat drdz = mesh->vgeo[base+mesh->Np*RZID];
          dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
          dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
          dfloat dsdz = mesh->vgeo[base+mesh->Np*SZID];
          dfloat dtdx = mesh->vgeo[base+mesh->Np*TXID];
          dfloat dtdy = mesh->vgeo[base+mesh->Np*TYID];
          dfloat dtdz = mesh->vgeo[base+mesh->Np*TZID];
          dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
          
          int idn = n*mesh->Np+i;
          int idm = m*mesh->Np+i;
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
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
            dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
            dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
            dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
            dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];

            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP    = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh->vmapP[idM]/mesh->Np;
            dlong baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat drdzP = mesh->vgeo[baseP+mesh->Np*RZID];
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
            dfloat dsdzP = mesh->vgeo[baseP+mesh->Np*SZID];
            dfloat dtdxP = mesh->vgeo[baseP+mesh->Np*TXID];
            dfloat dtdyP = mesh->vgeo[baseP+mesh->Np*TYID];
            dfloat dtdzP = mesh->vgeo[baseP+mesh->Np*TZID];
            
            // grab surface geometric factors
            dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat nz = mesh->sgeo[base+NZID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];
            
            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM; 
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

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
            
            dfloat penalty = adaptive->tau*hinv;     

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP<0) {
              int qSgn, gradqSgn;
              int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
              int bcType = adaptive->BCType[bc];          //find its type (Dirichlet/Neumann)
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
          if(fabs(AnmP)>tol){
            //#pragma omp critical
            {
              // remote info
              dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
              (*A)[nnz].row = globalIds[eM*mesh->Np + n];
              (*A)[nnz].col = globalIds[eP*mesh->Np + m];
              (*A)[nnz].val = AnmP;
              (*A)[nnz].ownerRank = rankM;
              ++nnz;
            }
          } 
        }
        if(fabs(Anm)>tol){
          //#pragma omp critical
          {
            // local block
            (*A)[nnz].row = globalIds[eM*mesh->Np+n];
            (*A)[nnz].col = globalIds[eM*mesh->Np+m];
            (*A)[nnz].val = Anm;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
  }

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

  free(globalIds);
  free(B);  free(Br); free(Bs); free(Bt);
}
