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

/*Build a graph from mesh connectivity info*/
graph_t::graph_t(platform_t &_platform,
                 const dlong _Nelements,
                 const int _dim,
                 const int _Nverts,
                 const int _Nfaces,
                 const int _NfaceVerts,
                 const memory<int>& faceVertices,
                 const memory<hlong>& EToV,
                 const memory<dfloat>& EX,
                 const memory<dfloat>& EY,
                 const memory<dfloat>& EZ,
                 comm_t _comm):
  platform(_platform),
  Nverts(_Nelements),
  Nelements(_Nelements),
  dim(_dim),
  Nfaces(_Nfaces),
  NelementVerts(_Nverts),
  NfaceVerts(_NfaceVerts) {

  gcomm = _comm.Dup();
  grank = gcomm.rank();
  gsize = gcomm.size();

  comm  = _comm.Dup();
  rank = comm.rank();
  size = comm.size();

  for (int n=0;n<Nfaces*NfaceVerts;++n)
    faceVerts[n] = faceVertices[n];

  /*Global number of elements*/
  NVertsGlobal=static_cast<hlong>(Nverts);
  comm.Allreduce(NVertsGlobal);

  /*Get global element count offsets*/
  hlong localNverts=static_cast<hlong>(Nverts);
  comm.Scan(localNverts, VoffsetU);
  VoffsetL = VoffsetU-Nverts;

  gNVertsGlobal = NVertsGlobal;
  gVoffsetL = VoffsetL;
  gVoffsetU = VoffsetU;

  /*Create array of packed element data*/
  elements.malloc(Nelements);

  if (dim==2) {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        elements[e].EX[v] = EX[v+e*NelementVerts];
        elements[e].EY[v] = EY[v+e*NelementVerts];

        elements[e].V[v] = EToV[v+e*NelementVerts];
      }
      for (int f=0;f<Nfaces;++f) {
        elements[e].E[f] = -1;
        elements[e].F[f] = -1;
      }
    }
  } else {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        elements[e].EX[v] = EX[v+e*NelementVerts];
        elements[e].EY[v] = EY[v+e*NelementVerts];
        elements[e].EZ[v] = EZ[v+e*NelementVerts];

        elements[e].V[v] = EToV[v+e*NelementVerts];
      }
      for (int f=0;f<Nfaces;++f) {
        elements[e].E[f] = -1;
        elements[e].F[f] = -1;
      }
    }
  }
}

/*Globally divide graph into two pieces according to a bipartition*/
void graph_t::Split(const memory<int>& partition) {

  /*Count how much of each partition we have locally*/
  dlong Nverts0=0;
  dlong Nverts1=0;
  for (dlong n=0;n<Nverts;++n) {
    if (partition[n]==0) Nverts0++;
    else                 Nverts1++;
  }

  hlong globalNverts0=static_cast<hlong>(Nverts0);
  hlong globalNverts1=static_cast<hlong>(Nverts1);
  comm.Allreduce(globalNverts0);
  comm.Allreduce(globalNverts1);

  /*Get offsets of partitions on each rank*/
  memory<hlong> starts0(size+1);
  memory<hlong> starts1(size+1);
  starts0[0]=0;
  starts1[0]=0;
  hlong localNverts0 = static_cast<hlong>(Nverts0);
  hlong localNverts1 = static_cast<hlong>(Nverts1);
  comm.Allgather(localNverts0, starts0+1);
  comm.Allgather(localNverts1, starts1+1);

  for(int r=0;r<size;++r) {
    starts0[r+1] += starts0[r];
    starts1[r+1] += starts1[r];
  }

  /*Determine number of ranks to hold left and right partitions*/
  const int size0 = (size+1)/2;
  const int size1 = size-size0;

  const hlong chunk0 = globalNverts0/size0;
  const hlong chunk1 = globalNverts1/size1;

  const int remainder0 = static_cast<int>(globalNverts0 - chunk0*size0);
  const int remainder1 = static_cast<int>(globalNverts1 - chunk1*size1);

  memory<int> Nsend0(size,0);
  memory<int> Nsend1(size,0);
  memory<int> Nrecv0(size);
  memory<int> Nrecv1(size);
  memory<int> sendOffsets0(size);
  memory<int> sendOffsets1(size);
  memory<int> recvOffsets0(size);
  memory<int> recvOffsets1(size);

  memory<hlong> newIds(Nverts+Nhalo);

  /*Determine new ids and send counts*/
  dlong cnt0=0;
  dlong cnt1=0;
  for(dlong e=0;e<Nverts;++e){
    if (partition[e]==0) {
      // new global element index
      const hlong ep = starts0[rank]+cnt0++;
      newIds[e] = ep;

      // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
      int r;
      if(ep<remainder0*(chunk0+1))
        r = ep/(chunk0+1);
      else
        r = remainder0 + ((ep-remainder0*(chunk0+1))/chunk0);

      ++Nsend0[r];
    } else {
      // new global element index
      const hlong ep = starts1[rank]+cnt1++;
      newIds[e] = ep;

      // 0, chunk+1, 2*(chunk+1) ..., remainder*(chunk+1), remainder*(chunk+1) + chunk
      int r;
      if(ep<remainder1*(chunk1+1))
        r = ep/(chunk1+1);
      else
        r = remainder1 + ((ep-remainder1*(chunk1+1))/chunk1);

      ++Nsend1[r+size0];
    }
  }

  starts0.free();
  starts1.free();

  if (L[0].Nglobal) {
    /*If we have connected the elements, share the newIds*/
    L[0].A.halo.Exchange(newIds, 1);

    /*Then update the connectivity*/
    dlong cnt=0;
    for(dlong e=0;e<Nverts;++e){
      const int part = partition[e];
      for (int f=0;f<Nfaces;++f) {
        const hlong gE = elements[e].E[f];
        if (gE!=-1) {
          dlong eN;
          if (gE>=VoffsetL && gE<VoffsetU) { /*local neighbor*/
            eN = static_cast<dlong>(gE-VoffsetL);
          } else { /*halo neighbor*/
            eN = colIds[cnt++]; /*Get the local id in the halo (we make this when building the Laplacian)*/
          }

          const int partN = partition[eN];
          if (partN==part) { /*If both elements are in the same partition*/
            elements[e].E[f] = newIds[eN]; /*Re index*/
          } else {
            elements[e].E[f] = -1;/*else break connections across the partitions*/
          }
        }
      }
    }
  }
  newIds.free();

  // find send offsets
  sendOffsets0[0]=0;
  sendOffsets1[0]=0;
  for(int r=1;r<size;++r) {
    sendOffsets0[r] = sendOffsets0[r-1] + Nsend0[r-1];
    sendOffsets1[r] = sendOffsets1[r-1] + Nsend1[r-1];
  }
  int NsendTotal0=0;
  int NsendTotal1=0;
  for(int r=0;r<size;++r) {
    NsendTotal0 += Nsend0[r];
    NsendTotal1 += Nsend1[r];
  }

  // exchange counts
  comm.Alltoall(Nsend0, Nrecv0);
  comm.Alltoall(Nsend1, Nrecv1);

  // find recv offsets
  recvOffsets0[0]=0;
  recvOffsets1[0]=0;
  for(int r=1;r<size;++r) {
    recvOffsets0[r] = recvOffsets0[r-1] + Nrecv0[r-1];
    recvOffsets1[r] = recvOffsets1[r-1] + Nrecv1[r-1];
  }

  // count incoming clusters
  dlong newNverts = 0;

  if (rank<size0) {
    for(int r=0;r<size;++r) {
      newNverts += Nrecv0[r];
    }
  } else {
    for(int r=0;r<size;++r) {
      newNverts += Nrecv1[r];
    }
  }

  /*make send buffers*/
  memory<element_t> sendElements0(NsendTotal0);
  memory<element_t> sendElements1(NsendTotal1);

  cnt0=0;
  cnt1=0;
  for(dlong e=0;e<Nverts;++e){
    if (partition[e]==0) {
      sendElements0[cnt0++] = elements[e];
    } else {
      sendElements1[cnt1++] = elements[e];
    }
  }

  /*make new list*/
  Nverts = newNverts;
  Nelements = newNverts;
  elements.malloc(Nverts);

  memory<element_t> null;

  // exchange elements
  if (rank<size0) {
    comm.Alltoallv(sendElements0, Nsend0, sendOffsets0,
                        elements, Nrecv0, recvOffsets0);
    comm.Alltoallv(sendElements1, Nsend1, sendOffsets1,
                            null, Nrecv1, recvOffsets1);
  } else {
    comm.Alltoallv(sendElements0, Nsend0, sendOffsets0,
                            null, Nrecv0, recvOffsets0);
    comm.Alltoallv(sendElements1, Nsend1, sendOffsets1,
                        elements, Nrecv1, recvOffsets1);
  }

  comm_t newComm = comm.Split(rank<size0, rank);
  comm.Free();
  comm = newComm;

  rank = comm.rank();
  size = comm.size();

  /*Global number of elements*/
  NVertsGlobal=static_cast<hlong>(Nverts);
  comm.Allreduce(NVertsGlobal);

  /*Get global element count offsets*/
  hlong localNverts=static_cast<hlong>(Nverts);
  comm.Scan(localNverts, VoffsetU);
  VoffsetL = VoffsetU-Nverts;
}

void graph_t::Report() {

  /* Min,Avg,Max Element counts*/
  hlong globalNverts = static_cast<hlong>(Nverts);
  gcomm.Allreduce(globalNverts);
  dfloat avgNverts = static_cast<dfloat>(globalNverts)/gsize;

  dlong minNverts=Nverts;
  dlong maxNverts=Nverts;
  gcomm.Allreduce(minNverts, Comm::Min);
  gcomm.Allreduce(maxNverts, Comm::Max);


  dlong cut=0.0;
  for (dlong n=0;n<Nverts;++n) {
    for (int f=0;f<Nfaces;++f) {
      const hlong eN = elements[n].E[f];
      if (eN!=-1) {
        if ((eN<gVoffsetL) || (eN>=gVoffsetU) ) {
          cut++;
        }
      }
    }
  }

  hlong gCut = static_cast<hlong>(cut);
  gcomm.Allreduce(gCut);
  hlong avgCut = gCut/gsize;

  dlong minCut=cut;
  dlong maxCut=cut;
  gcomm.Allreduce(minCut, Comm::Min);
  gcomm.Allreduce(maxCut, Comm::Max);

  if(grank==0) {
    printf("--------------------------------------ParAdogs Report------------------------------------------\n");
    printf("-----------------------------------------------------------------------------------------------\n");
    printf("   Nranks   |    Elements   |   Per Rank Elements   |   Halo Faces   |   Per Rank Halo Faces  |\n");
    printf("            |               |       (min,avg,max)   |                |         (min,avg,max)  |\n");
    printf("-----------------------------------------------------------------------------------------------\n");
    printf(      "%9d   | %11lld   |       %13lld   | %12lld   |         %13lld  |\n",
            gsize,
            static_cast<long long int>(globalNverts),
            static_cast<long long int>(minNverts),
            static_cast<long long int>(gCut),
            static_cast<long long int>(minCut));
    printf("            |               |       %13lld   |                |         %13lld  |\n",
            static_cast<long long int>(avgNverts),
            static_cast<long long int>(avgCut));
    printf("            |               |       %13lld   |                |         %13lld  |\n",
            static_cast<long long int>(maxNverts),
            static_cast<long long int>(maxCut));
    printf("-----------------------------------------------------------------------------------------------\n");
  }
}

void graph_t::ExtractMesh(dlong &Nelements_,
                          memory<hlong>& EToV,
                          memory<hlong>& EToE,
                          memory<int>& EToF,
                          memory<dfloat>& EX,
                          memory<dfloat>& EY,
                          memory<dfloat>& EZ) {

  /*Destroy any exiting mesh data and create new data from current graph*/
  Nelements_ = Nelements;

  EToV.malloc(Nelements*NelementVerts);
  EToE.malloc(Nelements*NelementVerts);
  EToF.malloc(Nelements*NelementVerts);

  EX.malloc(Nelements*NelementVerts);
  EY.malloc(Nelements*NelementVerts);
  if (dim==3)
    EZ.malloc(Nelements*NelementVerts);

  if (dim==2) {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        EToV[v+e*NelementVerts] = elements[e].V[v];
        EX[v+e*NelementVerts] = elements[e].EX[v];
        EY[v+e*NelementVerts] = elements[e].EY[v];
      }
      for (int f=0;f<Nfaces;++f) {
        EToE[f+e*Nfaces] = elements[e].E[f];
        EToF[f+e*Nfaces] = elements[e].F[f];
      }
    }
  } else {
    for (dlong e=0;e<Nelements;++e) {
      for (int v=0;v<NelementVerts;++v) {
        EToV[v+e*NelementVerts] = elements[e].V[v];
        EX[v+e*NelementVerts] = elements[e].EX[v];
        EY[v+e*NelementVerts] = elements[e].EY[v];
        EZ[v+e*NelementVerts] = elements[e].EZ[v];
      }
      for (int f=0;f<Nfaces;++f) {
        EToE[f+e*Nfaces] = elements[e].E[f];
        EToF[f+e*Nfaces] = elements[e].F[f];
      }
    }
  }
}

} //namespace paradogs

} //namespace libp
