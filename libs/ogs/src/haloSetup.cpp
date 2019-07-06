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

#include "ogs.hpp"
#include "ogsKernels.hpp"

halo_t *halo_t::Setup(dlong N, hlong *ids, MPI_Comm &comm,
                      int verbose, occa::device &device){

  halo_t *halo = new halo_t(comm, device);

  halo->N = N;
  halo->Nlocal = 0;
  halo->Nsend = 0;
  halo->Nhalo = 0;

  halo->haloSendIds = NULL;

  halo->gshHalo = NULL;

  //Keep track of how many gs handles we've created, and
  // build kernels if this is the first
  if (!ogs::Nrefs) ogs::initKernels(comm, device);
  ogs::Nrefs++;

  int rank, size;
  MPI_Comm_rank(halo->comm, &rank);
  MPI_Comm_size(halo->comm, &size);

  hlong *flagIds = (hlong *) calloc(N,sizeof(hlong));
  for (dlong i=0;i<N;i++) flagIds[i] = abs(ids[i]);

  //make a host gs handle (calls gslib)
  void *gsHandle = ogs::gsSetup(comm, N, flagIds, 0, 0);

  //use the host gs to find what nodes are local to this rank
  int *minRank = (int *) calloc(N,sizeof(int));
  int *maxRank = (int *) calloc(N,sizeof(int));
  for (dlong i=0;i<N;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  ogs::gsGatherScatter(minRank, 1, 1, 0, ogs_int, ogs_min, gsHandle); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogs::gsGatherScatter(maxRank, 1, 1, 0, ogs_int, ogs_max, gsHandle); //maxRank[n] contains the largest rank taking part in the gather of node n

  //discard the large gs handle
  ogs::gsFree(gsHandle); free(flagIds);

  //count local and halo nodes
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if (ids[i]>0) { //owned node
      halo->Nlocal++;
      if ((minRank[i]!=rank)||(maxRank[i]!=rank))
        halo->Nsend++; //need to send this node out
    } else { //node in halo
      halo->Nhalo++;
    }
  }

  //record what ids we need to pack into the outgoing buffer
  // also make the list of ids participating in the gsop
  halo->haloSendIds = (dlong*) calloc(halo->Nsend+1,sizeof(dlong)); //extra entry so the occa buffer will actually exist
  hlong *nonSymIds = (hlong *) calloc(halo->Nsend+halo->Nhalo,sizeof(hlong));
  dlong cnt = 0;
  dlong cnt2 = 0;
  for (dlong i=0;i<N;i++) {
    if (ids[i]==0) continue;

    if (ids[i]>0) { //owned node
      if ((minRank[i]!=rank)||(maxRank[i]!=rank)) {
        halo->haloSendIds[cnt++] = i;
        nonSymIds[cnt2++] = ids[i];
      }
    } else {
      nonSymIds[cnt2++] = ids[i];
    }
  }
  halo->o_haloSendIds = device.malloc((halo->Nsend+1)*sizeof(dlong), halo->haloSendIds);

  //make a host gs handle
  halo->gshHalo = ogs::gsSetup(comm, halo->Nsend+halo->Nhalo, nonSymIds, 0,0);

  free(nonSymIds);
  free(minRank); free(maxRank);

  halo->hostBuf = NULL;
  halo->haloBuf = NULL;
  halo->hostBufSize = 0;

  return halo;
}


void halo_t::Free() {

  if (Nsend) {
    free(haloSendIds);
    o_haloSendIds.free();
    Nsend = 0;
  }

  ogs::gsFree(gshHalo);

  ogs::Nrefs--;
  if (!ogs::Nrefs) ogs::freeKernels();
}

void halo_t::reallocHostBuffer(size_t Nbytes) {
  if (Nsend+Nhalo) {
    if (hostBufSize < (Nsend+Nhalo)*Nbytes) {
      if (hostBufSize) free(hostBuf);
      hostBuf = (void *) malloc((Nsend+Nhalo)*Nbytes);
      hostBufSize = (Nsend+Nhalo)*Nbytes;
    }
  }
}

void halo_t::reallocOccaBuffer(size_t Nbytes) {
  if (Nsend+Nhalo) {
    if (o_haloBuf.size() < (Nsend+Nhalo)*Nbytes) {
      if (o_haloBuf.size()) o_haloBuf.free();
      haloBuf = occaHostMallocPinned(device, (Nsend+Nhalo)*Nbytes, NULL,
                                     o_haloBuf, h_haloBuf);
    }
  }
}