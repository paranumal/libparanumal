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

#include "mesh.hpp"

// send data from partition boundary elements
// and receive data to ghost elements
void mesh_t::HaloExchange(size_t Nbytes,         // message size per element
                          void *sourceBuffer,
                          void *_sendBuffer,    // temporary buffer
                          void *_recvBuffer){

  // count outgoing and incoming meshes
  int tag = 999;

  // copy data from outgoing elements into temporary send buffer
  for(int i=0;i<totalHaloPairs;++i){
    // outgoing element
    int e = haloElementList[i];
    // copy element e data to _sendBuffer
    memcpy(((char*)_sendBuffer)+i*Nbytes, ((char*)sourceBuffer)+e*Nbytes, Nbytes);
  }

  // initiate immediate send  and receives to each other process as needed
  int offset = 0, message = 0;
  for(int rr=0;rr<size;++rr){
    if(rr!=rank){
      size_t count = NhaloPairs[rr]*Nbytes;
      if(count){
        //      printf("rank %d sending %d bytes to rank %d\n", rank, count, rr);
        MPI_Irecv(((char*)_recvBuffer)+offset, count, MPI_CHAR, rr, tag,
                  comm, haloRecvRequests+message);

        MPI_Isend(((char*)_sendBuffer)+offset, count, MPI_CHAR, rr, tag,
                  comm, haloSendRequests+message);
        offset += count;
        ++message;
      }
    }
  }
  //  printf("NhaloMessages = %d\n", NhaloMessages);

  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Waitall(NhaloMessages, haloRecvRequests, recvStatus);
  MPI_Waitall(NhaloMessages, haloSendRequests, sendStatus);
}


// start halo exchange (for q)
void mesh_t::HaloExchangeStart(size_t Nbytes,       // message size per element
                             void *_sendBuffer,    // temporary buffer
                             void *_recvBuffer){

  if(totalHaloPairs>0){
    // count outgoing and incoming meshes
    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(int rr=0;rr<size;++rr){
      if(rr!=rank){
        size_t count = NhaloPairs[rr]*Nbytes;
        if(count){
          MPI_Irecv(((char*)_recvBuffer)+offset, count, MPI_CHAR, rr, tag,
                    comm, (haloRecvRequests)+message);

          MPI_Isend(((char*)_sendBuffer)+offset, count, MPI_CHAR, rr, tag,
                    comm, (haloSendRequests)+message);
          offset += count;
          ++message;
        }
      }
    }
  }
}

void mesh_t::HaloExchangeFinish(){

  if(totalHaloPairs>0){
    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Waitall(NhaloMessages, haloRecvRequests, recvStatus);
    MPI_Waitall(NhaloMessages, haloSendRequests, sendStatus);
  }
}

// start halo exchange (for o_q)
void mesh_t::HaloExchangeStart(size_t Nentries,   // number of dfloat entries
                               occa::memory& o_q,
                               occa::stream& defaultStream){

  if(totalHaloPairs>0){
    //required buffer size
    size_t Nbytes = Nentries*sizeof(dfloat)*totalHaloPairs;

    if (o_haloBuffer.size() < Nbytes) {
      if (o_haloBuffer.size()) {
        o_haloBuffer.free();
        o_sendBuffer.free();
        o_recvBuffer.free();
      }
      o_haloBuffer = device.malloc(Nbytes);
      sendBuffer = (dfloat*) occaHostMallocPinned(device, Nbytes, NULL, o_sendBuffer, h_sendBuffer);
      recvBuffer = (dfloat*) occaHostMallocPinned(device, Nbytes, NULL, o_recvBuffer, h_recvBuffer);
    }

    device.setStream(defaultStream);
    haloExtractKernel(totalHaloPairs, Nentries, o_haloElementList, o_q, o_haloBuffer);
    device.finish();
  }
}


void mesh_t::HaloExchangeInterm(size_t Nentries,   // number of dfloat entries
                               occa::stream& defaultStream,
                               occa::stream& dataStream){

  if(totalHaloPairs>0){
    // count outgoing and incoming meshes
    int tag = 999;

    //required buffer size
    size_t Nbytes = Nentries*sizeof(dfloat)*totalHaloPairs;

    device.setStream(dataStream);
    o_haloBuffer.copyTo(sendBuffer, Nbytes, 0, "async: true");
    device.finish();

    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(int rr=0;rr<size;++rr){
      if(rr!=rank){
        size_t count = NhaloPairs[rr]*Nentries;
        if(count){
          MPI_Irecv(recvBuffer+offset, count, MPI_DFLOAT, rr, tag,
                    comm, haloRecvRequests+message);

          MPI_Isend(sendBuffer+offset, count, MPI_DFLOAT, rr, tag,
                    comm, haloSendRequests+message);
          offset += count;
          ++message;
        }
      }
    }
    device.setStream(defaultStream);
  }
}

void mesh_t::HaloExchangeFinish(size_t Nentries,
                                occa::memory& o_q,
                                size_t offset){

  if(totalHaloPairs>0){
    MPI_Waitall(NhaloMessages, haloRecvRequests, recvStatus);
    MPI_Waitall(NhaloMessages, haloSendRequests, sendStatus);

    size_t Nbytes = Nentries*sizeof(dfloat)*totalHaloPairs;

    o_q.copyFrom(recvBuffer, Nbytes, offset*sizeof(dfloat));
  }
}
