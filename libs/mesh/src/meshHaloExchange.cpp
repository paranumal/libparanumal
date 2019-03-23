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
                          void *sendBuffer,    // temporary buffer
                          void *recvBuffer){

  // count outgoing and incoming meshes
  int tag = 999;

  // copy data from outgoing elements into temporary send buffer
  for(int i=0;i<totalHaloPairs;++i){
    // outgoing element
    int e = haloElementList[i];
    // copy element e data to sendBuffer
    memcpy(((char*)sendBuffer)+i*Nbytes, ((char*)sourceBuffer)+e*Nbytes, Nbytes);
  }

  // initiate immediate send  and receives to each other process as needed
  int offset = 0, message = 0;
  for(int r=0;r<size;++r){
    if(r!=rank){
      size_t count = NhaloPairs[r]*Nbytes;
      if(count){
        //      printf("rank %d sending %d bytes to rank %d\n", rank, count, r);
        MPI_Irecv(((char*)recvBuffer)+offset, count, MPI_CHAR, r, tag,
                  comm, (MPI_Request*)haloRecvRequests+message);

        MPI_Isend(((char*)sendBuffer)+offset, count, MPI_CHAR, r, tag,
                  comm, (MPI_Request*)haloSendRequests+message);
        offset += count;
        ++message;
      }
    }
  }
  //  printf("NhaloMessages = %d\n", NhaloMessages);

  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Status *sendStatus = (MPI_Status*) calloc(NhaloMessages, sizeof(MPI_Status));
  MPI_Status *recvStatus = (MPI_Status*) calloc(NhaloMessages, sizeof(MPI_Status));

  MPI_Waitall(NhaloMessages, (MPI_Request*)haloRecvRequests, recvStatus);
  MPI_Waitall(NhaloMessages, (MPI_Request*)haloSendRequests, sendStatus);

  free(recvStatus);
  free(sendStatus);
}


// start halo exchange (for q)
void mesh_t::HaloExchangeStart(size_t Nbytes,       // message size per element
                             void *sendBuffer,    // temporary buffer
                             void *recvBuffer){

  if(totalHaloPairs>0){
    // count outgoing and incoming meshes
    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(int r=0;r<size;++r){
      if(r!=rank){
        size_t count = NhaloPairs[r]*Nbytes;
        if(count){
          MPI_Irecv(((char*)recvBuffer)+offset, count, MPI_CHAR, r, tag,
                    comm, ((MPI_Request*)haloRecvRequests)+message);

          MPI_Isend(((char*)sendBuffer)+offset, count, MPI_CHAR, r, tag,
                    comm, ((MPI_Request*)haloSendRequests)+message);
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
    MPI_Status *sendStatus = (MPI_Status*) calloc(NhaloMessages, sizeof(MPI_Status));
    MPI_Status *recvStatus = (MPI_Status*) calloc(NhaloMessages, sizeof(MPI_Status));

    MPI_Waitall(NhaloMessages, (MPI_Request*)haloRecvRequests, recvStatus);
    MPI_Waitall(NhaloMessages, (MPI_Request*)haloSendRequests, sendStatus);

    free(recvStatus);
    free(sendStatus);
  }
}
