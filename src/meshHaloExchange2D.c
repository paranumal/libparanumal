#include <stdio.h>
#include "mpi.h"
#include "mesh2D.h"

// send data from partition boundary elements
// and receive data to ghost elements
void meshHaloExchange2D(mesh2D *mesh,
			size_t Nbytes,       // message size per element
			void *sourceBuffer,  
			void *sendBuffer,    // temporary buffer
			void *recvBuffer){

  if(mesh->totalHaloPairs>0){
    // MPI info
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // count outgoing and incoming meshes
    iint tag = 999;

    // copy data from outgoing elements into temporary send buffer
    for(iint i=0;i<mesh->totalHaloPairs;++i){
      // outgoing element
      iint e = mesh->haloElementList[i];
      memcpy(((char*)sendBuffer)+i*Nbytes, ((char*)sourceBuffer)+e*Nbytes, Nbytes);
    }

    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(iint r=0;r<size;++r){
      if(r!=rank){
	size_t count = mesh->NhaloPairs[r]*Nbytes;
	if(count){
	  MPI_Irecv(((char*)recvBuffer)+offset, count, MPI_CHAR, r, tag,
		    MPI_COMM_WORLD, (MPI_Request*)mesh->haloRecvRequests+message);
	
	  MPI_Isend(((char*)sendBuffer)+offset, count, MPI_CHAR, r, tag,
		    MPI_COMM_WORLD, (MPI_Request*)mesh->haloSendRequests+message);
	  offset += count;
	  ++message;
	}
      }
    }

    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Status *sendStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
    MPI_Status *recvStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
  
    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloRecvRequests, recvStatus);
    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloSendRequests, sendStatus);
  
    free(recvStatus);
    free(sendStatus);
  }
}      


// start halo exchange (for q)
void meshHaloExchangeStart2D(mesh2D *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // temporary buffer
			     void *recvBuffer){

  if(mesh->totalHaloPairs>0){
    // MPI info
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // count outgoing and incoming meshes
    iint tag = 999;
    
    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(iint r=0;r<size;++r){
      if(r!=rank){
	size_t count = mesh->NhaloPairs[r]*Nbytes;
	if(count){
	  MPI_Irecv(((char*)recvBuffer)+offset, count, MPI_CHAR, r, tag,
		    MPI_COMM_WORLD, (MPI_Request*)mesh->haloRecvRequests+message);
	
	  MPI_Isend(((char*)sendBuffer)+offset, count, MPI_CHAR, r, tag,
		    MPI_COMM_WORLD, (MPI_Request*)mesh->haloSendRequests+message);
	  offset += count;
	  ++message;
	}
      }
    }
  }  
}

void meshHaloExchangeFinish2D(mesh2D *mesh){

  if(mesh->totalHaloPairs>0){
    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Status *sendStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
    MPI_Status *recvStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
  
    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloRecvRequests, recvStatus);
    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloSendRequests, sendStatus);

    free(recvStatus);
    free(sendStatus);
  }
}      

