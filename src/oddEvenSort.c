#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define mymin(a,b)  ((a<b)?(a):(b))

void oddEvenSort(void *base, size_t Nlocal, size_t sz, int (*compareFn)(const void *, const void *)){

  MPI_Status status;
  char *cbase = (char*) base;
  
  int rank, size, tag = 999;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // check to see if any process has a single entry
  int minNlocal;
  MPI_Allreduce(&Nlocal, &minNlocal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if(minNlocal==1){
    if(rank==0) printf("Exiting oddEvenSort: at least one process has less than two entries to sort\n");
    MPI_Finalize();
    exit(0);
  }
  
  // set up
  // figure out how many to receive from below (A) and how many to send above (C)
  size_t A=0, B=0, C=0;

  C = mymin(1,Nlocal/2);  // round up
  B = Nlocal-C;

  int it, Nit = size;
  
  // send C to process above
  if(rank+1<size)
    MPI_Send(&C, 1, MPI_INT, rank+1, tag, MPI_COMM_WORLD);
  
  // recv A from process below
  if(rank-1>=0)
    MPI_Recv(&A, 1, MPI_INT, rank-1, tag, MPI_COMM_WORLD, &status);
  
  // storage for combined buffer
  char *tmp = (char*) calloc(sz*(A+B+C), sizeof(char));
  memcpy(tmp+sz*A, cbase, (B+C)*sz); 

  for(it=0;it<Nit;++it){
    // S0: sort locally  [ B C ]
    qsort(tmp+sz*A, (B+C), sz, compareFn);
    
    // S1: send new C up
    if(rank+1<size)
      MPI_Send(tmp+sz*(A+B), sz*C, MPI_CHAR, rank+1, tag, MPI_COMM_WORLD);

    // recv new A from below
    if(rank-1>=0)
      MPI_Recv(tmp, sz*A, MPI_CHAR, rank-1, tag, MPI_COMM_WORLD, &status);
    
    // S2: sort A+B data (could use merge sort here)
    qsort(tmp, (A+B), sz, compareFn);
    
    // S3: send A data back down
    if(rank-1>=0)
      MPI_Send(tmp, sz*A, MPI_CHAR, rank-1, tag, MPI_COMM_WORLD);
    if(rank+1<size)
      MPI_Recv(tmp+sz*(A+B), sz*C, MPI_CHAR, rank+1, tag, MPI_COMM_WORLD, &status);
  }

  // copy back to original array
  memcpy(base, tmp+sz*A, sz*(B+C));
  free(tmp);
}

int compare(const void *a, const void *b){

  int ia = *((int*) a);
  int ib = *((int*) b);

  if(ia<ib) return -1;
  if(ia>ib) return +1;
  return 0;
}

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  int rank, size, tag = 999;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int Nglobal = atoi(argv[1]);
  int Nchunk  = Nglobal/size;
  int Nremainder = Nglobal-size*Nchunk;
  int Nlocal = Nchunk + (rank<Nremainder);

  
  int *data = (int*) calloc(Nlocal, sizeof(int));
  int *origData = (int*) calloc(Nlocal, sizeof(int));
  int n;
  srand48(rank);
  for(n=0;n<Nlocal;++n){
    data[n] = 1000*drand48();
    origData[n] = data[n];
  }

  oddEvenSort(data, Nlocal, sizeof(int), compare);

  char fileName[BUFSIZ];
  sprintf(fileName, "sortedData%02d.dat", rank);
  FILE *fp = fopen(fileName, "w");

  for(n=0;n<Nlocal;++n){
    fprintf(fp, " %d => %d \n", origData[n], data[n]);
  }

  fclose(fp);
	   
  
  MPI_Finalize();
}
