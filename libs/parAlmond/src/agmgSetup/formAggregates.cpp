/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {

static void formAggregatesDefault(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts);

static void formAggregatesLPSCN(parCSR *A, parCSR *C,
                         hlong* FineToCoarse,
                         hlong* globalAggStarts,
                         setupAide options);

void formAggregates(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts,
                     setupAide options)
{
  if (options.compareArgs("PARALMOND AGGREGATION STRATEGY", "DEFAULT")) {
    formAggregatesDefault(A, C, FineToCoarse, globalAggStarts);
  } else if (options.compareArgs("PARALMOND AGGREGATION STRATEGY", "LPSCN")) {
    formAggregatesLPSCN(A, C, FineToCoarse, globalAggStarts, options);
  } else {
    printf("WARNING:  Missing or bad value for option PARALMOND AGGREGATION STRATEGY.  Using default.\n");
    formAggregatesDefault(A, C, FineToCoarse, globalAggStarts);
  }
}

/*****************************************************************************/
// Default aggregation algorithm
//
// TODO:  Does this have a name?

static void formAggregatesDefault(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;
  const dlong offdNNZ = C->offd->nnz;

  dfloat *rands = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)    calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  hlong *globalRowStarts = A->globalRowStarts;

  for(dlong i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  // add the number of non-zeros in each column
  int *colCnt = (int *) calloc(M,sizeof(int));
  for(dlong i=0; i<diagNNZ; i++)
    colCnt[C->diag->cols[i]]++;

  for(dlong i=0; i<offdNNZ; i++)
    colCnt[C->offd->cols[i]]++;

  //gs for total column counts
  ogsGatherScatter(colCnt, ogsInt, ogsAdd, A->ogs);

  //add random pertubation
  for(int i=0;i<N;++i)
    rands[i] += colCnt[i];

  //gs to fill halo region
  ogsGatherScatter(rands, ogsDfloat, ogsAdd, A->ogs);

  hlong done = 0;
  while(!done){
    // first neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      hlong imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
          const dlong col = C->diag->cols[jj];
          if (col==i) continue;
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
          const dlong col = C->offd->cols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    for (dlong n=N;n<M;n++) {
      Tr[n] = 0.;
      Ts[n] = 0;
      Ti[n] = 0;
    }
    ogsGatherScatter(Tr, ogsDfloat, ogsAdd, A->ogs);
    ogsGatherScatter(Ts, ogsInt,    ogsAdd, A->ogs);
    ogsGatherScatter(Ti, ogsHlong,  ogsAdd, A->ogs);

    // second neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == (i + globalRowStarts[rank])))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    //share results
    for (dlong n=N;n<M;n++) states[n] = 0;
    ogsGatherScatter(states, ogsInt, ogsAdd, A->ogs);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = std::count(states, states+N, 0);
    MPI_Allreduce(&cnt,&done,1,MPI_HLONG, MPI_SUM,A->comm);
    done = (done == 0) ? 1 : 0;
  }

  dlong numAggs = 0;
  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));

  // count the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++) {
    if(states[i] == 1) {
      FineToCoarse[i] = globalAggStarts[rank] + numAggs++;
    } else {
      FineToCoarse[i] = -1;
    }
  }
  for(dlong i=N; i<M; i++) FineToCoarse[i] = 0;

  //share the initial aggregate flags
  ogsGatherScatter(FineToCoarse, ogsHlong, ogsAdd, A->ogs);

  // form the aggregates
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    hlong  imax = i + globalRowStarts[rank];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->colMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) {
    FineToCoarse[n] = 0;
    Tr[n] = 0.;
    Ts[n] = 0;
    Ti[n] = 0;
    Tc[n] = 0;
  }
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tr,     ogsDfloat, ogsAdd, A->ogs);
  ogsGatherScatter(Ts,     ogsInt,    ogsAdd, A->ogs);
  ogsGatherScatter(Ti,     ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tc,     ogsHlong,  ogsAdd, A->ogs);

  // second neighbours
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    //local entries
    for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
      const dlong col = C->diag->cols[jj];
      if (col==i) continue;
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
      const dlong col = C->offd->cols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);

  delete C;
}

/*****************************************************************************/
// Alan's "locally partial strong connected nodes" (LPSCN) algorithm.

typedef struct{
	dlong  index;
	dlong  Nnbs;
	dlong LNN;
} nbs_t;


int compareNBSmaxLPSCN(const void *a, const void *b)
{
	nbs_t *pa = (nbs_t *)a;
	nbs_t *pb = (nbs_t *)b;

	if (pa->Nnbs + pa->LNN < pb->Nnbs + pb->LNN)	return +1;
	if (pa->Nnbs + pa->LNN > pb->Nnbs + pb->LNN)	return -1;

	if (pa->index < pa->index )	return +1;
	if (pa->index > pa->index )	return -1;

	return 0;
}


int compareNBSminLPSCN(const void *a, const void *b)
{
	nbs_t *pa = (nbs_t *)a;
	nbs_t *pb = (nbs_t *)b;

	if (pa->Nnbs + pa->LNN < pb->Nnbs + pb->LNN)	return -1;
	if (pa->Nnbs + pa->LNN > pb->Nnbs + pb->LNN)	return +1;

	if (pa->index < pa->index )	return +1;
	if (pa->index > pa->index )	return -1;

	return 0;
}

static void formAggregatesLPSCN(parCSR *A, parCSR *C,
                         hlong* FineToCoarse,
                         hlong* globalAggStarts,
                         setupAide options)
{
  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;

  //  printf("\n-----------------------------------------\n");
  //printf("rank = %d (LQHSSN) =>  N = %d \t M = %d    ",rank,N,M);
  //printf("\n-----------------------------------------\n");

  dlong   *states = (dlong *)    calloc(M, sizeof(dlong));    //  M > N

  for(dlong i=0; i<N; i++)   // initialize states to -1
    states[i] = -1;

  hlong *globalRowStarts = A->globalRowStarts;

  // construct the local neigbors
  nbs_t *V = (nbs_t *) calloc(N,sizeof(nbs_t));

  for(dlong i=0; i<N; i++){
	V[i].index    = i;
	V[i].Nnbs     = C->diag->rowStarts[i+1] - C->diag->rowStarts[i];
	dlong dist     = C->diag->cols[C->diag->rowStarts[i+1]-1] - C->diag->cols[C->diag->rowStarts[i]];
	if (dist  > 0 )
	  V[i].LNN   =  V[i].Nnbs/dist;
	else
	  V[i].LNN = 0;
  }


  int MySort = 1;
  //options.getArgs("SORT",MySort);

  // sort the fake nbs
  if (options.compareArgs("PARALMOND LPSCN ORDERING", "MAX")) {
    if (MySort>0)		qsort(V,N,sizeof(nbs_t),compareNBSmaxLPSCN);
  } else if (options.compareArgs("PARALMOND LPSCN ORDERING", "MIN")) {
    if (MySort>0)		qsort(V,N,sizeof(nbs_t),compareNBSminLPSCN);
  } else {
    // TODO:  Calling exit() within a library is not a good way to handle errors...
    printf("ERROR:  Bad value for option PARALMOND LPSCN ORDERING.\n");
    exit(-1);
  }

  dlong R_num = 0; 	// number of non-isolated nodes to be used

  for (dlong i=0;i<N;i++)
	if(V[i].Nnbs>1)		R_num++;


  dlong R_nodes[N];   // N is number of local elements so at most N aggregates
  dlong R_pos[N];     //

  dlong k = 0;


  // get the sorted indices
  for (dlong i=0;i<N;i++){
	if(V[i].Nnbs>1){
		R_nodes[k] = V[i].index;	// R_nodes list of nodes to be used
		k++;
	}
  }

  dlong Agg_num = 0;

  // First aggregates
  //#pragma omp parallel for
  for(dlong i=0; i<R_num; i++){
	if (states[R_nodes[i]] == -1){
		int ok = 0;
		// verify that all NBS are free
		for(dlong j=C->diag->rowStarts[R_nodes[i]]; j<C->diag->rowStarts[R_nodes[i]+1];j++){
			if (states[C->diag->cols[j]]>-1){
				ok=1;
				break;
			}
		}
		// construct the aggregate
		if (ok == 0){
		        for(dlong j=C->diag->rowStarts[R_nodes[i]]; j<C->diag->rowStarts[R_nodes[i]+1];j++)
		                states[C->diag->cols[j]] = Agg_num;

			Agg_num++;
		}
	 }
  }

  R_num=0;   // reset the number of nodes to be used


  for (dlong i=0;i<N;i++){  // update list of  non-agreggate nodes
	if (states[V[i].index]==-1){      // cambie k por R_num
                R_nodes[R_num] = V[i].index;  // update the list of nodes
		R_num++;
	}
  }

  dlong *psudoAgg = (dlong *) calloc(M,sizeof(dlong));   // what is different bwtween dlong and hlong?

  // count the number of nodes at each aggregate
  for (dlong i=0;i<N;i++)
	if (states[V[i].index]>-1)		psudoAgg[states[V[i].index]]++;

  // #pragma omp parallel for
  for(dlong i=0; i<R_num; i++){
	if (states[R_nodes[i]] == -1){ 	 // sanity check
		if (C->diag->rowStarts[R_nodes[i]+1] - C->diag->rowStarts[R_nodes[i]] > 1){  	 // at most one neigbor
			dlong Agg_max;
			dlong posIndx = 0;
			dlong MoreAgg[Agg_num];
			for (dlong j=0;j<Agg_num;j++)
				MoreAgg[j]=0;
			// count the number of strong conections with each aggregate
			for(dlong j=C->diag->rowStarts[R_nodes[i]] ; j<C->diag->rowStarts[R_nodes[i]+1] ; j++){  // index 0 is itself
				if (states[C->diag->cols[j]] > -1){
					MoreAgg[states[ C->diag->cols[j]  ]]++;
				}
			}
			// look for the agregates with more strong connections & less nodes
			Agg_max = -1;
			for (dlong j=0;j<Agg_num;j++){
				if (Agg_max <= MoreAgg[j]){
					if (j == 0){
						Agg_max = MoreAgg[j];
						posIndx = j;
					}
					else if (Agg_max < MoreAgg[j]){
						Agg_max = MoreAgg[j];
						posIndx = j;
					}
					else if (psudoAgg[posIndx] > psudoAgg[j]){
						Agg_max = MoreAgg[j];
						posIndx = j;
					}
				}
			}
			states[R_nodes[i]] = posIndx;
			psudoAgg[posIndx]++;
		}
		else{  // no neighbors (isolated node)

		  // if (V[R_nodes[i]].offNnbs==0){
		    states[R_nodes[i]] = Agg_num;  // becomes a new aggregate
		    psudoAgg[Agg_num]++;
		    Agg_num++;
		    //
		    //}
		   // MULTIGRID becomes more slow if its added to an aggregate
		   /* dlong Min = psudoAgg[0];
		    for (dlong k = 1 ; k < Agg_num ; k++){
		      if (psudoAgg[k] < Min )     Min = psudoAgg[k];
		    }

		    states[R_nodes[i]] = Min;
		   */

		}
	}
  }

  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));
  //globalAggStarts = (hlong *) calloc(size+1,sizeof(hlong));

  // count the coarse nodes/aggregates in each rannk
  MPI_Allgather(&Agg_num,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++){
    //    if (states[i] >= 0)
      FineToCoarse[i] = globalAggStarts[rank] + states[i];
      //else
      // FineToCoarse[i] = -1;
  }


  // share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
	ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);


/*

// print info of aggregates

  dlong *AggNum = (dlong *) calloc(Agg_num,sizeof(dlong));

  for (dlong i = 0; i<M;i++)
	AggNum[FineToCoarse[i]]++;

  printf("\n rank %d elementos por Agregates ==========\n",rank);
  for (dlong i=0;i<Agg_num;i++)
  	printf("Agg %d  total = %d\n",i,AggNum[i]);
  printf("\n=======================================\n");
	*/

  free(states);
  free(V);
  free(psudoAgg);

  delete C;
}

} //namespace parAlmond
