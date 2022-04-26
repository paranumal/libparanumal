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
#include "parAdogs/parAdogsMatrix.hpp"
#include "parAdogs/parAdogsPartition.hpp"
#include <random>

namespace libp {

namespace paradogs {

extern std::mt19937 RNG;

/*Create a vertex matching using distance-2 aggregation*/
void parCSR::Aggregate(dlong& Nc,
                       const dfloat theta,
                       memory<hlong>& FineToCoarse) {

  /*Create rng*/
  std::uniform_real_distribution<> distrib(-0.25, 0.25);

  parCSR strong(Nrows, Ncols, platform, comm);
  strong.diag.rowStarts.malloc(Nrows+1);

  #pragma omp parallel for
  for(dlong i=0; i<Nrows+1; i++) {
    strong.diag.rowStarts[i]=0;
  }

  #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){
    int strong_per_row = 0;

    const dfloat Aii = diagA[i];

    //local entries
    dlong Jstart = diag.rowStarts[i];
    dlong Jend   = diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = diag.cols[jj];
      if (col==i) continue;

      const dfloat Ajj = std::abs(diagA[col]);

      if(std::abs(diag.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        strong_per_row++;
    }
    //non-local entries
    Jstart = offd.rowStarts[i];
    Jend   = offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = offd.cols[jj];
      const dfloat Ajj = std::abs(diagA[col]);

      if(std::abs(offd.vals[jj]) > theta*(sqrt(Aii*Ajj)))
        strong_per_row++;
    }

    strong.diag.rowStarts[i+1] = strong_per_row;
  }

  // cumulative sum
  for(dlong i=1; i<Nrows+1 ; i++) {
    strong.diag.rowStarts[i] += strong.diag.rowStarts[i-1];
  }
  strong.diag.nnz = strong.diag.rowStarts[Nrows];
  strong.diag.cols.malloc(strong.diag.nnz);
  strong.diag.vals.malloc(strong.diag.nnz);

  // fill in the columns for strong connections
  // #pragma omp parallel for
  for(dlong i=0; i<Nrows; i++){
    const dfloat Aii = diagA[i];

    dlong counter = strong.diag.rowStarts[i];

    //local entries
    dlong Jstart = diag.rowStarts[i];
    dlong Jend   = diag.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = diag.cols[jj];
      if (col==i) continue;

      const dfloat Ajj = std::abs(diagA[col]);

      if(std::abs(diag.vals[jj]) > theta*(sqrt(Aii*Ajj))) {
        strong.diag.cols[counter] = col;
        strong.diag.vals[counter++] = std::abs(diag.vals[jj]) + distrib(paradogs::RNG);
      }
    }
    //non-local entries
    Jstart = offd.rowStarts[i];
    Jend   = offd.rowStarts[i+1];
    for(dlong jj= Jstart; jj<Jend; jj++){
      const dlong col = offd.cols[jj];

      const dfloat Ajj = std::abs(diagA[col]);

      if(std::abs(offd.vals[jj]) > theta*(sqrt(Aii*Ajj))) {
        strong.diag.cols[counter] = col;
        strong.diag.vals[counter++] = std::abs(offd.vals[jj]) + distrib(paradogs::RNG);
      }
    }
  }

  memory<float> rand(Ncols);
  memory<int>   Ts(Ncols);
  memory<float> Tr(Ncols);
  memory<hlong> Tn(Ncols);

  /*Initialize state array*/
  /*  0 - Undecided */
  /* -1 - Not MIS */
  /*  1 - MIS */
  memory<int>   state(Ncols, 0);

  /*Use vertex degree with random noise to break ties*/
  // #pragma omp parallel for
  for (dlong n=0;n<Nrows;++n) {
    rand[n] = strong.diag.rowStarts[n+1]
              - strong.diag.rowStarts[n]
              + distrib(paradogs::RNG);
  }

  //fill halo region
  halo.Exchange(rand, 1);

  do {
    // first neighbours
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++){
      int    smax = state[n];

      if (smax==1) continue;

      float  rmax = rand[n];
      hlong  nmax = colMap[n];

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k  = strong.diag.cols[j];
        const int   sk = state[k];
        const float rk = rand[k];
        const hlong nk = colMap[k];
        if ((sk>smax)              || /*If neighbor is MIS node*/
           ((sk==smax)&&(rk>rmax)) || /*Else if it has a bigger weight*/
           ((sk==smax)&&(rk==rmax)&&(nk>nmax))) { /*Rare, but just in case, break tie with index number*/
          smax = sk;
          rmax = rk;
          nmax = nk;
        }
      }
      Ts[n] = smax;
      Tr[n] = rmax;
      Tn[n] = nmax;
    }

    //share results
    halo.Exchange(Ts, 1);
    halo.Exchange(Tr, 1);
    halo.Exchange(Tn, 1);

    // second neighbours
    #pragma omp parallel for
    for(dlong n=0; n<Nrows; n++){
      if (state[n]!=0) continue;

      int   smax = Ts[n];
      float rmax = Tr[n];
      hlong nmax = Tn[n];

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k = strong.diag.cols[j];
        const int   sk = Ts[k];
        const float rk = Tr[k];
        const dlong nk = Tn[k];
        if ((sk>smax)              || /*If neighbor is MIS node*/
           ((sk==smax)&&(rk>rmax)) || /*Else if it has a bigger weight*/
           ((sk==smax)&&(rk==rmax)&&(nk>nmax))) { /*Rare, but just in case, break tie with index number*/
          smax = sk;
          rmax = rk;
          nmax = nk;
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if(nmax == colMap[n]) state[n] = 1;

      // if there is an MIS node within distance 2, I am removed
      if(smax>0) state[n] = -1;
    }

    //share results
    halo.Exchange(state, 1);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = 0;
    for (dlong n=0;n<Nrows;n++) if (state[n]==0) cnt++;
    comm.Allreduce(cnt);

    if (cnt==0) break;

  } while(true);

  rand.free();
  Tr.free();
  Tn.free();


  // count the coarse nodes/aggregates
  Nc=0;
  for(dlong i=0; i<Nrows; i++)
    if(state[i] == 1) Nc++;

  /*Get global offsets*/
  hlong localNc=static_cast<hlong>(Nc);
  hlong NcOffsetL=0, NcOffsetU=0;
  comm.Scan(localNc, NcOffsetU);
  NcOffsetL = NcOffsetU-Nc;

  /*Initialize Matching array*/
  Nc=0;
  for(dlong i=0; i<Nrows; i++) {
    if(state[i] == 1) {
      Ts[i] = 1;
      FineToCoarse[i] = NcOffsetL+Nc++;
    } else {
      Ts[i] = -1;
      FineToCoarse[i] = -1;
    }
  }

  //share the initial aggregate flags
  halo.Exchange(Ts, 1);
  halo.Exchange(FineToCoarse, 1);

  // first neighbours
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++){
    if (FineToCoarse[n]==-1) {
      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k  = strong.diag.cols[j];
        const int   sk = FineToCoarse[k];

        /*If this node is an MIS node, join the aggregate*/
        if (state[k]==1) {
          FineToCoarse[n] = sk;
          Ts[n] = 1;
          break;
        }
      }
    }
  }

  halo.Exchange(Ts, 1);
  halo.Exchange(FineToCoarse, 1);

  // second neighbours
  #pragma omp parallel for
  for(dlong n=0; n<Nrows; n++){
    if (FineToCoarse[n]==-1) { //If we're still undecided
      hlong cmax = -1;
      float rmax = -1.0;
      hlong kmax = -1;

      for(dlong j=strong.diag.rowStarts[n];j<strong.diag.rowStarts[n+1];j++){
        const dlong k = strong.diag.cols[j];
        const int   sk = Ts[k];
        const hlong nk = colMap[k];
        if (sk!=-1) { /*If the neighbor is in the neighborhood of an MIS node*/
          // const float rk = rand[k];
          const float rk = strong.diag.vals[j];
          if( (rk>rmax)            || /*If edge is strongest*/
             ((rk==rmax)&&(nk>kmax))) { /*Rare, but just in case, break tie with index number*/
            cmax = FineToCoarse[k];
            rmax = rk;
            kmax = nk;
          }
        }
      }
      FineToCoarse[n] = cmax;
    }
  }

  //share results
  halo.Exchange(FineToCoarse, 1);
}

} //namespace paradogs

} //namespace libp
