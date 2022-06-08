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
#include "parAdogs/parAdogsPartition.hpp"

namespace libp {

namespace paradogs {

/****************************************/
/* Multilevel Spectral Bipartition      */
/****************************************/
void graph_t::SpectralBipartition(const dfloat targetFraction[2]) {

  /*Create multilevel heirarchy*/
  MultigridSetup();

  /*Compute Fiedler vector */
  memory<dfloat>& Fiedler = FiedlerVector();

  /*Use Fiedler vector to bipartion graph*/
  const hlong K = std::ceil(targetFraction[0]*NVertsGlobal);
  const dfloat pivot = ParallelPivot(Nverts, Fiedler, K, comm);

  memory<int> partition(L[0].A.Ncols);

  for (dlong n=0;n<Nverts;++n) {
    if (Fiedler[n]<=pivot) {
      partition[n] = 0;
    } else {
      partition[n] = 1;
    }
  }

  /*Fill halo region of partition vector*/
  L[0].A.halo.Exchange(partition, 1);

  /*Split the graph according to this partitioning*/
  Split(partition);

  /*Clear the coarse levels*/
  MultigridDestroy();
}

} //namespace paradogs

} //namespace libp

