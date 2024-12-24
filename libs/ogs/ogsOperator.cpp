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

#include <limits>
#include "ogs.hpp"
#include "ogs/ogsUtils.hpp"
#include "ogs/ogsOperator.hpp"
#include "primitives.hpp"

namespace libp {

namespace ogs {

template<typename T>
struct Op_Add {
  inline const T init() const { return T{0}; }
  inline void operator()(T& gv, const T v) const { gv += v; }
};
template<typename T>
struct Op_Mul {
  inline const T init() const { return T{1}; }
  inline void operator()(T& gv, const T v) const { gv *= v; }
};
template<typename T>
struct Op_Max {
  inline const T init() const { return -std::numeric_limits<T>::max(); }
  inline void operator()(T& gv, const T v) const { gv = (v>gv) ? v : gv; }
};
template<typename T>
struct Op_Min {
  inline const T init() const {return  std::numeric_limits<T>::max(); }
  inline void operator()(T& gv, const T v) const { gv = (v<gv) ? v : gv; }
};

/********************************
 * Gather Operation
 ********************************/
template <template<typename> class U,
          template<typename> class V,
          template<typename> class Op,
          typename T>
void ogsOperator_t::Gather(U<T> gv,
                           const V<T> v,
                           const int K,
                           const Transpose trans) {

  dlong Nrows;
  dlong *__restrict__ rowStarts, *__restrict__ colIds;
  if (trans==NoTrans) {
    Nrows = NrowsN;
    rowStarts = rowStartsN.ptr();
    colIds = colIdsN.ptr();
  } else {
    Nrows = NrowsT;
    rowStarts = rowStartsT.ptr();
    colIds = colIdsT.ptr();
  }

  const T*__restrict__ v_ptr  = v.ptr();
  T*__restrict__ gv_ptr = gv.ptr();

  const Op<T> op;

  if (K==1) {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong start = rowStarts[n];
      const dlong end   = rowStarts[n+1];

      T val = op.init();
      for(dlong g=start;g<end;++g){
        op(val, v_ptr[colIds[g]]);
      }
      gv_ptr[n] = val;
    }
  } else {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong start = rowStarts[n];
      const dlong end   = rowStarts[n+1];

      for (int k=0;k<K;++k) {
        T val = op.init();
        for(dlong g=start;g<end;++g){
          op(val, v_ptr[k+colIds[g]*K]);
        }
        gv_ptr[k+n*K] = val;
      }
    }
  }
}

template <template<typename> class U,
          template<typename> class V,
          typename T>
void ogsOperator_t::Gather(U<T> gv,
                           const V<T> v,
                           const int k,
                           const Op op,
                           const Transpose trans) {
  switch (op){
    case Add:
      Gather<U, V, Op_Add, T>(gv, v, k, trans); break;
    case Mul:
      Gather<U, V, Op_Mul, T>(gv, v, k, trans); break;
    case Max:
      Gather<U, V, Op_Max, T>(gv, v, k, trans); break;
    case Min:
      Gather<U, V, Op_Min, T>(gv, v, k, trans); break;
  }
}

template
void ogsOperator_t::Gather(memory<float> gv, const memory<float> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(memory<double> gv, const memory<double> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(memory<int> gv, const memory<int> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(memory<long long int> gv, const memory<long long int> v,
                           const int k, const Op op, const Transpose trans);

template
void ogsOperator_t::Gather(pinnedMemory<float> gv, const memory<float> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<double> gv, const memory<double> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<int> gv, const memory<int> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<long long int> gv, const memory<long long int> v,
                           const int k, const Op op, const Transpose trans);

template
void ogsOperator_t::Gather(pinnedMemory<float> gv, const pinnedMemory<float> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<double> gv, const pinnedMemory<double> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<int> gv, const pinnedMemory<int> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(pinnedMemory<long long int> gv, const pinnedMemory<long long int> v,
                           const int k, const Op op, const Transpose trans);


template<typename T>
void ogsOperator_t::Gather(deviceMemory<T> o_gv,
                           deviceMemory<T> o_v,
                           const int k,
                           const Op op,
                           const Transpose trans) {
  constexpr Type type = ogsType<T>::get();
  InitializeKernels(platform, type, op);

  if (trans==NoTrans) {
    if (gBlocking.NrowBlocksN)
      gatherKernel[type][op](gBlocking.NrowBlocksN,
                              k,
                              gBlocking.o_blockRowStartsN,
                              o_rowStartsN,
                              o_colIdsN,
                              o_v,
                              o_gv);
  } else {
    if (gBlocking.NrowBlocksT)
      gatherKernel[type][op](gBlocking.NrowBlocksT,
                              k,
                              gBlocking.o_blockRowStartsT,
                              o_rowStartsT,
                              o_colIdsT,
                              o_v,
                              o_gv);
  }
}

template
void ogsOperator_t::Gather(deviceMemory<float> gv, const deviceMemory<float> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(deviceMemory<double> gv, const deviceMemory<double> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(deviceMemory<int> gv, const deviceMemory<int> v,
                           const int k, const Op op, const Transpose trans);
template
void ogsOperator_t::Gather(deviceMemory<long long int> gv, const deviceMemory<long long int> v,
                           const int k, const Op op, const Transpose trans);


/********************************
 * Scatter Operation
 ********************************/
template <template<typename> class U,
          template<typename> class V,
          typename T>
void ogsOperator_t::Scatter(U<T> v, const V<T> gv,
                            const int K, const Transpose trans) {

  dlong Nrows;
  dlong *__restrict__ rowStarts, *__restrict__ colIds;
  if (trans==Trans) {
    Nrows = NrowsN;
    rowStarts = rowStartsN.ptr();
    colIds = colIdsN.ptr();
  } else {
    Nrows = NrowsT;
    rowStarts = rowStartsT.ptr();
    colIds = colIdsT.ptr();
  }

  T*__restrict__ v_ptr  = v.ptr();
  const T*__restrict__ gv_ptr = gv.ptr();

  if (K==1) {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong start = rowStarts[n];
      const dlong end   = rowStarts[n+1];

      for(dlong g=start;g<end;++g){
        v_ptr[colIds[g]] = gv_ptr[n];
      }
    }
  } else {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong start = rowStarts[n];
      const dlong end   = rowStarts[n+1];

      for(dlong g=start;g<end;++g){
        for (int k=0;k<K;++k) {
          v_ptr[k+colIds[g]*K] = gv_ptr[k+n*K];
        }
      }
    }
  }
}

template
void ogsOperator_t::Scatter(memory<float> v, const memory<float> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<double> v, const memory<double> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<int> v, const memory<int> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<long long int> v, const memory<long long int> gv,
                            const int K, const Transpose trans);

template
void ogsOperator_t::Scatter(memory<float> v, const pinnedMemory<float> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<double> v, const pinnedMemory<double> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<int> v, const pinnedMemory<int> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(memory<long long int> v, const pinnedMemory<long long int> gv,
                            const int K, const Transpose trans);

template<typename T>
void ogsOperator_t::Scatter(deviceMemory<T> o_v,
                            deviceMemory<T> o_gv,
                            const int k,
                            const Transpose trans) {
  constexpr Type type = ogsType<T>::get();
  InitializeKernels(platform, type, Add);

  if (trans==Trans) {
    if (sBlocking.NrowBlocksN)
      scatterKernel[type](sBlocking.NrowBlocksN,
                          k,
                          sBlocking.o_blockRowStartsN,
                          o_rowStartsN,
                          o_colIdsN,
                          o_gv,
                          o_v);
  } else {
    if (sBlocking.NrowBlocksT)
      scatterKernel[type](sBlocking.NrowBlocksT,
                          k,
                          sBlocking.o_blockRowStartsT,
                          o_rowStartsT,
                          o_colIdsT,
                          o_gv,
                          o_v);
  }
}

template
void ogsOperator_t::Scatter(deviceMemory<float> v, const deviceMemory<float> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(deviceMemory<double> v, const deviceMemory<double> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(deviceMemory<int> v, const deviceMemory<int> gv,
                            const int K, const Transpose trans);
template
void ogsOperator_t::Scatter(deviceMemory<long long int> v, const deviceMemory<long long int> gv,
                            const int K, const Transpose trans);

/********************************
 * GatherScatter Operation
 ********************************/
template <template<typename> class U,
          template<typename> class Op,
          typename T>
void ogsOperator_t::GatherScatter(U<T> v, const int K,
                                  const Transpose trans) {

  dlong Nrows;
  dlong *__restrict__ gRowStarts, *__restrict__ gColIds;
  dlong *__restrict__ sRowStarts, *__restrict__ sColIds;

  if (trans==Trans) {
    Nrows = NrowsN;
    gRowStarts = rowStartsT.ptr();
    gColIds    = colIdsT.ptr();
    sRowStarts = rowStartsN.ptr();
    sColIds    = colIdsN.ptr();
  } else if (trans==Sym) {
    Nrows = NrowsT;
    gRowStarts = rowStartsT.ptr();
    gColIds    = colIdsT.ptr();
    sRowStarts = rowStartsT.ptr();
    sColIds    = colIdsT.ptr();
  } else {
    Nrows = NrowsT;
    gRowStarts = rowStartsN.ptr();
    gColIds    = colIdsN.ptr();
    sRowStarts = rowStartsT.ptr();
    sColIds    = colIdsT.ptr();
  }

  T*__restrict__ v_ptr = v.ptr();

  const Op<T> op;

  if (K==1) {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong gstart = gRowStarts[n];
      const dlong gend   = gRowStarts[n+1];
      const dlong sstart = sRowStarts[n];
      const dlong send   = sRowStarts[n+1];

      T val = op.init();
      for(dlong g=gstart;g<gend;++g){
        op(val, v_ptr[gColIds[g]]);
      }
      for(dlong s=sstart;s<send;++s){
        v_ptr[sColIds[s]] = val;
      }
    }
  } else {
    #pragma omp parallel for
    for(dlong n=0;n<Nrows;++n){
      const dlong gstart = gRowStarts[n];
      const dlong gend   = gRowStarts[n+1];
      const dlong sstart = sRowStarts[n];
      const dlong send   = sRowStarts[n+1];

      for (int k=0;k<K;++k) {
        T val = op.init();
        for(dlong g=gstart;g<gend;++g){
          op(val, v_ptr[k+gColIds[g]*K]);
        }
        for(dlong s=sstart;s<send;++s){
          v_ptr[k+sColIds[s]*K] = val;
        }
      }
    }
  }
}

template <template<typename> class U,
          typename T>
void ogsOperator_t::GatherScatter(U<T> v,
                                  const int k,
                                  const Op op,
                                  const Transpose trans) {
  switch (op){
    case Add:
      GatherScatter<U, Op_Add, T>(v, k, trans); break;
    case Mul:
      GatherScatter<U, Op_Mul, T>(v, k, trans); break;
    case Max:
      GatherScatter<U, Op_Max, T>(v, k, trans); break;
    case Min:
      GatherScatter<U, Op_Min, T>(v, k, trans); break;
  }
}

template
void ogsOperator_t::GatherScatter(memory<float> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(memory<double> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(memory<int> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(memory<long long int> v,const int k,
                                  const Op op, const Transpose trans);

template<typename T>
void ogsOperator_t::GatherScatter(deviceMemory<T> o_v,
                                  const int k,
                                  const Op op,
                                  const Transpose trans) {
  constexpr Type type = ogsType<T>::get();
  InitializeKernels(platform, type, Add);

  if (trans==Trans) {
    if (gsBlocking.NrowBlocksT)
      gatherScatterKernel[type][Add](gsBlocking.NrowBlocksT,
                                     k,
                                     gsBlocking.o_blockRowStartsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_rowStartsN,
                                     o_colIdsN,
                                     o_v);
  } else if (trans==Sym) {
    if (gsBlocking.NrowBlocksT)
      gatherScatterKernel[type][Add](gsBlocking.NrowBlocksT,
                                     k,
                                     gsBlocking.o_blockRowStartsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_v);
  } else {
    if (gsBlocking.NrowBlocksT)
      gatherScatterKernel[type][Add](gsBlocking.NrowBlocksT,
                                     k,
                                     gsBlocking.o_blockRowStartsT,
                                     o_rowStartsN,
                                     o_colIdsN,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_v);
  }
}

template
void ogsOperator_t::GatherScatter(deviceMemory<float> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(deviceMemory<double> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(deviceMemory<int> v,const int k,
                                  const Op op, const Transpose trans);
template
void ogsOperator_t::GatherScatter(deviceMemory<long long int> v,const int k,
                                  const Op op, const Transpose trans);

/*
Binary search for the first entry between v[start] and v[end]
which is >= val. Returns end if no such index exist
*/
static dlong upperBound(dlong first,
                        dlong last,
                        const dlong *v,
                        const dlong val) {

  dlong count = last - first;

  while (count > 0) {
    const dlong step = count / 2;
    const dlong mid = first + step;

    if (v[mid] < val) {
      first = mid + 1;
      count -= step + 1;
    } else {
      count = step;
    }
  }

  return first;
}

static void blockRows(const dlong Nrows,
                      const int NodesPerBlock,
                      const memory<dlong> rowStarts,
                      dlong& Nblocks,
                      memory<dlong>& blockStarts) {

  if (!Nrows) return;

  //Check for a pathalogically big row. We can't currently run this
  memory<dlong> rowSizes(Nrows+1);
  prim::adjacentDifference(Nrows+1, rowStarts, rowSizes);
  dlong maxRowSize = prim::max(Nrows, rowSizes);
  rowSizes.free();

  LIBP_ABORT("Multiplicity of a global node in ogsOperator_t::setupRowBlocks is larger than requested blocking factor: " << NodesPerBlock,
             maxRowSize > NodesPerBlock);

  // We're going to resursively bisect the list of rows into blocks,
  //  so we need the scratch space to be some power of 2.
  //  Worst case is every block as only one row,
  //  so scratch space is at most Nrows blocks
  dlong maxNblocks = 1;
  while (maxNblocks < Nrows) { maxNblocks *= 2; }

  memory<dlong> blockStartsOld(maxNblocks+1);
  memory<dlong> blockStartsNew(maxNblocks+1);
  memory<dlong> blockSizes(maxNblocks);

  Nblocks = 1;
  blockStartsOld[0] = 0;
  blockStartsOld[1] = Nrows;
  blockStartsNew[0] = 0;

  blockSizes[0] = rowStarts[Nrows];
  dlong maxSize = blockSizes[0];

  while (maxSize > NodesPerBlock) {
    blockStartsNew[2*Nblocks] = Nrows;
    /*Recursively bisect the list of rows until the max block size is < NodesPerBlock*/
    #pragma omp parallel for
    for (dlong n=0;n<Nblocks;++n) {
      const dlong start = blockStartsOld[n];
      const dlong end   = blockStartsOld[n+1];
      const dlong rowBlockSize = rowStarts[end] - rowStarts[start];

      if (rowBlockSize > NodesPerBlock) {
        //Find the index ~middle of this block
        const dlong midSize = rowStarts[start] + (rowBlockSize + 1)/2;
        dlong mid = upperBound(start, end, rowStarts.ptr(), midSize);
        if (mid == end) --mid; // need at least one row in the right block
        blockStartsNew[2*n] = start;
        blockStartsNew[2*n+1] = mid;
        blockSizes[2*n] = rowStarts[mid] - rowStarts[start];
        blockSizes[2*n+1] = rowStarts[end] - rowStarts[mid];
      } else {
        blockStartsNew[2*n] = start;
        blockStartsNew[2*n+1] = end;
        blockSizes[2*n] = rowBlockSize;
        blockSizes[2*n+1] = 0;
      }
    }

    //swap blockStarts arrays
    blockStartsOld.swap(blockStartsNew);

    //Check if we're done bisecting
    Nblocks *= 2;
    maxSize = prim::max(Nblocks, blockSizes);
  }

  dlong Nunique=0;
  prim::unique(Nblocks+1, blockStartsOld, Nunique, blockStarts);
  Nblocks = Nunique-1;
}

//divide the list of colIds into roughly equal sized blocks so that each
// threadblock loads approximately an equal amount of data
void ogsOperator_t::createBlocking(const int NodesPerBlock,
                                   ogsOperator_t::rowBlocking_t& blocking) {
  blockRows(NrowsT,
            NodesPerBlock,
            rowStartsT,
            blocking.NrowBlocksT,
            blocking.blockRowStartsT);
  blocking.o_blockRowStartsT = platform.malloc(blocking.blockRowStartsT);

  if (kind==Signed) {
    blockRows(NrowsN,
              NodesPerBlock,
              rowStartsN,
              blocking.NrowBlocksN,
              blocking.blockRowStartsN);
    blocking.o_blockRowStartsN = platform.malloc(blocking.blockRowStartsN);
  } else {
    blocking.NrowBlocksN = blocking.NrowBlocksT;
    blocking.blockRowStartsN = blocking.blockRowStartsT;
    blocking.o_blockRowStartsN = blocking.o_blockRowStartsT;
  }
}

//Make gather operator using nodes list. List of non-zeros must be sorted by row index
ogsOperator_t::ogsOperator_t(platform_t &platform_,
                             Kind kind_,
                             const dlong NrowsN_,
                             const dlong NrowsT_,
                             const dlong Ncols_,
                             const dlong Nids,
                             memory<hlong> baseIds,
                             memory<dlong> rows,
                             memory<dlong> cols):
  platform(platform_),
  Ncols(Ncols_),
  NrowsN(NrowsN_),
  NrowsT(NrowsT_),
  kind(kind_)
{
  nnzT = Nids;
  rowStartsT.malloc(NrowsT+1);
  prim::runLengthEncodeConsecutive(nnzT, rows, NrowsT, rowStartsT);

  colIdsT = cols;

  o_rowStartsT = platform.malloc(rowStartsT);
  o_colIdsT = platform.malloc(colIdsT);

  if (kind == Signed) {
    memory<int> flags(Nids);

    #pragma omp parallel for
    for (dlong n=0;n<Nids;++n) {
      flags[n] = (baseIds[n]>0) ? 1 : 0;
    }

    nnzN = prim::count(Nids, flags, 1);
    memory<dlong> idsN(nnzN);
    prim::select(Nids, flags, 1, idsN);

    memory<dlong> rowsN(nnzN);
    prim::transformGather(nnzN, idsN, rows, rowsN);

    rowStartsN.malloc(NrowsN+1);
    prim::runLengthEncodeConsecutive(nnzN, rowsN, NrowsN, rowStartsN);

    colIdsN.malloc(nnzN);
    prim::transformGather(nnzN, idsN, cols, colIdsN);

    o_rowStartsN = platform.malloc(rowStartsN);
    o_colIdsN = platform.malloc(colIdsN);
  } else {
    nnzN = nnzT;
    rowStartsN = rowStartsT;
    colIdsN = colIdsT;
    o_rowStartsN = o_rowStartsT;
    o_colIdsN = o_colIdsT;
  }

  //divide the list of colIds into roughly equal sized blocks so that each
  // threadblock loads approximately an equal amount of data
  createBlocking(gNodesPerBlock, gBlocking);

  if (gNodesPerBlock==sNodesPerBlock) {
    sBlocking = gBlocking;
  } else {
    createBlocking(sNodesPerBlock, sBlocking);
  }

  if (gsNodesPerBlock==gNodesPerBlock) {
    gsBlocking = gBlocking;
  } else if (gsNodesPerBlock==sNodesPerBlock) {
    gsBlocking = sBlocking;
  } else {
    createBlocking(gsNodesPerBlock, gsBlocking);
  }
}


void ogsOperator_t::Free() {
  rowStartsT.free();
  colIdsT.free();
  rowStartsN.free();
  colIdsN.free();

  o_rowStartsT.free();
  o_colIdsT.free();
  o_rowStartsN.free();
  o_colIdsN.free();

  gBlocking.blockRowStartsT.free();
  gBlocking.blockRowStartsN.free();
  gBlocking.o_blockRowStartsN.free();
  gBlocking.o_blockRowStartsT.free();

  sBlocking.blockRowStartsT.free();
  sBlocking.blockRowStartsN.free();
  sBlocking.o_blockRowStartsN.free();
  sBlocking.o_blockRowStartsT.free();

  gsBlocking.blockRowStartsT.free();
  gsBlocking.blockRowStartsN.free();
  gsBlocking.o_blockRowStartsN.free();
  gsBlocking.o_blockRowStartsT.free();

  nnzN=0;
  nnzT=0;
  NrowsN=0;
  NrowsT=0;
  Ncols=0;

  gBlocking.NrowBlocksN=0;
  gBlocking.NrowBlocksT=0;
  sBlocking.NrowBlocksN=0;
  sBlocking.NrowBlocksT=0;
  gsBlocking.NrowBlocksN=0;
  gsBlocking.NrowBlocksT=0;
}


template <template<typename> class U,
          template<typename> class V,
          typename T>
void extract(const dlong N,
             const int K,
             const memory<dlong> ids,
             const U<T> q,
             V<T> gatherq) {

  const T*__restrict__ q_ptr = q.ptr();
  T*__restrict__ gatherq_ptr = gatherq.ptr();

  if (K==1) {
    for(dlong n=0;n<N;++n){
      const dlong gid = ids[n];
      gatherq_ptr[n] = q_ptr[gid];
    }
  } else {
    for(dlong n=0;n<N;++n){
      const dlong gid = ids[n];
      for (int k=0;k<K;++k) {
        gatherq_ptr[k+n*K] = q_ptr[k+gid*K];
      }
    }
  }
}

template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const memory<float> q, memory<float> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const memory<double> q, memory<double> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const memory<int> q, memory<int> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const memory<long long int> q, memory<long long int> gatherq);

template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const pinnedMemory<float> q, pinnedMemory<float> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const pinnedMemory<double> q, pinnedMemory<double> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const pinnedMemory<int> q, pinnedMemory<int> gatherq);
template void extract(const dlong N, const int K, const memory<dlong> ids,
                      const pinnedMemory<long long int> q, pinnedMemory<long long int> gatherq);

} //namespace ogs

} //namespace libp
