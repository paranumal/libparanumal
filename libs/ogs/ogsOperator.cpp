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

namespace libp {

namespace ogs {

template<typename T>
struct Op_Add {
  inline const T init(){ return T{0}; }
  inline void operator()(T& gv, const T v) { gv += v; }
};
template<typename T>
struct Op_Mul {
  inline const T init(){ return T{1}; }
  inline void operator()(T& gv, const T v) { gv *= v; }
};
template<typename T>
struct Op_Max {
  inline const T init(){ return -std::numeric_limits<T>::max(); }
  inline void operator()(T& gv, const T v) { gv = (v>gv) ? v : gv; }
};
template<typename T>
struct Op_Min {
  inline const T init() {return  std::numeric_limits<T>::max(); }
  inline void operator()(T& gv, const T v) { gv = (v<gv) ? v : gv; }
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
  dlong *rowStarts, *colIds;
  if (trans==NoTrans) {
    Nrows = NrowsN;
    rowStarts = rowStartsN.ptr();
    colIds = colIdsN.ptr();
  } else {
    Nrows = NrowsT;
    rowStarts = rowStartsT.ptr();
    colIds = colIdsT.ptr();
  }

  const T* v_ptr  = v.ptr();
  T* gv_ptr = gv.ptr();

  #pragma omp parallel for
  for(dlong n=0;n<Nrows;++n){
    const dlong start = rowStarts[n];
    const dlong end   = rowStarts[n+1];

    for (int k=0;k<K;++k) {
      T val = Op<T>().init();
      for(dlong g=start;g<end;++g){
        Op<T>()(val, v_ptr[k+colIds[g]*K]);
      }
      gv_ptr[k+n*K] = val;
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
    if (NrowBlocksN)
      gatherKernel[type][op](NrowBlocksN,
                              k,
                              o_blockRowStartsN,
                              o_rowStartsN,
                              o_colIdsN,
                              o_v,
                              o_gv);
  } else {
    if (NrowBlocksT)
      gatherKernel[type][op](NrowBlocksT,
                              k,
                              o_blockRowStartsT,
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
  dlong *rowStarts, *colIds;
  if (trans==Trans) {
    Nrows = NrowsN;
    rowStarts = rowStartsN.ptr();
    colIds = colIdsN.ptr();
  } else {
    Nrows = NrowsT;
    rowStarts = rowStartsT.ptr();
    colIds = colIdsT.ptr();
  }

  T* v_ptr  = v.ptr();
  const T* gv_ptr = gv.ptr();

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
    if (NrowBlocksN)
      scatterKernel[type](NrowBlocksN,
                          k,
                          o_blockRowStartsN,
                          o_rowStartsN,
                          o_colIdsN,
                          o_gv,
                          o_v);
  } else {
    if (NrowBlocksT)
      scatterKernel[type](NrowBlocksT,
                          k,
                          o_blockRowStartsT,
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
  dlong *gRowStarts, *gColIds;
  dlong *sRowStarts, *sColIds;

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

  T* v_ptr = v.ptr();

  #pragma omp parallel for
  for(dlong n=0;n<Nrows;++n){
    const dlong gstart = gRowStarts[n];
    const dlong gend   = gRowStarts[n+1];
    const dlong sstart = sRowStarts[n];
    const dlong send   = sRowStarts[n+1];

    for (int k=0;k<K;++k) {
      T val = Op<T>().init();
      for(dlong g=gstart;g<gend;++g){
        Op<T>()(val, v_ptr[k+gColIds[g]*K]);
      }
      for(dlong s=sstart;s<send;++s){
        v_ptr[k+sColIds[s]*K] = val;
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
    if (NrowBlocksT)
      gatherScatterKernel[type][Add](NrowBlocksT,
                                     k,
                                     o_blockRowStartsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_rowStartsN,
                                     o_colIdsN,
                                     o_v);
  } else if (trans==Sym) {
    if (NrowBlocksT)
      gatherScatterKernel[type][Add](NrowBlocksT,
                                     k,
                                     o_blockRowStartsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_rowStartsT,
                                     o_colIdsT,
                                     o_v);
  } else {
    if (NrowBlocksT)
      gatherScatterKernel[type][Add](NrowBlocksT,
                                     k,
                                     o_blockRowStartsT,
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

void ogsOperator_t::setupRowBlocks() {

  dlong blockSumN=0, blockSumT=0;
  NrowBlocksN=0, NrowBlocksT=0;

  if (NrowsN) NrowBlocksN++;
  if (NrowsT) NrowBlocksT++;

  for (dlong i=0;i<NrowsT;i++) {
    const dlong rowSizeN  = rowStartsN[i+1]-rowStartsN[i];
    const dlong rowSizeT  = rowStartsT[i+1]-rowStartsT[i];

    //this row is pathalogically big. We can't currently run this
    LIBP_ABORT("Multiplicity of global node id: " << i
               << " in ogsOperator_t::setupRowBlocks is too large.",
               rowSizeN > gatherNodesPerBlock);
    LIBP_ABORT("Multiplicity of global node id: " << i
               << " in ogsOperator_t::setupRowBlocks is too large.",
               rowSizeT > gatherNodesPerBlock);

    if (blockSumN+rowSizeN > gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      NrowBlocksN++; //count the previous block
      blockSumN=rowSizeN; //start a new row block
    } else {
      blockSumN+=rowSizeN; //add this row to the block
    }

    if (blockSumT+rowSizeT > gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      NrowBlocksT++; //count the previous block
      blockSumT=rowSizeT; //start a new row block
    } else {
      blockSumT+=rowSizeT; //add this row to the block
    }
  }

  blockRowStartsN.calloc(NrowBlocksN+1);
  blockRowStartsT.calloc(NrowBlocksT+1);

  blockSumN=0, blockSumT=0;
  NrowBlocksN=0, NrowBlocksT=0;
  if (NrowsN) NrowBlocksN++;
  if (NrowsT) NrowBlocksT++;

  for (dlong i=0;i<NrowsT;i++) {
    const dlong rowSizeN  = rowStartsN[i+1]-rowStartsN[i];
    const dlong rowSizeT  = rowStartsT[i+1]-rowStartsT[i];

    if (blockSumN+rowSizeN > gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      blockRowStartsN[NrowBlocksN++] = i; //mark the previous block
      blockSumN=rowSizeN; //start a new row block
    } else {
      blockSumN+=rowSizeN; //add this row to the block
    }
    if (blockSumT+rowSizeT > gatherNodesPerBlock) { //adding this row will exceed the nnz per block
      blockRowStartsT[NrowBlocksT++] = i; //mark the previous block
      blockSumT=rowSizeT; //start a new row block
    } else {
      blockSumT+=rowSizeT; //add this row to the block
    }
  }
  blockRowStartsN[NrowBlocksN] = NrowsT;
  blockRowStartsT[NrowBlocksT] = NrowsT;

  o_blockRowStartsN = platform.malloc(blockRowStartsN);
  o_blockRowStartsT = platform.malloc(blockRowStartsT);
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

  blockRowStartsT.free();
  blockRowStartsN.free();
  o_blockRowStartsN.free();
  o_blockRowStartsT.free();

  nnzN=0;
  nnzT=0;
  NrowsN=0;
  NrowsT=0;
  Ncols=0;
  NrowBlocksN=0;
  NrowBlocksT=0;
}


template <template<typename> class U,
          template<typename> class V,
          typename T>
void extract(const dlong N,
             const int K,
             const memory<dlong> ids,
             const U<T> q,
             V<T> gatherq) {

  const T* q_ptr = q.ptr();
  T* gatherq_ptr = gatherq.ptr();

  for(dlong n=0;n<N;++n){
    const dlong gid = ids[n];

    for (int k=0;k<K;++k) {
      gatherq_ptr[k+n*K] = q_ptr[k+gid*K];
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
