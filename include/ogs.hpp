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

/*
  OCCA Gather/Scatter Library

  The code

    MPI_Comm comm;
  	dlong N;
    hlong id[N];    // the hlong and dlong types are defined in "types.h"
    int verbose;
    ...
    ogs_t *ogs = ogs_t::Setup(N, id, comm, verbose);

  defines a partition of the set of (processor, local index) pairs,
    (p,i) \in S_j  iff   abs(id[i]) == j  on processor p
  That is, all (p,i) pairs are grouped together (in group S_j) that have the
    same id (=j).
  S_0 is treated specially --- it is ignored completely
    (i.e., when id[i] == 0, local index i does not participate in any
    gather/scatter operation

  When a partition of pairs, is shared between MPI processes, one specific
  member of the parition is chosen in an arbitrary but consistent way. This member
  is considered the 'owning' member and is used to define the nonsymmetric gather
  and scatter operations.

  When "ogs" is no longer needed, free it with

    ogsFree(ogs);

  A basic gatherScatter operation is, e.g.,

    occa::memory o_v;
    ...
    ogsGatherScatter(o_v, ogsDouble, ogsAdd, ogs);

  This gs call has the effect,

    o_v[i] <--  \sum_{ (p,j) \in S_{id[i]} } o_v_(p) [j]

  where o_v_(p) [j] means o_v[j] on proc p. In other words, every o_v[i] is replaced
  by the sum of all o_v[j]'s with the same id, given by id[i]. This accomplishes
  "direct stiffness summation" corresponding to the action of QQ^T, where
  "Q" is a boolean matrix that copies from a global vector (indexed by id)
  to the local vectors indexed by (p,i) pairs.

  Summation on doubles is not the only operation and datatype supported. Support
  includes the operations
    ogsAdd, ogsMul, ogsMax, ogsMin
  and datatypes
    ogsDfloat, ogsDouble, ogsFloat, ogsInt, ogsLong, ogsDlong, ogsHlong.

  For the nonsymmetric behavior, the operations are

    ogsGather (o_Gv, o_v, gsDouble, gsAdd, ogs);
    ogsScatter(o_Sv, o_v, gsDouble, gsAdd, ogs);

  A version for vectors (contiguously packed) is, e.g.,

    occa::memory o_v;
    ogsGatherScatterVec(o_v,k, ogsDouble,ogsAdd, transpose, ogs);

  which is like "gs" operating on the datatype double[k],
  with summation here being vector summation. Number of messages sent
  is independent of k.

  For combining the communication for "gs" on multiple arrays:

    occa::memory o_v1, o_v2, ..., o_vk;

    ogsGatherScatterMany(o_v, k, stride, ogsDouble, op, ogs);

  when the arrays o_v1, o_v2, ..., o_vk are packed in o_v as

    o_v1 = o_v + 0*stride, o_v2 = o_v + 1*stride, ...

  This call is equivalent to

    ogsGatherScatter(o_v1, gsDouble, op, ogs);
    ogsGatherScatter(o_v2, gsDouble, op, ogs);
    ...
    ogsGatherScatter(o_vk, gsDouble, op, ogs);

  except that all communication is done together.

  TODO: Comment on Start and Finish ordering

*/

#ifndef OGS_HPP
#define OGS_HPP

#include <math.h>
#include <stdlib.h>
#include <occa.hpp>
#include <mpi.h>

#include "types.h"
#include "utils.hpp"
#include "core.hpp"

//ogs defs
#include "../libs/ogs/include/ogsDefs.h"

/* type enum */
#define LIST OGS_FOR_EACH_TYPE(ITEM) ogs_type_n
#define ITEM(T) ogs_##T,
typedef enum { LIST } ogs_type;
#undef ITEM
#undef LIST

/* operation enum */
#define LIST OGS_FOR_EACH_OP(T,ITEM) ogs_op_n
#define ITEM(T,op) ogs_##op,
typedef enum { LIST } ogs_op;
#undef ITEM
#undef LIST

// OCCA+gslib gather scatter
class ogs_t {
public:
  MPI_Comm& comm;
  occa::device& device;

  dlong         N;
  dlong         Ngather;        //  total number of gather nodes
  dlong         Nlocal;         //  number of local nodes
  dlong         NlocalGather;   //  number of local gathered nodes
  dlong         Nhalo;          //  number of halo nodes
  dlong         NhaloGather;    //  number of gathered nodes on halo
  dlong         NownedHalo;     //  number of owned halo nodes

  dlong         *localGatherOffsets;
  dlong         *localGatherIds;
  occa::memory o_localGatherOffsets;
  occa::memory o_localGatherIds;

  dlong         *haloGatherOffsets;
  dlong         *haloGatherIds;
  occa::memory o_haloGatherOffsets;
  occa::memory o_haloGatherIds;

  void         *gshSym;       // gslib gather
  void         *gshNonSym;    // gslib gather

  void* hostBuf;
  size_t hostBufSize;

  void* haloBuf;
  occa::memory o_haloBuf;
  occa::memory h_haloBuf;

  ogs_t(MPI_Comm& _comm, occa::device& _device):
    comm(_comm), device(_device) {};

  void Free();

  static ogs_t *Setup(dlong N, hlong *ids, MPI_Comm &comm,
                      int verbose, occa::device& device);

  // Host buffer versions
  void GatherScatter    (void  *v,
                         const ogs_type type, const ogs_op op);
  void GatherScatterVec (void  *v, const int k,
                         const ogs_type type, const ogs_op op);
  void GatherScatterMany(void  *v, const int k, const dlong stride,
                         const ogs_type type, const ogs_op op);

  void Gather    (void  *gv, void  *v,
                  const ogs_type type, const ogs_op op);
  void GatherVec (void  *gv, void  *v, const int k,
                  const ogs_type type, const ogs_op op);
  void GatherMany(void  *gv, void  *v, const int k,
                  const dlong gstride, const dlong stride,
                  const ogs_type type, const ogs_op op);

  void Scatter    (void  *v, void  *gv,
                   const ogs_type type, const ogs_op op);
  void ScatterVec (void  *v, void  *gv, const int k,
                   const ogs_type type, const ogs_op op);
  void ScatterMany(void  *v, void  *gv, const int k,
                   const dlong stride, const dlong gstride,
                   const ogs_type type, const ogs_op op);

  // Synchronous device buffer versions
  void GatherScatter    (occa::memory&  o_v,
                         const ogs_type type, const ogs_op op);
  void GatherScatterVec (occa::memory&  o_v, const int k,
                         const ogs_type type, const ogs_op op);
  void GatherScatterMany(occa::memory&  o_v, const int k,
                         const dlong stride,
                         const ogs_type type, const ogs_op op);

  void Gather    (occa::memory&  o_gv, occa::memory&  o_v,
                  const ogs_type type, const ogs_op op);
  void GatherVec (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                  const ogs_type type, const ogs_op op);
  void GatherMany(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                  const dlong gstride, const dlong stride,
                  const ogs_type type, const ogs_op op);

  void Scatter    (occa::memory&  o_v, occa::memory&  o_gv,
                   const ogs_type type, const ogs_op op);
  void ScatterVec (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                   const ogs_type type, const ogs_op op);
  void ScatterMany(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                   const dlong stride, const dlong gstride,
                   const ogs_type type, const ogs_op op);

  // Asynchronous device buffer versions
  void GatherScatterStart     (occa::memory&  o_v,
                               const ogs_type type, const ogs_op op);
  void GatherScatterFinish    (occa::memory&  o_v,
                               const ogs_type type, const ogs_op op);
  void GatherScatterVecStart  (occa::memory&  o_v, const int k,
                               const ogs_type type, const ogs_op op);
  void GatherScatterVecFinish (occa::memory&  o_v, const int k,
                               const ogs_type type, const ogs_op op);
  void GatherScatterManyStart (occa::memory&  o_v, const int k, const dlong stride,
                               const ogs_type type, const ogs_op op);
  void GatherScatterManyFinish(occa::memory&  o_v, const int k, const dlong stride,
                               const ogs_type type, const ogs_op op);

  void GatherStart     (occa::memory&  o_gv, occa::memory&  o_v,
                        const ogs_type type, const ogs_op op);
  void GatherFinish    (occa::memory&  o_gv, occa::memory&  o_v,
                        const ogs_type type, const ogs_op op);
  void GatherVecStart  (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const ogs_type type, const ogs_op op);
  void GatherVecFinish (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const ogs_type type, const ogs_op op);
  void GatherManyStart (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const dlong gstride, const dlong stride,
                        const ogs_type type, const ogs_op op);
  void GatherManyFinish(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const dlong gstride, const dlong stride,
                        const ogs_type type, const ogs_op op);

  void ScatterStart     (occa::memory&  o_v, occa::memory&  o_gv,
                         const ogs_type type, const ogs_op op);
  void ScatterFinish    (occa::memory&  o_v, occa::memory&  o_gv,
                         const ogs_type type, const ogs_op op);
  void ScatterVecStart  (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const ogs_type type, const ogs_op op);
  void ScatterVecFinish (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const ogs_type type, const ogs_op op);
  void ScatterManyStart (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const dlong stride, const dlong gstride,
                         const ogs_type type, const ogs_op op);
  void ScatterManyFinish(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const dlong stride, const dlong gstride,
                         const ogs_type type, const ogs_op op);

  void reallocHostBuffer(size_t Nbytes);
  void reallocOccaBuffer(size_t Nbytes);
};

// OCCA halo exchange
class halo_t {
public:
  MPI_Comm& comm;
  occa::device& device;

  dlong         N;
  dlong         Nlocal;         //  number of local nodes
  dlong         Nsend;          //  number of outgong nodes
  dlong         Nhalo;          //  number of halo nodes recieved

  dlong         *haloSendIds;
  occa::memory o_haloSendIds;

  void         *gshHalo;       // gslib gather

  void* hostBuf;
  size_t hostBufSize;

  void* haloBuf;
  occa::memory o_haloBuf;
  occa::memory h_haloBuf;

  halo_t(MPI_Comm& _comm, occa::device& _device):
    comm(_comm), device(_device) {};

  void Free();

  static halo_t *Setup(dlong N, hlong *ids, MPI_Comm &comm,
                       int verbose, occa::device& device);

  // Synchronous Host buffer version
  void Exchange(void  *v, const int k, const ogs_type type);

  // Asynchronous Host buffer version
  void ExchangeStart (void* v, const int k, const ogs_type type);
  void ExchangeFinish(void* v, const int k, const ogs_type type);

  // Synchronous device buffer version
  void Exchange      (occa::memory &o_v, const int k, const ogs_type type);

  // Asynchronous device buffer version
  void ExchangeStart (occa::memory &o_v, const int k, const ogs_type type);
  void ExchangeFinish(occa::memory &o_v, const int k, const ogs_type type);

  void reallocHostBuffer(size_t Nbytes);
  void reallocOccaBuffer(size_t Nbytes);
};

#endif
