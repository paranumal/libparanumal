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
    occa::device device
    ...
    ogs_t *ogs = ogs_t::Setup(N, id, comm, verbose, device);

  defines a partition of the set of (processor, local index) pairs,
    (p,i) \in S_j  iff   abs(id[i]) == j  on processor p
  That is, all (p,i) pairs are grouped together (in group S_j) that have the
    same id (=j).
  S_0 is treated specially --- it is ignored completely
    (i.e., when id[i] == 0, local index i does not participate in any
    gather/scatter operation
  If id[i] on proc p is negative then the pair (p,i) is "flagged". This
  determines the non-symmetric behavior. For the simpler, symmetric case,
  all id's should be positive.

  When "ogs" is no longer needed, free it with

    ogsFree(ogs);

  A basic gatherScatter operation is, e.g.,

    occa::memory o_v;
    ...
    ogs->GatherScatter(o_v, ogs_double, ogs_add, ogs_sym);

  This gs call has the effect,

    o_v[i] <--  \sum_{ (p,j) \in S_{id[i]} } o_v_(p) [j]

  where o_v_(p) [j] means o_v[j] on proc p. In other words, every o_v[i] is replaced
  by the sum of all o_v[j]'s with the same id, given by id[i]. This accomplishes
  "direct stiffness summation" corresponding to the action of QQ^T, where
  "Q" is a boolean matrix that copies from a global vector (indexed by id)
  to the local vectors indexed by (p,i) pairs.

  Summation on doubles is not the only operation and datatype supported. Support
  includes the operations
    ogs_add, ogs_mul, ogs_max, ogs_min
  and datatypes
    ogs_dfloat, ogs_double, ogs_float, ogs_int, ogs_longlong, ogs_dlong, ogs_hlong.

  For the nonsymmetric behavior, the "transpose" parameter is important:

    ogs->GatherScatter(o_v, ogs_double, ogs_add, [ogs_notrans/ogs_trans/ogs_sym]);

  When transpose == ogs_notrans, any "flagged" (p,i) pairs (id[i] negative on p)
  do not participate in the sum, but *do* still receive the sum on output.
  As a special case, when only one (p,i) pair is unflagged per group this
  corresponds to the rectangular "Q" matrix referred to above.

  When transpose == ogs_trans, the "flagged" (p,i) pairs *do* participate in the sum,
  but do *not* get set on output. In the special case of only one unflagged
  (p,i) pair, this corresponds to the transpose of "Q" referred to above.

  When transpose == ogs_sym, all ids are considered "unflagged". That is,
  the "flagged" (p,i) pairs *do* participate in the sum, and *do* get set
  on output.

  An additional nonsymmetric operation is

    ogs->Gather(o_Gv, o_v, ogs_double, ogs_add, ogs_notrans);

  this has the effect of "assembling" the vector o_Gv. That is

    o_Gv[gid[j]] <--  \sum_{ (p,j) \in S_{id[i]} } o_v_(p) [j]

  for some ordering gid. As with the GatherScatter operation, when
  transpose == ogs_notrans, any "flagged" (p,i) pairs (id[i] negative on p)
  do not participate in the sum, whereas when transpose == ogs_trans the "flagged"
  (p,i) pairs *do* participate in the sum. Using transpose == ogs_sym is not
  supported (the symmetrized version of this operation is just GatherScatter).

  The reverse of this operation is

    ogs->Scatter(o_v, o_Gv, ogs_double, ogs_add, ogs_notrans);

  which has the effect of scattering in the assembled entries in o_Gv back to the
  orginal ordering. When transpose == ogs_notrans, "flagged" (p,i) pairs (id[i]
  negative on p) recieve their corresponding entry from o_Gv, and when
  transpose == ogs_trans the "flagged" (p,i) pairs do *not* recieve an entry.
  Using transpose == ogs_sym is not supported.

  A versions for vectors (contiguously packed) is, e.g.,

    occa::memory o_v;
    ogs->GatherScatterVec(o_v, k, ogs_double, ogs_add, ogs_sym);

  which is like "GatherScatter" operating on the datatype double[k],
  with summation here being vector summation. Number of messages sent
  is independent of k.

  For combining the communication for "GatherScatter" on multiple arrays:

    occa::memory o_v1, o_v2, ..., o_vk;

    ogs->GatherScatterMany(o_v, k, stride, ogs_double, op, trans);

  when the arrays o_v1, o_v2, ..., o_vk are packed in o_v as

    o_v1 = o_v + 0*stride, o_v2 = o_v + 1*stride, ...

  This call is equivalent to

    ogs->GatherScatter(o_v1, ogs_double, op, trans);
    ogs->GatherScatter(o_v2, ogs_double, op, trans);
    ...
    ogs->GatherScatter(o_vk, ogs_double, op, trans);

  except that all communication is done together.

  Asynchronous versions of the various GatherScatter functions are provided by

    ogs->GatherScatterStart(o_v, ogs_double, ogs_add, ogs_sym);
    ...
    ogs->GatherScatterFinish(o_v, ogs_double, ogs_add, ogs_sym);

  MPI communication is not initiated in GatherScatterStart, rather some initial
  message packing and host<->device transfers are queued. The user can then queue
  their own local work to the device which overlapps with this work before
  calling GatherScatterFinish to maximize the amount of communication hiding.



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

/* transpose switch */
typedef enum { ogs_sym, ogs_notrans, ogs_trans } ogs_transpose;

// OCCA+gslib gather scatter
class ogs_t {
public:
  MPI_Comm& comm;
  occa::device& device;

  dlong         N;
  dlong         Ngather;        //  total number of gather nodes

  dlong         Nlocal;         //  number of local nodes
  dlong         NlocalGather;   //  number of local gathered nodes
  dlong         NlocalScatter;  //  number of local scattered nodes

  dlong         Nhalo;          //  number of halo nodes
  dlong         NhaloGather;    //  number of gathered nodes on halo
  dlong         NhaloScatter;   //  number of scattered nodes on halo

  dlong         NlocalFused;
  dlong         NlocalFusedSym;

  dlong         *localGatherOffsets;
  dlong         *localScatterOffsets;
  dlong         *localGatherIds;
  dlong         *localScatterIds;
  occa::memory o_localGatherOffsets;
  occa::memory o_localScatterOffsets;
  occa::memory o_localGatherIds;
  occa::memory o_localScatterIds;

  //shortend offset and ids lists for fused kernels
  dlong         *localFusedGatherOffsets;
  dlong         *localFusedScatterOffsets;
  dlong         *localFusedOffsets;
  dlong         *localFusedGatherIds;
  dlong         *localFusedScatterIds;
  dlong         *localFusedIds;
  occa::memory o_localFusedGatherOffsets;
  occa::memory o_localFusedScatterOffsets;
  occa::memory o_localFusedOffsets;
  occa::memory o_localFusedGatherIds;
  occa::memory o_localFusedScatterIds;
  occa::memory o_localFusedIds;

  dlong         *haloGatherOffsets;
  dlong         *haloScatterOffsets;
  dlong         *haloGatherIds;
  dlong         *haloScatterIds;
  occa::memory o_haloGatherOffsets;
  occa::memory o_haloScatterOffsets;
  occa::memory o_haloGatherIds;
  occa::memory o_haloScatterIds;

  void *gsh;       // gslib handle
  void *gshSym;    // Symmetrized gslib handle (all ids made positive)

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

  static void Unique(hlong *ids, dlong _N, MPI_Comm _comm);

  // Host buffer versions
  void GatherScatter    (void  *v,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterVec (void  *v, const int k,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterMany(void  *v, const int k, const dlong stride,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void Gather    (void  *gv, void  *v,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherVec (void  *gv, void  *v, const int k,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherMany(void  *gv, void  *v, const int k,
                  const dlong gstride, const dlong stride,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void Scatter    (void  *v, void  *gv,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterVec (void  *v, void  *gv, const int k,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterMany(void  *v, void  *gv, const int k,
                   const dlong stride, const dlong gstride,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);

  // Synchronous device buffer versions
  void GatherScatter    (occa::memory&  o_v,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterVec (occa::memory&  o_v, const int k,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterMany(occa::memory&  o_v, const int k,
                         const dlong stride,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void Gather    (occa::memory&  o_gv, occa::memory&  o_v,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherVec (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherMany(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                  const dlong gstride, const dlong stride,
                  const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void Scatter    (occa::memory&  o_v, occa::memory&  o_gv,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterVec (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterMany(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                   const dlong stride, const dlong gstride,
                   const ogs_type type, const ogs_op op, const ogs_transpose trans);

  // Asynchronous device buffer versions
  void GatherScatterStart     (occa::memory&  o_v,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterFinish    (occa::memory&  o_v,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterVecStart  (occa::memory&  o_v, const int k,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterVecFinish (occa::memory&  o_v, const int k,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterManyStart (occa::memory&  o_v, const int k, const dlong stride,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherScatterManyFinish(occa::memory&  o_v, const int k, const dlong stride,
                               const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void GatherStart     (occa::memory&  o_gv, occa::memory&  o_v,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherFinish    (occa::memory&  o_gv, occa::memory&  o_v,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherVecStart  (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherVecFinish (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherManyStart (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const dlong gstride, const dlong stride,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void GatherManyFinish(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                        const dlong gstride, const dlong stride,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void ScatterStart     (occa::memory&  o_v, occa::memory&  o_gv,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterFinish    (occa::memory&  o_v, occa::memory&  o_gv,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterVecStart  (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterVecFinish (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterManyStart (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const dlong stride, const dlong gstride,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);
  void ScatterManyFinish(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                         const dlong stride, const dlong gstride,
                         const ogs_type type, const ogs_op op, const ogs_transpose trans);

  void reallocHostBuffer(size_t Nbytes);
  void reallocOccaBuffer(size_t Nbytes);
};

// OCCA halo exchange (thin wrapper of an ogs_t object)
class halo_t {
public:
  ogs_t* ogs;

  void Free() { if (ogs) { ogs->Free(); ogs=NULL; } }

  static halo_t *Setup(dlong N, hlong *ids, MPI_Comm &comm,
                       int verbose, occa::device& device) {
    halo_t *halo = new halo_t();
    halo->ogs = ogs_t::Setup(N, ids, comm, verbose, device);
    return halo;
  }

  // Synchronous Host buffer version
  void Exchange(void  *v, const int k, const ogs_type type) {
    ogs->GatherScatterVec(v, k, type, ogs_add, ogs_notrans);
  }

  // Synchronous device buffer version
  void Exchange(occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVec(o_v, k, type, ogs_add, ogs_notrans);
  }

  // Asynchronous device buffer version
  void ExchangeStart (occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVecStart(o_v, k, type, ogs_add, ogs_notrans);
  }
  void ExchangeFinish(occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVecFinish(o_v, k, type, ogs_add, ogs_notrans);
  }

  // Synchronous Host buffer version
  void Combine(void  *v, const int k, const ogs_type type) {
    ogs->GatherScatterVec(v, k, type, ogs_add, ogs_sym);
  }

  // Synchronous device buffer version
  void Combine(occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVec(o_v, k, type, ogs_add, ogs_sym);
  }

  // Asynchronous device buffer version
  void CombineStart (occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVecStart(o_v, k, type, ogs_add, ogs_sym);
  }
  void CombineFinish(occa::memory &o_v, const int k, const ogs_type type) {
    ogs->GatherScatterVecFinish(o_v, k, type, ogs_add, ogs_sym);
  }
};

#endif
