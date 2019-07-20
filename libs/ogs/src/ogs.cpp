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

#include "ogs.hpp"
#include "ogsKernels.hpp"

// Host buffer versions
void ogs_t::GatherScatter    (void  *v,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGatherScatter(v, 1, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterVec (void  *v, const int k,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGatherScatter(v, k, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterMany(void  *v, const int k, const dlong stride,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGatherScatter(v, 1, k, stride, type, op, trans, *this); }

void ogs_t::Gather    (void  *gv, void  *v,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGather(gv, v, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherVec (void  *gv, void  *v, const int k,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGather(gv, v, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherMany(void  *gv, void  *v, const int k,
                       const dlong gstride, const dlong stride,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostGather(gv, v, 1, k, gstride, stride, type, op, trans, *this); }

void ogs_t::Scatter    (void  *v, void  *gv,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostScatter(v, gv, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterVec (void  *v, void  *gv, const int k,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostScatter(v, gv, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterMany(void  *v, void  *gv, const int k,
                        const dlong stride, const dlong gstride,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::hostScatter(v, gv, 1, k, stride, gstride, type, op, trans, *this); }


// Synchronous device buffer versions
void ogs_t::GatherScatter    (occa::memory&  o_v,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherScatterStart (o_v, 1, 1, 0, type, op, trans, *this);
  ogs::occaGatherScatterFinish(o_v, 1, 1, 0, type, op, trans, *this);
}

void ogs_t::GatherScatterVec (occa::memory&  o_v, const int k,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherScatterStart (o_v, k, 1, 0, type, op, trans, *this);
  ogs::occaGatherScatterFinish(o_v, k, 1, 0, type, op, trans, *this);
}

void ogs_t::GatherScatterMany(occa::memory&  o_v, const int k,
                              const dlong stride,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherScatterStart (o_v, 1, k, stride, type, op, trans, *this);
  ogs::occaGatherScatterFinish(o_v, 1, k, stride, type, op, trans, *this);
}

void ogs_t::Gather    (occa::memory&  o_gv, occa::memory&  o_v,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherStart (o_gv, o_v, 1, 1, 0, 0, type, op, trans, *this);
  ogs::occaGatherFinish(o_gv, o_v, 1, 1, 0, 0, type, op, trans, *this);
}

void ogs_t::GatherVec (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherStart (o_gv, o_v, k, 1, 0, 0, type, op, trans, *this);
  ogs::occaGatherFinish(o_gv, o_v, k, 1, 0, 0, type, op, trans, *this);
}

void ogs_t::GatherMany(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                       const dlong gstride, const dlong stride,
                       const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaGatherStart (o_gv, o_v, 1, k, gstride, stride, type, op, trans, *this);
  ogs::occaGatherFinish(o_gv, o_v, 1, k, gstride, stride, type, op, trans, *this);
}

void ogs_t::Scatter    (occa::memory&  o_v, occa::memory&  o_gv,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaScatterStart (o_v, o_gv, 1, 1, 0, 0, type, op, trans, *this);
  ogs::occaScatterFinish(o_v, o_gv, 1, 1, 0, 0, type, op, trans, *this);
}

void ogs_t::ScatterVec (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaScatterStart (o_v, o_gv, k, 1, 0, 0, type, op, trans, *this);
  ogs::occaScatterFinish(o_v, o_gv, k, 1, 0, 0, type, op, trans, *this);
}

void ogs_t::ScatterMany(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                        const dlong stride, const dlong gstride,
                        const ogs_type type, const ogs_op op, const ogs_transpose trans) {
  ogs::occaScatterStart (o_v, o_gv, 1, k, stride, gstride, type, op, trans, *this);
  ogs::occaScatterFinish(o_v, o_gv, 1, k, stride, gstride, type, op, trans, *this);
}

// Asynchronous device buffer versions
void ogs_t::GatherScatterStart     (occa::memory&  o_v,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterStart (o_v, 1, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterFinish    (occa::memory&  o_v,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterFinish(o_v, 1, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterVecStart  (occa::memory&  o_v, const int k,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterStart (o_v, k, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterVecFinish (occa::memory&  o_v, const int k,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterFinish(o_v, k, 1, 0, type, op, trans, *this); }

void ogs_t::GatherScatterManyStart (occa::memory&  o_v, const int k, const dlong stride,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterStart (o_v, 1, k, stride, type, op, trans, *this); }

void ogs_t::GatherScatterManyFinish(occa::memory&  o_v, const int k, const dlong stride,
                                    const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherScatterFinish(o_v, 1, k, stride, type, op, trans, *this); }

void ogs_t::GatherStart     (occa::memory&  o_gv, occa::memory&  o_v,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherStart (o_gv, o_v, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherFinish    (occa::memory&  o_gv, occa::memory&  o_v,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherFinish(o_gv, o_v, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherVecStart  (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherStart (o_gv, o_v, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherVecFinish (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherFinish(o_gv, o_v, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::GatherManyStart (occa::memory&  o_gv, occa::memory&  o_v, const int k,
                             const dlong gstride, const dlong stride,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherStart (o_gv, o_v, 1, k, gstride, stride, type, op, trans, *this); }

void ogs_t::GatherManyFinish(occa::memory&  o_gv, occa::memory&  o_v, const int k,
                             const dlong gstride, const dlong stride,
                             const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaGatherFinish(o_gv, o_v, 1, k, gstride, stride, type, op, trans, *this); }

void ogs_t::ScatterStart     (occa::memory&  o_v, occa::memory&  o_gv,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterStart (o_v, o_gv, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterFinish    (occa::memory&  o_v, occa::memory&  o_gv,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterFinish(o_v, o_gv, 1, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterVecStart  (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterStart (o_v, o_gv, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterVecFinish (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterFinish(o_v, o_gv, k, 1, 0, 0, type, op, trans, *this); }

void ogs_t::ScatterManyStart (occa::memory&  o_v, occa::memory&  o_gv, const int k,
                              const dlong stride, const dlong gstride,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterStart (o_v, o_gv, 1, k, stride, gstride, type, op, trans, *this); }

void ogs_t::ScatterManyFinish(occa::memory&  o_v, occa::memory&  o_gv, const int k,
                              const dlong stride, const dlong gstride,
                              const ogs_type type, const ogs_op op, const ogs_transpose trans)
{ ogs::occaScatterFinish(o_v, o_gv, 1, k, stride, gstride, type, op, trans, *this); }

void ogs_t::Unique(hlong *ids, dlong _N, MPI_Comm _comm) {
  ogs::gsUnique(ids, _N, _comm);
}