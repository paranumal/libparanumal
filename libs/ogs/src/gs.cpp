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

extern "C" {
#include "gslib.h"
}

namespace ogs {

OGS_DEFINE_TYPE_SIZES()
OGS_GS_DEFINE_TYPE_MAP()
OGS_GS_DEFINE_OP_MAP()

// MPI based gather scatter using libgs
void gsGatherScatter(void* v,
                     const dlong Nentries,
                     const dlong Nvectors,
                     const dlong stride,
                     const ogs_type type,
                     const ogs_op op,
                     void *gshSym){

  const gs_op  gsop  = ogs_gs_op_map[op];
  const gs_dom gsdom = ogs_gs_type_map[type];

  //call libgs (symmetric behaviour)
  if (Nentries==1 && Nvectors==1)
    gs(v, gsdom, gsop, 0, (gs_data*)gshSym, 0);
  else if (Nvectors==1)
    gs_vec(v, Nentries, gsdom, gsop, 0, (gs_data*)gshSym, 0);
  else if (Nentries==1) {
    const size_t Nbytes = ogs_type_size[type];
    void* V[Nvectors];
    for (int i=0;i<Nvectors;i++)
      V[i] = (char*)v + i*stride*Nbytes;

    gs_many(V, Nvectors, gsdom, gsop, 0, (gs_data*)gshSym, 0);
  }
}

// MPI based gather using libgs
void gsGather(void* v,
              const dlong Nentries,
              const dlong Nvectors,
              const dlong stride,
              const ogs_type type,
              const ogs_op op,
              void *gshNonSym){

  const gs_op  gsop  = ogs_gs_op_map[op];
  const gs_dom gsdom = ogs_gs_type_map[type];

  //call libgs (non-symmetric behaviour)
  if (Nentries==1 && Nvectors==1)
    gs(v, gsdom, gsop, 1, (gs_data*)gshNonSym, 0);
  else if (Nvectors==1)
    gs_vec(v, Nentries, gsdom, gsop, 1, (gs_data*)gshNonSym, 0);
  else if (Nentries==1) {
    const size_t Nbytes = ogs_type_size[type];
    void* V[Nvectors];
    for (int i=0;i<Nvectors;i++)
      V[i] = (char*)v + i*stride*Nbytes;

    gs_many(V, Nvectors, gsdom, gsop, 1, (gs_data*)gshNonSym, 0);
  }
}

// MPI based scatter using libgs
void gsScatter(void* v,
               const dlong Nentries,
               const dlong Nvectors,
               const dlong stride,
               const ogs_type type,
               const ogs_op op,
               void *gshNonSym){

  const gs_op  gsop  = ogs_gs_op_map[op];
  const gs_dom gsdom = ogs_gs_type_map[type];

  //call libgs (non-symmetric behaviour)
  if (Nentries==1 && Nvectors==1)
    gs(v, gsdom, gsop, 0, (gs_data*)gshNonSym, 0);
  else if (Nvectors==1)
    gs_vec(v, Nentries, gsdom, gsop, 0, (gs_data*)gshNonSym, 0);
  else if (Nentries==1) {
    const size_t Nbytes = ogs_type_size[type];
    void* V[Nvectors];
    for (int i=0;i<Nvectors;i++)
      V[i] = (char*)v + i*stride*Nbytes;

    gs_many(V, Nvectors, gsdom, gsop, 0, (gs_data*)gshNonSym, 0);
  }
}

//Setup a gslib struct
void *gsSetup(MPI_Comm meshComm,
              dlong NuniqueBases,
              hlong *gatherGlobalNodes,
              int unique, int verbose){

  /* gslib stuff */
  comm_ext world;
  struct comm com;

  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)meshComm; // MPI_COMM_WORLD;

  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);

  dlong n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = (slong) gatherGlobalNodes[n];
  }

  struct gs_data *gsh = gs_setup(id, NuniqueBases, &com, unique, gs_auto, verbose); // gs_auto, gs_crystal_router, gs_pw

  free(id);

  return gsh;
}

void gsUnique(hlong *gatherGlobalNodes,
              dlong NuniqueBases,
              MPI_Comm meshComm){

  /* gslib stuff */
  comm_ext world;
  struct comm com;

  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)meshComm; // MPI_COMM_WORLD;

  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);

  dlong n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = (slong) gatherGlobalNodes[n];
  }

  gs_unique(id, NuniqueBases, &com);

  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    gatherGlobalNodes[n] = (hlong) id[n];
  }

  free(id);
}

void gsFree(void* gs) {
  gs_free((gs_data*)gs);
}

} //namespace ogs