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

#include "mesh.hpp"

/* Set up trace halo infomation for inter-processor MPI
   exchange of trace nodes */

// Setup assumes field to be exchanged is Nelements*Nfields**Nfaces*Nfp in size
// with Nfp being the fastest running index (hence each field entry is strided
// Nfaces*Nfp apart)
halo_t** mesh_t::MultiRateHaloTraceSetup(int Nfields){

  hlong *globalOffsets = (hlong *) calloc(size+1,sizeof(hlong));
  hlong localNelements = (hlong) Nelements;

  //gather number of elements on each rank
  MPI_Allgather(&localNelements, 1, MPI_HLONG, globalOffsets+1, 1, MPI_HLONG, comm);

  for(int rr=0;rr<size;++rr)
    globalOffsets[rr+1] = globalOffsets[rr]+globalOffsets[rr+1];

  hlong globalOffset = globalOffsets[rank];
  free(globalOffsets);

  //populate a global numbering system which has the Nfields stride
  hlong *globalids = (hlong *) calloc((Nelements+totalHaloPairs)
                                       *Np*Nfields,sizeof(hlong));
  for (dlong e=0;e<Nelements;e++) {
    for (int k=0;k<Nfields;k++) {
      for (int n=0;n<Np;n++) {
        dlong id = e*Np*Nfields + k*Np + n;
        globalids[id] = (e+globalOffset)*Np*Nfields + k*Np + n + 1;
      }
    }
  }

  //make a trace array populated with the global ids
  hlong *traceIds = (hlong *) calloc((Nelements+totalHaloPairs)
                                       *Nfp*Nfaces*Nfields,sizeof(hlong));

  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Nfp*Nfaces;n++) {
      const dlong vid = e*Nfp*Nfaces + n;
      const dlong idM = vmapM[vid];

      const dlong eM = idM/(Nfp*Nfaces);
      const int fidM = idM%(Nfp*Nfaces);

      const dlong id  = eM*Nfields*Np + fidM;
      const dlong tid = e*Nfields*Nfp*Nfaces + n;

      for (int k=0;k<Nfields;k++) {
        traceIds[tid+k*Nfp*Nfaces] = globalids[id+k*Np];
      }
    }
  }

  //exchange full Nfp*Nfaces*Nfields per element global trace ids
  halo->Exchange(traceIds, Nfp*Nfaces*Nfields, ogs_hlong);

  //the halo region is filled, but there are duplicate IDs in the local section
  // bad news for the halo exchange, so remove them
  for (dlong e=0;e<Nelements;e++) {
    for (int n=0;n<Nfp*Nfaces;n++) {
      const dlong vid = e*Nfp*Nfaces + n;
      const dlong idM = vmapM[vid];

      const dlong eM = idM/(Nfp*Nfaces);
      const int fidM = idM%(Nfp*Nfaces);

      const dlong id  = eM*Nfields*Np + fidM;
      const dlong tid = e*Nfields*Nfp*Nfaces + n;

      for (int k=0;k<Nfields;k++) {
        traceIds[tid+k*Nfp*Nfaces] = globalids[id+k*Np];
        globalids[id+k*Np] = 0; //zero out the id after copying it to avoid duplicate
                                // ids in the halo exchange (corners and edges)
      }
    }
  }

  //make array of halo exchangers
  halo_t** mrTraceHalo = (halo_t **) malloc(mrNlevels*sizeof(halo_t*));

  //make a global trace id array to be used for exchange on each multirate level
  hlong *mrTraceIds = (hlong *) calloc((Nelements+totalHaloPairs)
                                       *Nfp*Nfaces*Nfields,sizeof(hlong));
  memcpy(mrTraceIds, traceIds, Nelements*Nfp*Nfaces*Nfields*sizeof(hlong)); //copy local part

  //for each multirate level
  for (int lev=0;lev<mrNlevels;lev++) {

    //fill the trace ids we need
    for (dlong m=0;m<mrNelements[lev];m++) { //for all elements in multirate level

      const dlong e = mrElements[lev][m]; //element id

      for (int f=0;f<Nfaces;f++) {
        for (int n=0;n<Nfp;n++) {
          dlong id  = e*Nfp*Nfaces + f*Nfp + n;
          dlong idP = mapP[id];

          dlong eP = idP/(Nfp*Nfaces);
          int fidP = idP%(Nfp*Nfaces);
          if (eP >= Nelements) { //neighbor is in halo
            dlong iid = eP*Nfp*Nfaces*Nfields + fidP;
            for (int k=0;k<Nfields;k++) {
              mrTraceIds[iid+k*Nfp*Nfaces] = -traceIds[iid+k*Nfp*Nfaces]; //flag trace ids
            }
          }
        }
      }
    }

    int verbose = 0;
    mrTraceHalo[lev] = halo_t::Setup((Nelements+totalHaloPairs)*Nfp*Nfaces*Nfields,
                                      mrTraceIds, comm, verbose, device);

    //no need to zero out mrTraceIds for next multirate level
    // the next level set includes the lower level elements
  }

  free(globalids);
  free(traceIds);
  free(mrTraceIds);

  return mrTraceHalo;
}
