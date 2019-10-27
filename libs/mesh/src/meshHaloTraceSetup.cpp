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

// Setup assumes field to be exchanged is Nelements*Nfields*Np in size
// with Np being the fastest running index (hence each field entry is strided
// Np apart)
halo_t* mesh_t::HaloTraceSetup(int Nfields){

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

  //exchange full Np*Nfields per element global ids
  halo->Exchange(globalids, Np*Nfields, ogs_hlong);

  //flag the trace ids we need
  for (dlong e=0;e<Nelements;e++) {
    for (int f=0;f<Nfaces;f++) {
      for (int n=0;n<Nfp;n++) {
        dlong id  = e*Nfp*Nfaces + f*Nfp + n;
        dlong idP = vmapP[id];

        dlong eP = idP/Np;
        int nP   = idP%Np;
        if (eP >= Nelements) { //neighbor is in halo
          dlong iid = eP*Np*Nfields + nP;
          for (int k=0;k<Nfields;k++) {
            globalids[iid+k*Np] *= -1; //flag trace ids
          }
        }
      }
    }
  }
  //set the remaining globalids to zero so they are ignored
  for (dlong e=Nelements;e<Nelements+totalHaloPairs;e++) {
    for (int n=0;n<Np*Nfields;n++) {
      if (globalids[e*Np*Nfields + n]>0) globalids[e*Np*Nfields + n] = 0;
    }
  }

  int verbose = 0;
  halo_t* traceHalo = halo_t::Setup((Nelements+totalHaloPairs)*Np*Nfields,
                                    globalids, comm, verbose, device);

  free(globalids);

  return traceHalo;
}
