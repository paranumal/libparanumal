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

#include "mesh.hpp"

namespace libp {

/* Set up trace halo infomation for inter-processor MPI
   exchange of trace nodes */

// Setup assumes field to be exchanged is Nelements*Nfields*Np in size
// with Np being the fastest running index (hence each field entry is strided
// Np apart)
ogs::halo_t mesh_t::HaloTraceSetup(int Nfields){

  hlong localNelements = Nelements;
  hlong globalOffset = Nelements;
  comm.Scan(localNelements, globalOffset);
  globalOffset -= localNelements;

  //populate a global numbering system which has the Nfields stride
  memory<hlong> globalids((Nelements+totalHaloPairs)*Np*Nfields);
  for (dlong e=0;e<Nelements;e++) {
    for (int k=0;k<Nfields;k++) {
      for (int n=0;n<Np;n++) {
        dlong id = e*Np*Nfields + k*Np + n;
        globalids[id] = (e+globalOffset)*Np*Nfields + k*Np + n + 1;
      }
    }
  }

  //exchange full Np*Nfields per element global ids
  halo.Exchange(globalids, Np*Nfields);

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
  ogs::halo_t traceHalo;
  traceHalo.Setup((Nelements+totalHaloPairs)*Np*Nfields,
                  globalids, comm,
                  ogs::Pairwise, verbose, platform);

  return traceHalo;
}

} //namespace libp
