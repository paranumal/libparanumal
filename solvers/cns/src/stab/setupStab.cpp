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

#include "cns.hpp"

void cns_t::setupStab(){
  // needed gradients of all fields 
  Ngrads = Nfields*mesh.dim;
  // setup trace halo exchange 
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  // Define Global Variables
  props["defines/" "p_Ngrads"]    = Ngrads;
  props["defines/" "s_Nverts"]    = int(mesh.Nverts);
  props["defines/" "s_DGDG_TYPE"] = int(0); 
  props["defines/" "s_FVFV_TYPE"] = int(1); 
  props["defines/" "s_DGFV_TYPE"] = int(2); 

  settings.getSetting("DETECTOR TYPE", detectType);
  settings.getSetting("STAB TYPE", stabType);

  
  switch(detectType){
    case Stab::NODETECT:   break;
    case Stab::ALL:        break; 
    case Stab::KLOCKNER:   setupKlocknerDetector(); break;
    case Stab::PERSSON:    break;
    case Stab::DUCROS:     break; 
    default: LIBP_FORCE_ABORT("Unknown detector type: " << detectType);
  }

  switch (stabType) {
    case Stab::NOSTAB:   break; 
    case Stab::FILTER:   break;
    case Stab::LIMITER:  break;
    case Stab::ARTDIFF:  setupArtificialDiffusion(); break;
    case Stab::SUBCELL:  break;
    default:  LIBP_FORCE_ABORT("Unknown stab type: " << stabType);
  } 
}


void cns_t::applyStab(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T){
switch (stabType) {
    case Stab::NOSTAB:   break; 
    case Stab::FILTER:   break;
    case Stab::LIMITER:  break;
    case Stab::ARTDIFF:  applyArtificialDiffusion(o_Q, o_gradQ, T); break;
    case Stab::SUBCELL:  break;
    default:  LIBP_FORCE_ABORT("Unknown stab type: " << stabType);
  } 
}


void cns_t::applyDetect(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_gradQ, const dfloat T){
 switch(detectType){
    case Stab::NODETECT:   break;
    case Stab::ALL:        break; 
    case Stab::KLOCKNER:   applyKlocknerDetector(o_Q, o_gradQ, T); break;
    case Stab::PERSSON:    break;
    case Stab::DUCROS:     break; 
    default: LIBP_FORCE_ABORT("Unknown detector type: " << detectType);
  }

}




