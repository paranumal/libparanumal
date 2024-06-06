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

#include "stab.hpp"

namespace libp {

void stab_t::setTypes(const Stab::Solver _Solver, 
                      const Stab::Detector _Detector, 
                      const Stab::Type _Type) {

// Set Solver Dependent Parameters
  if (_Solver==Stab::HJS) {
    solver = Stab::HJS;
    // Local Hamilton Jacobi Solver with 2 fields i.e. qp and qm;  
    Ndfields = 2;        // number of fields to be detected
    Nsfields = 2;        // number of fields to be stabilized
  } else if (_Solver==Stab::CNS) {
    solver = Stab::CNS;
    Ndfields  = 1; 
    Nsfields  = mesh.dim==2 ? 4:5; // non-isothermal flow only
  } else if (_Solver==Stab::INS) {
    LIBP_FORCE_ABORT("Unknown solver type: " << solver);
  } else {
    LIBP_FORCE_ABORT("Unknown solver type: " << solver);
  }

// Set Detector Type
if (_Detector==Stab::NODETECT) {
    detector= Stab::NODETECT;
 }else if (_Detector==Stab::ALL) {
    detector = Stab::ALL;
 }else if (_Detector==Stab::KLOCKNER) {
    detector = Stab::KLOCKNER;
 }else if (_Detector==Stab::PERSSON) {
    detector = Stab::PERSSON;
 }else if (_Detector==Stab::DUCROS) {
    detector = Stab::DUCROS;
 }else {
    LIBP_FORCE_ABORT("Unknown detector type: " << detector);
  }


// Set Stabilization Type
if (_Type==Stab::FILTER) {
    type = Stab::FILTER;
 }else if (_Type==Stab::LIMITER) {
     type = Stab::LIMITER;
 }else if (_Type==Stab::ARTDIFF) {
     type = Stab::ARTDIFF;
 }else if (_Type==Stab::SUBCELL) {
     type = Stab::SUBCELL;
 }else if (_Type==Stab::NOSTAB) {
     type = Stab::NOSTAB;
 }else {
    LIBP_FORCE_ABORT("Unknown solver type: " << type);
  }

props["defines/" "s_Ndfields"]= Ndfields;
props["defines/" "s_Nsfields"]= Nsfields;

}


} //namespace libp
