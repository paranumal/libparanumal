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

#ifndef PRECON_HPP
#define PRECON_HPP

#include "core.hpp"
#include "operator.hpp"

namespace libp {

/*Abstracted Preconditioner Object*/
class precon_t: public operator_t {
 public:
  void Operator(deviceMemory<dfloat> &o_r, deviceMemory<dfloat> &o_Mr) {
    assertInitialized();
    precon->Operator(o_r, o_Mr);
  }

  /*Generic setup. Create a Precon object and wrap it in a shared_ptr*/
  template<class Precon, class... Args>
  void Setup(Args&& ... args) {
    precon = std::make_shared<Precon>(args...);
  }

 private:
  std::shared_ptr<operator_t> precon=nullptr;

  void assertInitialized() {
    LIBP_ABORT("Precon not initialized",
               precon==nullptr);
  }
};

//Identity operator
class IdentityPrecon: public operator_t {
private:
  dlong N;

public:
  IdentityPrecon(dlong _N): N(_N) {}

  void Operator(deviceMemory<dfloat> &o_r, deviceMemory<dfloat> &o_Mr){
    o_Mr.copyFrom(o_r, N); //identity
  }
};

} //namespace libp

#endif
