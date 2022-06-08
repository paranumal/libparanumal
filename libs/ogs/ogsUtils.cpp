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

#include <limits>
#include "ogs.hpp"
#include "ogs/ogsOperator.hpp"
#include "ogs/ogsExchange.hpp"
#include "ogs/ogsUtils.hpp"

namespace libp {

namespace ogs {

stream_t ogsBase_t::dataStream;

kernel_t ogsOperator_t::gatherScatterKernel[4][4];
kernel_t ogsOperator_t::gatherKernel[4][4];
kernel_t ogsOperator_t::scatterKernel[4];

kernel_t ogsExchange_t::extractKernel[4];


void InitializeKernels(platform_t& platform, const Type type, const Op op) {

  //check if the gather kernel is initialized
  if (!ogsOperator_t::gatherKernel[type][op].isInitialized()) {

    properties_t kernelInfo = platform.props();

    kernelInfo["defines/p_blockSize"] = ogsOperator_t::blockSize;
    kernelInfo["defines/p_gatherNodesPerBlock"] = ogsOperator_t::gatherNodesPerBlock;

    switch (type) {
      case Float:  kernelInfo["defines/T"] =  "float"; break;
      case Double: kernelInfo["defines/T"] =  "double"; break;
      case Int32:  kernelInfo["defines/T"] =  "int32_t"; break;
      case Int64:  kernelInfo["defines/T"] =  "int64_t"; break;
    }

    switch (type) {
      case Float:
        switch (op) {
          case Add: kernelInfo["defines/OGS_OP_INIT"] =  float{0}; break;
          case Mul: kernelInfo["defines/OGS_OP_INIT"] =  float{1}; break;
          case Min: kernelInfo["defines/OGS_OP_INIT"] =  std::numeric_limits<float>::max(); break;
          case Max: kernelInfo["defines/OGS_OP_INIT"] = -std::numeric_limits<float>::max(); break;
        }
        break;
      case Double:
        switch (op) {
          case Add: kernelInfo["defines/OGS_OP_INIT"] =  double{0}; break;
          case Mul: kernelInfo["defines/OGS_OP_INIT"] =  double{1}; break;
          case Min: kernelInfo["defines/OGS_OP_INIT"] =  std::numeric_limits<double>::max(); break;
          case Max: kernelInfo["defines/OGS_OP_INIT"] = -std::numeric_limits<double>::max(); break;
        }
        break;
      case Int32:
        switch (op) {
          case Add: kernelInfo["defines/OGS_OP_INIT"] =  int32_t{0}; break;
          case Mul: kernelInfo["defines/OGS_OP_INIT"] =  int32_t{1}; break;
          case Min: kernelInfo["defines/OGS_OP_INIT"] =  std::numeric_limits<int32_t>::max(); break;
          case Max: kernelInfo["defines/OGS_OP_INIT"] = -std::numeric_limits<int32_t>::max(); break;
        }
        break;
      case Int64:
        switch (op) {
          case Add: kernelInfo["defines/OGS_OP_INIT"] =  int64_t{0}; break;
          case Mul: kernelInfo["defines/OGS_OP_INIT"] =  int64_t{1}; break;
          case Min: kernelInfo["defines/OGS_OP_INIT"] =  std::numeric_limits<int64_t>::max(); break;
          case Max: kernelInfo["defines/OGS_OP_INIT"] = -std::numeric_limits<int64_t>::max(); break;
        }
        break;
    }

    switch (op) {
      case Add: kernelInfo["defines/OGS_OP(a,b)"] = "a+=b"; break;
      case Mul: kernelInfo["defines/OGS_OP(a,b)"] = "a*=b"; break;
      case Min: kernelInfo["defines/OGS_OP(a,b)"] = "if(b<a) a=b"; break;
      case Max: kernelInfo["defines/OGS_OP(a,b)"] = "if(b>a) a=b"; break;
    }

    ogsOperator_t::gatherScatterKernel[type][op] = platform.buildKernel(OGS_DIR "/okl/ogsKernels.okl",
                                                         "gatherScatter",
                                                         kernelInfo);


    ogsOperator_t::gatherKernel[type][op] = platform.buildKernel(OGS_DIR "/okl/ogsKernels.okl",
                                                "gather",
                                                kernelInfo);

    if (!ogsOperator_t::scatterKernel[type].isInitialized()) {
      ogsOperator_t::scatterKernel[type] = platform.buildKernel(OGS_DIR "/okl/ogsKernels.okl",
                                                 "scatter",
                                                 kernelInfo);

      ogsExchange_t::extractKernel[type] = platform.buildKernel(OGS_DIR "/okl/ogsKernels.okl",
                                                "extract", kernelInfo);\
    }
  }
}

} //namespace ogs

} //namespace libp
