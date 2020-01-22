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

#ifndef UTILS_HPP
#define UTILS_HPP

#include <occa.hpp>
#include <mpi.h>
#include <string>
#include "types.h"

//error codes
#define LIBP_SUCCESS 0
#define LIBP_ERROR -1

#ifndef __PRETTY_FUNCTION__
#  define __PRETTY_FUNCTION__ __FUNCTION__
#endif

#define LIBP_ABORT2(filename, function, line, message)              \
  {                                                                 \
    std::string banner = "---[ Error ]";                            \
    std::cerr << '\n'                                               \
       << std::string(74, '=') << '\n'                              \
       << banner << std::string(74 - banner.size(), '-') << '\n'    \
       << "    File     : " << filename << '\n'                     \
       << "    Line     : " << line     << '\n'                     \
       << "    Function : " << function << '\n'                     \
       << "    Message  : " << message  << '\n'                     \
       << std::string(74, '=') << '\n';                             \
    MPI_Abort(MPI_COMM_WORLD,LIBP_ERROR);                           \
  }
#define LIBP_ABORT(message) LIBP_ABORT2(__FILE__, __PRETTY_FUNCTION__, __LINE__, message)

#define LIBP_WARNING(message)                                       \
  {                                                                 \
    std::string banner = "---[ Warning ]";                          \
    std::cerr << '\n'                                               \
       << std::string(74, '=') << '\n'                              \
       << banner << std::string(74 - banner.size(), '-') << '\n'    \
       << "     " << message  << '\n'                               \
       << std::string(74, '=') << '\n';                             \
  }

#define mymax(a,b) (((a)>(b))?(a):(b))
#define mymin(a,b) (((a)<(b))?(a):(b))

// block size for reduction (hard coded)
#define BLOCKSIZE 256



#endif
