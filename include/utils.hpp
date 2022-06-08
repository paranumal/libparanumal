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

#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <cstring>
#include <ostream>
#include <iostream>
#include <cstddef>
#include <memory>
#include <algorithm>
#include <typeinfo>
#include <cmath>
#include <occa.hpp>
#include "types.h"

namespace libp {

using properties_t = occa::json;
using device_t = occa::device;
using kernel_t = occa::kernel;
using stream_t = occa::stream;

//error codes
#define LIBP_SUCCESS 0
#define LIBP_ERROR -1

#ifndef __PRETTY_FUNCTION__
#  define __PRETTY_FUNCTION__ __FUNCTION__
#endif

#define LIBP_TEMPLATE_CHECK(checkFunction, expr, filename, function, line, message) \
  do {                                                                  \
    const bool isErr = (bool) (expr);                                   \
    if (isErr) {                                                        \
      std::stringstream _check_ss;                                      \
      _check_ss << message;                                             \
      checkFunction(filename, function, line, _check_ss.str());         \
    }                                                                   \
  } while (false)

#define LIBP_ABORT3(expr, filename, function, line, message) LIBP_TEMPLATE_CHECK(libp::abort, expr, filename, function, line, message)
#define LIBP_ABORT2(expr, filename, function, line, message) LIBP_ABORT3(expr, filename, function, line, message)
#define LIBP_ABORT(message, expr)                            LIBP_ABORT2(expr, __FILE__, __PRETTY_FUNCTION__, __LINE__, message)

#define LIBP_WARNING3(expr, filename, function, line, message) LIBP_TEMPLATE_CHECK(libp::warn, expr, filename, function, line, message)
#define LIBP_WARNING2(expr, filename, function, line, message) LIBP_WARNING3(expr, filename, function, line, message)
#define LIBP_WARNING(message, expr)                            LIBP_WARNING2(expr, __FILE__, __PRETTY_FUNCTION__, __LINE__, message)

#define LIBP_FORCE_ABORT(message)   LIBP_ABORT(message, true)
#define LIBP_FORCE_WARNING(message) LIBP_WARNING(message, true)

class exception : public std::exception {
 public:
  const std::string header;
  const std::string filename;
  const std::string function;
  const std::string message;
  const int line;

  std::string exceptionMessage;

  exception(const std::string &header_,
            const std::string &filename_,
            const std::string &function_,
            const int line_,
            const std::string &message_);
  ~exception() throw();

  const char* what() const throw();
  std::string toString() const;
  std::string location() const;
};

std::ostream& operator << (std::ostream& out,
                           const exception &exc);

void abort(const std::string &filename,
           const std::string &function,
           const int line,
           const std::string &message);

void warn(const std::string &filename,
          const std::string &function,
          const int line,
          const std::string &message);

} //namespace libp

#endif
