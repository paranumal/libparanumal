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

#include "utils.hpp"

namespace libp {


exception::exception(const std::string &header_,
                      const std::string &filename_,
                      const std::string &function_,
                      const int line_,
                      const std::string &message_) :
    header(header_),
    filename(filename_),
    function(function_),
    message(message_),
    line(line_),
    exceptionMessage(toString()) {}

exception::~exception() throw() {}

const char* exception::what() const throw() {
  return exceptionMessage.c_str();
}

std::string exception::toString() const {
  std::stringstream ss;
  std::string banner = "---[ " + header + " ]";
  ss << '\n'
     << banner << std::string(80 - banner.size(), '-') << '\n'
     << location()
     << "    Message  : " << message << '\n'
     << std::string(80, '=') << '\n';
  return ss.str();
}

std::string exception::location() const {
  std::stringstream ss;
  ss << "    File     : " << filename << '\n'
     << "    Line     : " << line     << '\n'
     << "    Function : " << function << '\n';
  return ss.str();
}

std::ostream& operator << (std::ostream& out,
                           const exception &exc) {
  out << exc.toString() << std::flush;
  return out;
}

void abort(const std::string &filename,
           const std::string &function,
           const int line,
           const std::string &message) {
  throw exception("Error", filename, function, line, message);
}

void warn(const std::string &filename,
          const std::string &function,
          const int line,
          const std::string &message) {
  exception exp("Warning", filename, function, line, message);
  std::cout << exp;
}

} //namespace libp
