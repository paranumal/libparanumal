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

#ifndef LIBP_TIMER_HPP
#define LIBP_TIMER_HPP

#include "core.hpp"
#include "comm.hpp"
#include "platform.hpp"
#include <chrono>

namespace libp {

using timePoint_t = std::chrono::time_point<std::chrono::high_resolution_clock>;

/* Host time*/
timePoint_t Time();

/* Host time after global sync*/
timePoint_t GlobalTime(comm_t comm);

/* Host time after platform sync*/
timePoint_t PlatformTime(platform_t &platform);

/* Host time after platform sync*/
timePoint_t GlobalPlatformTime(platform_t &platform);

/*Time between time points, in seconds*/
double ElapsedTime(const timePoint_t start, const timePoint_t end);

} //namespace libp

#endif
