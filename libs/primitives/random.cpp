/*

The MIT License (MIT)

Copyright (c) 2017-2023 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "primitives.hpp"
#include <limits>
#include "omp.h"

namespace libp {

namespace prim {

/*
The  pseudo-random  generator uses the linear congruential algorithm:
X(n+1) = (a * X(n) + c) mod m  as  described  in the  Art of Computer
Programming, Knuth 1973, Vol. 2.

Not the best RNG, but we care about it being reasonably fast and, more
importantly, we want prim::random to generate the same sequence of
random numbers, whether or not it was computed in parallel
*/

#define LCG_A 6364136223846793005ULL
#define LCG_C 1ULL

#define SEED 1ULL

struct RandCoeff {
  uint64_t a;
  uint64_t c;

  static RandCoeff default_vals() { return {LCG_A, LCG_C}; }

  RandCoeff operator*(const RandCoeff& rhs) const {
    return {a * rhs.a, a * rhs.c + c};
  }

  void operator*=(const RandCoeff& rhs) {
    c = a * rhs.c + c;
    a = a * rhs.a;
  }
};

RandCoeff pow(RandCoeff base, uint32_t n) {
  RandCoeff result{1, 0};
  while(n != 0) {
    if(n & 1) result *= base;
    n >>= 1;
    base *= base;
  }
  return result;
}

struct RandState {
  uint64_t x;

  static RandState initialize(const uint64_t seed = SEED,
                              RandCoeff coef = RandCoeff::default_vals()) {
    return coef * RandState{std::numeric_limits<uint64_t>::max() ^ seed};
  }

  void operator*=(RandCoeff coef) {
    x = coef.a * x + coef.c;
  }

  template<typename T>
  T get();

  friend RandState operator*(RandCoeff coef, RandState stat) {
    return {coef.a * stat.x + coef.c};
  }
};

template<>
int RandState::get() {
  return std::abs(static_cast<int>(x));
}

template<>
long long int RandState::get() {
  return std::abs(static_cast<long long int>(x));
}

template<>
float RandState::get() {
  return x*(1.0/std::numeric_limits<uint64_t>::max());
}

template<>
double RandState::get() {
  return x*(1.0/std::numeric_limits<uint64_t>::max());
}

static RandState state = RandState::initialize();


void seedRNG(const uint64_t seed) {
  state = RandState::initialize(seed);
}

template <typename T>
void random(const dlong N, memory<T> v) {

  constexpr int blockSize = 512;
  const dlong Nblocks = (N+blockSize-1)/blockSize;

  // compute increments
  RandCoeff step1 = RandCoeff::default_vals();
  RandCoeff stepBlock = pow(step1, blockSize);

  #pragma omp parallel
  {
#if !defined(LIBP_DEBUG)
    const int thread = omp_get_thread_num();
    const int Nthreads = omp_get_num_threads();
#else
    const int thread = 0;
    const int Nthreads = 1;
#endif

    // Get the rng state
    RandState rng = state;

    for (dlong b=0;b<Nblocks;++b) {
      if ((b % Nthreads) == thread) { //check if my thread does this block
        RandState curr_rng = rng;

        const int M = std::min(N - b*blockSize, blockSize);
        for (int n=0;n<M;++n) {
          v[b*blockSize + n] = curr_rng.get<T>();
          curr_rng *= step1;
        }
      }

      // Shift rng down a block
      rng *= stepBlock;
    }
  }

  // Update the global rng state N steps
  state *= pow(step1, N);
}

template void random(const dlong N, memory<int> v);
template void random(const dlong N, memory<long long int> v);
template void random(const dlong N, memory<float> v);
template void random(const dlong N, memory<double> v);

} //namespace prim

} //namespace libp
