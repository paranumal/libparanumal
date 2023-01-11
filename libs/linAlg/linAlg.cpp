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

#include "linAlg.hpp"
#include "platform.hpp"

namespace libp {
  
  /*********************/
  /* vector operations */
  /*********************/
  
  // o_a[n] = alpha
  template<>
  void linAlg_t::set(const dlong N, const float alpha, deviceMemory<float> o_a) {
    setKernelFloat(N, alpha, o_a);
  }
  
    template<> void linAlg_t::set(const dlong N, const double alpha, deviceMemory<double> o_a) {
    setKernelDouble(N, alpha, o_a);
  }
  
  // o_a[n] += alpha
    template<> void linAlg_t::add(const dlong N, const float alpha, deviceMemory<float> o_a) {
    addKernelFloat(N, alpha, o_a);
  }
    template<> void linAlg_t::add(const dlong N, const double alpha, deviceMemory<double> o_a) {
    addKernelDouble(N, alpha, o_a);
  }
  
  // o_a[n] *= alpha
    template<> void linAlg_t::scale(const dlong N, const float alpha, deviceMemory<float> o_a)  {
    scaleKernelFloat(N, alpha, o_a);
  }
    template<> void linAlg_t::scale(const dlong N, const double alpha, deviceMemory<double> o_a)  {
    scaleKernelDouble(N, alpha, o_a);
  }
  
  // o_y[n] = beta*o_y[n] + alpha*o_x[n]
    template<> void linAlg_t::axpy(const dlong N, const float alpha, deviceMemory<float> o_x,
		      const float beta,  deviceMemory<float> o_y) {
    axpyKernelFloat(N, alpha, o_x, beta, o_y);
  }
  template<> void linAlg_t::axpy(const dlong N, const double alpha, deviceMemory<double> o_x,
		      const double beta,  deviceMemory<double> o_y) {
    axpyKernelDouble(N, alpha, o_x, beta, o_y);
  }
  
  // o_z[n] = beta*o_y[n] + alpha*o_x[n]
    template<> void linAlg_t::zaxpy(const dlong N, const float alpha, deviceMemory<float> o_x,
			  const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
    zaxpyKernelFloat(N, alpha, o_x, beta, o_y, o_z);
  }
    template<> void linAlg_t::zaxpy(const dlong N, const double alpha, deviceMemory<double> o_x,
		       const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
    zaxpyKernelDouble(N, alpha, o_x, beta, o_y, o_z);
  }

  // o_x[n] = alpha*o_a[n]*o_x[n]
    template<> void linAlg_t::amx(const dlong N, const float alpha,
		     deviceMemory<float> o_a, deviceMemory<float> o_x) {
    amxKernelFloat(N, alpha, o_a, o_x);
  }
    template<> void linAlg_t::amx(const dlong N, const double alpha,
		     deviceMemory<double> o_a, deviceMemory<double> o_x) {
    amxKernelDouble(N, alpha, o_a, o_x);
  }


  // o_y[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
    template<> void linAlg_t::amxpy(const dlong N, const float alpha,
			  deviceMemory<float> o_a, deviceMemory<float> o_x,
			  const float beta, deviceMemory<float> o_y) {
    amxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y);
  }
    template<> void linAlg_t::amxpy(const dlong N, const double alpha,
			  deviceMemory<double> o_a, deviceMemory<double> o_x,
			  const double beta, deviceMemory<double> o_y) {
    amxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y);
  }

  // o_z[n] = alpha*o_a[n]*o_x[n] + beta*o_y[n]
    template<> void linAlg_t::zamxpy(const dlong N, const float alpha,
			   deviceMemory<float> o_a, deviceMemory<float> o_x,
			   const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
    zamxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y, o_z);
  }
    template<> void linAlg_t::zamxpy(const dlong N, const double alpha,
			   deviceMemory<double> o_a, deviceMemory<double> o_x,
			   const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
    zamxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y, o_z);
  }

  // o_x[n] = alpha*o_x[n]/o_a[n]
    template<> void linAlg_t::adx(const dlong N, const float alpha,
			deviceMemory<float> o_a, deviceMemory<float> o_x) {
    adxKernelFloat(N, alpha, o_a, o_x);
  }
    template<> void linAlg_t::adx(const dlong N, const double alpha,
			deviceMemory<double> o_a, deviceMemory<double> o_x) {
    adxKernelDouble(N, alpha, o_a, o_x);
  }

  // o_y[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]

    template<> void linAlg_t::adxpy(const dlong N, const float alpha,
			  deviceMemory<float> o_a, deviceMemory<float> o_x,
			  const float beta, deviceMemory<float> o_y) {
    adxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y);
  }
    template<> void linAlg_t::adxpy(const dlong N, const double alpha,
			  deviceMemory<double> o_a, deviceMemory<double> o_x,
			  const double beta, deviceMemory<double> o_y) {
      adxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y);
  }


  // o_z[n] = alpha*o_x[n]/o_a[n] + beta*o_y[n]
    template<> void linAlg_t::zadxpy(const dlong N, const float alpha,
			   deviceMemory<float> o_a, deviceMemory<float> o_x,
			   const float beta, deviceMemory<float> o_y, deviceMemory<float> o_z) {
    zadxpyKernelFloat(N, alpha, o_a, o_x, beta, o_y, o_z);
  }
    template<> void linAlg_t::zadxpy(const dlong N, const double alpha,
			   deviceMemory<double> o_a, deviceMemory<double> o_x,
			   const double beta, deviceMemory<double> o_y, deviceMemory<double> o_z) {
    zadxpyKernelDouble(N, alpha, o_a, o_x, beta, o_y, o_z);
  }

  // \min o_a
  template<>
  float linAlg_t::min(const dlong N, deviceMemory<float> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    //pinned scratch buffer
    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    minKernelFloat(Nblock, N, o_a, o_scratch);

    float globalmin;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalmin = h_scratch[0];
    } else {
      globalmin = std::numeric_limits<float>::max();
    }

    comm.Allreduce(globalmin, Comm::Min);

    return globalmin;
  }

  template<>
  double linAlg_t::min(const dlong N, deviceMemory<double> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    //pinned scratch buffer
    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    minKernelDouble(Nblock, N, o_a, o_scratch);

    double globalmin;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalmin = h_scratch[0];
    } else {
      globalmin = std::numeric_limits<double>::max();
    }

    comm.Allreduce(globalmin, Comm::Min);

    return globalmin;
  }


  // \max o_a
  template<>  
  float linAlg_t::max(const dlong N, deviceMemory<float> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    maxKernelFloat(Nblock, N, o_a, o_scratch);

    float globalmax;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalmax = h_scratch[0];
    } else {
      globalmax = -std::numeric_limits<float>::max();
    }

    comm.Allreduce(globalmax, Comm::Max);

    return globalmax;
  }

  template<>
  double linAlg_t::max(const dlong N, deviceMemory<double> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    maxKernelDouble(Nblock, N, o_a, o_scratch);

    double globalmax;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalmax = h_scratch[0];
    } else {
      globalmax = -std::numeric_limits<double>::max();
    }

    comm.Allreduce(globalmax, Comm::Max);

    return globalmax;
  }


  // \sum o_a
  template<>
  float linAlg_t::sum(const dlong N, deviceMemory<float> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    sumKernelFloat(Nblock, N, o_a, o_scratch);

    float globalsum;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalsum = h_scratch[0];
    } else {
      globalsum = 0.0;
    }

    comm.Allreduce(globalsum, Comm::Sum);

    return globalsum;
  }

  template<>
  double linAlg_t::sum(const dlong N, deviceMemory<double> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    sumKernelDouble(Nblock, N, o_a, o_scratch);

    double globalsum;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalsum = h_scratch[0];
    } else {
      globalsum = 0.0;
    }

    comm.Allreduce(globalsum, Comm::Sum);

    return globalsum;
  }

  // ||o_a||_2

  template<>
  float linAlg_t::norm2(const dlong N, deviceMemory<float> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    norm2KernelFloat(Nblock, N, o_a, o_scratch);

    float globalnorm;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalnorm = h_scratch[0];
    } else {
      globalnorm = 0.0;
    }

    comm.Allreduce(globalnorm, Comm::Sum);

    return sqrt(globalnorm);
  }

  template<>
  double linAlg_t::norm2(const dlong N, deviceMemory<double> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    norm2KernelDouble(Nblock, N, o_a, o_scratch);

    double globalnorm;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalnorm = h_scratch[0];
    } else {
      globalnorm = 0.0;
    }

    comm.Allreduce(globalnorm, Comm::Sum);

    return sqrt(globalnorm);
  }

  // o_x.o_y

  template<>
  float linAlg_t::innerProd(const dlong N, deviceMemory<float> o_x, deviceMemory<float> o_y,
                           comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    innerProdKernelFloat(Nblock, N, o_x, o_y, o_scratch);

    float globaldot;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globaldot = h_scratch[0];
    } else {
      globaldot = 0.0;
    }

    comm.Allreduce(globaldot, Comm::Sum);

    return globaldot;
  }

  template<>
  double linAlg_t::innerProd(const dlong N, deviceMemory<double> o_x, deviceMemory<double> o_y,
                           comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    innerProdKernelDouble(Nblock, N, o_x, o_y, o_scratch);

    double globaldot;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globaldot = h_scratch[0];
    } else {
      globaldot = 0.0;
    }

    comm.Allreduce(globaldot, Comm::Sum);

    return globaldot;
  }

  
  // o_w.o_x.o_y
  template<>
  float linAlg_t::weightedInnerProd(const dlong N, deviceMemory<float> o_w,
				    deviceMemory<float> o_x, deviceMemory<float> o_y,
				    comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    weightedInnerProdKernelFloat(Nblock, N, o_w, o_x, o_y, o_scratch);

    float globaldot;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globaldot = h_scratch[0];
    } else {
      globaldot = 0.0;
    }

    comm.Allreduce(globaldot, Comm::Sum);

    return globaldot;
  }

  template<>
  double linAlg_t::weightedInnerProd(const dlong N, deviceMemory<double> o_w,
                                   deviceMemory<double> o_x, deviceMemory<double> o_y,
                                   comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    weightedInnerProdKernelDouble(Nblock, N, o_w, o_x, o_y, o_scratch);

    double globaldot;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globaldot = h_scratch[0];
    } else {
      globaldot = 0.0;
    }

    comm.Allreduce(globaldot, Comm::Sum);

    return globaldot;
  }


  
  // ||o_a||_w2
  template<>
  float linAlg_t::weightedNorm2(const dlong N, deviceMemory<float> o_w,
                               deviceMemory<float> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<float> h_scratch = platform->hostReserve<float>(1);
    deviceMemory<float> o_scratch = platform->reserve<float>(blocksize);

    weightedNorm2KernelFloat(Nblock, N, o_w, o_a, o_scratch);

    float globalnorm;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalnorm = h_scratch[0];
    } else {
      globalnorm = 0.0;
    }

    comm.Allreduce(globalnorm, Comm::Sum);

    return sqrt(globalnorm);
  }

  template<>
  double linAlg_t::weightedNorm2(const dlong N, deviceMemory<double> o_w,
				 deviceMemory<double> o_a, comm_t comm) {
    int Nblock = (N+blocksize-1)/blocksize;
    Nblock = (Nblock>blocksize) ? blocksize : Nblock; //limit to blocksize entries

    pinnedMemory<double> h_scratch = platform->hostReserve<double>(1);
    deviceMemory<double> o_scratch = platform->reserve<double>(blocksize);

    weightedNorm2KernelDouble(Nblock, N, o_w, o_a, o_scratch);

    double globalnorm;
    if (Nblock>0) {
      h_scratch.copyFrom(o_scratch, 1, 0, properties_t("async", true));
      platform->finish();
      globalnorm = h_scratch[0];
    } else {
      globalnorm = 0.0;
    }

    comm.Allreduce(globalnorm, Comm::Sum);

    return sqrt(globalnorm);
  }
} //namespace libp
