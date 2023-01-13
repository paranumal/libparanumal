/*

  The MIT License (MIT)

  Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Anthony Austin

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

#ifndef INITIALGUESS_HPP
#define INITIALGUESS_HPP

#include "core.hpp"
#include "platform.hpp"
#include "solver.hpp"

namespace libp {

  namespace InitialGuess {

    void AddSettings(settings_t& settings, const std::string prefix = "");
  
    // Abstract base class for different initial guess strategies.
    class initialGuessStrategy_t {
    protected:
      platform_t platform;
      settings_t settings;
      comm_t   comm;

      dlong Ntotal;     // Degrees of freedom

    public:
      initialGuessStrategy_t(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm):
	platform(_platform), settings(_settings), comm(_comm), Ntotal(_N) {}
  
      virtual void FormInitialGuess(deviceMemory<float>& o_x, deviceMemory<float>& o_rhs)  { printf("AROEOHAER\n"); }
      virtual void Update(operator_t& linearOperator, deviceMemory<float>& o_x, deviceMemory<float>& o_rhs) { printf("AROEOHAER\n"); }
  
      virtual void FormInitialGuess(deviceMemory<double>& o_x, deviceMemory<double>& o_rhs) { printf("AROEOHAER\n"); }
      virtual void Update(operator_t& linearOperator, deviceMemory<double>& o_x, deviceMemory<double>& o_rhs) { printf("AROEOHAER\n"); }
    };
  
    // Default initial guess strategy:  use whatever the last solution was (starting at the zero vector)
    template <typename T>
    class Last : public initialGuessStrategy_t {
    private:
      deviceMemory<T> o_xLast;
    
    public:
      Last(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

      Last<T>& operator=(const Last<T> &m)=default;
      
      void FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
      void Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
    };

    template class Last<double>;
    template class Last<float>;
  
  
    // Zero initial guess strategy:  use a zero initial guess.
    template <typename T>  class Zero : public initialGuessStrategy_t {
    public:
      Zero(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);
      void FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
      void Update(operator_t& linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs); 
    };

    template class Zero<double>;
    template class Zero<float>;

    // Initial guess strategies based on RHS projection.
    template<typename T>
    class Projection : public initialGuessStrategy_t {
    protected:

      using initialGuessStrategy_t::platform;
      using initialGuessStrategy_t::comm;
      using initialGuessStrategy_t::settings;
      using initialGuessStrategy_t::Ntotal;      
      
      int curDim;           // Current dimension of the initial guess space
      int maxDim;           // Maximum dimension of the initial guess space

      deviceMemory<T> o_Btilde;  //  space (orthogonalized)
      deviceMemory<T> o_Xtilde;  // Solution space corresponding to  space

      kernel_t igBasisInnerProductsKernel;
      kernel_t igReconstructKernel;
      kernel_t igUpdateKernel;

      void igBasisInnerProducts(deviceMemory<T>& o_x,
				deviceMemory<T>& o_Q,
				deviceMemory<T>& o_alphas,
				pinnedMemory<T>& alphas);

      void igReconstruct(const T a,
			 deviceMemory<T>& o_u,
			 const T b,
			 deviceMemory<T>& o_alphas,
			 deviceMemory<T>& o_Q,
			 deviceMemory<T>& o_unew);
  
    public:
      Projection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);
    
      void FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
      virtual void Update(operator_t& linearOperator, deviceMemory<float>& o_x, deviceMemory<float>& o_rhs) {printf("Update fail\n"); }
      virtual void Update(operator_t& linearOperator, deviceMemory<double>& o_x, deviceMemory<double>& o_rhs) {printf("Update fail\n"); }
    };

    template class Projection<double>;
    template class Projection<float>;
  
    // "Classic" initial guess strategy from Fischer's 1998 paper.
    template <typename T>
    class ClassicProjection : public Projection<T> {
    public:

      using Projection<T>::curDim;
      using Projection<T>::maxDim;
      using Projection<T>::o_Btilde;
      using Projection<T>::o_Xtilde;
      
      ClassicProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);
    
      void Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
    };


    template class ClassicProjection<double>;
    template class ClassicProjection<float>;

  
    // Rolling QR update for projection history space a la Christensen's thesis.
    template <typename T>
    class RollingQRProjection : public Projection<T> {
    private:

      using Projection<T>::curDim;
      using Projection<T>::maxDim;
      using Projection<T>::o_Btilde;
      using Projection<T>::o_Xtilde;
      
      using Projection<T>::initialGuessStrategy_t::platform;
      using Projection<T>::initialGuessStrategy_t::comm;
      using Projection<T>::initialGuessStrategy_t::settings;
      using Projection<T>::initialGuessStrategy_t::Ntotal;      

      
      memory<T>   R;   // R factor in QR decomposition (row major)

      pinnedMemory<T> h_c, h_s;
      deviceMemory<T> o_c, o_s;

      kernel_t igDropQRFirstColumnKernel;

      void givensRotation(T a, T b, T& c, T& s);

    public:
      RollingQRProjection(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

      void Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
    };

    template class RollingQRProjection<double>;
    template class RollingQRProjection<float>;

#if 0
    // Extrapolation initial guess strategy.
    template <typename T>
    class Extrap : public initialGuessStrategy_t {
    private:
      int Nhistory;
      int ExtrapDegree;
      int shift;
      int entry;
      deviceMemory<T> o_xh;

      pinnedMemory<T> h_coeffs;
      deviceMemory<T> o_coeffs;

      int Nsparse;
      pinnedMemory<int> h_sparseIds;
      deviceMemory<int> o_sparseIds;
      pinnedMemory<T> h_sparseCoeffs;
      deviceMemory<T> o_sparseCoeffs;

      kernel_t igExtrapKernel;
      kernel_t igExtrapSparseKernel;

      void extrapCoeffs(int m, int M, memory<T> c);

    public:
      Extrap(dlong _N, platform_t& _platform, settings_t& _settings, comm_t _comm);

      void FormInitialGuess(deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
      void Update(operator_t &linearOperator, deviceMemory<T>& o_x, deviceMemory<T>& o_rhs);
    };

    template class Extrap<double>;
    template class Extrap<float>;
#endif

  } //namespace InitialGuess

} //namespace libp

#endif /* INITIALGUESS_HPP */
