
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

#include "platform.hpp"
#include "linAlgMatrix.hpp"

namespace libp {

/**************************/
/* host matrix operations */
/**************************/


//class linAlgMatrix_t<double>;
//class linAlgMatrix_t<float>;

//std::ostream & operator << (std::ostream &out, const linAlgMatrix_t<double> &matrix);
//std::ostream & operator << (std::ostream &out, const linAlgMatrix_t<float> &matrix);


void meh(){
}

void eig(linAlgMatrix_t<float> &A,
         linAlgMatrix_t<float> &B,
         linAlgMatrix_t<float> &VR,
         linAlgMatrix_t<float> &alphaR,
         linAlgMatrix_t<float> &alphaI,
         linAlgMatrix_t<float> &beta){
  
  if(A.rows()!=A.cols()) std::cout << "in linAlgMatrix eig: matrix A not square" << std::endl;
  if(B.rows()!=B.cols()) std::cout << "in linAlgMatrix eig: matrix B not square" << std::endl;
  
  int n = A.rows();
  char JOBVL  = 'V'; // compute left eigenvectors  because A and B are assumed row major
  char JOBVR  = 'N';
  int LDA = n;
  int LDB = n;
  int LDVL = n;
  int LDVR = n;
  int LWORK = 8*n;
  
  linAlgMatrix_t<float> VRT(n,n);

  linAlgMatrix_t<float> copyA = A;
  linAlgMatrix_t<float> copyB = B;

  
  float* Aptr = const_cast<float*>(copyA.data.ptr());
  float* Bptr = const_cast<float*>(copyB.data.ptr());
  float* VRTptr = const_cast<float*>(VRT.data.ptr());
  float* alphaRptr = const_cast<float*>(alphaR.data.ptr());
  float* alphaIptr = const_cast<float*>(alphaI.data.ptr());
  float* betaptr   = const_cast<float*>(beta.data.ptr());     
  
  memory<float> WORK(LWORK);
  
  int INFO = -999;
  
  sggev_ (&JOBVL, &JOBVR, &n,
           Aptr, &LDA,
           Bptr, &LDB,
           alphaRptr, alphaIptr,
           betaptr,
           VRTptr,  &LDVR, 
           nullptr, &LDVL,
           WORK.ptr(), &LWORK,
           &INFO);
  
  LIBP_ABORT("sggev_ reports info = " << INFO, INFO);
  VR = VRT.transpose();
}


void eig(linAlgMatrix_t<double> &A,
         linAlgMatrix_t<double> &B,
         linAlgMatrix_t<double> &VR,
         linAlgMatrix_t<double> &alphaR,
         linAlgMatrix_t<double> &alphaI,
         linAlgMatrix_t<double> &beta){
  
  if(A.rows()!=A.cols()) std::cout << "in linAlgMatrix eig: matrix A not square" << std::endl;
  if(B.rows()!=B.cols()) std::cout << "in linAlgMatrix eig: matrix B not square" << std::endl;
  
  int n = A.rows();
  char JOBVL  = 'V'; // compute left eigenvectors  because A and B are assumed row major
  char JOBVR  = 'N';
  int LDA = n;
  int LDB = n;
  int LDVL = n;
  int LDVR = n;
  int LWORK = 8*n;
  
  linAlgMatrix_t<double> VRT(n,n);

  linAlgMatrix_t<double> copyA = A;
  linAlgMatrix_t<double> copyB = B;
  
  double* Aptr = const_cast<double*>(copyA.data.ptr());
  double* Bptr = const_cast<double*>(copyB.data.ptr());
  double* VRTptr = const_cast<double*>(VRT.data.ptr());
  double* alphaRptr = const_cast<double*>(alphaR.data.ptr());
  double* alphaIptr = const_cast<double*>(alphaI.data.ptr());
  double* betaptr   = const_cast<double*>(beta.data.ptr());     
  
  memory<double> WORK(LWORK);
  
  int INFO = -999;

#if 0
  dggev_ (&JOBVL, &JOBVR, &n,
	  Aptr, &LDA,
	  Bptr, &LDB,
	  alphaRptr, alphaIptr,
	  betaptr,
	  VRTptr,  &LDVR, 
	  nullptr, &LDVL,
	  WORK.ptr(), &LWORK,
	  &INFO);

   VR = VRT.transpose();
  LIBP_ABORT("dggev_ reports info = " << INFO, INFO);
#else
  int ITYPE = 1;
  char JOBZ = 'V';
  char UPLO = 'U';
  dsygv_(&ITYPE, &JOBZ, &UPLO, &n, Aptr, &LDA, Bptr, &LDB, alphaRptr, WORK.ptr(), &LWORK, &INFO);

  VR = copyA.transpose();
  for(int i=0;i<n;++i){
    alphaIptr[i] = 0;
    betaptr[i] =1;
  }

  LIBP_ABORT("dsygv_ reports info = " << INFO, INFO);
#endif  
  


#if 0
  // A*x = B*x*lambda
  linAlgMatrix_t<double> AVR = A*VR;
  linAlgMatrix_t<double> BVR = B*VR;

  // check for symm
  std::cout << "symm eig error" << std::endl;
  for(int r=1;r<=A.rows();++r){
    for(int c=1;c<=A.cols();++c){
      double res = AVR(r,c) - (alphaR(c)/beta(c))*BVR(r,c);
      std::cout << res << " " ;
    }
    std::cout << std::endl;
  }
  std::cout << "end symm error" << std::endl;
#endif
}

void eig(linAlgMatrix_t<double> &A,
         linAlgMatrix_t<double> &VR,
         linAlgMatrix_t<double> &lambdaR,
         linAlgMatrix_t<double> &lambdaI){
  
  if(A.rows()!=A.cols()) std::cout << "in linAlgMatrix eig: matrix A not square" << std::endl;
  
  int n = A.rows();
  char JOBVL  = 'V'; // compute left eigenvectors  because A and B are assumed row major
  char JOBVR  = 'N';
  int LDA = n;
  int LDVL = n;
  int LDVR = n;
  int LWORK = 8*n;
  
  linAlgMatrix_t<double> VRT(n,n);

  linAlgMatrix_t<double> copyA = A;
  
  double* Aptr = const_cast<double*>(copyA.data.ptr());
  double* VRTptr = const_cast<double*>(VRT.data.ptr());
  double* lambdaRptr = const_cast<double*>(lambdaR.data.ptr());
  double* lambdaIptr = const_cast<double*>(lambdaI.data.ptr());
  
  memory<double> WORK(LWORK);
  
  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &n,
	  Aptr, &LDA,
	  lambdaRptr, lambdaIptr,
	  VRTptr,  &LDVR, 
	  nullptr, &LDVL,
	  WORK.ptr(), &LWORK,
	  &INFO);

  LIBP_ABORT("dggev_ reports info = " << INFO, INFO);
  VR = VRT.transpose();

#if 0
  // A*x = B*x*lambda
  linAlgMatrix_t<double> AVR = A*VR;
  linAlgMatrix_t<double> BVR = B*VR;

  // check for symm
  std::cout << "symm eig error" << std::endl;
  for(int r=1;r<=A.rows();++r){
    for(int c=1;c<=A.cols();++c){
      double res = AVR(r,c) - (alphaR(c)/beta(c))*BVR(r,c);
      std::cout << res << " " ;
    }
    std::cout << std::endl;
  }
  std::cout << "end symm error" << std::endl;
#endif
}


  
  template<typename T>
  std::ostream & operator <<  (std::ostream &out, const linAlgMatrix_t<T> &matrix){
    
    out << " [" << std::endl;
    for(int r=1;r<=matrix.rows();++r){
      for(int c=1;c<=matrix.cols();++c){
	out << " " << matrix(r,c);
      }
      out << std::endl;
    }
    out << "]" << std::endl;
    return out;
  }


}

