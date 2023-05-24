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

#ifndef LINALGMATRIX_HPP
#define LINALGMATRIX_HPP 1

#include <iostream>
#include "linAlg.hpp"

extern "C" {
  void sggev_(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA, float *B, int *LDB,
               float *ALPHAR, float *ALPHAI, float *BETA, float *VL, int *LDVL, float *VR, int *LDVR, 
               float *WORK, int *LWORK, int *INFO );
  void dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *B, int *LDB,
               double *ALPHAR, double *ALPHAI, double *BETA, double *VL, int *LDVL, double *VR, int *LDVR, 
               double *WORK, int *LWORK, int *INFO );
  void sgeev_(char *JOBVL, char *JOBVR, int *N, float *A, int *LDA,
               float *ALPHAR, float *ALPHAI,  float *VL, int *LDVL, float *VR, int *LDVR, 
               float *WORK, int *LWORK, int *INFO );
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
               double *ALPHAR, double *ALPHAI, double *VL, int *LDVL, double *VR, int *LDVR, 
               double *WORK, int *LWORK, int *INFO );

  void dsygv_(int *ITYPE, char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *B, int *LDB,
	      double *ALPHAR, 
	      double *WORK, int *LWORK, int *INFO );

  void sgesv_(int *, int *, float *, int *, int *, float *, int *, int *);
  void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);  
}

// this assumes 1-indexing


namespace libp {
template <typename T>
class linAlgMatrix_t  {
private:
   int Nrows;
   int Ncols;

   const int index(int row, int col) const{
     return (row-1)*Ncols + col-1;
   }
public:

   memory<T> data;
   
   // default constructor
   linAlgMatrix_t() : data(1){
     
     Nrows = 1;
     Ncols = 1;

     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         (*this)(r,c) = 0;
       }
     }
   }


   
   // sized constructor
   linAlgMatrix_t(const size_t _Nrows,
                  const size_t _Ncols) : data(_Nrows*_Ncols){
     
     Nrows = _Nrows;
     Ncols = _Ncols;

     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         (*this)(r,c) = 0;
       }
     }
   }

   // copy constructor
   linAlgMatrix_t(const linAlgMatrix_t &A) : data(A.Nrows*A.Ncols){
     
     Nrows = A.Nrows;
     Ncols = A.Ncols;

     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         (*this)(r,c) = A(r,c);
       }
     }
   }

   // sized constructor
   void import(const size_t _Nrows,
               const size_t _Ncols,
               const T* ptr){
     
     reshape(_Nrows, _Ncols);
     
     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         (*this)(r,c) = ptr[(r-1)*Ncols+(c-1)];
       }
     }
   }
   
   
   // resize matrix, while retaining entries, and padding with zero
   void reshape(int _Nrows, int _Ncols){
     
     linAlgMatrix_t<T> tmp(_Nrows,_Ncols);

     for(int r=1;r<=std::min(_Nrows,Nrows);++r){
       for(int c=1;c<=std::min(_Ncols,Ncols);++c){
         tmp(r,c) = (*this)(r,c);
       }
     }

     data.realloc(_Nrows*_Ncols);
     Nrows = _Nrows;
     Ncols = _Ncols;

     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
	 (*this)(r,c) = tmp(r,c);
       }
     }

   }

   // resize matrix, while retaining entries, and padding with zero
   void calloc(int _Nrows, int _Ncols){

     Nrows = _Nrows;
     Ncols = _Ncols;
     
     data.realloc(_Nrows*_Ncols);
     for(int r=1;r<=_Nrows;++r){
       for(int c=1;c<=_Ncols;++c){
         (*this)(r,c) = 0;
       }
     }


   }

   
   
   T &operator()(const int row, const int col){
     if(row<=0 || col<=0 || row>Nrows || col>Ncols){
       std::cout << "Bounds error in linAlgMatrix_t::operator()" << std::endl;
       exit(-1);
     }
     return data[index(row,col)];
   }
   
   const T operator()(const int row, const int col) const{
     if(row<=0 || col<=0 || row>Nrows || col>Ncols){
       std::cout << "Bounds error in linAlgMatrix_t::operator()" << std::endl;
       exit(-1);
     }
     return data[index(row,col)];
   }

   T &operator()(const int ind){
     if(ind<=0 || ind>Nrows*Ncols){
       std::cout << "Bounds error in linAlgMatrix_t::operator()" << std::endl;
       exit(-1);
     }
     return data[ind-1];
   }
   
   const T operator()(const int ind) const {
     if(ind<=0 || ind>Nrows*Ncols){
       std::cout << "Bounds error in linAlgMatrix_t::operator()" << std::endl;
       exit(-1);
     }
     return data[ind-1];
   }


   linAlgMatrix_t<T>& operator=(const linAlgMatrix_t<T> &other){

     if(this==&other){
       return *this;
     }

     if(Nrows!=other.rows() || Ncols!=other.cols()){
       data.calloc(other.rows()*other.cols());
       Nrows = other.rows();
       Ncols = other.cols();
     }
     
     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         (*this)(r,c) = other(r,c);
       }
     }

     return *this;

   }

   
   linAlgMatrix_t<T> transpose(){
     linAlgMatrix_t<T> AT(Ncols, Nrows);
     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=Ncols;++c){
         AT(c,r) = (*this)(r,c);
       }
     }
     // check this
     return AT;
   }

   linAlgMatrix_t<T> inverse(){
     if(Nrows!=Ncols){
       std::cout << "ERROR in linAlgMatrix_t::inverse - matrix not square" << std::endl;
     }

     // use copy constructor
     linAlgMatrix_t Ainv(*this);

     linAlg_t::matrixInverse(Nrows, Ainv.data);

     return Ainv;
   }

   const int rows() const{
     return Nrows;
   }

   const int cols() const{
     return Ncols;
   }

   // C = this/B
   linAlgMatrix_t<T> rightDivide(linAlgMatrix_t<T> &B){
     linAlgMatrix_t<T> C(Nrows, B.cols());
     linAlg_t::matrixRightSolve(Nrows, Ncols, data,
                                B.rows(), B.cols(), B.data, C.data);
     return C;
   }

   // hack this for the moment
  linAlgMatrix_t<T> leftDivide(linAlgMatrix_t<T> &B){

    dlong N = Nrows;
    dlong NRHS = B.cols();
    
    linAlgMatrix_t<int> IPIV(N,1);

    linAlgMatrix_t<T> Acm = this->transpose();
    linAlgMatrix_t<T> Bcm = B.transpose();

    
    int INFO = 0;
    if(sizeof(T)==sizeof(double)){
      
      dgesv_(&N, &NRHS, (double*)Acm.ptr(), &N, IPIV.ptr(), (double*)Bcm.ptr(), &N, &INFO);
    }else{
      sgesv_(&N, &NRHS, (float*)Acm.ptr(), &N, IPIV.ptr(), (float*)Bcm.ptr(), &N, &INFO);
     }

     linAlgMatrix_t<T> C = Bcm.transpose();
     return C;
   }
   
   void eigenVectors(linAlgMatrix_t<T> &VR,
                     linAlgMatrix_t<T> &WR,
                     linAlgMatrix_t<T> &WI){

     if(Nrows!=Ncols){
       std::cout << "ERROR in linAlgMatrix_t::eigenVectors - matrix not square" << std::endl;
     }

     
     // assume sized already ?
     linAlg_t::matrixEigenVectors(Nrows,
                                  data,
                                  VR.data,
                                  WR.data,
                                  WI.dat);
   }

   friend void eig(linAlgMatrix_t<float> &A,
                   linAlgMatrix_t<float> &B,
                   linAlgMatrix_t<float> &VR,
                   linAlgMatrix_t<float> &alphaR,
                   linAlgMatrix_t<float> &alphaI,
                   linAlgMatrix_t<float> &beta);
   friend void eig(linAlgMatrix_t<double> &A,
                   linAlgMatrix_t<double> &B,
                   linAlgMatrix_t<double> &VR,
                   linAlgMatrix_t<double> &alphaR,
                   linAlgMatrix_t<double> &alphaI,
                   linAlgMatrix_t<double> &beta);
   friend void eig(linAlgMatrix_t<double> &A,
                   linAlgMatrix_t<double> &VR,
                   linAlgMatrix_t<double> &alphaR,
                   linAlgMatrix_t<double> &alphaI);



   

   T* ptr(){
     return data.ptr();
   }
   
//   template <typename U>
//   friend std::ostream & operator << (std::ostream &out, const linAlgMatrix_t<U> &matrix);

   // could replace with (d,s)gemm
   linAlgMatrix_t<T> operator * (const linAlgMatrix_t<T> &B){

     linAlgMatrix_t<T> AB(Nrows, B.cols());

#if 0
     for(int r=1;r<=Nrows;++r){
       for(int c=1;c<=B.cols();++c){
         T val = (T)0;
         for(int n=1;n<=Ncols;++n){
           val += (*this)(r,n)*B(n,c);
         }
         AB(r,c) = val;
       }
     }
#else
     linAlg_t::matrixMultiply(Nrows, Ncols, data, B.rows(), B.cols(), B.data, AB.data);
#endif
     return AB;
   }

   linAlgMatrix_t<T> operator + (const linAlgMatrix_t<T> &B){

     if(Nrows*Ncols!=B.rows()*B.cols()){
       std::cout << "linAlgMatrix_t::operator + different number of entries" << std::endl;
     }

     linAlgMatrix_t<T> AplusB(Nrows, Ncols);
     
     for(int n=1;n<=Nrows*Ncols;++n){
       AplusB(n) = (*this)(n)+B(n);
     }
     return AplusB;
   }


  template <typename U>
  friend std::ostream & operator << (std::ostream &out, const linAlgMatrix_t<U> &matrix);

  
};




}

#endif
