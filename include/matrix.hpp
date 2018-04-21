#ifndef __MATRIX
#define __MATRIX

#include "mesh.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <execinfo.h>
//#include "blaslapack.hpp"

static int clarrayrequest = 0;

template <class T>
class matrix {

private:
  int Nrows;
  int Ncolumns;
  T *data;

  void allocate();

public:

  matrix();

  matrix(int nr, int nc);

  // copy constructor
  matrix(const matrix <T> &A);

  // assignment operator
  matrix <T> operator= (const matrix <T> &A);

  void operator += (const matrix <T> &A);

  // assignment operator from double **array (assumes initialized)
  matrix <T> operator= (const T *A);

  // assignment operator (all values set to d)
  matrix <T> operator= (const double &d);

  ~matrix();

  int nrows() const;

  int ncolumns() const;

  T *c_array();

  void resize(int nr, int nc);

  void resize(int nr);

  matrix <T> transpose();

  // column major serial access - 1-indexed
  T operator[] (int r) const ;

  // column major serial access - 1-indexed
  T & operator[] (int r) ;

  T operator() (int r, int c) const ;

  T & operator() (int r, int c) ;

  T operator() (int r) const;

  T & operator() (int r) ;

  // permute
  matrix <T> operator[] (const matrix <int> &ind)  ;

  // in-place sort in ascending order
  void slow_sort(int index);

  // member function to randomize entries
  void randomize();

  // sort columns using user supplied comparison function
  void sort( int (*compare)(const void *, const void *) );

  void symeig(matrix <T> &d, matrix <T> &v);

  void eig(matrix <T> &WR, matrix <T> &WI, matrix <T> &VL, matrix <T> &VR);

  int byteCount() const;

  int entryCount() const;

  int size();

  int size() const;

  long double frobenius() const;

  T maxentry() const;

  T minentry() const;

  occa::memory occaCopy(occa::device &helper);

  matrix <T> inverse();

};

using namespace std;

template <class T>
ostream & operator << (ostream &os, matrix <T> & A);

template <class T>
matrix <T> operator* (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator+ (const matrix <T> & A, const matrix <T> &B);

template <class T>
matrix <T> operator- (const matrix <T> & A, const matrix <T> &B);

// general left matrix inverse not implemented
template <class T>
matrix <T> operator| (const matrix <T> & A, const matrix <T> &B);

// general left matrix inverse not implemented
//template <>
//matrix <double> operator| (const matrix <double> & A, const matrix <double> &B);

// general left matrix inverse not implemented
//template <>
//matrix <float> operator| (const matrix <float> & A, const matrix <float> &B);



// use #define to create a "compiler substitution rule"
#define imatrix matrix<int>
#define dmatrix matrix<double>

//int get_clarrayrequest();

template<class T>
void occadebug(const char *mess, matrix<T> &a, occa::memory &c_a);

template<class T>
void occadebugsum(const char *mess, matrix<T> &a, occa::memory &c_a);

template<class T>
void occadebug(const char *mess, matrix<T> &a, occa::memory &c_a, int rmin, int rmax, int cmin, int cmax);

#include "matrix.tpp"

#endif
