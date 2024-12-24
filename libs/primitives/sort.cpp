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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
using __gnu_parallel::stable_sort;
#else
using std::sort;
using std::stable_sort;
#endif

namespace libp {

namespace prim {


template <typename T>
class ZipRef {
  std::pair<T*, dlong*> ptr;

 public:
  ZipRef() = delete;
  ZipRef(const ZipRef& z) = default;
  ZipRef(ZipRef&& z) = default;
  ZipRef(T* const val, dlong* const id): ptr{val, id} {}

  ZipRef& operator=(const ZipRef& z)             { return copy_assign(z); }
  ZipRef& operator=(const std::pair<T, dlong>& val) { return val_assign(val); }

  ZipRef& copy_assign(const ZipRef& z) {
    *(std::get<0>(ptr)) = *(std::get<0>(z.ptr));
    *(std::get<1>(ptr)) = *(std::get<1>(z.ptr));
    return *this;
  }

  ZipRef& val_assign(const std::pair<T, dlong>& t) {
    *(std::get<0>(ptr)) = std::get<0>(t);
    *(std::get<1>(ptr)) = std::get<1>(t);
    return *this;
  }

  T val() const {return *(std::get<0>(ptr));}

  operator std::pair<T, dlong>() const { return std::pair(*(std::get<0>(ptr)), *(std::get<1>(ptr))); }

  void swap(const ZipRef& o) const {
    std::swap(*std::get<0>(ptr), *std::get<0>(o.ptr));
    std::swap(*std::get<1>(ptr), *std::get<1>(o.ptr));
  }

  /*Comparisons*/
  bool operator==(const ZipRef &o) const { return val() == o.val(); }
  inline friend bool operator==(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() == std::get<0>(t); }
  inline friend bool operator==(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) == r.val(); }
  bool operator<=(const ZipRef &o) const { return val() <= o.val(); }
  inline friend bool operator<=(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() <= std::get<0>(t); }
  inline friend bool operator<=(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) <= r.val(); }
  bool operator>=(const ZipRef &o) const { return val() >= o.val(); }
  inline friend bool operator>=(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() >= std::get<0>(t); }
  inline friend bool operator>=(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) >= r.val(); }
  bool operator!=(const ZipRef &o) const { return val() != o.val(); }
  inline friend bool operator!=(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() != std::get<0>(t); }
  inline friend bool operator!=(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) != r.val(); }
  bool operator< (const ZipRef &o) const { return val() <  o.val(); }
  inline friend bool operator<(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() <  std::get<0>(t); }
  inline friend bool operator<(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) <  r.val(); }
  bool operator> (const ZipRef &o) const { return val() >  o.val(); }
  inline friend bool operator>(const ZipRef& r, const std::pair<T, dlong>& t) { return r.val() >  std::get<0>(t); }
  inline friend bool operator>(const std::pair<T, dlong>& t, const ZipRef& r) { return std::get<0>(t) >  r.val(); }
};

template<typename T>
class ZipIter {
  std::pair<T*, dlong*> it;

 public:
  using iterator_category = typename std::iterator_traits<T*>::iterator_category;
  using difference_type   = typename std::iterator_traits<T*>::difference_type;
  using value_type        = std::pair<typename std::iterator_traits<T*>::value_type, dlong>;
  using pointer           = std::pair<typename std::iterator_traits<T*>::pointer, dlong*>;
  using reference         = ZipRef<typename std::iterator_traits<T*>::value_type>;

  ZipIter() = default;
  ZipIter(const ZipIter &rhs) = default;
  ZipIter(ZipIter&& rhs) = default;
  ZipIter(T* val, dlong* id): it(std::move(val), std::move(id)) {}

  ZipIter& operator=(const ZipIter& rhs) = default;
  ZipIter& operator=(ZipIter&& rhs) = default;

  ZipIter& operator+=(const difference_type d) { std::get<0>(it)+=d; std::get<1>(it)+=d; return *this; }
  ZipIter& operator-=(const difference_type d) { std::get<0>(it)-=d; std::get<1>(it)-=d; return *this; }

  reference operator* () const { return reference(std::get<0>(it), std::get<1>(it)); }
  pointer   operator->() const { return pointer(std::get<0>(it), std::get<1>(it)); }
  reference operator[](difference_type rhs) const {return *(operator+(rhs));}

  ZipIter& operator++() { return operator+=(1); }
  ZipIter& operator--() { return operator+=(-1); }
  ZipIter operator++(dlong) {ZipIter tmp(*this); operator++(); return tmp;}
  ZipIter operator--(dlong) {ZipIter tmp(*this); operator--(); return tmp;}

  difference_type operator-(const ZipIter& rhs) const {return std::get<0>(it)-std::get<0>(rhs.it);}
  ZipIter operator+(const difference_type d) const {ZipIter tmp(*this); tmp += d; return tmp;}
  ZipIter operator-(const difference_type d) const {ZipIter tmp(*this); tmp -= d; return tmp;}
  inline friend ZipIter operator+(const difference_type d, const ZipIter& z) {return z+d;}
  inline friend ZipIter operator-(const difference_type d, const ZipIter& z) {return z-d;}

  bool operator==(const ZipIter& rhs) const {return it == rhs.it;}
  bool operator<=(const ZipIter& rhs) const {return it <= rhs.it;}
  bool operator>=(const ZipIter& rhs) const {return it >= rhs.it;}
  bool operator!=(const ZipIter& rhs) const {return it != rhs.it;}
  bool operator< (const ZipIter& rhs) const {return it <  rhs.it;}
  bool operator> (const ZipIter& rhs) const {return it >  rhs.it;}
};

template<typename T>
class Zip {
  std::pair<T*, dlong*> zip;
  size_t len;

 public:
  Zip() = delete;
  Zip(const Zip& z) = default;
  Zip(Zip&& z) = default;
  Zip(memory<T>& val, memory<dlong>& ids, size_t len_) : zip {val.ptr(), ids.ptr()}, len(len_) {}

  ZipIter<T> begin() {return ZipIter(std::get<0>(zip), std::get<1>(zip));}
  ZipIter<T> end() {return ZipIter(std::get<0>(zip)+len, std::get<1>(zip)+len);}
};

template<typename T>
void swap(const ZipRef<T>& a, const ZipRef<T>& b) {
  a.swap(b);
}


template<typename T>
void sort(const dlong N, memory<T> v) {
	if (N<=0) return;
	::sort(v.ptr(), v.ptr() + N);
}

template void sort(const dlong N, memory<int> v);
template void sort(const dlong N, memory<long long int> v);
template void sort(const dlong N, memory<float> v);
template void sort(const dlong N, memory<double> v);

template<typename T>
void sort(const dlong N, memory<T> v, memory<dlong> sortIds) {

	if (N<=0) return;

	/*Write initial ordering*/
	#pragma omp parallel for
	for (dlong n=0; n<N; ++n) { sortIds[n] = n; }

	/*Zip vals and ids together*/
	Zip<T> zip(v, sortIds, N);

	/*Simultaneous sort*/
	::sort(zip.begin(), zip.end());
}

template void sort(const dlong N, memory<int> v          , memory<dlong> sortIds);
template void sort(const dlong N, memory<long long int> v, memory<dlong> sortIds);
template void sort(const dlong N, memory<float> v        , memory<dlong> sortIds);
template void sort(const dlong N, memory<double> v       , memory<dlong> sortIds);

template<typename T>
void stableSort(const dlong N, memory<T> v) {
  if (N<=0) return;
  ::stable_sort(v.ptr(), v.ptr() + N);
}

template void stableSort(const dlong N, memory<int> v);
template void stableSort(const dlong N, memory<long long int> v);
template void stableSort(const dlong N, memory<float> v);
template void stableSort(const dlong N, memory<double> v);

template<typename T>
void stableSort(const dlong N, memory<T> v, memory<dlong> sortIds) {

  if (N<=0) return;

  /*Write initial ordering*/
  #pragma omp parallel for
  for (dlong n=0; n<N; ++n) { sortIds[n] = n; }

  /*Zip vals and ids together*/
  Zip<T> zip(v, sortIds, N);

  /*Simultaneous sort*/
  ::stable_sort(zip.begin(), zip.end());
}

template void stableSort(const dlong N, memory<int> v          , memory<dlong> sortIds);
template void stableSort(const dlong N, memory<long long int> v, memory<dlong> sortIds);
template void stableSort(const dlong N, memory<float> v        , memory<dlong> sortIds);
template void stableSort(const dlong N, memory<double> v       , memory<dlong> sortIds);

} //namespace prim

} //namespace libp
