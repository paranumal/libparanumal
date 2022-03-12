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

#ifndef LIBP_MEMORY_HPP
#define LIBP_MEMORY_HPP

#include "utils.hpp"

namespace libp {

template<typename T>
class memory {
  template <typename U> friend class memory;

 private:
  using size_t = std::size_t;
  using ptrdiff_t = std::ptrdiff_t;

  std::shared_ptr<T[]> shrdPtr;
  size_t lngth;
  size_t offset;

 public:
  memory() :
    lngth{0},
    offset{0} {}

  memory(const size_t lngth_) :
    shrdPtr(new T[lngth_]),
    lngth{lngth_},
    offset{0} {}

  memory(const size_t lngth_,
         const T val) :
    shrdPtr(new T[lngth_]),
    lngth{lngth_},
    offset{0} {
    #pragma omp parallel for
    for (size_t i=0;i<lngth;++i) {
      shrdPtr[i] = val;
    }
  }

  /*Conversion constructor*/
  template<typename U>
  memory(const memory<U> &m):
    shrdPtr{std::reinterpret_pointer_cast<T[]>(m.shrdPtr)},
    lngth{m.lngth*sizeof(T)/sizeof(U)},
    offset{m.offset*sizeof(T)/sizeof(U)} {
    // Check that this conversion made sense
    LIBP_ABORT("libp::memory type conversion failed. Trying to convert "
                << m.lngth << " " << sizeof(T) << "-byte words to "
                << lngth << " " << sizeof(U) << "-byte words.",
                lngth*sizeof(U) != m.lngth*sizeof(T));

    LIBP_ABORT("libp::memory type conversion failed. Source memory has offset at "
               << m.lngth << " " << sizeof(T) << "-byte words, destination memory would have offset at"
               << lngth << " " << sizeof(U) << "-byte words.",
               offset*sizeof(U) != m.offset*sizeof(T));
  }

  memory(const memory<T> &m)=default;
  memory& operator = (const memory<T> &m)=default;
  ~memory()=default;

  void malloc(const size_t lngth_) {
    *this = memory<T>(lngth_);
  }

  void calloc(const size_t lngth_) {
    *this = memory<T>(lngth_, T{0});
  }

  void realloc(const size_t lngth_) {
    memory<T> m(lngth_);
    const ptrdiff_t cnt = std::min(lngth, lngth_);
    m.copyFrom(*this, cnt);
    *this = m;
  }

  memory& swap(memory<T> &m) {
    std::swap(shrdPtr, m.shrdPtr);
    std::swap(lngth, m.lngth);
    std::swap(offset, m.offset);
    return *this;
  }

  T* ptr() {
    return shrdPtr.get()+offset;
  }
  const T* ptr() const {
    return shrdPtr.get()+offset;
  }

  T* begin() {return ptr();}
  T* end() {return ptr() + length();}

  size_t length() const {
    return lngth;
  }

  size_t size() const {
    return lngth*sizeof(T);
  }

  size_t use_count() const {
    return shrdPtr.use_count();
  }

  T& operator[](const ptrdiff_t idx) const {
    return shrdPtr[idx+offset];
  }

  bool operator == (const memory<T> &other) const {
    return (shrdPtr==other.shrdPtr && offset==other.offset);
  }
  bool operator != (const memory<T> &other) const {
    return (shrdPtr!=other.shrdPtr || offset!=other.offset);
  }

  memory<T> operator + (const ptrdiff_t offset_) const {
    return slice(offset_);
  }
  memory<T>& operator += (const ptrdiff_t offset_) {
    *this = slice(offset_);
    return *this;
  }

  memory<T> slice(const ptrdiff_t offset_,
                  const ptrdiff_t count = -1) const {
    memory<T> m(*this);
    m.offset = offset + offset_;
    m.lngth = (count==-1)
                ? (lngth - offset_)
                : count;
    return m;
  }

  /*Copy from raw ptr*/
  void copyFrom(const T* src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) {

    const ptrdiff_t cnt = (count==-1) ? lngth : count;

    LIBP_ABORT("libp::memory::copyFrom Cannot have negative count ("
               << cnt << ")",
               cnt < 0);
    LIBP_ABORT("libp::memory::copyFrom Cannot have negative offset ("
               << offset_ << ")",
               offset_ < 0);
    LIBP_ABORT("libp::memory::copyFrom Destination memory has size [" << lngth << "],"
               << " trying to access [" << offset_ << ", " << offset_+static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt)+offset_ > lngth);

    std::copy(src,
              src+cnt,
              ptr()+offset_);
  }

  /*Copy from memory*/
  void copyFrom(const memory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset_ = 0) {
    const ptrdiff_t cnt = (count==-1) ? lngth : count;

    LIBP_ABORT("libp::memory::copyFrom Cannot have negative count ("
               << cnt << ")",
               cnt < 0);
    LIBP_ABORT("libp::memory::copyFrom Cannot have negative offset ("
               << offset_ << ")",
               offset_ < 0);
    LIBP_ABORT("libp::memory::copyFrom Source memory has size [" << src.length() << "],"
               << " trying to access [0, " << static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt) > src.length());
    LIBP_ABORT("libp::memory::copyFrom Destination memory has size [" << lngth << "],"
               << " trying to access [" << offset_ << ", " << offset_+static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt)+offset_ > lngth);

    std::copy(src.ptr(),
              src.ptr()+cnt,
              ptr()+offset_);
  }

  /*Copy to raw pointer*/
  void copyTo(T *dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset_ = 0) const {
    const ptrdiff_t cnt = (count==-1) ? lngth : count;

    LIBP_ABORT("libp::memory::copyTo Cannot have negative count ("
               << cnt << ")",
               cnt < 0);
    LIBP_ABORT("libp::memory::copyTo Cannot have negative offset ("
               << offset_ << ")",
               offset_ < 0);
    LIBP_ABORT("libp::memory::copyTo Source memory has size [" << lngth << "],"
               << " trying to access [" << offset_ << ", " << offset_+static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt)+offset_ > lngth);

    std::copy(ptr()+offset_,
              ptr()+offset_+cnt,
              dest);
  }

  /*Copy to memory*/
  void copyTo(memory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset_ = 0) const {
    const ptrdiff_t cnt = (count==-1) ? lngth : count;

    LIBP_ABORT("libp::memory::copyTo Cannot have negative count ("
               << cnt << ")",
               cnt < 0);
    LIBP_ABORT("libp::memory::copyTo Cannot have negative offset ("
               << offset_ << ")",
               offset_ < 0);
    LIBP_ABORT("libp::memory::copyTo Destination memory has size [" << dest.length() << "],"
               << " trying to access [0, " << cnt << "]",
               static_cast<size_t>(cnt) > dest.length());
    LIBP_ABORT("libp::memory::copyTo Source memory has size [" << lngth << "],"
               << " trying to access [" << offset_ << ", " << offset_+static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt)+offset_ > lngth);

    std::copy(ptr()+offset_,
              ptr()+offset_+cnt,
              dest.ptr());
  }

  memory<T> clone() const {
    memory<T> m(lngth);
    m.copyFrom(*this);
    return m;
  }

  void free() {
    shrdPtr = nullptr;
    lngth=0;
    offset=0;
  }
};

template <typename T>
std::ostream& operator << (std::ostream &out,
                         const memory<T> &memory) {
  out << "memory - "
      << "type: " << typeid(T).name() << ", "
      << "ptr : " << memory.ptr() << ", "
      << "length : " << memory.length() << ", "
      << "use_count : " << memory.use_count();
  return out;
}

/*Extern declare common instantiations for faster compilation*/
extern template class memory<int>;
extern template class memory<long long int>;
extern template class memory<float>;
extern template class memory<double>;

/*libp::deviceMemory is a wrapper around occa::memory*/
template<typename T>
class deviceMemory: public occa::memory {
 public:
  deviceMemory() = default;
  deviceMemory(const deviceMemory<T> &m)=default;
  deviceMemory(occa::memory m):
    occa::memory(m)
  {
    if (isInitialized())
      occa::memory::setDtype(occa::dtype::get<T>());
  }

  /*Conversion constructor*/
  template<typename U>
  deviceMemory(const deviceMemory<U> &m):
    occa::memory(m)
  {
    if (isInitialized())
      occa::memory::setDtype(occa::dtype::get<T>());
  }

  deviceMemory<T>& operator = (const deviceMemory<T> &m)=default;
  ~deviceMemory()=default;

  T* ptr() {
    return static_cast<T*>(occa::memory::ptr());
  }
  const T* ptr() const {
    return static_cast<const T*>(occa::memory::ptr());
  }

  T& operator[](const ptrdiff_t idx) {
    return ptr()[idx];
  }

  deviceMemory<T> operator + (const ptrdiff_t offset) const {
    if (isInitialized())
      return deviceMemory<T>(occa::memory::operator+(offset));
    else
      return deviceMemory<T>();
  }

  deviceMemory<T>& operator += (const ptrdiff_t offset) {
    *this = deviceMemory<T>(occa::memory::slice(offset));
    return *this;
  }

  /*Copy from libp::memory*/
  void copyFrom(const libp::memory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    LIBP_ABORT("libp::memory::copyFrom Source memory has size [" << src.length() << "],"
               << " trying to access [0, " << static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt) > src.length());

    occa::memory::copyFrom(src.ptr(),
                           cnt*sizeof(T),
                           offset*sizeof(T),
                           props);
  }

  void copyFrom(const libp::memory<T> src,
                const properties_t &props) {

    if (length()==0) return;

    LIBP_ABORT("libp::memory::copyFrom Source memory has size [" << src.length() << "],"
               << " trying to access [0, " << length() << "]",
               length() > src.length());

    occa::memory::copyFrom(src.ptr(),
                           length()*sizeof(T),
                           0,
                           props);
  }

  /*Copy from libp::deviceMemory*/
  void copyFrom(const deviceMemory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyFrom(src,
                           cnt*sizeof(T),
                           offset*sizeof(T),
                           0,
                           props);
  }

  void copyFrom(const deviceMemory<T> src,
                const properties_t &props) {

    if (length()==0) return;

    occa::memory::copyFrom(src,
                           length()*sizeof(T),
                           0,
                           0,
                           props);
  }

  /*Copy to libp::memory*/
  void copyTo(libp::memory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset = 0,
              const properties_t &props = properties_t()) const {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    LIBP_ABORT("libp::memory::copyTo Destination memory has size [" << dest.length() << "],"
               << " trying to access [0, " << static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt) > dest.length());

    occa::memory::copyTo(dest.ptr(),
                         cnt*sizeof(T),
                         offset*sizeof(T),
                         props);
  }

  void copyTo(libp::memory<T> dest,
              const properties_t &props) const {

    if (length()==0) return;

    LIBP_ABORT("libp::memory::copyTo Destination memory has size [" << dest.length() << "],"
               << " trying to access [0, " << length() << "]",
               length() > dest.length());

    occa::memory::copyTo(dest.ptr(),
                         length()*sizeof(T),
                         0,
                         props);
  }

  /*Copy to libp::deviceMemory*/
  void copyTo(deviceMemory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset = 0,
              const properties_t &props = properties_t()) const {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyTo(dest,
                         cnt*sizeof(T),
                         0,
                         offset*sizeof(T),
                         props);
  }

  void copyTo(deviceMemory<T> dest,
              const properties_t &props) const {

    if (length()==0) return;

    occa::memory::copyTo(dest,
                         length()*sizeof(T),
                         0,
                         0,
                         props);
  }
};

/*Extern declare common instantiations for faster compilation*/
extern template class deviceMemory<int>;
extern template class deviceMemory<long long int>;
extern template class deviceMemory<float>;
extern template class deviceMemory<double>;

/*libp::pinnedMemory is another wrapper around occa::memory,
  but is allocated slightly differently*/
template<typename T>
class pinnedMemory: public occa::memory {
 public:
  pinnedMemory() = default;
  pinnedMemory(const pinnedMemory<T> &m)=default;
  pinnedMemory(occa::memory m):
    occa::memory(m)
  {
    if (isInitialized())
      occa::memory::setDtype(occa::dtype::get<T>());
  };

  /*Conversion constructor*/
  template<typename U>
  pinnedMemory(const pinnedMemory<U> &m):
    occa::memory(m)
  {
    if (isInitialized())
      occa::memory::setDtype(occa::dtype::get<T>());
  }

  pinnedMemory<T>& operator = (const pinnedMemory<T> &m)=default;
  ~pinnedMemory()=default;

  T* ptr() {
    return static_cast<T*>(occa::memory::ptr());
  }
  const T* ptr() const {
    return static_cast<const T*>(occa::memory::ptr());
  }

  T& operator[](const ptrdiff_t idx) {
    return ptr()[idx];
  }

  pinnedMemory<T> operator + (const ptrdiff_t offset) const {
    if (isInitialized())
      return pinnedMemory<T>(occa::memory::operator+(offset));
    else
      return pinnedMemory<T>();
  }

  pinnedMemory<T>& operator += (const ptrdiff_t offset) {
    *this = pinnedMemory<T>(occa::memory::slice(offset));
    return *this;
  }

  /*Copy from libp::memory*/
  void copyFrom(const libp::memory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    LIBP_ABORT("libp::memory::copyFrom Source memory has size [" << src.length() << "],"
               << " trying to access [0, " << static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt) > src.length());

    occa::memory::copyFrom(src.ptr(),
                           cnt*sizeof(T),
                           offset*sizeof(T),
                           props);
  }

  void copyFrom(const libp::memory<T> src,
                const properties_t &props) {

    if (length()==0) return;

    LIBP_ABORT("libp::memory::copyFrom Source memory has size [" << src.length() << "],"
               << " trying to access [0, " << length() << "]",
               length() > src.length());

    occa::memory::copyFrom(src.ptr(),
                           length()*sizeof(T),
                           0,
                           props);
  }

  /*Copy from libp::deviceMemory*/
  void copyFrom(const deviceMemory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyFrom(src,
                           cnt*sizeof(T),
                           offset*sizeof(T),
                           0,
                           props);
  }

  void copyFrom(const deviceMemory<T> src,
                const properties_t &props) {

    if (length()==0) return;

    occa::memory::copyFrom(src,
                           length()*sizeof(T),
                           0,
                           0,
                           props);
  }

  /*Copy from libp::pinnedMemory*/
  void copyFrom(const pinnedMemory<T> src,
                const ptrdiff_t count = -1,
                const ptrdiff_t offset = 0,
                const properties_t &props = properties_t()) {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyFrom(src,
                           cnt*sizeof(T),
                           offset*sizeof(T),
                           0,
                           props);
  }

  void copyFrom(const pinnedMemory<T> src,
                const properties_t &props) {

    if (length()==0) return;

    occa::memory::copyFrom(src,
                           length()*sizeof(T),
                           0,
                           0,
                           props);
  }

  /*Copy to libp::memory*/
  void copyTo(libp::memory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset = 0,
              const properties_t &props = properties_t()) const {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    LIBP_ABORT("libp::memory::copyTo Destination memory has size [" << dest.length() << "],"
               << " trying to access [0, " << static_cast<size_t>(cnt) << "]",
               static_cast<size_t>(cnt) > dest.length());

    occa::memory::copyTo(dest.ptr(),
                         cnt*sizeof(T),
                         offset*sizeof(T),
                         props);
  }

  void copyTo(libp::memory<T> dest,
              const properties_t &props) const {

    if (length()==0) return;

    LIBP_ABORT("libp::memory::copyTo Destination memory has size [" << dest.length() << "],"
               << " trying to access [0, " << length() << "]",
               length() > dest.length());

    occa::memory::copyTo(dest.ptr(),
                         length()*sizeof(T),
                         0,
                         props);
  }

  /*Copy to libp::deviceMemory*/
  void copyTo(deviceMemory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset = 0,
              const properties_t &props = properties_t()) const {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyTo(dest,
                         cnt*sizeof(T),
                         0,
                         offset*sizeof(T),
                         props);
  }

  void copyTo(deviceMemory<T> dest,
              const properties_t &props) const {

    if (length()==0) return;

    occa::memory::copyTo(dest,
                         length()*sizeof(T),
                         0,
                         0,
                         props);
  }

  /*Copy to libp::pinnedMemory*/
  void copyTo(pinnedMemory<T> dest,
              const ptrdiff_t count = -1,
              const ptrdiff_t offset = 0,
              const properties_t &props = properties_t()) const {
    const ptrdiff_t cnt = (count==-1) ? length() : count;

    if (cnt==0) return;

    occa::memory::copyTo(dest,
                         cnt*sizeof(T),
                         0,
                         offset*sizeof(T),
                         props);
  }

  void copyTo(pinnedMemory<T> dest,
              const properties_t &props) const {

    if (length()==0) return;

    occa::memory::copyTo(dest,
                         length()*sizeof(T),
                         0,
                         0,
                         props);
  }
};

} //namespace libp

#endif
