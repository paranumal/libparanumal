/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifndef MEMPOOL_HPP
#define MEMPOOL_HPP

#include <string>
#include <vector>
#include <ostream>
#include <sstream>
#include "utils.h"

#include <occa.hpp>

class memPool {
private:
  occa::device _device;
  occa::memory data;

  size_t _size;
  size_t _reserved;

public:
  memPool();
  memPool(occa::device device);
  ~memPool();

  size_t size();
  size_t reserved();

  occa::memory getMem(size_t offset);

  class buffer {
  private:
    memPool* _pool;

    occa::memory _data;
    bool _inPool;

    size_t _size;
    size_t _offset;

  public:
    buffer();
    buffer(memPool* pool, occa::memory data, size_t size=0, size_t offset=0, bool inPool=true);
    ~buffer() = default;

    size_t size();
    size_t offset();

    void reserve(memPool* pool, size_t size_);
    void release();

    occa::memory getMem();

    buffer operator+(const size_t offset_) const;
  };

  buffer reserve(size_t Nbytes);
  void release(size_t Nbytes);
  void increase(size_t Nbytes);

};

class pinnedMemPool {
private:
  occa::device _device;
  occa::memory data;

  size_t _size;
  size_t _reserved;

public:
  pinnedMemPool();
  pinnedMemPool(occa::device device);
  ~pinnedMemPool();

  size_t size();
  size_t reserved();

  void* getMem(size_t offset);

  class buffer {
  private:
    pinnedMemPool* _pool;

    size_t _size;
    size_t _offset;

  public:
    buffer();
    buffer(pinnedMemPool* pool, size_t size=0, size_t offset=0);
    ~buffer() = default;

    size_t size();
    size_t offset();

    void reserve(pinnedMemPool* pool, size_t size_);
    void release();

    void* getMem();

    buffer operator+(const size_t offset_) const;
  };

  buffer reserve(size_t Nbytes);
  void release(size_t Nbytes);

};

#endif
