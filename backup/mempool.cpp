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

#include "mempool.hpp"

using std::stringstream;

memPool::memPool():
  _size{0}, _reserved{0} {}

memPool::memPool(occa::device device):
  _size{0}, _reserved{0} {
  _device = device;
}

memPool::~memPool() {
  _size = 0;
  _reserved = 0;
  if (data.size()) data.free();
}

size_t memPool::size() { return _size; };

size_t memPool::reserved() { return _reserved; };

occa::memory memPool::getMem(size_t offset) {
  return data + offset;
}

memPool::buffer memPool::reserve(size_t Nbytes) {
  const size_t offset = _reserved;
  if (Nbytes > _size-_reserved) {//requested space doesn't exist in pool
    //make a new memory buffer (separate from contiguous memory pool)
    occa::memory data_ = _device.malloc(Nbytes);

    //return a wrapper for this buffer
    memPool::buffer buf(this, data_, Nbytes, offset, false);
    return buf;
  } else {//requested space exists in pool
    //return a wrapper for the contiguous memory pool ptr + offset
    memPool::buffer buf(this, data, Nbytes, offset, true);
    _reserved += Nbytes;
    return buf;
  }
}

void memPool::release(size_t Nbytes) {
  //we assume this buffer is at the end of the pool
  if (_reserved<Nbytes) {
    stringstream ss;
    ss << "Cannot release buffer, invalid size";
    LIBP_ABORT(ss.str());
  }
  _reserved -= Nbytes;
}

void memPool::increase(size_t Nbytes) {
  //increase the size of the contiguous memory region
  const size_t newSize = _size + Nbytes;
  occa::memory newData = _device.malloc(newSize);

  //copy currently reserved data region into new pool
  if (_reserved)
    newData.copyFrom(data,_reserved);

  if (data.size()) data.free();
  data = newData;
  _size = newSize;
}

memPool::buffer::buffer():
  _pool{NULL}, _size{0}, _offset{0} {}

memPool::buffer::buffer(memPool* pool, occa::memory data, size_t size, size_t offset, bool inPool):
  _pool{pool}, _size{size}, _offset{offset}, _data{data}, _inPool{inPool} {}

size_t memPool::buffer::size() { return _size; };

size_t memPool::buffer::offset(){ return _offset; };

void memPool::buffer::release() {
  if (!_inPool) {//buffer is not from pool
    if (_pool) _pool->increase(_size);
  } else {
    if (_pool) _pool->release(_size);
  }
  _size = 0;
  _offset = 0;
}

occa::memory memPool::buffer::getMem() {
  if (!_inPool) {
    return _data;
  } else {
    if (_pool)
      return _pool->getMem(_offset);
    else {
      stringstream ss;
      ss << "Memory buffer not bound to a pool";
      LIBP_ABORT(ss.str());
    }
  }
}

memPool::buffer memPool::buffer::operator+(const size_t offset) const {
  memPool::buffer buf(this->_pool, this->_data, 0, 0);

  if (offset > _size) {
    stringstream ss;
    ss << "Invalid offset when shifting mempool buffer";
    LIBP_ABORT(ss.str());
  }
  buf._size = this->_size - offset;
  buf._offset = this->_offset + offset;
  buf._data = this->_data + offset;

  return buf;
}




pinnedMemPool::pinnedMemPool():
  _size{0}, _reserved{0} {}

pinnedMemPool::pinnedMemPool(occa::device device):
  _size{0}, _reserved{0} {
  _device = device;
}

pinnedMemPool::~pinnedMemPool() {
  _size = 0;
  _reserved = 0;
  if (data.size()) data.free();
}

size_t pinnedMemPool::size() { return _size; };

size_t pinnedMemPool::reserved() { return _reserved; };

void* pinnedMemPool::getMem(size_t offset) {
  return ((char*)data.ptr() + offset);
}

pinnedMemPool::buffer pinnedMemPool::reserve(size_t Nbytes) {
  size_t offset = _reserved;
  if (Nbytes) {
    if (_reserved + Nbytes > _size) {
      size_t newSize = _reserved + Nbytes;

      occa::properties props;
      props["mapped"] = true;
      occa::memory data_ = _device.malloc(newSize,props);

      if (_size) {
        memcpy(data_.ptr(), data.ptr(), _size);
        data.free();
      }

      data = data_;
      _size = newSize;
    }
  }
  _reserved += Nbytes;
  pinnedMemPool::buffer buf(this, Nbytes, offset);
  return buf;
}

void pinnedMemPool::release(size_t Nbytes) {
  //we assume this buffer is at the end of the pool
  if (_reserved<Nbytes) {
    stringstream ss;
    ss << "Cannot release buffer, invalid size";
    LIBP_ABORT(ss.str());
  }
  _reserved -= Nbytes;
}

pinnedMemPool::buffer::buffer():
  _pool{NULL}, _size{0}, _offset{0} {}

pinnedMemPool::buffer::buffer(pinnedMemPool* pool, size_t size, size_t offset):
  _pool{pool}, _size{size}, _offset{offset} {}

size_t pinnedMemPool::buffer::size() { return _size; };

size_t pinnedMemPool::buffer::offset(){ return _offset; };

void pinnedMemPool::buffer::release() {
  if (_pool)
    _pool->release(_size);
  _size = 0;
  _offset = 0;
}

void* pinnedMemPool::buffer::getMem() {
  if (_pool)
    return _pool->getMem(_offset);
  else {
    stringstream ss;
    ss << "Memory buffer not bound to a pool";
    LIBP_ABORT(ss.str());
  }
}

pinnedMemPool::buffer pinnedMemPool::buffer::operator+(const size_t offset) const {
  pinnedMemPool::buffer buf(this->_pool, 0, 0);

  if (offset > _size) {
    stringstream ss;
    ss << "Invalid offset when shifting mempool buffer";
    LIBP_ABORT(ss.str());
  }
  buf._size = this->_size - offset;
  buf._offset = this->_offset + offset;

  return buf;
}
