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

#include "comm.hpp"

namespace libp {

/*Static MPI_Init and MPI_Finalize*/
void comm_t::Init(int &argc, char** &argv) { MPI_Init(&argc, &argv); }
void comm_t::Finalize() { MPI_Finalize(); }

/*Static handle to MPI_COMM_WORLD*/
comm_t comm_t::world() {
  comm_t c;
  c.comm_ptr = std::make_shared<MPI_Comm>();
  *(c.comm_ptr) = MPI_COMM_WORLD;
  MPI_Comm_rank(c.comm(), &(c._rank));
  MPI_Comm_size(c.comm(), &(c._size));
  return c;
}

/*MPI_Comm_dup and free*/
comm_t comm_t::Dup() const {
  comm_t c;
  /*Make a new comm shared_ptr, which will call MPI_Comm_free when destroyed*/
  c.comm_ptr = std::shared_ptr<MPI_Comm>(new MPI_Comm,
                                        [](MPI_Comm *comm) {
                                          if (*comm != MPI_COMM_NULL)
                                            MPI_Comm_free(comm);
                                          delete comm;
                                        });
  MPI_Comm_dup(comm(), c.comm_ptr.get());
  MPI_Comm_rank(c.comm(), &(c._rank));
  MPI_Comm_size(c.comm(), &(c._size));
  c.setGpuAware(gpuAware());
  return c;
}
void comm_t::Free() {
  comm_ptr = nullptr;
  _rank=0;
  _size=0;
}
/*Split*/
comm_t comm_t::Split(const int color, const int key) const {
  comm_t c;
  /*Make a new comm shared_ptr, which will call MPI_Comm_free when destroyed*/
  c.comm_ptr = std::shared_ptr<MPI_Comm>(new MPI_Comm,
                                        [](MPI_Comm *comm) {
                                          if (*comm != MPI_COMM_NULL)
                                            MPI_Comm_free(comm);
                                          delete comm;
                                        });

  MPI_Comm_split(comm(), color, key, c.comm_ptr.get());
  MPI_Comm_rank(c.comm(), &(c._rank));
  MPI_Comm_size(c.comm(), &(c._size));
  return c;
}

/*Rank and size getters*/
const int comm_t::rank() const {
  return _rank;
}
const int comm_t::size() const {
  return _size;
}

MPI_Comm comm_t::comm() const {
  if (comm_ptr == nullptr) {
    return MPI_COMM_NULL;
  } else {
    return *comm_ptr;
  }
}

/*GPU-aware setter*/
void comm_t::setGpuAware(const bool flag) {
  _gpuAware = flag;
}

/*GPU-aware getter*/
const bool comm_t::gpuAware() const {
  return _gpuAware;
}

void comm_t::Wait(request_t &request) const {
  MPI_Wait(&request, MPI_STATUS_IGNORE);
}

void comm_t::Waitall(const int count, memory<request_t> &requests) const {
  MPI_Waitall(count, requests.ptr(), MPI_STATUSES_IGNORE);
}

void comm_t::Waitall(const int count, request_t* requests) const {
  MPI_Waitall(count, requests, MPI_STATUSES_IGNORE);
}

void comm_t::Barrier() const {
  MPI_Barrier(comm());
}

} //namespace libp
