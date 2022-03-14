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

#ifndef LIBP_COMM_HPP
#define LIBP_COMM_HPP

#include <mpi.h>
#include "core.hpp"

namespace libp {

#define MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME

/*Generic data type*/
template<typename T>
struct mpiType {
  static MPI_Datatype getMpiType() {
    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(T), MPI_CHAR, &type);
    MPI_Type_commit(&type);
    return type;
  }
  static void freeMpiType(MPI_Datatype type) {
    MPI_Type_free(&type);
  }
  static constexpr bool isMpiType() { return false; }
};

/*Pre-defined MPI datatypes*/
#define TYPE(T, MPI_T)                               \
template<> struct mpiType<T> {                       \
  static MPI_Datatype getMpiType() { return MPI_T; } \
  static void freeMpiType(MPI_Datatype type) { }     \
  static constexpr bool isMpiType() { return true; } \
}

TYPE(char,   MPI_CHAR);
TYPE(int,    MPI_INT);
TYPE(long long int, MPI_LONG_LONG_INT);
TYPE(float,  MPI_FLOAT);
TYPE(double, MPI_DOUBLE);
#undef TYPE

/*Communicator class*/
class comm_t {

 private:
  std::shared_ptr<MPI_Comm> comm_ptr;
  int _rank=0;
  int _size=0;

 public:
  comm_t() = default;
  comm_t(const comm_t &c) = default;
  comm_t& operator = (const comm_t &c)=default;

  /*Static MPI_Init and MPI_Finalize*/
  static void Init(int &argc, char** &argv);
  static void Finalize();

  /*Static handle to MPI_COMM_WORLD*/
  static comm_t world();

  /*MPI_Comm_dup and MPI_Comm_delete*/
  comm_t Dup() const;
  comm_t Split(const int color, const int key) const;
  void Free();

  /*Rank and size getters*/
  const int rank() const;
  const int size() const;

  /*MPI_Comm getter*/
  MPI_Comm comm() const;

  using request_t = MPI_Request;

  /*Predefined ops*/
  using op_t = MPI_Op;
  static constexpr op_t Max  = MPI_MAX;
  static constexpr op_t Min  = MPI_MIN;
  static constexpr op_t Sum  = MPI_SUM;
  static constexpr op_t Prod = MPI_PROD;
  static constexpr op_t And  = MPI_LAND;
  static constexpr op_t Or   = MPI_LOR;
  static constexpr op_t Xor  = MPI_LXOR;

  /*libp::memory send*/
  template <template<typename> class mem, typename T>
  void Send(mem<T> m,
            const int dest,
            const int count=-1,
            const int tag=0) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    MPI_Send(m.ptr(), cnt, type, dest, tag, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory recv*/
  template <template<typename> class mem, typename T>
  void Recv(mem<T> m,
            const int source,
            const int count=-1,
            const int tag=0) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    MPI_Recv(m.ptr(), cnt, type, source, tag, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar send*/
  template <typename T>
  void Send(T& val,
            const int dest,
            const int tag=0) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Send(&val, 1, type, dest, tag, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar recv*/
  template <typename T>
  void Recv(T& val,
            const int source,
            const int tag=0) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Recv(&val, 1, type, source, tag, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory non-blocking send*/
  template <template<typename> class mem, typename T>
  void Isend(mem<T> m,
             const int dest,
             const int count,
             const int tag,
             request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Isend(m.ptr(), count, type, dest, tag, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory non-blocking recv*/
  template <template<typename> class mem, typename T>
  void Irecv(mem<T> m,
             const int source,
             const int count,
             const int tag,
             request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Irecv(m.ptr(), count, type, source, tag, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*scalar non-blocking send*/
  template <typename T>
  void Isend(T& val,
             const int dest,
             const int tag,
             request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Isend(&val, 1, type, dest, tag, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*scalar non-blocking recv*/
  template <typename T>
  void Irecv(T& val,
             const int source,
             const int tag,
             request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Irecv(&val, 1, type, source, tag, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory broadcast*/
  template <template<typename> class mem, typename T>
  void Bcast(mem<T> m,
             const int root,
             const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    MPI_Bcast(m.ptr(), cnt, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar broadcast*/
  template <typename T>
  void Bcast(T& val,
             const int root) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Bcast(&val, 1, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory reduce*/
  template <template<typename> class mem, typename T>
  void Reduce(const mem<T> snd,
                    mem<T> rcv,
              const int root,
              const op_t op = Sum,
              const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(snd.length()) : count;
    MPI_Reduce(snd.ptr(), rcv.ptr(), cnt, type, op, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory in-place reduce*/
  template <template<typename> class mem, typename T>
  void Reduce(mem<T> m,
              const int root,
              const op_t op = Sum,
              const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    if (_rank==root) {
      MPI_Reduce(MPI_IN_PLACE, m.ptr(), cnt, type, op, root, comm());
    } else {
      MPI_Reduce(m.ptr(), nullptr, cnt, type, op, root, comm());
    }
    mpiType<T>::freeMpiType(type);
  }

  /*scalar reduce*/
  template <typename T>
  void Reduce(const T& snd,
                    T& rcv,
              const int root,
              const op_t op = Sum) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Reduce(&snd, &rcv, 1, type, op, root, comm());
    mpiType<T>::freeMpiType(type);
  }
  template <typename T>
  void Reduce(T& val,
              const int root,
              const op_t op = Sum) const {
    T rcv=val;
    Reduce(val, rcv, root, op);
    if (rank()==root) val=rcv;
  }

  /*libp::memory allreduce*/
  template <template<typename> class mem, typename T>
  void Allreduce(const mem<T> snd,
                       mem<T> rcv,
                 const op_t op = Sum,
                 const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(snd.length()) : count;
    MPI_Allreduce(snd.ptr(), rcv.ptr(), cnt, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory in-place allreduce*/
  template <template<typename> class mem, typename T>
  void Allreduce(mem<T> m,
                 const op_t op = Sum,
                 const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    MPI_Allreduce(MPI_IN_PLACE, m.ptr(), cnt, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar allreduce*/
  template <typename T>
  void Allreduce(const T& snd,
                       T& rcv,
                 const op_t op = Sum) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Allreduce(&snd, &rcv, 1, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }
  template <typename T>
  void Allreduce(T& val,
                 const op_t op = Sum) const {
    T rcv=val;
    Allreduce(val, rcv, op);
    val = rcv;
  }

  /*libp::memory non-blocking allreduce*/
  template <template<typename> class mem, typename T>
  void Iallreduce(const mem<T> snd,
                        mem<T> rcv,
                  const op_t op,
                  const int count,
                  request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Iallreduce(snd.ptr(), rcv.ptr(), count, type, op, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory non-blocking in-place allreduce*/
  template <template<typename> class mem, typename T>
  void Iallreduce(mem<T> m,
                  const int root,
                  const op_t op,
                  const int count,
                  request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Iallreduce(MPI_IN_PLACE, m.ptr(), count, type, op, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*scalar non-blocking allreduce*/
  template <template<typename> class mem, typename T>
  void Iallreduce(const T& snd,
                        T& rcv,
                  const op_t op,
                  request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Iallreduce(&snd, &rcv, 1, type, op, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }
  /*scalar non-blocking in-place allreduce*/
  template <template<typename> class mem, typename T>
  void Iallreduce(T& val,
                  const op_t op,
                  request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Iallreduce(MPI_IN_PLACE, &val, 1, type, op, comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory scan*/
  template <template<typename> class mem, typename T>
  void Scan(const mem<T> snd,
                  mem<T> rcv,
            const op_t op = Sum,
            const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(snd.length()) : count;
    MPI_Scan(snd.ptr(), rcv.ptr(), cnt, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory in-place scan*/
  template <template<typename> class mem, typename T>
  void Scan(mem<T> m,
            const op_t op = Sum,
            const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(m.length()) : count;
    MPI_Scan(MPI_IN_PLACE, m.ptr(), cnt, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar scan*/
  template <typename T>
  void Scan(const T& snd,
                  T& rcv,
            const op_t op = Sum) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Scan(&snd, &rcv, 1, type, op, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory gather*/
  template <template<typename> class mem, typename T>
  void Gather(const mem<T> snd,
                    mem<T> rcv,
              const int root,
              const int sendCount=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (sendCount==-1) ? static_cast<int>(snd.length()) : sendCount;
    MPI_Gather(snd.ptr(), cnt, type,
               rcv.ptr(), cnt, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory gatherv*/
  template <template<typename> class mem, typename T>
  void Gatherv(const mem<T> snd,
               const int sendcount,
                     mem<T> rcv,
               const memory<int> recvCounts,
               const memory<int> recvOffsets,
               const int root) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Gatherv(snd.ptr(), sendcount, type,
                rcv.ptr(), recvCounts.ptr(), recvOffsets.ptr(), type,
                root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar gather*/
  template <template<typename> class mem, typename T>
  void Gather(const T& snd,
                    mem<T> rcv,
              const int root) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Gather(&snd,      1, type,
               rcv.ptr(), 1, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory scatter*/
  template <template<typename> class mem, typename T>
  void Scatter(const mem<T> snd,
                     mem<T> rcv,
               const int root,
               const int count=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (count==-1) ? static_cast<int>(rcv.length()) : count;
    MPI_Scatter(snd.ptr(), cnt, type,
                rcv.ptr(), cnt, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory scatterv*/
  template <template<typename> class mem, typename T>
  void Scatterv(const mem<T> snd,
                const memory<int> sendCounts,
                const memory<int> sendOffsets,
                      mem<T> rcv,
                const int recvcount,
                const int root) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Scatterv(snd.ptr(), sendCounts.ptr(), sendOffsets.ptr(), type,
                 rcv.ptr(), recvcount, type,
                 root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar scatter*/
  template <template<typename> class mem, typename T>
  void Scatter(T& rcv,
               const mem<T> snd,
               const int root) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Scatter(snd.ptr,   1, type,
                &rcv,      1, type, root, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory allgather*/
  template <template<typename> class mem, typename T>
  void Allgather(const mem<T> snd,
                       mem<T> rcv,
                 const int sendCount=-1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    const int cnt = (sendCount==-1) ? static_cast<int>(snd.length()) : sendCount;
    MPI_Allgather(snd.ptr(), cnt, type,
                  rcv.ptr(), cnt, type, comm());
    mpiType<T>::freeMpiType(type);
  }
  template <template<typename> class mem, typename T>
  void Allgather(mem<T> m,
                 const int cnt) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Allgather(MPI_IN_PLACE, cnt, type,
                  m.ptr(),      cnt, type, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory allgatherv*/
  template <template<typename> class mem, typename T>
  void Allgatherv(const mem<T> snd,
                  const int sendcount,
                        mem<T> rcv,
                  const memory<int> recvCounts,
                  const memory<int> recvOffsets) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Allgatherv(snd.ptr(), sendcount, type,
                   rcv.ptr(), recvCounts.ptr(), recvOffsets.ptr(), type,
                   comm());
    mpiType<T>::freeMpiType(type);
  }

  /*scalar allgather*/
  template <template<typename> class mem, typename T>
  void Allgather(const T& snd,
                       mem<T> rcv) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Allgather(&snd,      1, type,
                  rcv.ptr(), 1, type, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory alltoall*/
  template <template<typename> class mem, typename T>
  void Alltoall(const mem<T> snd,
                      mem<T> rcv,
                const int cnt=1) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Alltoall(snd.ptr(), cnt, type,
                 rcv.ptr(), cnt, type, comm());
    mpiType<T>::freeMpiType(type);
  }

  /*libp::memory alltoallv*/
  template <template<typename> class mem, typename T>
  void Alltoallv(const mem<T> snd,
                 const memory<int> sendCounts,
                 const memory<int> sendOffsets,
                       mem<T> rcv,
                 const memory<int> recvCounts,
                 const memory<int> recvOffsets) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Alltoallv(snd.ptr(), sendCounts.ptr(), sendOffsets.ptr(), type,
                  rcv.ptr(), recvCounts.ptr(), recvOffsets.ptr(), type,
                  comm());
    mpiType<T>::freeMpiType(type);
  }

  template <template<typename> class mem, typename T>
  void Ialltoallv(const mem<T> snd,
                  const memory<int> sendCounts,
                  const memory<int> sendOffsets,
                        mem<T> rcv,
                  const memory<int> recvCounts,
                  const memory<int> recvOffsets,
                  request_t &request) const {
    MPI_Datatype type = mpiType<T>::getMpiType();
    MPI_Ialltoallv(snd.ptr(), sendCounts.ptr(), sendOffsets.ptr(), type,
                  rcv.ptr(), recvCounts.ptr(), recvOffsets.ptr(), type,
                  comm(), &request);
    mpiType<T>::freeMpiType(type);
  }

  void Wait(request_t &request) const;
  void Waitall(const int count, memory<request_t> &requests) const;
  void Barrier() const;

  static void GetProcessorName(char* name, int &namelen) {
    MPI_Get_processor_name(name,&namelen);
  }
};

} //namespace libp

#endif
