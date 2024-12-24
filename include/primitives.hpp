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

#ifndef PRIMITIVES_HPP
#define PRIMITIVES_HPP

#include "core.hpp"
#include "memory.hpp"

namespace libp {

namespace prim {

/// min: Parallel primitive returning minimum value in input array
///
/// Equivalent to the following code
/// min_v = std::numeric_limits<T>::max();
/// for(std::size_t i = 0; i < N; ++i) {
///   min_v = std::min(min_v, v[n]);
/// }
/// return min_v;
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 3, 3, 7, 3]
///
/// n = prim::min(8, input);
/// // n: 1
template<typename T>
T min(const dlong N, const memory<T> v);


/// max: Parallel primitive returning maximum value in input array
///
/// Equivalent to the following code
/// max_v = -std::numeric_limits<T>::max();
/// for(std::size_t i = 0; i < N; ++i) {
///   max_v = std::max(max_v, v[n]);
/// }
/// return max_v;
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 3, 3, 7, 3]
///
/// n = prim::max(8, input);
/// // n: 7
template<typename T>
T max(const dlong N, const memory<T> v);


/// sum: Parallel primitive returning sum of all values in input array
///
/// Equivalent to the following code
/// sum_v = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   sum_v += v[n];
/// }
/// return sum_v;
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 3, 3, 7, 3]
///
/// n = prim::sum(8, input);
/// // n: 26
template<typename T>
T sum(const dlong N, const memory<T> v);


/// abs: Parallel primitive returning absolute value of an input array
///
/// Equivalent to the following code
/// for(std::size_t i = 0; i < N; ++i) {
///   absv[n] = std::abs(v[n]);
/// }
///
/// Example
/// memory<int> input; // e.g., [-1, 2, -3, -4, 3, -3, -7, -3]
/// memory<int> output; // empty array of 8 elements
///
/// n = prim::abs(8, input, output);
/// // output: [1, 2, 3, 4, 3, 3, 7, 3]
template<typename T>
void abs(const dlong N, const memory<T> v, memory<T> absv);


/// set: Parallel primitive assigning an input array to a value
///
/// Equivalent to the following code
/// for(std::size_t i = 0; i < N; ++i) {
///   v[n] = val;
/// }
///
/// Example
/// memory<int> output; // empty array of 8 elements
///
/// n = prim::set(8, 2, output);
/// // output: [2, 2, 2, 2, 2, 2, 2, 2]
template<typename T>
void set(const dlong N, const T val, memory<T> v);


/// range: Parallel primitive creating an ordered range in an
/// output array
///
/// Equivalent to the following code
/// for(std::size_t i = 0; i < N; ++i) {
///   v[n] = start + step * i;
/// }
///
/// Example
/// memory<int> output; // empty array of 8 elements
///
/// n = prim::range(8, 2, 4 output);
/// // output: [2, 6, 10, 14, 18, 22, 26, 30]
template<typename T>
void range(const dlong N, const T start, const T step, memory<T> v);


/// count: Parallel primitive for counting occurances of given value in
/// input array
///
/// Equivalent to the following code
/// count = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   count += (input[i] == value) ? 1 : 0;
/// }
/// return count;
///
/// Example
/// memory<int> input; // e.g., [1, 2, 1, 4, 1, 1, 7, 1]
///
/// n = prim::count(8, input, 1);
/// // n: 5
template<typename T>
dlong count(const dlong N, const memory<T> v, const T& value);


/// exclusiveScan: Parallel primitive for in-place exclusive prefix sum.
///
/// Equivalent to the following code
/// T sum = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   T v = v[i];
///   v[i] = sum;
///   sum += v;
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
///
/// prim::exclusiveScan(8, input);
/// // input: [0, 1, 3, 6, 10, 15, 21, 28]
template<typename T>
void exclusiveScan(const dlong N, memory<T> v);


/// exclusiveScan: Parallel primitive for exclusive prefix sum.
///
/// Equivalent to the following code
/// T sum = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   output[i] = sum;
///   sum += input[i];
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
/// memory<int> output; // empty array of 8 elements
///
/// prim::exclusiveScan(8, input, output);
/// // output: [0, 1, 3, 6, 10, 15, 21, 28]
template<typename T>
void exclusiveScan(const dlong N, const memory<T> v, memory<T> w);


/// inclusiveScan: Parallel primitive for in-place inclusive prefix sum.
///
/// Equivalent to the following code
/// T sum = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   sum += v[i];
///   v[i] = sum;
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
///
/// prim::inclusiveScan(8, input);
/// // input: [1, 3, 6, 10, 15, 21, 28, 36]
template<typename T>
void inclusiveScan(const dlong N, memory<T> v);


/// inclusiveScan: Parallel primitive for inclusive prefix sum.
///
/// Equivalent to the following code
/// T sum = 0;
/// for(std::size_t i = 0; i < N; ++i) {
///   sum += input[i];
///   output[i] = sum;
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
/// memory<int> output; // empty array of 8 elements
///
/// prim::inclusiveScan(8, input, output);
/// // output: [1, 3, 6, 10, 15, 21, 28, 36]
template<typename T>
void inclusiveScan(const dlong N, const memory<T> v, memory<T> w);


/// sort: Parallel increasing sort. Sorts an input array based on
/// the < comparison
///
/// Example
/// memory<double> input; // e.g., [0.6, 0.3, 0.65, 0.4, 0.2, 0.08, 1, 0.7]
///
/// prim::sort(8, input);
/// // input: [0.08, 0.2, 0.3, 0.4, 0.6, 0.65, 0.7, 1]
template<typename T>
void sort(const dlong N, memory<T> v);


/// sort: Parallel increasing sort with ids. Sorts an input array based on
/// the < comparison and returns an array of the original locations
/// for each entry.
///
/// Example
/// memory<double> input; // e.g., [0.6, 0.3, 0.65, 0.4, 0.2, 0.08, 1, 0.7]
/// memory<int> ids; // empty array of 8 elements
///
/// prim::sort(8, input, ids);
/// // input: [0.08, 0.2, 0.3, 0.4, 0.6, 0.65, 0.7, 1]
/// // ids: [5, 4, 1, 3, 0, 2, 7, 6]
template<typename T>
void sort(const dlong N, memory<T> v, memory<dlong> sortIds);


/// stableSort: Parallel stable increasing sort. Sorts an input array based on
/// the < comparison, preserving the relative order of equal elements.
///
/// Example
/// memory<double> input; // e.g., [0.6, 0.3, 0.65, 0.4, 0.2, 0.08, 1, 0.7]
///
/// prim::sort(8, input);
/// // input: [0.08, 0.2, 0.3, 0.4, 0.6, 0.65, 0.7, 1]
template<typename T>
void stableSort(const dlong N, memory<T> v);


/// stableSort: Parallel stable increasing sort with ids. Sorts an input array based on
/// the < comparison and returns an array of the original locations
/// for each entry. Preserves the relative order of equal elements.
///
/// Example
/// memory<double> input; // e.g., [0.6, 0.3, 0.65, 0.4, 0.2, 0.08, 1, 0.7]
/// memory<int> ids; // empty array of 8 elements
///
/// prim::sort(8, input, ids);
/// // input: [0.08, 0.2, 0.3, 0.4, 0.6, 0.65, 0.7, 1]
/// // ids: [5, 4, 1, 3, 0, 2, 7, 6]
template<typename T>
void stableSort(const dlong N, memory<T> v, memory<dlong> sortIds);


/// transformGather: Parallel primitive for indirectly reading
/// an input array
///
/// Equivalent to the following code
/// for(std::size_t i = 0; i < N; ++i) {
///   output[n] = input[ids[n]];
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
/// memory<int> ids;   // e.g., [2, 2, 7, 6, 1, 3, 3, 3]
/// memory<int> output; // empty array of 8 elements
///
/// prim::transformGather(8, ids, input, output);
/// // output: [3, 3, 8, 7, 2, 4, 4, 4]
template<typename T>
void transformGather(const dlong N,
                     const memory<dlong> ids,
                     const memory<T> v,
                           memory<T> w);


/// transformScatter: Parallel primitive for indirectly writing to
/// an output array. Note that the array of ids to be written
/// must not contain duplicates or data races may occur.
///
/// Equivalent to the following code
/// for(std::size_t i = 0; i < N; ++i) {
///   output[ids[n]] = input[n];
/// }
///
/// Example
/// memory<int> input; // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
/// memory<int> ids;   // e.g., [1, 0, 3, 2, 5, 4, 7, 6]
/// memory<int> output; // empty array of 8 elements
///
/// prim::transformScatter(8, ids, input, output);
/// // output: [2, 1, 4, 3, 6, 5, 8, 7]
template<typename T>
void transformScatter(const dlong N,
                      const memory<dlong> ids,
                      const memory<T> v,
                            memory<T> w);


/// adjacentDifference: Parallel primitive for applying a difference operation across pairs of
/// consecutive elements. Writes the output to the position of the left item.
///
/// Copies the first item to the output then performs the difference operator with each pair
/// of neighboring elements and writes its result to the location of the second element.
/// Equivalent to the following code
/// output[0] = input[0];
/// for(std::size_t i = 1; i < N; ++i) {
///     output[i] = input[i] - input[i - 1];
/// }
///
/// Example
/// memory<int> input; // e.g., [8, 7, 6, 5, 4, 3, 2, 1]
/// memory<int> output; // empty array of 8 elements
///
/// prim::adjacentDifference(8, input, output);
/// // output: [8, -1, -1, -1, -1, -1, -1, -1]
template<typename T>
void adjacentDifference(const dlong N, const memory<T> v, memory<T> diff);


/// adjacentDifferenceFlag: Parallel primitive for flagging differences across pairs of
/// consecutive elements. Writes the output to the position of the left item.
///
/// Writes a 1 to the first output location then performs the difference operator with each pair
/// of neighboring elements and writes a 1 to the location of the second element if a non-zero differnce exists.
/// Equivalent to the following code
/// output[0] = 1;
/// for(std::size_t i = 1; i < N; ++i) {
///     output[i] = (input[i] - input[i - 1]) ? 1 : 0;
/// }
///
/// Example
/// memory<int> input; // e.g., [8, 7, 7, 5, 5, 3, 3, 3]
/// memory<int> output; // empty array of 8 elements
///
/// prim::adjacentDifferenceFlag(8, input, output);
/// // output: [1, 1, 0, 1, 0, 1, 0, 0]
template<typename T>
void adjacentDifferenceFlag(const dlong N, const memory<T> v, memory<dlong> flag);


/// select: Parallel select primitive. Performs a selection of entries of an input
/// array based on an input flag. The output array must be large enough to hold all
/// locations.
///
/// Example
/// memory<int> input; // e.g., [1, 2, 2, 1, 3, 2, 4, 4]
/// memory<int> output(prim::count(8, input, 2)); // empty array of 3 elements
///
/// prim::select(8, input, 2, output);
/// // output: [1, 2, 5]
template<typename T>
void select(const dlong N, const memory<T> v, const T& val, memory<dlong> ids);


/// unique: Parallel unique primitive. From given input range, eliminates
/// all but the first element from every consecutive group of equivalent elements
/// and copies them into output.
///
/// Example
/// memory<int> input; // e.g., [1, 4, 2, 4, 4, 7, 7, 7]
/// memory<int> output; // uninitialized array
///
/// prim::select(8, input, Nunique, output);
/// // Nunique: 5
/// // output: [1, 4, 2, 4, 7]
template<typename T>
void unique(const dlong N,
            const memory<T> v,
                  dlong& Nunique,
                  memory<T>& v_unique);


/// runLengthEncode: Parallel run-length encoding. Encodes the size of groups of
/// consequtive values of the input array.
///
/// Example
/// memory<int> input; // e.g., [1, 1, 1, 2, 10, 10, 10, 88]
/// memory<int> output; // uninitialized array
///
/// prim::runLengthEncode(8, input, Ngroups, output);
/// // Ngroups: 4
/// // output: [0, 3, 4, 7, 8]
template<typename T>
void runLengthEncode(const dlong N,
                     const memory<T> v,
                     dlong& Ngroups,
                     memory<dlong>& offset);


/// runLengthEncodeConsecutive: Parallel run-length encoding.
/// Encodes the size of groups of cconsequtive values of the input array,
/// assuming precisely Ngroups of values exist from the range [0, Ngroups-1]
/// in order
///
/// Example
/// memory<int> input; // e.g., [1, 1, 1, 2, 4, 5, 5, 5]
/// memory<int> output; // empty array of 9 elements
///
/// prim::runLengthEncodeConsecutive(8, input, 6, output);
/// // output: [0, 0, 3, 4, 4, 5, 8]
template<typename T>
void runLengthEncodeConsecutive(const dlong N,
                                const memory<T> v,
                                dlong Ngroups,
                                memory<dlong> offset);


/// random: Parallel random number generation.
/// Creates an array of randomly generated numbers. The generation is
/// the same regardless of whether execution was parallel or serial.
/// For floating-point types, the randomly generated values will be
/// between 0. and 1., while for integer types the generated values
/// will be between 0 and the maximum signed value of the type.
template <typename T>
void random(const dlong N, memory<T> v);


/// seedRNG: Seed the random number generator used in prim::random
void seedRNG(const uint64_t seed);

}

} //namespace libp

#endif
