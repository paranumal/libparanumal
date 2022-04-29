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

#ifndef TYPES_HPP
#define TYPES_HPP

// precision of AMG storage
#if 0
#define pfloat float
#else
#define pfloat double
#endif


//float data type
#if 0
#define dfloat float
#define dfloatFormat "%f"
#define dfloatString "float"
#else
#define dfloat double
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

//host index data type
#if 0
#define hlong int
#define hlongFormat "%d"
#define hlongString "int"
#else
#define hlong long long int
#define hlongFormat "%lld"
#define hlongString "long long int"
#endif

//device index data type
#if 1
#define dlong int
#define dlongFormat "%d"
#define dlongString "int"
#else
#define dlong long long int
#define dlongFormat "%lld"
#define dlongString "long long int"
#endif

#endif
