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

#ifndef TYPES_H
#define TYPES_H

/* 
  Define the integer types used throughout the code,
  controlled by preprocessor macros.
  
  The integer type sint/uint (signed/unsigned) is used
  most frequently, e.g., for indexing into local arrays,
  and for processor ids. It can be one of
  
    macro             sint/uint type
    
    (default)         int
    USE_LONG          long
    USE_LONG_LONG     long long
    
  The slong/ulong type is used in relatively few places
  for global identifiers and indices. It can be one of

    macro             slong/ulong type
    
    (default)         int
    GLOBAL_LONG       long
    GLOBAL_LONG_LONG  long long

  Since the long long type is not ISO C90, it is never
  used unless explicitly asked for.
*/

#if defined(USE_LONG_LONG) || defined(GLOBAL_LONG_LONG)
typedef long long long_long;
#  define WHEN_LONG_LONG(x) x
#  if !defined(LLONG_MAX)
#    if defined(LONG_LONG_MAX)
#      define LLONG_MAX LONG_LONG_MAX
#    else
#      define LLONG_MAX 9223372036854775807
#    endif
#  endif
#  if !defined(LLONG_MIN)
#    if defined(LONG_LONG_MIN)
#      define LLONG_MIN LONG_LONG_MIN
#    else
#      define LLONG_MIN -9223372036854775807
#    endif
#  endif
#else
#  define WHEN_LONG_LONG(x)
#endif

#if !defined(USE_LONG) && !defined(USE_LONG_LONG)
#  define TYPE_LOCAL(i,l,ll) i
#elif defined(USE_LONG)
#  define TYPE_LOCAL(i,l,ll) l
#elif defined(USE_LONG_LONG)
#  define TYPE_LOCAL(i,l,ll) ll
#endif

#if !defined(GLOBAL_LONG) && !defined(GLOBAL_LONG_LONG)
#  define TYPE_GLOBAL(i,l,ll) i
#elif defined(GLOBAL_LONG)
#  define TYPE_GLOBAL(i,l,ll) l
#else
#  define TYPE_GLOBAL(i,l,ll) ll
#endif

/* local integer type: for quantities O(N/P) */
#define sint   signed TYPE_LOCAL(int,long,long long)
#define uint unsigned TYPE_LOCAL(int,long,long long)
#define iabs TYPE_LOCAL(abs,labs,llabs)

/* global integer type: for quantities O(N) */
#define slong   signed TYPE_GLOBAL(int,long,long long)
#define ulong unsigned TYPE_GLOBAL(int,long,long long)
#define iabsl TYPE_GLOBAL(abs,labs,llabs)

#endif

