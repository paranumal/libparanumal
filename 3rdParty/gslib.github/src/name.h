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

#ifndef NAME_H
#define NAME_H

/* establishes some macros to establish
   * the FORTRAN naming convention
     default      gs_setup, etc.
     -DUPCASE     GS_SETUP, etc.
     -DUNDERSCORE gs_setup_, etc.
   * a prefix for all external (non-FORTRAN) function names
     for example, -DPREFIX=jl_   transforms fail -> jl_fail
   * a prefix for all external FORTRAN function names     
     for example, -DFPREFIX=jlf_ transforms gs_setup_ -> jlf_gs_setup_
*/

/* the following macro functions like a##b,
   but will expand a and/or b if they are themselves macros */
#define TOKEN_PASTE_(a,b) a##b
#define TOKEN_PASTE(a,b) TOKEN_PASTE_(a,b)

#ifdef PREFIX
#  define PREFIXED_NAME(x) TOKEN_PASTE(PREFIX,x)
#else
#  define PREFIXED_NAME(x) x
#endif

#ifdef FPREFIX
#  define FPREFIXED_NAME(x) TOKEN_PASTE(FPREFIX,x)
#else
#  define FPREFIXED_NAME(x) x
#endif

#if defined(UPCASE)
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(up)
#  define FORTRAN_UNPREFIXED(low,up) up
#elif defined(UNDERSCORE)
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(TOKEN_PASTE(low,_))
#  define FORTRAN_UNPREFIXED(low,up) TOKEN_PASTE(low,_)
#else
#  define FORTRAN_NAME(low,up) FPREFIXED_NAME(low)
#  define FORTRAN_UNPREFIXED(low,up) low
#endif

#endif

