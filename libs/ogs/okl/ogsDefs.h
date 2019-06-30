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

/* the supported domains */
#define OGS_FOR_EACH_DOMAIN(macro) \
  macro(double   ) \
  macro(float    ) \
  macro(int      ) \
  macro(long     ) \
  macro(long long int)

/* the supported ops */
#define OGS_FOR_EACH_OP(T,macro) \
  macro(T,add) \
  macro(T,mul) \
  macro(T,min) \
  macro(T,max)

#define OGS_DO_add(a,b) a+=b
#define OGS_DO_mul(a,b) a*=b
#define OGS_DO_min(a,b) if(b<a) a=b
#define OGS_DO_max(a,b) if(b>a) a=b
