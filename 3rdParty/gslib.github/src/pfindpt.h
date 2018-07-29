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

#ifndef PFINDPT_H
#define PFINDPT_H

/* requires "types.h", "tuple_list.h", and, 
   when MPI is defined, "crystal.h" */
#if !defined(TYPES_H) || !defined(TUPLE_LIST_H) || \
    ( defined(MPI) && !defined(CRYSTAL_H) )
#warning "pfindpt.h" requires "types.h", "tuple_list.h", and "crystal.h"
#endif

typedef struct pfindpt_data_ pfindpt_data;

#ifndef MPI
#  define crystal_data void
#endif

pfindpt_data *pfindpt_setup(unsigned ndim, const real *const*xw,
                            const unsigned *n, uint nel,
                            uint max_hash_size, real bbox_tol,
                            crystal_data *crystal);

#ifndef MPI
#  undef crystal_data
#endif
                            
void pfindpt_free(pfindpt_data *p);
void pfindpt_transfer(pfindpt_data *p, tuple_list *list, int dynamic);
void pfindpt(pfindpt_data *p, tuple_list *list, int guess);
void pfindpt_weights(pfindpt_data *p, const real *r);
real pfindpt_eval(pfindpt_data *p, const real *u);

#endif

