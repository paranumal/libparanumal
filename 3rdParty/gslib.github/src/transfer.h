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

#ifndef TRANSFER_H
#define TRANSFER_H

#ifdef MPI

#if !defined(TUPLE_LIST_H) || !defined(CRYSTAL_H)
#warning "transfer.h" requires "tuple_list.h" and "crystal.h"
#endif

/*------------------------------------------------------------------------------
  
  Transfer
 
  Treats one integer (not long) member of the tuple list as a target proc;
  Sends out tuples accordingly, using the crystal router.
  Target proc member overwritten with source proc.
  
  dynamic: non-zero if the tuple list should grow to accomodate arrivals
  tl:      the tuple list
  pf:      which tuple member specifies target proc
  crystal: an initialized crystal router structure (cf. crystal.h)

  ----------------------------------------------------------------------------*/
void transfer(int dynamic, tuple_list* tl,
              unsigned pf, crystal_data *crystal);

#endif

#endif

