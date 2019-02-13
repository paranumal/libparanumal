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

#ifndef ADAPTIVE_H
#define ADAPTIVE_H 1

// {{{ Dimension
#ifdef APP_DIMENSION
#define DIM APP_DIMENSION
#else
#define DIM 3
#endif
// }}}

// {{{ Headers
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wcast-align"
#endif
#if DIM == 2
#include <p4est.h>
#include <p4est_bits.h>
#include <p4est_connectivity.h>
#include <p4est_extended.h>
#include <p4est_ghost.h>
#include <p4est_iterate.h>
#include <p4est_lnodes.h>
#include <p4est_mesh.h>
#include <p4est_nodes.h>
#include <p4est_vtk.h>
#define P4EST_EDGES 0
#elif DIM == 3
#include <p4est_to_p8est.h>
#include <p8est.h>
#include <p8est_bits.h>
#include <p8est_connectivity.h>
#include <p8est_extended.h>
#include <p8est_ghost.h>
#include <p8est_iterate.h>
#include <p8est_lnodes.h>
#include <p8est_mesh.h>
#include <p8est_nodes.h>
#include <p8est_tets_hexes.h>
#include <p8est_vtk.h>
#define P4EST_EDGES 12
#else
#error "bad dimension"
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include "occa.hpp"
#include "setupAide.hpp"
#include "asd.h"
#include "number_types.h"
#include "solver_info.h"
#include "app.h"
#include "myoccautil.h"
#include "connectivity.h"
// }}}

#endif
