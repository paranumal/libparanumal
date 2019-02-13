#ifndef OPERATORS_H
#define OPERATORS_H 1

#include "adaptive.h"

void get_operators(int N, occa::device &device, occa::memory &o_r,
                   occa::memory &o_w, occa::memory &o_D, occa::memory &o_Ib,
                   occa::memory &o_It, occa::memory &o_Pb, occa::memory &o_Pt);

#endif
