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

#ifndef MESH_DEFINES3D_H
#define MESH_DEFINES3D_H 1

/* offsets for geometric factors */
#define RXID 0
#define RYID 1
#define SXID 2
#define SYID 3
#define  JID 4
#define JWID 5
#define IJWID 6
#define RZID 7
#define SZID 8
#define TXID 9
#define TYID 10
#define TZID 11

/* offsets for second order geometric factors */
#define G00ID 0
#define G01ID 1
#define G11ID 2
#define GWJID 3
#define G12ID 4
#define G02ID 5
#define G22ID 6


/* offsets for nx, ny, sJ, 1/J */
#define NXID 0
#define NYID 1
#define SJID 2
#define IJID 3
#define IHID 4
#define WSJID 5
#define WIJID 6
#define NZID 7
#define STXID 8
#define STYID 9
#define STZID 10
#define SBXID 11
#define SBYID 12
#define SBZID 13
#define SURXID 14
#define SURYID 15
#define SURZID 16

// //offsets for boltzmann PML variables
// #define QXID1 0
// #define QXID2 1
// #define QXID3 2
// #define QXID4 3
// #define QXID5 4
// #define QXID6 5
// #define QXID8 6
// //
// #define QYID1 7
// #define QYID2 8
// #define QYID3 9
// #define QYID4 10
// #define QYID5 11
// #define QYID7 12
// #define QYID9 13
// //
// #define QZID1 14
// #define QZID2 15
// #define QZID3 16
// #define QZID4 17
// #define QZID6 18
// #define QZID7 19
// #define QZID10  20

#endif

