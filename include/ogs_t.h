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

// OCCA+gslib gather scatter
typedef struct {
  
  dlong         Ngather;     //  total number of gather nodes
  dlong         NtotalGather;     //  total number of gather nodes
  dlong         NnonHaloGather;       //  number of local gathered nodes 
  dlong         NhaloGather;          //  number of gathered nodes on halo

  dlong       *nonHaloGatherOffsets;
  int         *nonHaloGatherHaloFlags;
  int         *nonHaloGatherBaseRanks;
  dlong       *nonHaloGatherLocalIds;
  hlong       *nonHaloGatherBaseIds;

  dlong       *haloGatherOffsets;
  int         *haloGatherHaloFlags;
  int         *haloGatherBaseRanks;
  dlong       *haloGatherLocalIds;
  hlong       *haloGatherBaseIds;

  dlong    *ownedHaloGatherIds;

  dfloat * haloGatherTmp;
  occa::memory o_nonHaloGatherOffsets;  //  start of local bases
  occa::memory o_nonHaloGatherLocalIds; //  base connected nodes
  occa::memory o_nonHaloGatherTmp;      //  DEVICE gather buffer

  occa::memory o_haloGatherOffsets;  //  start of local bases
  occa::memory o_haloGatherLocalIds; //  base connected nodes
  occa::memory o_haloGatherTmp;      //  DEVICE gather buffer
  
  occa::memory o_ownedHaloGatherIds;

  void         *haloGsh;       // gslib gather 
  dlong         Nhalo;            //  number of halo nodes
  dlong         NownedHalo;       //  number of owned halo nodes
  
  //degree vectors
  dfloat *invDegree, *gatherInvDegree;
  occa::memory o_invDegree;
  occa::memory o_gatherInvDegree;

}ogs_t;
