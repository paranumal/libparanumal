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

#ifndef OGSBASE_HPP
#define OGSBASE_HPP

#include "ogs.hpp"

namespace libp {

namespace ogs {

//forward declarations
class ogsOperator_t;
class ogsExchange_t;

class ogsBase_t {
public:
  platform_t platform;
  comm_t comm;

  dlong         N=0;
  dlong         Ngather=0;        //  total number of local positive gather nodes

  dlong         NlocalT=0;        //  number of local gather nodes
  dlong         NhaloT=0;         //  number of halo gather nodes
  dlong         NlocalP=0;        //  number of positive local gather nodes
  dlong         NhaloP=0;         //  number of positive halo gather nodes

  hlong         NgatherGlobal=0;  //  global number of positive gather nodes

  Kind kind;
  bool unique=false;
  bool gather_defined=false;

  static stream_t dataStream;

  ogsBase_t()=default;
  virtual ~ogsBase_t()=default;

  virtual void Setup(const dlong _N,
                      memory<hlong> ids,
                      comm_t _comm,
                      const Kind _kind,
                      const Method method,
                      const bool _unique,
                      const bool verbose,
                      platform_t& _platform);
  void Free();

protected:
  std::shared_ptr<ogsOperator_t> gatherLocal;
  std::shared_ptr<ogsOperator_t> gatherHalo;
  std::shared_ptr<ogsExchange_t> exchange;

  void AssertGatherDefined();

private:
  void FindSharedGroups(const dlong Nids,
                        memory<hlong> baseIds,
                        memory<int> sharedFlag,
                        const int verbose);
  void ConstructGatherOperators(const dlong Nlocal,
                                const memory<dlong> localRowIds,
                                const memory<dlong> localColIds,
                                const memory<hlong> localBaseIds,
                                const dlong Nhalo,
                                const memory<dlong> haloRowIds,
                                const memory<dlong> haloColIds,
                                const memory<hlong> haloBaseIds);

  void ConstructSharedNodes(const memory<hlong> haloBaseIds,
                            dlong &Nshared,
                            memory<int>&   sharedRemoteRanks,
                            memory<dlong>& sharedLocalRows,
                            memory<dlong>& sharedRemoteRows,
                            memory<hlong>& sharedLocalBaseIds,
                            memory<hlong>& sharedRemoteBaseIds);

  ogsExchange_t* AutoSetup(const dlong Nshared,
                           const memory<int>   sharedRemoteRanks,
                           const memory<dlong> sharedLocalRows,
                           const memory<dlong> sharedRemoteRows,
                           const memory<hlong> sharedLocalBaseIds,
                           const memory<hlong> sharedRemoteBaseIds,
                           const memory<hlong> haloBaseIds,
                           ogsOperator_t& gatherHalo,
                           comm_t _comm,
                           platform_t &_platform,
                           const int verbose);
};

} //namespace ogs

} //namespace libp

#endif
