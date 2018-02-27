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
