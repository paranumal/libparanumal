// OCCA+gslib gather scatter
typedef struct {
  
  iint         Ngather;     //  total number of gather nodes
  iint         NtotalGather;     //  total number of gather nodes
  iint         NnonHaloGather;       //  number of local gathered nodes 
  iint         NhaloGather;          //  number of gathered nodes on halo

  iint         *nonHaloGatherOffsets;
  iint         *nonHaloGatherHaloFlags;
  iint         *nonHaloGatherBaseRanks;
  iint         *nonHaloGatherLocalIds;
  iint         *nonHaloGatherBaseIds;

  iint         *haloGatherOffsets;
  iint         *haloGatherHaloFlags;
  iint         *haloGatherBaseRanks;
  iint         *haloGatherLocalIds;
  iint         *haloGatherBaseIds;

  iint *ownedHaloGatherIds;

  dfloat * haloGatherTmp;
  occa::memory o_nonHaloGatherOffsets;  //  start of local bases
  occa::memory o_nonHaloGatherLocalIds; //  base connected nodes
  occa::memory o_nonHaloGatherTmp;      //  DEVICE gather buffer

  occa::memory o_haloGatherOffsets;  //  start of local bases
  occa::memory o_haloGatherLocalIds; //  base connected nodes
  occa::memory o_haloGatherTmp;      //  DEVICE gather buffer
  
  occa::memory o_ownedHaloGatherIds;

  void         *haloGsh;       // gslib gather 
  iint         Nhalo;            //  number of halo nodes
  iint         NownedHalo;       //  number of owned halo nodes
  
  //degree vectors
  dfloat *invDegree, *gatherInvDegree;
  occa::memory o_invDegree;
  occa::memory o_gatherInvDegree;

}ogs_t;
