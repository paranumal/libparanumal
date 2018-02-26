// OCCA+gslib gather scatter
typedef struct {
  
  int         Ngather;     //  total number of gather nodes
  int         NtotalGather;     //  total number of gather nodes
  int         NnonHaloGather;       //  number of local gathered nodes 
  int         NhaloGather;          //  number of gathered nodes on halo

  int         *nonHaloGatherOffsets;
  int         *nonHaloGatherHaloFlags;
  int         *nonHaloGatherBaseRanks;
  int         *nonHaloGatherLocalIds;
  int         *nonHaloGatherBaseIds;

  int         *haloGatherOffsets;
  int         *haloGatherHaloFlags;
  int         *haloGatherBaseRanks;
  int         *haloGatherLocalIds;
  int         *haloGatherBaseIds;

  int *ownedHaloGatherIds;

  dfloat * haloGatherTmp;
  occa::memory o_nonHaloGatherOffsets;  //  start of local bases
  occa::memory o_nonHaloGatherLocalIds; //  base connected nodes
  occa::memory o_nonHaloGatherTmp;      //  DEVICE gather buffer

  occa::memory o_haloGatherOffsets;  //  start of local bases
  occa::memory o_haloGatherLocalIds; //  base connected nodes
  occa::memory o_haloGatherTmp;      //  DEVICE gather buffer
  
  occa::memory o_ownedHaloGatherIds;

  void         *haloGsh;       // gslib gather 
  int         Nhalo;            //  number of halo nodes
  int         NownedHalo;       //  number of owned halo nodes
  
  //degree vectors
  dfloat *invDegree, *gatherInvDegree;
  occa::memory o_invDegree;
  occa::memory o_gatherInvDegree;

}ogs_t;
