// OCCA+gslib gather scatter
typedef struct {
  
  iint         Ngather;          //  number of owned gather nodes
  iint         NtotalGather;     //  total number of gather nodes

  iint         *gatherOffsets;
  iint         *gatherHaloFlags;
  iint         *gatherBaseRanks;
  iint         *gatherLocalIds;
  iint         *gatherBaseIds;

  dfloat * gatherTmp;
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  
  void         *haloGsh;       // gslib gather 
  iint         Nhalo;            //  number of halo nodes
  iint         NownedHalo;       //  number of owned halo nodes
  occa::memory o_haloLocalIds;   //  list of halo nodes to
  occa::memory o_haloTmp;        //  temporary halo buffer
  dfloat       *haloTmp;         //  temporary HOST halo buffer

  iint *haloLocalIds;
  iint *haloGlobalIds;

  //degree vectors
  dfloat *invDegree, *gatherInvDegree;
  occa::memory o_invDegree;
  occa::memory o_gatherInvDegree;

}ogs_t;
