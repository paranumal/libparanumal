// OCCA+gslib gather scatter
typedef struct {
  
  iint         Ngather;          //  number of gather nodes

  iint         *gatherOffsets;
  iint         *gatherHaloFlags;
  iint         *gatherBaseRanks;
  iint         *gatherLocalIds;
  iint         *gatherBaseIds;
  
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  void         *gatherGsh;       // gslib gather 

  iint         Nscatter;
  occa::memory o_scatterOffsets; //  start of local bases
  occa::memory o_scatterLocalIds;//  base connected nodes
  
  iint         Nhalo;            //  number of halo nodes
  occa::memory o_haloLocalIds;   //  list of halo nodes to
  occa::memory o_haloTmp;        //  temporary halo buffer
  void         *haloTmp;         //  temporary HOST halo buffer

}ogs_t;
