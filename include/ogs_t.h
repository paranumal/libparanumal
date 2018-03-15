// OCCA+gslib gather scatter
typedef struct {
  
  int         Ngather;          //  number of gather nodes

  int         *gatherOffsets;
  int         *gatherHaloFlags;
  int         *gatherBaseRanks;
  int         *gatherLocalIds;
  int         *gatherBaseIds;

  char * gatherTmp;
  
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  void         *gatherGsh;       // gslib gather 

  int         Nscatter;
  occa::memory o_scatterOffsets; //  start of local bases
  occa::memory o_scatterLocalIds;//  base connected nodes
  
  int         Nhalo;            //  number of halo nodes
  occa::memory o_haloLocalIds;   //  list of halo nodes to
  occa::memory o_haloTmp;        //  temporary halo buffer
  void         *haloTmp;         //  temporary HOST halo buffer

}ogs_t;
