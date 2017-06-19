
typedef struct csr_t {

  iint Nrows;
  iint Ncols;

  iint NlocalCols;

  //local
  iint diagNNZ;
  iint   *diagRowStarts;
  iint   *diagCols;
  dfloat *diagCoefs;

  //non-local
  iint offdNNZ;
  iint   *offdRowStarts;
  iint   *offdCols;
  dfloat *offdCoefs;

  iint   *colMap;

  // MPI halo exchange info
  iint  NHalo;
  iint  NrecvTotal;  // number of elements to be sent in halo exchange
  iint  NsendTotal;
  iint  totalHaloPairs;
  iint *haloElementList; // sorted list of elements to be sent in halo exchange
  iint *NsendPairs;      // number of elements worth of data to send
  iint *NrecvPairs;      // number of elements worth of data to recv
  iint  NsendMessages;   // number of messages to send
  iint  NrecvMessages;   // number of messages to recv
  dfloat *sendBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} csr;


typedef struct ell_t {

  iint Nrows;
  iint Ncols;
  iint nnzPerRow;
  iint strideLength;
  iint actualNNZ;

  occa::memory o_cols;
  occa::memory o_coefs;

} ell;

typedef struct coo_t {

  iint Nrows;
  iint Ncols;
  iint nnz;

  // device memory
  occa::memory o_offsets;
  occa::memory o_cols;
  occa::memory o_coefs;

} coo;

typedef struct hyb_t {

  iint Nrows;
  iint Ncols;

  iint NlocalCols;

  coo *C;
  ell *E;

  occa::memory o_diagInv;
  occa::memory o_temp1;

  // MPI halo exchange info
  iint  NHalo;
  iint *colMap;
  iint  NrecvTotal;  // number of elements to be sent in halo exchange
  iint  NsendTotal;
  iint *haloElementList; // sorted list of elements to be sent in halo exchange
  occa::memory o_haloElementList;
  iint *NsendPairs;      // number of elements worth of data to send
  iint *NrecvPairs;      // number of elements worth of data to recv
  iint  NsendMessages;   // number of messages to send
  iint  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} hyb;


typedef struct dcsr_t {

  iint Nrows;
  iint Ncols;

  iint NlocalCols;

  //local
  iint diagNNZ;
  occa::memory o_diagRowStarts;
  occa::memory o_diagCols;
  occa::memory o_diagCoefs;

  //non-local
  iint offdNNZ;
  occa::memory o_offdRowStarts;
  occa::memory o_offdCols;
  occa::memory o_offdCoefs;

  occa::memory o_diagInv;
  occa::memory o_temp1;

  // MPI halo exchange info
  iint  NHalo;
  iint  NrecvTotal;  // number of elements to be sent in halo exchange
  iint  NsendTotal;
  iint  totalHaloPairs;
  iint *haloElementList; // sorted list of elements to be sent in halo exchange
  iint *NsendPairs;      // number of elements worth of data to send
  iint *NrecvPairs;      // number of elements worth of data to recv
  iint  NsendMessages;   // number of messages to send
  iint  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;

  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} dcsr;

typedef struct agmgLevel_t {
  iint Nrows;
  iint Ncols;

  csr *A;
  csr *P;
  csr *R;

  iint *globalRowStarts; //global partitioning of fine level
  iint *globalAggStarts; //global partitioning of coarse level

  hyb  *deviceA;
  dcsr  *dcsrP;
  hyb  *deviceR;

  dfloat *nullA;

  dfloat *rhs, *res, *x;

  dfloat *ckp1, *vkp1, *wkp1;

  occa::memory o_rhs, o_res, o_x;
  occa::memory o_ckp1, o_vkp1, o_wkp1;

  dfloat *smoother_params;

  dfloat threshold;
  iint numAggregates;
  SmoothType stype;

} agmgLevel;
