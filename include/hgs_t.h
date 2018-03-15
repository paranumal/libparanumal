// Holmes' gather + scatter
typedef struct {
  
  int         Ngather;          //  number of gather nodes
  int         Nlocal;           //  number of local nodes before gather      

  int         *sortIds;         
  int         *gatherOffsets;   
  int         *gatherLocalIds;
  dfloat       *invDegree;

  occa::memory o_sortIds;        //  sorting for exchagne + gather
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  occa::memory o_invDegree;      //  inverse degree weighting
  
  int         NsendTotal;            //  number of nodes to send
  int         NrecvTotal;            //  number of nodes to recv
  int         haloOffset; 
  int        *Nsend;
  int        *Nrecv;
  dfloat      *sendBuffer;
  dfloat      *recvBuffer;
  void        *haloSendRequests;
  void        *haloRecvRequests;

}hgs_t;
