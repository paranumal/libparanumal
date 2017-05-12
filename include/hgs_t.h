// Holmes' gather + scatter
typedef struct {
  
  iint         Ngather;          //  number of gather nodes
  iint         Nlocal;           //  number of local nodes before gather      

  iint         *sortIds;         
  iint         *gatherOffsets;   
  iint         *gatherLocalIds;
  dfloat       *invDegree;

  occa::memory o_sortIds;        //  sorting for exchagne + gather
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  occa::memory o_invDegree;      //  inverse degree weighting
  
  iint         NsendTotal;            //  number of nodes to send
  iint         NrecvTotal;            //  number of nodes to recv
  iint         haloOffset; 
  iint        *Nsend;
  iint        *Nrecv;
  dfloat      *sendBuffer;
  dfloat      *recvBuffer;
  void        *haloSendRequests;
  void        *haloRecvRequests;

}hgs_t;
