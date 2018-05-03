typedef struct{

  hlong row;
  hlong col;
  int ownerRank;
  dfloat val;

}nonZero_t;

typedef struct {

  long long int preconBytes;

  ogs_t *ogs;
  ogs_t *FEMogs;

  dfloat *zP;
  occa::memory o_zP;

  occa::memory o_Gr;
  occa::memory o_Gz;
  occa::memory o_Sr;

  occa::memory o_vmapPP;
  occa::memory o_faceNodesP;

  occa::memory o_oasForward;
  occa::memory o_oasBack;
  occa::memory o_oasDiagInvOp;

  occa::memory o_oasForwardDg;
  occa::memory o_oasBackDg;
  occa::memory o_oasDiagInvOpDg;
  occa::memory o_invDegreeDGP;

  occa::memory o_oasForwardDgT;
  occa::memory o_oasBackDgT;

  occa::kernel restrictKernel;

  occa::kernel coarsenKernel;
  occa::kernel prolongateKernel;

  occa::kernel overlappingPatchKernel;
  occa::kernel exactPatchSolverKernel;
  occa::kernel approxPatchSolverKernel;
  occa::kernel exactFacePatchSolverKernel;
  occa::kernel approxFacePatchSolverKernel;
  occa::kernel exactBlockJacobiSolverKernel;
  occa::kernel approxBlockJacobiSolverKernel;
  occa::kernel patchGatherKernel;
  occa::kernel facePatchGatherKernel;
  occa::kernel CGLocalPatchKernel;

  occa::memory o_rFEM;
  occa::memory o_zFEM;
  occa::memory o_GrFEM;
  occa::memory o_GzFEM;

  occa::kernel SEMFEMInterpKernel;
  occa::kernel SEMFEMAnterpKernel;

  ogs_t *ogsP, *ogsDg;

  occa::memory o_diagA;
  occa::memory o_invDiagA;
  occa::memory o_invAP;
  occa::memory o_invDegreeAP;
  occa::memory o_patchesIndex;

  // coarse grid basis for preconditioning
  occa::memory o_V1, o_Vr1, o_Vs1, o_Vt1;
  occa::memory o_r1, o_z1;
  dfloat *r1, *z1;

  void *xxt;

  occa::memory o_coarseInvDegree;

  int coarseNp;
  hlong coarseTotal;
  hlong *coarseOffsets;
  dfloat *B, *tmp2;
  occa::memory *o_B, o_tmp2;
  void *xxt2;
  parAlmond_t *parAlmond;

  // block Jacobi precon
  occa::memory o_invMM;
  occa::kernel blockJacobiKernel;
  occa::kernel partialblockJacobiKernel;

  //dummy almond level to store the OAS smoothing op
  agmgLevel *OASLevel;
  void **OASsmoothArgs;

  //SEMFEM variables
  mesh_t *femMesh;

} precon_t;


//Multigrid function callbacks
void AxTri2D        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void coarsenTri2D   (void **args, occa::memory &o_x, occa::memory &o_Rx);
void prolongateTri2D(void **args, occa::memory &o_x, occa::memory &o_Px);
void ellipticGather (void **args, occa::memory &o_x, occa::memory &o_Gx);
void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx);
void ellipticMultigridSmooth         (void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);
void ellipticMultigridSmoothChebyshev(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

//smoother ops
void LocalPatch  (void **args, occa::memory &o_r, occa::memory &o_Sr);
void dampedJacobi(void **args, occa::memory &o_r, occa::memory &o_Sr);