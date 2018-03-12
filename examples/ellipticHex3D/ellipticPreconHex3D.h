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
  occa::memory o_invDegreeP;

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
  mesh3D *femMesh;

} precon_t;

extern "C"
{
  void dgetrf_ (int *, int *, double *, int *, int *, int *);
  void dgetri_ (int *, double *, int *, int *, double *, int *, int *);
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );
}

//Multigrid function callbacks
void AxHex3D        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void coarsenHex3D   (void **args, occa::memory &o_x, occa::memory &o_Rx);
void prolongateHex3D(void **args, occa::memory &o_x, occa::memory &o_Px);
void ellipticGather (void **args, occa::memory &o_x, occa::memory &o_Gx);
void ellipticScatter(void **args, occa::memory &o_x, occa::memory &o_Sx);
void smoothHex3D    (void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);
void smoothChebyshevHex3D(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

//smoother ops
void overlappingPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr);
void FullPatchIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void FacePatchIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void LocalPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr);
void dampedJacobi  (void **args, occa::memory &o_r, occa::memory &o_Sr);