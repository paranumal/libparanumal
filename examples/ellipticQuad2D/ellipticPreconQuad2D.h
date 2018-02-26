typedef struct{

  int row;
  int col;
  int ownerRank;
  dfloat val;

}nonZero_t;

typedef struct {

  long long int preconBytes;

  dfloat *zP;
  occa::memory o_zP;

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

  ogs_t *ogsP, *ogsDg;
  hgs_t *hgsP, *hgsDg;

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
  int coarseTotal;
  int *coarseOffsets;
  dfloat *B, *tmp2;
  occa::memory *o_B, o_tmp2;
  void *xxt2;
  parAlmond_t *parAlmond;

  // block Jacobi precon
  occa::memory o_invMM;
  occa::kernel blockJacobiKernel;

  //dummy almond level to store the OAS smoothing op
  agmgLevel *OASLevel;
  void **OASsmoothArgs;

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

void ellipticBuildIpdgQuad2D(mesh2D *mesh, dfloat tau, dfloat lambda, int *BCType, nonZero_t **A,
                              int *nnzA, int *globalStarts, const char *options);

void ellipticBuildContinuousQuad2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, int *nnz,
                              hgs_t **hgs, int *globalStarts, const char* options);

void ellipticBuildPatchesIpdgQuad2D(mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   int *BCType, nonZero_t **A, int *nnzA,
                                   hgs_t **hgs, int *globalStarts,
                                   int *Npataches, int **patchesIndex, dfloat **patchesInvA, dfloat **localA,
                                   const char *options);

void ellipticCoarsePreconditionerSetupQuad2D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **V1, nonZero_t **A, int *nnzA,
                                   hgs_t **hgs, int *globalStarts, const char *options);

void ellipticBuildJacobiIpdgQuad2D(mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options);

void ellipticBuildFullPatchesIpdgQuad2D(mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   int *Npataches, int **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

void ellipticBuildFacePatchesIpdgQuad2D(mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   int *Npataches, int **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

void ellipticBuildLocalPatchesIpdgQuad2D(mesh2D *mesh, int basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   int *Npataches, int **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

//Multigrid function callbacks
void AxQuad2D        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void coarsenQuad2D   (void **args, occa::memory &o_x, occa::memory &o_Rx);
void prolongateQuad2D(void **args, occa::memory &o_x, occa::memory &o_Px);
void smoothQuad2D    (void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);
void smoothChebyshevQuad2D(void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

//smoother ops
void overlappingPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr);
void FullPatchIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void FacePatchIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void LocalPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr);
void dampedJacobi  (void **args, occa::memory &o_r, occa::memory &o_Sr);