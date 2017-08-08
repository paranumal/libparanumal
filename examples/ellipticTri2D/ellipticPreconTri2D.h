typedef struct{

  iint row;
  iint col;
  iint ownerRank;
  dfloat val;

}nonZero_t;

typedef struct {

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
  occa::kernel exactBlockJacobiSolverKernel;
  occa::kernel approxBlockJacobiSolverKernel;
  occa::kernel patchGatherKernel;

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

  iint coarseNp;
  iint coarseTotal;
  iint *coarseOffsets;
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

void ellipticBuildIpdgTri2D(mesh2D *mesh, dfloat tau, dfloat lambda, iint *BCType, nonZero_t **A, iint *nnzA,
                              hgs_t **hgs, iint *globalStarts, const char *options);

void ellipticBuildContinuousTri2D(mesh2D *mesh, dfloat lambda, nonZero_t **A, iint *nnz,
                              hgs_t **hgs, iint *globalStarts, const char* options);

void ellipticBuildPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA, dfloat **localA,
                                   const char *options);

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **V1, nonZero_t **A, iint *nnzA,
                                   hgs_t **hgs, iint *globalStarts, const char *options);

void ellipticBuildJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options);


void ellipticBuildExactPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   dfloat **patchesInvA, const char *options);

void ellipticBuildApproxPatchesIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);
void ellipticBuildExactBlockJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   dfloat **patchesInvA, const char *options);

void ellipticBuildApproxBlockJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);


void ellipticBuildJacobiIpdgTri2D(mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options);

//Multigrid function callbacks
void AxTri2D        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void coarsenTri2D   (void **args, occa::memory &o_x, occa::memory &o_Rx);
void prolongateTri2D(void **args, occa::memory &o_x, occa::memory &o_Px);
void smoothTri2D    (void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);
void smoothChebyshevTri2D    (void **args, occa::memory &o_r, occa::memory &o_x, bool xIsZero);

//smoother ops
void overlappingPatchIpdg(void **args, occa::memory &o_r, occa::memory &o_Sr);
void exactFullPatchIpdg  (void **args, occa::memory &o_r, occa::memory &o_Sr);
void approxFullPatchIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void exactBlockJacobiIpdg  (void **args, occa::memory &o_r, occa::memory &o_Sr);
void approxBlockJacobiIpdg (void **args, occa::memory &o_r, occa::memory &o_Sr);
void dampedJacobi        (void **args, occa::memory &o_r, occa::memory &o_Sr);