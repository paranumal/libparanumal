#include "adaptive.h"


void asd_jacobi_gauss_quadrature(ldouble alpha, ldouble beta,
                                        int N, ldouble *x,
				 ldouble *w);

// {{{ Level

// {{{ Kernel Info
static void level_kernelinfo(occa::properties &info, occa::device &device, int N, int NqGJ,
    int periodic_brick)
{
  const int Nq = N + 1;
  const int Nqk = Nq;
  const int Np = Nq * Nq * Nq;
  const int Nfp = Nq * Nq;
  const int Nfaces = 2 * DIM;

 // FIXME
  // occaKernelInfoAddParserFlag(info, "automate-add-barriers", "disabled");

  info["defines"].asObject();
  
  info["defines/ONE"] = (dfloat_t) 1.0;
  info["defines/TWO"] = (dfloat_t) 2.0;
  info["defines/HALF"] = (dfloat_t) 0.5;
  info["defines/ZERO"] = (dfloat_t) 0.0;

  info["defines/dlong"] = occa_iint_name;

  const char *const dfloatStringLocal =
    (sizeof(double) == sizeof(dfloat_t)) ? "double" : "float";

  info["defines/dfloat"] = dfloatStringLocal;
  if (sizeof(double) == sizeof(dfloat_t))
    info["defines/p_DFLOAT_DOUBLE"] = 1;
  else
    info["defines/p_DFLOAT_FLOAT"] = 1;

  info["defines/p_DFLOAT_MAX"] = DFLOAT_MAX;

  info["defines/p_KblkV"] = KERNEL_KBLKV;
  info["defines/p_KblkS"] = KERNEL_KBLKS;
  info["defines/p_Nt"] = KERNEL_NT;

  info["defines/p_REDUCE_LDIM"] = KERNEL_REDUCE_LDIM;

  info["defines/p_DIM"] = DIM;

  ASD_ASSERT(sizeof(p4est_qcoord_t) <= sizeof(iint_t));
  info["defines/p_P4EST_ROOT_LEN"] = P4EST_ROOT_LEN;
  info["defines/p_P4EST_HALF"] = P4EST_HALF;
  info["defines/p_P4EST_FACES"] = P4EST_FACES;
  info["defines/p_P4EST_EDGES"] = P4EST_EDGES;

  info["defines/p_NX"] = NX;

  if (periodic_brick)
    info["defines/p_PERIODIC_BRICK"] = 1;
#if 0
  info["defines/p_FIELD_UX"] = FIELD_UX;
  info["defines/p_FIELD_UY"] = FIELD_UY;
  info["defines/p_FIELD_UZ"] = FIELD_UZ;
#endif
  info["defines/p_NFIELDS"] = NFIELDS;

  info["defines/p_VGEO_X"] = VGEO_X;
  info["defines/p_VGEO_Y"] = VGEO_Y;
  info["defines/p_VGEO_RX"] = VGEO_RX;
  info["defines/p_VGEO_SX"] = VGEO_SX;
  info["defines/p_VGEO_RY"] = VGEO_RY;
  info["defines/p_VGEO_SY"] = VGEO_SY;
  info["defines/p_VGEO_J"] = VGEO_J;
  info["defines/p_NVGEO"] = NVGEO;

  info["defines/p_VGEO_Z"] = VGEO_Z;
  info["defines/p_VGEO_TX"] = VGEO_TX;
  info["defines/p_VGEO_TY"] = VGEO_TY;
  info["defines/p_VGEO_RZ"] = VGEO_RZ;
  info["defines/p_VGEO_SZ"] = VGEO_SZ;
  info["defines/p_VGEO_TZ"] = VGEO_TZ;

  info["defines/p_SGEO_NX"] = SGEO_NX;
  info["defines/p_SGEO_NY"] = SGEO_NY;
  info["defines/p_SGEO_NZ"] = SGEO_NZ;

  info["defines/p_SGEO_SJ"] = SGEO_SJ;
  info["defines/p_NSGEO"] = NSGEO;

  info["defines/p_GGEO_RR"] = GGEO_RR;
  info["defines/p_GGEO_RS"] = GGEO_RS;
  info["defines/p_GGEO_RT"] = GGEO_RT;
  info["defines/p_GGEO_SS"] = GGEO_SS;
  info["defines/p_GGEO_ST"] = GGEO_ST;
  info["defines/p_GGEO_TT"] = GGEO_TT;
  info["defines/p_GGEO_JW"] = GGEO_JW;
  info["defines/p_NGGEO"] = NGGEO;
  
  info["defines/p_BC_SKIP"] = BC_SKIP;
  info["defines/p_BC_NONE"] = BC_NONE;
  info["defines/p_BC_DEFAULT"] = BC_DEFAULT;

  info["defines/p_M_PI"] = M_PI;
  info["defines/p_M_PI_4"] = M_PI_4;

  info["defines/conn_mapping"] = "conn_mapping_identity";

  info["defines/p_NCORNERS"] = P4EST_CHILDREN;
  info["defines/p_NHANG"] = P4EST_HALF;

  info["defines/p_N"] = N;
  info["defines/p_Nq"] = Nq;
  info["defines/p_Nqk"] = Nqk;
  info["defines/p_Nq2"] = Nq * Nq;
  info["defines/p_Nq3"] = Nq * Nq * Nq;
  info["defines/p_Np"] = Np;
  info["defines/p_Nfaces"] = Nfaces;
  info["defines/p_Nfp"] = Nfp;

  info["defines/p_NqGJ"] = NqGJ;
  info["defines/p_NqGJ2"] = NqGJ*NqGJ;
  info["defines/p_NpGJ"] = NqGJ*NqGJ*NqGJ;
  
  info["defines/p_NelementsblkV"] = 1;
  info["defines/p_NelementsblkS"] = 1;

  int blockSize = KERNEL_REDUCE_LDIM;
  
  info["defines/p_blockSize"] = blockSize;
  info["defines/p_NthreadsUpdatePCG"] = blockSize;
  info["defines/p_NwarpsUpdatePCG"] = (int) (blockSize/32); // WARNING: CUDA SPECIFIC
  
  info["includes"] = (char*)strdup(DADAPTIVE "/okl/adaptiveOcca.h");

  std::cout << info << std::endl;
  
}
// }}}

static void level_get_mesh_constants(level_t *lvl, mesh_t *mesh)
{
  lvl->N = mesh->N;
  lvl->Nq = mesh->Nq;
  lvl->Nfaces = mesh->Nfaces;
  lvl->Nqk = mesh->Nqk;
  lvl->Nfp = mesh->Nfp;
  lvl->Np = mesh->Np;

  lvl->Nmortar = mesh->Nmortar;
  lvl->Ncontinuous = mesh->Ncontinuous;
  lvl->Ncindices = mesh->Ncindices;

  lvl->Klocal = mesh->Klocal;
  lvl->Kghost = mesh->Kghost;
  lvl->Kuniqmirror = mesh->Kuniqmirror;
  lvl->Kmirror = mesh->Kmirror;
  lvl->Ktotal = mesh->Ktotal;
  lvl->Kintra = mesh->Kintra;
}


static void level_get_mesh(level_t *lvl, mesh_t *mesh, p4est_t *pxest,
                           p4est_ghost_t *ghost, occa::device &device,
                           int brick)
{
  // paranoid checks
  ASD_ABORT_IF(mesh->Kintra > lvl->Kmax, "Kmax not big enough (Kintra)");
  ASD_ABORT_IF(mesh->Kmirror > lvl->Kmax, "Kmax not big enough (Kmirror)");
  ASD_ABORT_IF(mesh->Ktotal > lvl->Kmax, "Kmax not big enough (Ktotal)");
  ASD_ABORT_IF(mesh->Ncontinuous > lvl->Kmax * lvl->Np,
               "Kmax not big enough (Ncontinuous)");
  ASD_ABORT_IF(mesh->Ncindices > lvl->Kmax * lvl->Np,
               "Kmax not big enough (Ncindices)");
  ASD_ABORT_IF(mesh->Nmortar > lvl->Kmax * lvl->Nfaces,
               "Kmax not big enough (Nmortar)");

  lvl->o_IToE.copyFrom(mesh->IToE, sizeof(iint_t) * mesh->Kintra, 0);

  lvl->o_MToE.copyFrom(mesh->MToE, sizeof(iint_t) * mesh->Kmirror, 0);

  lvl->o_UMToE.copyFrom(mesh->UMToE, sizeof(iint_t) * mesh->Kuniqmirror, 0);
  lvl->o_GToE.copyFrom(mesh->GToE, sizeof(iint_t) * mesh->Kghost, 0);

  lvl->o_EToL.copyFrom(mesh->EToL, sizeof(iint_t) * mesh->Ktotal, 0);
  lvl->o_EToT.copyFrom(mesh->EToT, sizeof(iint_t) * mesh->Ktotal, 0);
  lvl->o_EToX.copyFrom(mesh->EToX, sizeof(iint_t) * mesh->Ktotal, 0);
  lvl->o_EToY.copyFrom(mesh->EToY, sizeof(iint_t) * mesh->Ktotal, 0);
  lvl->o_EToZ.copyFrom(mesh->EToZ, sizeof(iint_t) * mesh->Ktotal, 0);
  lvl->o_EToB.copyFrom(mesh->EToB, sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,0);
  lvl->o_EToE.copyFrom(mesh->EToE, sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,0);
  lvl->o_EToF.copyFrom(mesh->EToF, sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,0);
  lvl->o_EToO.copyFrom(mesh->EToO, sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,0);
  lvl->o_EToC.copyFrom(mesh->EToC, sizeof(iint_t) * mesh->Ktotal, 0);

  if (brick) {
    lvl->o_EToP.copyFrom(mesh->EToP, sizeof(iint_t) * mesh->Ktotal, 0);
  }

  lvl->o_CToD_starts.copyFrom(mesh->CToD_starts,   sizeof(iint_t) * (mesh->Ncontinuous + 1),0);
  lvl->o_CToD_indices.copyFrom(mesh->CToD_indices, sizeof(iint_t) * mesh->Ncindices, 0);

  lvl->o_MFToEM.copyFrom(mesh->MFToEM, sizeof(iint_t) * mesh->Nmortar, 0);
  lvl->o_MFToFM.copyFrom(mesh->MFToFM, sizeof(iint_t) * mesh->Nmortar, 0);
  lvl->o_MFToEP.copyFrom(mesh->MFToEP, sizeof(iint_t) * mesh->Nmortar * P4EST_HALF,  0);
  lvl->o_MFToFP.copyFrom(mesh->MFToFP, sizeof(iint_t) * mesh->Nmortar, 0);
  lvl->o_MFToOP.copyFrom(mesh->MFToOP, sizeof(iint_t) * mesh->Nmortar, 0);

  // {{{ mirror and ghost communication information
  lvl->Nn = 0;
  for (int r = 0; r < pxest->mpisize; ++r)
    if (ghost->proc_offsets[r + 1] - ghost->proc_offsets[r] ||
        ghost->mirror_proc_offsets[r + 1] - ghost->mirror_proc_offsets[r])
      lvl->NToR[lvl->Nn++] = r;
  // }}}

  // {{{
  /* Fill metric terms */
  lvl->compute_X(mesh->Ktotal,
		 lvl->o_EToL,
		 lvl->o_EToT,
		 lvl->o_EToX,
		 lvl->o_EToY,
		 lvl->o_EToZ,
		 lvl->o_tree_to_vertex,
		 lvl->o_tree_vertices,
		 lvl->o_r,
		 lvl->o_vgeo);

  lvl->compute_cubature_X(mesh->Ktotal,
			  lvl->o_EToL,
			  lvl->o_EToT,
			  lvl->o_EToX,
			  lvl->o_EToY,
			  lvl->o_EToZ,
			  lvl->o_tree_to_vertex,
			  lvl->o_tree_vertices,
			  lvl->o_rGJ,
			  lvl->o_vgeoGJ);

  
  if (1) { // prefs->mesh_continuous){
    
    // Compute the periodic shifts for the brick case (which is the only case we
    // support periodicity
    const dfloat_t px =
      pxest->connectivity->vertices[3 * (pxest->connectivity->num_vertices - 1) + 0] -
      pxest->connectivity->vertices[0];
    const dfloat_t py =
      pxest->connectivity->vertices[3 * (pxest->connectivity->num_vertices - 1) + 1] -
      pxest->connectivity->vertices[1];
    const dfloat_t pz =
      pxest->connectivity->vertices[3 * (pxest->connectivity->num_vertices - 1) + 2] -
      pxest->connectivity->vertices[2];

    lvl->coarse_X(mesh->Ncontinuous,
                  lvl->o_CToD_starts,
		  lvl->o_CToD_indices,
		  lvl->o_EToP,
                  px,
		  py,
		  pz,
		  lvl->o_vgeo);

    lvl->interp_X(mesh->Ktotal,
		  lvl->o_EToC,
		  lvl->o_Ib,
                  lvl->o_It,
		  lvl->o_vgeo);

  }
  
  lvl->compute_geo(mesh->Ktotal, lvl->o_D, lvl->o_w, lvl->o_vgeo, lvl->o_sgeo, lvl->o_ggeo);
  // }}}
  printf("calling compute_cubature_geo\n");
  lvl->compute_cubature_geo(mesh->Ktotal, lvl->o_DGJ, lvl->o_wGJ, lvl->o_vgeoGJ, lvl->o_sgeoGJ, lvl->o_ggeoGJ);
}

void occa_p4est_topidx_to_iint(occa::device &device, size_t N,
                               p4est_topidx_t *a, occa::memory &o_ia)
{

  iint_t *ia = (iint_t *)asd_malloc_aligned(N * sizeof(iint_t));
  for (size_t n = 0; n < N; ++n)
    ia[n] = (iint_t)a[n];

  o_ia = device.malloc(N * sizeof(iint_t), ia);
  asd_free_aligned(ia);

  return;
}

void occa_double_to_dfloat(occa::device &device, size_t N, double *a,
                           occa::memory &o_ia)
{

  dfloat_t *ia = (dfloat_t *)asd_malloc_aligned(N * sizeof(dfloat_t));
  for (size_t n = 0; n < N; ++n)
    ia[n] = (dfloat_t)a[n];

  o_ia = device.malloc(N * sizeof(dfloat_t), ia);
  asd_free_aligned(ia);

  return;
}

level_t *adaptiveLevelSetup(setupAide &options, p4est_t *pxest,
			    p4est_ghost_t *ghost, occa::device &device,
			    int *brick_n, int *brick_p, int *brick_TToC,
			    int N, double occa_kmax_mem_frac, MPI_Comm &comm)
{
  level_t *lvl = new level_t[1];

  ASD_ABORT_IF(sizeof(p4est_locidx_t) > sizeof(iint_t),
               "p4est_locidx_t not compatible with iint_t");

  ASD_ABORT_IF(sizeof(p4est_topidx_t) > sizeof(iint_t),
               "p4est_topidx_t not compatible with iint_t");

  // {{{ send connectivity to the device
  {
    p4est_connectivity_t *conn = pxest->connectivity;

    occa_double_to_dfloat(device, NX * conn->num_vertices, conn->vertices,
                          lvl->o_tree_vertices);

    occa_p4est_topidx_to_iint(device, P4EST_CHILDREN * conn->num_trees,
                              conn->tree_to_vertex, lvl->o_tree_to_vertex);
  }
  // }}}

  // {{{ Build Operators

  get_operators(N, device, lvl->o_r, lvl->o_w, lvl->o_D,
                lvl->o_Ib, lvl->o_It, lvl->o_Pb, lvl->o_Pt);
  // }}}

  mesh_t *mesh = adaptiveMeshSetup(pxest, ghost, brick_n, brick_p, brick_TToC, N);

  // {{{ Mesh Constants
  level_get_mesh_constants(lvl, mesh);
  // }}}

  printf("Klocal=%d\n", lvl->Klocal);
  
  // {{{ Compute Kmax
  // FIXME?  Right now we just use a fixed Kmax for all of the mesh arrays.
  // There are some that may get allocated too big.  We may want to move to
  // some sort of dynamic resizing in the future.
  size_t available_bytes =
      (size_t)(occa_kmax_mem_frac *
               (ldouble)(device.memorySize() - device.memoryAllocated()));

  const int NpcgWork = 10;
  lvl->NpcgWork = NpcgWork;

  const int Nfaces = lvl->Nfaces;
  const int Nfp = lvl->Nfp;
  const int Np = lvl->Np;

  // We are assuming brick for storage
  const int brick = 1;

  const size_t to_allocate_bytes_per_element =
      (sizeof(iint_t) *
           (11 + (8 + P4EST_HALF) * Nfaces + brick + 2 * Np ) +
       sizeof(dfloat_t) *
           ((3 * NFIELDS + NVGEO + 2 + NpcgWork) * Np + (NSGEO * Nfaces * Nfp)));

  const size_t uKmax = available_bytes / to_allocate_bytes_per_element;
  const iint_t Kmax = lvl->Kmax = (iint_t)ASD_MIN(uKmax, IINT_MAX);
  // }}}

  printf("available bytes = %llu, to_allocation = %llu, Kmax = %d uKmax = %llu\n",
	 available_bytes, to_allocate_bytes_per_element, Kmax, uKmax);
  

  
  // {{{ Allocate Mesh Indices
  lvl->o_IToE = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_MToE = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_UMToE = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_GToE = device.malloc(sizeof(iint_t) * Kmax, NULL);

  lvl->o_EToL = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToT = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToX = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToY = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToZ = device.malloc(sizeof(iint_t) * Kmax, NULL);

  lvl->o_EToB = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToE = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToF = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToO = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);

  lvl->o_EToC = device.malloc(sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToP =
      device.malloc(sizeof(iint_t) * (brick ? Kmax : 1), NULL);
  lvl->o_EToOff = device.malloc(sizeof(iint_t) * (Kmax + 1), NULL);

  lvl->o_CToD_starts =
      device.malloc(sizeof(iint_t) * (Kmax * Np + 1), NULL);
  lvl->o_CToD_indices = device.malloc(sizeof(iint_t) * Kmax * Np, NULL);

  lvl->o_MFToEM = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToFM = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToEP =
      device.malloc(sizeof(iint_t) * Kmax * Nfaces * P4EST_HALF, NULL);
  lvl->o_MFToFP = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToOP = device.malloc(sizeof(iint_t) * Kmax * Nfaces, NULL);
  // }}}

  // {{{ Allocate Volume Fields
  lvl->o_q =
      device.malloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->o_rhsq =
      device.malloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->o_q_buf =
      device.malloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->pin_q_send = device.mappedAlloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np);
  lvl->pin_q_recv = device.mappedAlloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np);

  lvl->q_send = (dfloat_t*)lvl->pin_q_send.getMappedPointer();
  lvl->q_recv = (dfloat_t*)lvl->pin_q_recv.getMappedPointer();

  lvl->q_send_buf = (dfloat_t *)asd_malloc_aligned(NFIELDS * sizeof(dfloat_t) * Kmax * Np);
  lvl->q_recv_buf = (dfloat_t *)asd_malloc_aligned(NFIELDS * sizeof(dfloat_t) * Kmax * Np);

  lvl->NToR = (int *)asd_malloc_aligned(sizeof(int) * pxest->mpisize);
  lvl->EToA = (int8_t *)asd_malloc_aligned(sizeof(int8_t) * Kmax);

  lvl->q_send_requests =
      (MPI_Request *)asd_malloc_aligned(sizeof(MPI_Request) * pxest->mpisize);
  lvl->q_recv_requests =
      (MPI_Request *)asd_malloc_aligned(sizeof(MPI_Request) * pxest->mpisize);

  for (int r = 0; r < pxest->mpisize; ++r)
  {
    lvl->q_send_requests[r] = MPI_REQUEST_NULL;
    lvl->q_recv_requests[r] = MPI_REQUEST_NULL;
  }

  lvl->q_send_statuses =
      (MPI_Status *)asd_malloc_aligned(sizeof(MPI_Status) * pxest->mpisize);
  lvl->q_recv_statuses =
      (MPI_Status *)asd_malloc_aligned(sizeof(MPI_Status) * pxest->mpisize);
  // }}}

  // {{{ Allocate Metric Terms
  lvl->o_vgeo =						      
    device.malloc(NVGEO * sizeof(dfloat_t) * Kmax * Np, NULL);
							      
  lvl->o_sgeo =						      
    device.malloc(NSGEO * sizeof(dfloat_t) * Kmax * Nfaces * Nfp, NULL);
							      
  lvl->o_ggeo =						      
    device.malloc(NGGEO * sizeof(dfloat_t) * Kmax * Np, NULL);  // }}}

  // {{{ reduction buffers
  {
    const int LDIM = KERNEL_REDUCE_LDIM;
    const iint_t n_reduce = Kmax * Np;
    iint_t n_groups = (n_reduce + LDIM - 1) / LDIM;
    n_groups = (n_groups + 8 - 1) / 8;

    lvl->o_red_buf[0] =
        device.malloc(sizeof(dfloat_t) * n_reduce, NULL);
    lvl->o_red_buf[1] =
        device.malloc(sizeof(dfloat_t) * n_groups, NULL);
  }

  iint_t Ntotal = lvl->Klocal*lvl->Np;

  int blockSize = KERNEL_REDUCE_LDIM;
  int NthreadsUpdatePCG = KERNEL_REDUCE_LDIM;

  iint_t NblocksUpdatePCG = ASD_MIN((Ntotal+NthreadsUpdatePCG-1)/NthreadsUpdatePCG, 160);
  lvl->NblocksUpdatePCG = NblocksUpdatePCG;
  lvl->tmpNormr = (dfloat_t*) calloc(NblocksUpdatePCG, sizeof(dfloat_t));
  lvl->o_tmpNormr = device.malloc(NblocksUpdatePCG*sizeof(dfloat_t),
				  lvl->tmpNormr);
  
  int Nblock  = ASD_MAX(1,(Ntotal+blockSize-1)/blockSize);
  int Nblock2 = ASD_MAX(1,(Nblock+blockSize-1)/blockSize);
  
  lvl->Nblock = Nblock;
  lvl->Nblock2 = Nblock2;

  // WARNING: USES NFIELDS
  lvl->tmp =(dfloat_t*)
    calloc(NFIELDS * Nblock, sizeof(dfloat_t));

  lvl->o_tmp =
    device.malloc(NFIELDS * sizeof(dfloat_t) * Nblock, NULL);
  lvl->o_tmp2 =
    device.malloc(NFIELDS * sizeof(dfloat_t) * Nblock, NULL);
  
  // }}}

  lvl->o_pcgWork = new occa::memory[lvl->NpcgWork];
  dfloat *zeros = (dfloat*) calloc(lvl->Kmax*lvl->Np, sizeof(dfloat_t));
  for(int n=0;n<lvl->NpcgWork;++n){
    lvl->o_pcgWork[n] = device.malloc(lvl->Kmax*lvl->Np*sizeof(dfloat_t), zeros);
  }
  free(zeros);

  // cubature info
  int NGJ = (lvl->Nq)-1;
  int NqGJ = lvl->Nq; 2*((NGJ+1)/2);
  ldouble *lr = (ldouble *)asd_malloc_aligned(lvl->Nq * sizeof(ldouble));
  ldouble *lw = (ldouble *)asd_malloc_aligned(lvl->Nq * sizeof(ldouble));
  ldouble *lV = (ldouble *)asd_malloc_aligned(lvl->Nq * lvl->Nq * sizeof(ldouble));
  
  ldouble *lrGJ = (ldouble *)asd_malloc_aligned(NqGJ * sizeof(ldouble));
  ldouble *lwGJ = (ldouble *)asd_malloc_aligned(NqGJ * sizeof(ldouble));
  ldouble *lVGJ = (ldouble *)asd_malloc_aligned(NqGJ * NqGJ * sizeof(ldouble));
  ldouble *lDGJ = (ldouble *)asd_malloc_aligned(NqGJ * NqGJ * sizeof(ldouble));
  ldouble *lIGJ = (ldouble *)asd_malloc_aligned(NqGJ * lvl->Nq * sizeof(ldouble));


  asd_jacobi_gauss_lobatto_quadrature(0, 0, lvl->Nq-1, lr, lw);
  asd_jacobi_p_vandermonde(0, 0, lvl->Nq-1, lvl->Nq, lr, lV);
  
  asd_jacobi_gauss_quadrature(0, 0, NGJ, lrGJ, lwGJ); // TWTWTWTWTW

  for(int n=0;n<NqGJ;++n){
    printf("rGJ[%d] = %Lf\n", n, lrGJ[n]);
  }

  for(int n=0;n<NqGJ;++n){
    printf("wGJ[%d] = %Lf\n", n, lwGJ[n]);
  }

  
  asd_jacobi_p_vandermonde    (0, 0, NqGJ-1, NqGJ, lrGJ, lVGJ);
  asd_jacobi_p_differentiation(0, 0, NqGJ-1, NqGJ, lrGJ, lVGJ, lDGJ);
  asd_jacobi_p_interpolation(0, 0, lvl->N, NqGJ, lrGJ, lV, lIGJ); // interpolate from GLL to GJ

  printf("lDGJ=[\n");
  for(int n=0;n<NqGJ;++n){
    for(int m=0;m<NqGJ;++m){
      printf("%LF ", lDGJ[n*NqGJ+m]);
    }
    printf(";\n");
  }
  printf("];\n");

  printf("lIGJ=[\n");
  for(int n=0;n<NqGJ;++n){
    for(int m=0;m<lvl->Nq;++m){
      printf("%LF ", lIGJ[n*lvl->Nq+m]);
    }
    printf(";\n");
  }
  printf("];\n");


  
  dfloat *r = (dfloat *)asd_malloc_aligned(lvl->Nq * sizeof(dfloat));
  dfloat *w = (dfloat *)asd_malloc_aligned(lvl->Nq * sizeof(dfloat));
  dfloat *V = (dfloat *)asd_malloc_aligned(lvl->Nq * lvl->Nq * sizeof(dfloat));
  
  dfloat *rGJ = (dfloat *)asd_malloc_aligned(NqGJ * sizeof(dfloat));
  dfloat *wGJ = (dfloat *)asd_malloc_aligned(NqGJ * sizeof(dfloat));
  dfloat *VGJ = (dfloat *)asd_malloc_aligned(NqGJ * NqGJ * sizeof(dfloat));
  dfloat *DGJ = (dfloat *)asd_malloc_aligned(NqGJ * NqGJ * sizeof(dfloat));
  dfloat *IGJ = (dfloat *)asd_malloc_aligned(NqGJ * lvl->Nq * sizeof(dfloat));

  for(int n=0;n<NqGJ;++n){
    rGJ[n] = lrGJ[n];
    wGJ[n] = lwGJ[n];
  }
  for(int n=0;n<NqGJ;++n){
    for(int m=0;m<NqGJ;++m){
      DGJ[n+m*NqGJ] = lDGJ[n*NqGJ+m];// switch to row major
      VGJ[n+m*NqGJ] = lVGJ[n*NqGJ+m];
    }
  }
  for(int n=0;n<lvl->Nq;++n){
    for(int m=0;m<NqGJ;++m){
      IGJ[m*lvl->Nq+n] = lIGJ[n*NqGJ+m]; // WHY ?
    }
  }

  lvl->NqGJ = NqGJ;
  lvl->o_rGJ = device.malloc(NqGJ*sizeof(dfloat), rGJ);
  lvl->o_wGJ = device.malloc(NqGJ*sizeof(dfloat), wGJ);
  lvl->o_DGJ = device.malloc(NqGJ*NqGJ*sizeof(dfloat), DGJ);
  lvl->o_VGJ = device.malloc(NqGJ*NqGJ*sizeof(dfloat), VGJ);
  lvl->o_IGJ = device.malloc(NqGJ*lvl->Nq*sizeof(dfloat), IGJ);

  lvl->o_vgeoGJ =						      
    device.malloc(NVGEO * sizeof(dfloat_t) * Kmax * NqGJ * NqGJ * NqGJ, NULL);
							      
  lvl->o_sgeoGJ =						      
    device.malloc(NSGEO * sizeof(dfloat_t) * Kmax * Nfaces * NqGJ * NqGJ, NULL);

  printf("populating ggeoGJ\n");
  lvl->o_ggeoGJ =						      
    device.malloc(NGGEO * sizeof(dfloat_t) * Kmax * NqGJ * NqGJ * NqGJ, NULL);  // }}}

  
  // {{{ Build Kernels
  occa::properties info;

  level_kernelinfo(info,
		   device,
		   N,
		   NqGJ,
		   (brick_p[0] || brick_p[1] || brick_p[2]));


  lvl->dotMultiplyKernel = device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
					      "dotMultiply",					      
					      info);
  

  lvl->compute_Ax = device.buildKernel(DADAPTIVE "/okl/adaptiveAxHex3D.okl",
				       "adaptiveAxHex3D",
				       info);
  
  lvl->compute_partial_Ax = device.buildKernel(DADAPTIVE "/okl/adaptiveAxHex3D.okl",
					       "adaptivePartialAxHex3D",
					       info);

  lvl->compute_cubature_Ax = device.buildKernel(DADAPTIVE "/okl/adaptiveAxHex3D.okl",
						"adaptiveCubatureAxHex3D",
						info);
  
  lvl->compute_diagonal_Jacobi = device.buildKernel(DADAPTIVE "/okl/adaptiveBuildJacobiDiagonalContinuousHex3D.okl",
					       "adaptiveBuildJacobiDiagonalContinuousHex3D",
					       info);

  
  lvl->updatePCGKernel = device.buildKernel(DADAPTIVE "/okl/adaptiveUpdatePCG.okl",
					       "adaptiveUpdatePCG",
					       info);
  
  lvl->compute_X = device.buildKernel(DADAPTIVE "/okl/adaptiveComputeX.okl",
				      "adaptiveComputeX",
				      info);

  lvl->interp_X = device.buildKernel(DADAPTIVE "/okl/adaptiveInterpX.okl",
				     "adaptiveInterpX",
				     info);

  lvl->compute_cubature_X = device.buildKernel(DADAPTIVE "/okl/adaptiveComputeX.okl",
				      "adaptiveCubatureComputeX",
				      info);
  
  lvl->coarse_X = device.buildKernel(DADAPTIVE "/okl/adaptiveCoarseX.okl",
				     "adaptiveCoarseX",
				     info);
  
  lvl->compute_geo = device.buildKernel(DADAPTIVE "/okl/adaptiveGeometricFactorsHex3D.okl",
					"adaptiveGeometricFactorsHex3D",
					info);

  lvl->compute_cubature_geo = device.buildKernel(DADAPTIVE "/okl/adaptiveGeometricFactorsHex3D.okl",
						 "adaptiveCubatureGeometricFactorsHex3D",
						 info);
  
  lvl->get_mirror_fields = device.buildKernel(DADAPTIVE "/okl/adaptiveGetMirrorFields.okl",
					      "adaptiveGetMirrorFields",
					      info);

  lvl->set_ghost_fields = device.buildKernel(DADAPTIVE "/okl/adaptiveSetGhostFields.okl",
					      "adaptiveSetGhostFields",
					      info);

  lvl->gather_scatter = device.buildKernel(DADAPTIVE "/okl/adaptiveGatherScatter.okl",
					  "adaptiveGatherScatter",
					  info);

  lvl->gather_noncon = device.buildKernel(DADAPTIVE "/okl/adaptiveGatherNoncon.okl",
					  "adaptiveGatherNoncon",
					  info);


  lvl->scatter_noncon = device.buildKernel(DADAPTIVE "/okl/adaptiveScatterNoncon.okl",
					  "adaptiveScatterNoncon",
					  info);

  lvl->zero_children = device.buildKernel(DADAPTIVE "/okl/adaptiveZeroChildren.okl",
					  "adaptiveZeroChildren",
					  info);

  lvl->filterKernel = device.buildKernel(DADAPTIVE "/okl/adaptiveRankOneProjectionHex3D.okl",
				       "adaptiveRankOneProjectionHex3D",
				       info);


  
#if 0

  lvl->reduce_min = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "reduce_min", info, OKL_LANG);

  lvl->reduce_sum = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "reduce_sum", info, OKL_LANG);

  lvl->coarsen_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "coarsen_fields", info, OKL_LANG);

  lvl->refine_and_fill_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "refine_and_fill_fields", info, OKL_LANG);

  lvl->mortar_advection = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "mortar_advection", info, OKL_LANG);

  lvl->update_advection = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "update_advection", info, OKL_LANG);

  lvl->zero_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "zero_fields", info, OKL_LANG);


  lvl->set_ghost_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "set_ghost_fields", info, OKL_LANG);

  occaKernelInfoFree(info);
  // }}}
#endif

  level_get_mesh(lvl, mesh, pxest, ghost, device, brick);

  // BUILD occa Gather Scatter
  // TW: what is actual order ?
  // TW: what is comm ?
  // TW: need proper conversion between (p4est_gloidx_t) and (hlong)
  hlong *mult = (hlong*) calloc(lvl->Klocal*lvl->Np, sizeof(hlong));
  for(int n=0;n<lvl->Klocal*lvl->Np;++n){
    ++(mult[mesh->DToC[n]]);
  }

  printf("DToC = [\n");
  for(int e=0;e<lvl->Klocal;++e){
    printf("e=%d: ", e);
    for(int n=0;n<lvl->Np;++n){
      printf("%d ", mesh->DToC[e*lvl->Np+n]);
    }
    printf("\n");
  }
  printf("];\n");
  

  int maxMult= 0;
  for(int n=0;n<lvl->Klocal*lvl->Np;++n){
    if(maxMult<mult[n]){
      maxMult = mult[n];
    }
  }
  printf("MAXIMUM MULTIPLICITY = %d\n", maxMult);
  
  hlong *shiftDToC = (hlong*) calloc(lvl->Klocal*lvl->Np, sizeof(hlong));
  for(int n=0;n<lvl->Klocal*lvl->Np;++n){
    shiftDToC[n] = mesh->DToC[n]+1;
    if(mesh->DToC[n]<=0) printf("mesh->DToC[%d]  = %llu\n", n, mesh->DToC[n]);
  }
  lvl->ogs = ogsSetup(lvl->Klocal*lvl->Np, (hlong*) shiftDToC, comm, 1, device); // 1 verbose

  free(shiftDToC);

  // TW: generate diagonal weighting
  dfloat *ones = (dfloat*) calloc(lvl->Klocal*lvl->Np, sizeof(dfloat));
  for(int n=0;n<lvl->Klocal*lvl->Np;++n)
    ones[n] = 1;

  lvl->o_invDegree = device.malloc(lvl->Klocal*lvl->Np*sizeof(dfloat), ones);

#if USE_GASPAR==1
  lvl->zero_children(lvl->Klocal, lvl->o_EToC, lvl->o_invDegree);
#endif
  
  adaptiveGatherScatter(lvl, lvl->o_invDegree);
  //ogsGatherScatter(lvl->o_invDegree, ogsDfloat, ogsAdd, lvl->ogs);

  lvl->o_invDegree.copyTo(ones);
  
  for(int n=0;n<lvl->Klocal*lvl->Np;++n){
    if(ones[n])
      ones[n] = 1./ones[n];
  }

  lvl->o_invDegree.copyFrom(ones);

#if USE_GASPAR==1
  lvl->zero_children(lvl->Klocal, lvl->o_EToC, lvl->o_invDegree);
#endif
  
#if USE_GASPAR==0
  lvl->o_invDegree.copyFrom(lvl->ogs->o_invDegree);
#endif
  
#if 0
  lvl->o_invDegree.copyTo(ones);
  
  for(int n=0;n<lvl->Klocal*lvl->Np;++n){
    printf("invDegree[%d] = %lf\n", n, ones[n]);
  }
#endif

  // build diagonal
  dfloat lambda= 1;
  options.getArgs("LAMBDA", lambda);
  
  lvl->o_invDiagA = device.malloc(lvl->Klocal*lvl->Np*sizeof(dfloat));

  lvl->compute_diagonal_Jacobi(lvl->Klocal, lvl->o_ggeo, lvl->o_D, lambda, lvl->o_invDiagA);

  lvl->gather_noncon(lvl->Klocal, lvl->o_EToC, lvl->o_Ib, lvl->o_It, lvl->o_invDiagA);

#ifdef PLUMADG_GATHER_SCATTER
  // FIXME for MPI
  lvl->gather_scatter(lvl->Ncontinuous, lvl->o_CToD_starts, lvl->o_CToD_indices, lvl->o_invDiagA);
#else
  ogsGatherScatter(lvl->o_invDiagA, ogsDfloat, ogsAdd, lvl->ogs);
#endif

#if USE_GASPAR==1
  lvl->scatter_noncon(lvl->Klocal, lvl->o_EToC, lvl->o_Ib, lvl->o_It, lvl->o_invDiagA);
#endif

  dfloat *tmp = (dfloat*) calloc(lvl->Np*lvl->Klocal, sizeof(dfloat));
  lvl->o_invDiagA.copyTo(tmp, lvl->Np*lvl->Klocal*sizeof(dfloat), 0);

  for(dlong n=0;n<lvl->Np*lvl->Klocal;++n){
    tmp[n] = 1./tmp[n];
  }

  lvl->o_invDiagA.copyFrom(tmp, lvl->Np*lvl->Klocal*sizeof(dfloat), 0);

  //  lvl->dotMultiplyKernel(lvl->Np*lvl->Klocal, lvl->o_invDegree, lvl->o_invDiagA, lvl->o_invDiagA);

#if 0
  // quick build bubble filter (N to N-1 leaving boundary invariant)
#if 0 
  dfloat filtU[6] = { -0.000000000000000e+00,
		     -2.547057573459081e-01,
		     2.103615014726592e-01,
		      -2.103615014726592e-01,
		      2.547057573459079e-01,
		      0.000000000000000e+00};

  dfloat filtV[6] = {
    4.119429204355501e-01,
    -9.815247311448987e-01,
    1.188430384123744e+00,
    -1.188430384123744e+00,
    9.815247311448986e-01,
    -4.119429204355499e-01,
  };
#else

  dfloat filtU[8] = {
    -0.000000000000000e+00,
    -2.870464195082011e-01,
    2.255970689586477e-01,
    -2.051627322522581e-01,
    2.051627322522581e-01,
    -2.255970689586476e-01,
    2.870464195082011e-01,
  0.000000000000000e+00};

  dfloat filtV[8] = {
    2.390457218668787e-01,
    -5.806261821771470e-01,
    7.387802839637823e-01,
    -8.123632632350666e-01,
    8.123632632350675e-01,
    -7.387802839637828e-01,
    5.806261821771467e-01,
    -2.390457218668786e-01};
#endif
  

  lvl->o_filtU = device.malloc(lvl->Nq*sizeof(dfloat), filtU);
  lvl->o_filtV = device.malloc(lvl->Nq*sizeof(dfloat), filtV);
#endif

  
  mesh_free(mesh);

  return lvl;
}

void level_free(level_t *lvl)
{
  lvl->o_tree_to_vertex.free();
  lvl->o_tree_vertices.free();

  asd_free_aligned(lvl->NToR);
  asd_free_aligned(lvl->EToA);

  lvl->o_IToE.free();
  lvl->o_MToE.free();
  lvl->o_UMToE.free();
  lvl->o_GToE.free();

  
  lvl->o_EToL.free();
  lvl->o_EToT.free();
  lvl->o_EToX.free();
  lvl->o_EToY.free();
  lvl->o_EToZ.free();
  lvl->o_EToB.free();
  lvl->o_EToE.free();
  lvl->o_EToF.free();
  lvl->o_EToO.free();
  lvl->o_EToC.free();
  lvl->o_EToP.free();
  lvl->o_EToOff.free();

  lvl->o_CToD_starts.free();
  lvl->o_CToD_indices.free();

  lvl->o_MFToEM.free();
  lvl->o_MFToFM.free();
  lvl->o_MFToEP.free();
  lvl->o_MFToFP.free();
  lvl->o_MFToOP.free();


  lvl->o_r.free();
  lvl->o_w.free();
  lvl->o_D.free();

  lvl->o_Ib.free();
  lvl->o_It.free();

  lvl->o_Pb.free();
  lvl->o_Pt.free();

  lvl->o_vgeo.free();
  lvl->o_sgeo.free();
  lvl->o_q.free();
  lvl->o_rhsq.free();

#if 0
  // leak for the moment
  lvl->pin_q_send.free();
  lvl->pin_q_recv.free();
#endif
  
  lvl->o_q_buf.free();
    
  lvl->o_red_buf[0].free();
  lvl->o_red_buf[1].free();

  
  asd_free_aligned(lvl->q_send_buf);
  asd_free_aligned(lvl->q_recv_buf);
  asd_free_aligned(lvl->q_send_requests);
  asd_free_aligned(lvl->q_recv_requests);
  asd_free_aligned(lvl->q_send_statuses);
  asd_free_aligned(lvl->q_recv_statuses);

  lvl->compute_X.free();
  lvl->interp_X.free();
  lvl->coarse_X.free();
  lvl->compute_geo.free();
  lvl->compute_ics.free();
  lvl->compute_dt.free();
  lvl->compute_energy.free();
  lvl->compute_error.free();
  lvl->coarsen_fields.free();
  lvl->refine_and_fill_fields.free();
  lvl->volume_advection.free();
  lvl->mortar_advection.free();
  lvl->update_advection.free();
  lvl->zero_fields.free();
  lvl->get_mirror_fields.free();
  lvl->set_ghost_fields.free();
  lvl->reduce_min.free();
  lvl->reduce_sum.free();

  lvl->gather_noncon.free();
  lvl->scatter_noncon.free();
  lvl->gather_scatter.free();

  delete [] lvl;
}
// }}}
