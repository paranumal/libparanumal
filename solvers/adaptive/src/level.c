#include "adaptive.h"

// {{{ Level

// {{{ Kernel Info
static void level_kernelinfo(occa::properties &info, occa::device &device, int N,
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
  info["defines/iint"] = occa_iint_name;

  const char *const dfloat =
      (sizeof(double) == sizeof(dfloat_t)) ? "double" : "float";

  info["defines/dfloat"] = dfloat;
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

  info["defines/p_FIELD_UX"] = FIELD_UX;
  info["defines/p_FIELD_UY"] = FIELD_UY;
  info["defines/p_FIELD_UZ"] = FIELD_UZ;
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
}
// }}}
#if 0
static void level_set_working_dims(level_t *lvl, prefs_t *prefs)
{
  const int Nq = lvl->Nq;
  const int Nqk = lvl->Nqk;

  const iint_t Nmortar = lvl->Nmortar;
  const iint_t Ncontinuous = lvl->Ncontinuous;

  const iint_t Klocal = lvl->Klocal;
  const iint_t Kghost = lvl->Kghost;
  const iint_t Kmirror = lvl->Kmirror;
  const iint_t Ktotal = lvl->Ktotal;

  const int KblkV = prefs->kernel_KblkV;
  int dim = 1;
  occaDim global = {1, 1, 1}, local = {1, 1, 1};
  dim = 3;
  local = (occaDim){Nq, Nq, Nqk * KblkV};

  global = (occaDim){(Ktotal + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->compute_X, dim, local, global);
  occaKernelSetWorkingDims(lvl->compute_geo, dim, local, global);
  occaKernelSetWorkingDims(lvl->compute_ics, dim, local, global);
  occaKernelSetWorkingDims(lvl->zero_fields, dim, local, global);

  global = (occaDim){(Klocal + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->compute_dt, dim, local, global);
  occaKernelSetWorkingDims(lvl->compute_energy, dim, local, global);
  occaKernelSetWorkingDims(lvl->compute_error, dim, local, global);
  occaKernelSetWorkingDims(lvl->update_advection, dim, local, global);

  global = (occaDim){(Kmirror + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->get_mirror_fields, dim, local, global);

  global = (occaDim){(Kghost + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->set_ghost_fields, dim, local, global);

  const int KblkS = prefs->kernel_KblkS;
  local = (occaDim){Nq, Nqk, KblkS};

  global = (occaDim){(Ktotal + KblkS - 1) / KblkS, 1, 1};
  occaKernelSetWorkingDims(lvl->interp_X, dim, local, global);

  global = (occaDim){(Nmortar + KblkS - 1) / KblkS, 1, 1};
  occaKernelSetWorkingDims(lvl->mortar_advection, dim, local, global);

  dim = 1;

  const int Nt = prefs->kernel_Nt;
  global = (occaDim){(Ncontinuous + Nt - 1) / Nt, 1, 1};
  local = (occaDim){Nt, 1, 1};
  occaKernelSetWorkingDims(lvl->coarse_X, dim, local, global);
}
#endif

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

#if 0

static void level_get_mesh(level_t *lvl, mesh_t *mesh, prefs_t *prefs,
                            p4est_t *pxest, p4est_ghost_t *ghost,
                            occaDevice device)
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

  device_async_ptr_to_mem(lvl->o_IToE, mesh->IToE,
                          sizeof(iint_t) * mesh->Kintra, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_MToE, mesh->MToE,
                          sizeof(iint_t) * mesh->Kmirror, occaNoOffset);

  device_async_ptr_to_mem(lvl->o_UMToE, mesh->UMToE,
                          sizeof(iint_t) * mesh->Kuniqmirror, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_GToE, mesh->GToE,
                          sizeof(iint_t) * mesh->Kghost, occaNoOffset);

  device_async_ptr_to_mem(lvl->o_EToL, mesh->EToL,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToT, mesh->EToT,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToX, mesh->EToX,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToY, mesh->EToY,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToZ, mesh->EToZ,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToB, mesh->EToB,
                          sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToE, mesh->EToE,
                          sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToF, mesh->EToF,
                          sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToO, mesh->EToO,
                          sizeof(iint_t) * mesh->Ktotal * mesh->Nfaces,
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_EToC, mesh->EToC,
                          sizeof(iint_t) * mesh->Ktotal, occaNoOffset);

  if (prefs->brick)
    device_async_ptr_to_mem(lvl->o_EToP, mesh->EToP,
                            sizeof(iint_t) * mesh->Ktotal, occaNoOffset);

  device_async_ptr_to_mem(lvl->o_CToD_starts, mesh->CToD_starts,
                          sizeof(iint_t) * (mesh->Ncontinuous + 1),
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_CToD_indices, mesh->CToD_indices,
                          sizeof(iint_t) * mesh->Ncindices, occaNoOffset);

  device_async_ptr_to_mem(lvl->o_MFToEM, mesh->MFToEM,
                          sizeof(iint_t) * mesh->Nmortar, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_MFToFM, mesh->MFToFM,
                          sizeof(iint_t) * mesh->Nmortar, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_MFToEP, mesh->MFToEP,
                          sizeof(iint_t) * mesh->Nmortar * P4EST_HALF,
                          occaNoOffset);
  device_async_ptr_to_mem(lvl->o_MFToFP, mesh->MFToFP,
                          sizeof(iint_t) * mesh->Nmortar, occaNoOffset);
  device_async_ptr_to_mem(lvl->o_MFToOP, mesh->MFToOP,
                          sizeof(iint_t) * mesh->Nmortar, occaNoOffset);

  // {{{ mirror and ghost communication information
  lvl->Nn = 0;
  for (int r = 0; r < pxest->mpisize; ++r)
    if (ghost->proc_offsets[r + 1] - ghost->proc_offsets[r] ||
        ghost->mirror_proc_offsets[r + 1] - ghost->mirror_proc_offsets[r])
      lvl->NToR[lvl->Nn++] = r;
  // }}}

  // {{{
  /* Fill metric terms */
  occaKernelRun(lvl->compute_X, occaIint(mesh->Ktotal), lvl->o_EToL,
                lvl->o_EToT, lvl->o_EToX, lvl->o_EToY, lvl->o_EToZ,
                lvl->o_tree_to_vertex, lvl->o_tree_vertices, lvl->o_r,
                lvl->o_vgeo);

  if (prefs->mesh_continuous)
  {
    if (prefs->brick &&
        (strcmp(prefs->conn_mapping, "conn_mapping_identity") != 0) &&
        (prefs->brick_p[0] || prefs->brick_p[1] || prefs->brick_p[2]))
      ASD_WARNING("Continuous node numbering for periodic levels assumes "
                  "conn_mapping does not change the periodic vertices.");

    // Compute the periodic shifts for the brick case (which is the only case we
    // support periodicity
    const dfloat_t px =
        pxest->connectivity
            ->vertices[3 * (pxest->connectivity->num_vertices - 1) + 0] -
        pxest->connectivity->vertices[0];
    const dfloat_t py =
        pxest->connectivity
            ->vertices[3 * (pxest->connectivity->num_vertices - 1) + 1] -
        pxest->connectivity->vertices[1];
    const dfloat_t pz =
        pxest->connectivity
            ->vertices[3 * (pxest->connectivity->num_vertices - 1) + 2] -
        pxest->connectivity->vertices[2];

    occaKernelRun(lvl->coarse_X, occaIint(mesh->Ncontinuous),
                  lvl->o_CToD_starts, lvl->o_CToD_indices, lvl->o_EToP,
                  occaDfloat(px), occaDfloat(py), occaDfloat(pz), lvl->o_vgeo);

    occaKernelRun(lvl->interp_X, occaIint(mesh->Ktotal), lvl->o_EToC, lvl->o_Ib,
                  lvl->o_It, lvl->o_vgeo);

    if (prefs->size > 1 && prefs->brick &&
        (prefs->brick_p[0] || prefs->brick_p[1] || prefs->brick_p[2]))
      ASD_ABORT("FIXME Continuous geometry for periodic levels requires "
                "communication of geometry.");
  }

  occaKernelRun(lvl->compute_geo, occaIint(mesh->Ktotal), lvl->o_D, lvl->o_vgeo,
                lvl->o_sgeo);
  // }}}
}
#endif

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

level_t *level_new(setupAide &options, p4est_t *pxest,
                   p4est_ghost_t *ghost, occa::device &device,
                   int *brick_n, int *brick_p, int *brick_TToC,
                   int N, double occa_kmax_mem_frac)
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

  mesh_t *mesh = mesh_new(pxest, ghost, brick_n,
      brick_p, brick_TToC, N);

  // {{{ Mesh Constants
  level_get_mesh_constants(lvl, mesh);
  // }}}

  // {{{ Compute Kmax
  // FIXME?  Right now we just use a fixed Kmax for all of the mesh arrays.
  // There are some that may get allocated too big.  We may want to move to
  // some sort of dynamic resizing in the future.
  size_t available_bytes =
      (size_t)(occa_kmax_mem_frac *
               (long double)(device.memorySize() - device.memoryAllocated()));

  const int Nfaces = lvl->Nfaces;
  const int Nfp = lvl->Nfp;
  const int Np = lvl->Np;

  // We are assuming brick for storage
  const int brick = 1;

  const size_t to_allocate_bytes_per_element =
      (sizeof(iint_t) *
           (11 + (8 + P4EST_HALF) * Nfaces + brick + 2 * Np) +
       sizeof(dfloat_t) *
           ((3 * NFIELDS + NVGEO + 2) * Np + (NSGEO * Nfaces * Nfp)));
  const size_t uKmax = available_bytes / to_allocate_bytes_per_element;
  const iint_t Kmax = lvl->Kmax = (iint_t)ASD_MIN(uKmax, IINT_MAX);
  // }}}

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

  lvl->pin_q_send = device.mappedAlloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);
  lvl->pin_q_recv = device.mappedAlloc(NFIELDS * sizeof(dfloat_t) * Kmax * Np,
      NULL);

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
  lvl->o_sgeo = device.malloc(
      NSGEO * sizeof(dfloat_t) * Kmax * Nfaces * Nfp, NULL);
  // }}}

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
  // }}}

  // {{{ Build Kernels
  occa::properties info;

  level_kernelinfo(info, device, N, (brick_p[0] || brick_p[1] ||
        brick_p[2]));

#if 0
  lvl->compute_X = occaDeviceBuildKernelFromString(device, prefs->kernels,
                                                   "compute_X", info, OKL_LANG);

  lvl->interp_X = occaDeviceBuildKernelFromString(device, prefs->kernels,
                                                  "interp_X", info, OKL_LANG);

  lvl->coarse_X = occaDeviceBuildKernelFromString(device, prefs->kernels,
                                                  "coarse_X", info, OKL_LANG);

  lvl->compute_geo = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "compute_geo", info, OKL_LANG);

  lvl->compute_ics = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "compute_ics", info, OKL_LANG);

  lvl->compute_dt = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "compute_dt", info, OKL_LANG);

  lvl->compute_energy = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "compute_energy", info, OKL_LANG);

  lvl->compute_error = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "compute_error", info, OKL_LANG);

  lvl->reduce_min = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "reduce_min", info, OKL_LANG);

  lvl->reduce_sum = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "reduce_sum", info, OKL_LANG);

  lvl->coarsen_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "coarsen_fields", info, OKL_LANG);

  lvl->refine_and_fill_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "refine_and_fill_fields", info, OKL_LANG);

  lvl->volume_advection = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "volume_advection", info, OKL_LANG);

  lvl->mortar_advection = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "mortar_advection", info, OKL_LANG);

  lvl->update_advection = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "update_advection", info, OKL_LANG);

  lvl->zero_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "zero_fields", info, OKL_LANG);

  lvl->get_mirror_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "get_mirror_fields", info, OKL_LANG);

  lvl->set_ghost_fields = occaDeviceBuildKernelFromString(
      device, prefs->kernels, "set_ghost_fields", info, OKL_LANG);

  occaKernelInfoFree(info);
  // }}}

  level_set_working_dims(lvl, prefs);
  level_get_mesh(lvl, mesh, prefs, pxest, ghost, device);
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

  lvl->o_q_buf.free();
  lvl->pin_q_send.free();
  lvl->pin_q_recv.free();

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

  delete [] lvl;
}
// }}}