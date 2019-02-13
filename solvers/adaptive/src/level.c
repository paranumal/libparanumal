#include "adaptive.h"

// {{{ Level

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
                   int N)
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

#if 0
  mesh_t *mesh = mesh_new(prefs, pxest, ghost);

  // {{{ Mesh Constants
  level_get_mesh_constants(lvl, mesh);
  // }}}

  // {{{ Compute Kmax
  // FIXME?  Right now we just use a fixed Kmax for all of the mesh arrays.
  // There are some that may get allocated too big.  We may want to move to
  // some sort of dynamic resizing in the future.
  size_t available_bytes =
      (size_t)(prefs->occa_kmax_mem_frac *
               (long double)(occaDeviceMemorySize(device) -
                             occaDeviceBytesAllocated(device)));

  const int Nfaces = lvl->Nfaces;
  const int Nfp = lvl->Nfp;
  const int Np = lvl->Np;

  const size_t to_allocate_bytes_per_element =
      (sizeof(iint_t) *
           (11 + (8 + P4EST_HALF) * Nfaces + prefs->brick + 2 * Np) +
       sizeof(dfloat_t) *
           ((3 * NFIELDS + NVGEO + 2) * Np + (NSGEO * Nfaces * Nfp)));
  const size_t uKmax = available_bytes / to_allocate_bytes_per_element;
  const iint_t Kmax = lvl->Kmax = (iint_t)ASD_MIN(uKmax, IINT_MAX);
  // }}}

  // {{{ Allocate Mesh Indices
  lvl->o_IToE = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_MToE = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_UMToE = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_GToE = device_malloc(device, sizeof(iint_t) * Kmax, NULL);

  lvl->o_EToL = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToT = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToX = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToY = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToZ = device_malloc(device, sizeof(iint_t) * Kmax, NULL);

  lvl->o_EToB = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToE = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToF = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_EToO = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);

  lvl->o_EToC = device_malloc(device, sizeof(iint_t) * Kmax, NULL);
  lvl->o_EToP =
      device_malloc(device, sizeof(iint_t) * (prefs->brick ? Kmax : 1), NULL);
  lvl->o_EToOff = device_malloc(device, sizeof(iint_t) * (Kmax + 1), NULL);

  lvl->o_CToD_starts =
      device_malloc(device, sizeof(iint_t) * (Kmax * Np + 1), NULL);
  lvl->o_CToD_indices = device_malloc(device, sizeof(iint_t) * Kmax * Np, NULL);

  lvl->o_MFToEM = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToFM = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToEP =
      device_malloc(device, sizeof(iint_t) * Kmax * Nfaces * P4EST_HALF, NULL);
  lvl->o_MFToFP = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  lvl->o_MFToOP = device_malloc(device, sizeof(iint_t) * Kmax * Nfaces, NULL);
  // }}}

  // {{{ Allocate Volume Fields
  lvl->o_q =
      device_malloc(device, NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);
  lvl->o_rhsq =
      device_malloc(device, NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->o_q_buf =
      device_malloc(device, NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->pin_q_send = occaDeviceMappedAlloc(
      device, NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);
  lvl->pin_q_recv = occaDeviceMappedAlloc(
      device, NFIELDS * sizeof(dfloat_t) * Kmax * Np, NULL);

  lvl->q_send = occaMemoryGetMappedPointer(lvl->pin_q_send);
  lvl->q_recv = occaMemoryGetMappedPointer(lvl->pin_q_recv);

  lvl->q_send_buf = asd_malloc_aligned(NFIELDS * sizeof(dfloat_t) * Kmax * Np);
  lvl->q_recv_buf = asd_malloc_aligned(NFIELDS * sizeof(dfloat_t) * Kmax * Np);

  lvl->NToR = asd_malloc_aligned(sizeof(int) * pxest->mpisize);
  lvl->EToA = asd_malloc_aligned(sizeof(int8_t) * Kmax);

  lvl->q_send_requests =
      asd_malloc_aligned(sizeof(MPI_Request) * pxest->mpisize);
  lvl->q_recv_requests =
      asd_malloc_aligned(sizeof(MPI_Request) * pxest->mpisize);

  for (int r = 0; r < pxest->mpisize; ++r)
  {
    lvl->q_send_requests[r] = MPI_REQUEST_NULL;
    lvl->q_recv_requests[r] = MPI_REQUEST_NULL;
  }

  lvl->q_send_statuses =
      asd_malloc_aligned(sizeof(MPI_Status) * pxest->mpisize);
  lvl->q_recv_statuses =
      asd_malloc_aligned(sizeof(MPI_Status) * pxest->mpisize);
  // }}}

  // {{{ Allocate Metric Terms
  lvl->o_vgeo =
      device_malloc(device, NVGEO * sizeof(dfloat_t) * Kmax * Np, NULL);
  lvl->o_sgeo = device_malloc(
      device, NSGEO * sizeof(dfloat_t) * Kmax * Nfaces * Nfp, NULL);
  // }}}

  // {{{ reduction buffers
  {
    const int LDIM = prefs->kernel_reduce_ldim;
    const iint_t n_reduce = Kmax * Np;
    iint_t n_groups = (n_reduce + LDIM - 1) / LDIM;
    n_groups = (n_groups + 8 - 1) / 8;

    lvl->o_red_buf[0] =
        device_malloc(device, sizeof(dfloat_t) * n_reduce, NULL);
    lvl->o_red_buf[1] =
        device_malloc(device, sizeof(dfloat_t) * n_groups, NULL);
  }
  // }}}

  // {{{ Build Kernels
  occaKernelInfo info = common_kernelinfo_new(prefs, device);

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

  mesh_free(mesh);
  asd_free(mesh);
#endif

  return lvl;
}

void level_free(level_t *lvl)
{
  lvl->o_tree_to_vertex.free();
  lvl->o_tree_vertices.free();

#if 0
  asd_free_aligned(lvl->NToR);
  asd_free_aligned(lvl->EToA);

  occaMemoryFree(lvl->o_IToE);
  occaMemoryFree(lvl->o_MToE);
  occaMemoryFree(lvl->o_UMToE);
  occaMemoryFree(lvl->o_GToE);

  occaMemoryFree(lvl->o_EToL);
  occaMemoryFree(lvl->o_EToT);
  occaMemoryFree(lvl->o_EToX);
  occaMemoryFree(lvl->o_EToY);
  occaMemoryFree(lvl->o_EToZ);
  occaMemoryFree(lvl->o_EToB);
  occaMemoryFree(lvl->o_EToE);
  occaMemoryFree(lvl->o_EToF);
  occaMemoryFree(lvl->o_EToO);
  occaMemoryFree(lvl->o_EToC);
  occaMemoryFree(lvl->o_EToP);
  occaMemoryFree(lvl->o_EToOff);

  occaMemoryFree(lvl->o_CToD_starts);
  occaMemoryFree(lvl->o_CToD_indices);

  occaMemoryFree(lvl->o_MFToEM);
  occaMemoryFree(lvl->o_MFToFM);
  occaMemoryFree(lvl->o_MFToEP);
  occaMemoryFree(lvl->o_MFToFP);
  occaMemoryFree(lvl->o_MFToOP);

  occaMemoryFree(lvl->o_r);
  occaMemoryFree(lvl->o_w);
  occaMemoryFree(lvl->o_D);

  occaMemoryFree(lvl->o_Ib);
  occaMemoryFree(lvl->o_It);

  occaMemoryFree(lvl->o_Pb);
  occaMemoryFree(lvl->o_Pt);

  occaMemoryFree(lvl->o_vgeo);
  occaMemoryFree(lvl->o_sgeo);
  occaMemoryFree(lvl->o_q);
  occaMemoryFree(lvl->o_rhsq);

  occaMemoryFree(lvl->o_q_buf);
  occaMemoryFree(lvl->pin_q_send);
  occaMemoryFree(lvl->pin_q_recv);

  occaMemoryFree(lvl->o_red_buf[0]);
  occaMemoryFree(lvl->o_red_buf[1]);

  asd_free_aligned(lvl->q_send_buf);
  asd_free_aligned(lvl->q_recv_buf);
  asd_free_aligned(lvl->q_send_requests);
  asd_free_aligned(lvl->q_recv_requests);
  asd_free_aligned(lvl->q_send_statuses);
  asd_free_aligned(lvl->q_recv_statuses);

  occaKernelFree(lvl->compute_X);
  occaKernelFree(lvl->interp_X);
  occaKernelFree(lvl->coarse_X);
  occaKernelFree(lvl->compute_geo);
  occaKernelFree(lvl->compute_ics);
  occaKernelFree(lvl->compute_dt);
  occaKernelFree(lvl->compute_energy);
  occaKernelFree(lvl->compute_error);
  occaKernelFree(lvl->coarsen_fields);
  occaKernelFree(lvl->refine_and_fill_fields);
  occaKernelFree(lvl->volume_advection);
  occaKernelFree(lvl->mortar_advection);
  occaKernelFree(lvl->update_advection);
  occaKernelFree(lvl->zero_fields);
  occaKernelFree(lvl->get_mirror_fields);
  occaKernelFree(lvl->set_ghost_fields);
  occaKernelFree(lvl->reduce_min);
  occaKernelFree(lvl->reduce_sum);
#endif

  delete [] lvl;
}
// }}}
