//
// plumadg
//




// {{{ OCCA
/* OCCA modes */
#define SERIAL 0
#define OPENMP 1
#define OPENCL 2
#define CUDA 3
const char *const occa_modes[] = {"SERIAL", "OpenMP", "OpenCL", "CUDA", NULL};

/* OCCA language */
#define OKL_LANG (1 << 0)
#define OFL_LANG (1 << 1)
#define NATIVE_LANG (1 << 2)

static int get_occa_mode(const char *info)
{
  int mode = -1;
  if (strstr(info, "Serial"))
    mode = SERIAL;
  else if (strstr(info, "OpenMP"))
    mode = OPENMP;
  else if (strstr(info, "OpenCL"))
    mode = OPENCL;
  else if (strstr(info, "CUDA"))
    mode = CUDA;

  return mode;
}

#define DEVICE_MEMFRAC 0.9 // fraction of device memory to use

static void device_async_ptr_to_mem(occaMemory dest, void *src, size_t bytes,
                                    size_t offset)
{
  if (bytes > 0)
    occaAsyncCopyPtrToMem(dest, src, bytes, offset);
}

static void device_async_mem_to_ptr(void *dest, occaMemory src, size_t bytes,
                                    size_t offset)
{
  if (bytes > 0)
    occaAsyncCopyMemToPtr(dest, src, bytes, offset);
}

static occaMemory device_malloc(occaDevice device, size_t bytecount, void *src)
{
  bytecount = ASD_MAX(1, bytecount);
  uintmax_t bytes = occaDeviceBytesAllocated(device);
  uintmax_t total_bytes = occaDeviceMemorySize(device);
  ASD_ABORT_IF(
      (double)(bytes + bytecount) > DEVICE_MEMFRAC * (double)total_bytes,
      "Over memory limit: \n"
      "      current: allocated %ju (%.2f GiB) out of %ju (%.2f GiB)\n"
      "      new val: allocated %ju (%.2f GiB) out of %ju (%.2f GiB)\n"
      "      (fudge factor is: %.2f)",
      bytes, ((double)bytes) / GiB, total_bytes, ((double)total_bytes) / GiB,
      bytes + bytecount, ((double)(bytes + bytecount)) / GiB, total_bytes,
      ((double)total_bytes) / GiB, DEVICE_MEMFRAC);
  if ((double)(bytes + bytecount) > 0.9 * DEVICE_MEMFRAC * (double)total_bytes)
    ASD_WARNING(
        "At 90%% of memory limit: \n"
        "      current: allocated %ju (%.2f GiB) out of %ju (%.2f GiB)\n"
        "      new val: allocated %ju (%.2f GiB) out of %ju (%.2f GiB)\n"
        "      (fudge factor is: %.2f)",
        bytes, ((double)bytes) / GiB, total_bytes, ((double)total_bytes) / GiB,
        bytes + bytecount, ((double)(bytes + bytecount)) / GiB, total_bytes,
        ((double)total_bytes) / GiB, DEVICE_MEMFRAC);

  return occaDeviceMalloc(device, bytecount, src);
}
// }}}


// {{{ Utilities
static void debug(MPI_Comm comm)
{
  int rank;
  ASD_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  /*
   * This snippet of code is used for parallel debugging.
   *
   * You then need to launch a fresh tmux session which is done by just
   * typing tmux at the command prompt. Next launch your code as usual
   *
   *     mpirun -np 4 ./debug_mpi
   *
   * You can run the following tmux commands to join and synchronise the
   * input to the running processes.
   *
   *     join-pane -s 2
   *     join-pane -s 3
   *     join-pane -s 4
   *     setw synchronize-panes
   */
  char command[ASD_BUFSIZ];
  snprintf(command, ASD_BUFSIZ, "tmux new-window -t %d 'lldb -p %d'", rank + 1,
           getpid());
  printf("command: %s\n", command);
  system(command);

  {
    int pause = 1;
    while (pause)
    {
      /* the code will wait in this loop until the debugger is attached */
    }
  }
}

// }}}

// {{{ Data <-> Device Utilities
static occaMemory occa_p4est_topidx_to_iint(occaDevice device, size_t N,
                                            p4est_topidx_t *a)
{

  iint_t *ia = asd_malloc_aligned(N * sizeof(iint_t));
  for (size_t n = 0; n < N; ++n)
    ia[n] = (iint_t)a[n];

  occaMemory o_ia = device_malloc(device, N * sizeof(iint_t), ia);
  asd_free_aligned(ia);
  return o_ia;
}

static occaMemory occa_double_to_dfloat(occaDevice device, size_t N, double *a)
{

  dfloat_t *ia = asd_malloc_aligned(N * sizeof(dfloat_t));
  for (size_t n = 0; n < N; ++n)
    ia[n] = (dfloat_t)a[n];

  occaMemory o_ia = device_malloc(device, N * sizeof(dfloat_t), ia);
  asd_free_aligned(ia);
  return o_ia;
}
// }}}

// {{{ Operators
static void get_operators(int N, occaDevice device, occaMemory *o_r,
                          occaMemory *o_w, occaMemory *o_D, occaMemory *o_Ib,
                          occaMemory *o_It, occaMemory *o_Pb, occaMemory *o_Pt)
{
  const int Nq = N + 1;

  long double *lr = asd_malloc_aligned(Nq * sizeof(long double));
  long double *lw = asd_malloc_aligned(Nq * sizeof(long double));
  long double *lV = asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lD = asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lM = asd_malloc_aligned(Nq * Nq * sizeof(long double));

  long double *lrb = asd_malloc_aligned(Nq * sizeof(long double));
  long double *lrt = asd_malloc_aligned(Nq * sizeof(long double));

  long double *lIb = asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lIt = asd_malloc_aligned(Nq * Nq * sizeof(long double));

  long double *lPb = asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lPt = asd_malloc_aligned(Nq * Nq * sizeof(long double));

  dfloat_t *I = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  for (int n = 0; n < Nq * Nq; ++n)
    I[n] = 0;
  for (int n = 0; n < Nq; ++n)
    I[n + n * Nq] = 1;

  asd_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
  asd_jacobi_p_vandermonde(0, 0, N, Nq, lr, lV);
  asd_jacobi_p_differentiation(0, 0, N, Nq, lr, lV, lD);
  asd_jacobi_p_mass(N, lV, lM);

  for (int n = 0; n < Nq; ++n)
  {
    lrb[n] = (lr[n] - 1) / 2;
    lrt[n] = (lr[n] + 1) / 2;
  }

  asd_jacobi_p_interpolation(0, 0, N, Nq, lrb, lV, lIb);
  asd_jacobi_p_interpolation(0, 0, N, Nq, lrt, lV, lIt);

  asd_jacobi_p_h_project(N, 0.5L, lV, lIb, lM, lPb);
  asd_jacobi_p_h_project(N, 0.5L, lV, lIt, lM, lPt);

  dfloat_t *r = asd_malloc_aligned(Nq * sizeof(dfloat_t));
  dfloat_t *w = asd_malloc_aligned(Nq * sizeof(dfloat_t));

  for (int n = 0; n < Nq; ++n)
  {
    r[n] = (dfloat_t)lr[n];
    w[n] = (dfloat_t)lw[n];
  }

  *o_r = device_malloc(device, Nq * sizeof(dfloat_t), r);
  *o_w = device_malloc(device, Nq * sizeof(dfloat_t), w);

  dfloat_t *DT = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
    for (int j = 0; j < Nq; ++j)
      DT[j * Nq + i] = (dfloat_t)lD[i * Nq + j];

  *o_D = device_malloc(device, Nq * Nq * sizeof(dfloat_t), DT);

  dfloat_t *Ib = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  dfloat_t *It = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
  {
    for (int j = 0; j < Nq; ++j)
    {
      Ib[j * Nq + i] = (dfloat_t)lIb[i * Nq + j];
      It[j * Nq + i] = (dfloat_t)lIt[i * Nq + j];
    }
  }

  *o_Ib = device_malloc(device, Nq * Nq * sizeof(dfloat_t), Ib);
  *o_It = device_malloc(device, Nq * Nq * sizeof(dfloat_t), It);

  dfloat_t *Pb = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  dfloat_t *Pt = asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
  {
    for (int j = 0; j < Nq; ++j)
    {
      Pb[j * Nq + i] = (dfloat_t)lPb[i * Nq + j];
      Pt[j * Nq + i] = (dfloat_t)lPt[i * Nq + j];
    }
  }

  *o_Pb = device_malloc(device, Nq * Nq * sizeof(dfloat_t), Pb);
  *o_Pt = device_malloc(device, Nq * Nq * sizeof(dfloat_t), Pt);

  asd_free_aligned(lr);
  asd_free_aligned(lw);
  asd_free_aligned(lV);
  asd_free_aligned(lD);
  asd_free_aligned(lM);
  asd_free_aligned(lrb);
  asd_free_aligned(lrt);
  asd_free_aligned(lIb);
  asd_free_aligned(lIt);
  asd_free_aligned(lPb);
  asd_free_aligned(lPt);

  asd_free_aligned(r);
  asd_free_aligned(w);
  asd_free_aligned(I);
  asd_free_aligned(DT);
  asd_free_aligned(Ib);
  asd_free_aligned(It);
  asd_free_aligned(Pb);
  asd_free_aligned(Pt);
}
// }}}

// {{{ Preferences
typedef struct prefs
{
  lua_State *L;

  MPI_Comm comm;
  int size;     // MPI comm size
  int rank;     // MPI comm rank
  int hostrank; // rank of process on a given host
  int loglevel; // asd log level for output

  char *occa_info;
  char *occa_flags;
  int occa_mode;
  double occa_kmax_mem_frac;

  int brick;
  char *conn_name;
  char *conn_mapping;
  int conn_vertex_mapping;

  dfloat_t shell_R1;
  dfloat_t shell_R2;

  int brick_n[3];
  int brick_p[3];
  p4est_topidx_t *brick_TToC; // tree id to cartesian coordinates

  int mesh_initial_refinement;
  int mesh_static_refinement; // number of levels of static refinement
  int mesh_start_level;       // refinement level of initial mesh
  int mesh_N;                 // polynomial order of initial mesh
  int mesh_continuous;        // make the mesh continuous

  char *kernel_include;       // Text to include with all of the kernels built
  char *kernel_bcs_advection; // Name of advection boundary condition function
  char *kernel_ics;           // Name of the initial condition function
  char *kernel_exact;         // Name of the exact solution function
  int kernel_KblkV; // Number of elements to process per workgroup in volume
  int kernel_KblkS; // Number of elements to process per workgroup in surface
  int kernel_Nt;    // Number of threads to use per workgroup for 1D kernels
  int kernel_reduce_ldim;     // local dimension for reductions
  int kernel_reduce_max_copy; // max number of elements to copy

  char *kernels;

  char *output_datadir; // directory for output data files
  char *output_prefix;  // prefix for output files

  int output_mfem; // write mfem files?
  int output_vtk;  // write vtk files?

  int output_vtk_binary;   // write binary vtk files?
  int output_vtk_compress; // write compressed vtk files?
} prefs_t;

static prefs_t *prefs_new(const char *filename, MPI_Comm comm, int loglevel)
{
  prefs_t *prefs = asd_malloc(sizeof(prefs_t));

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  ASD_ASSERT(lua_gettop(L) == 0);

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  int hostrank = asd_get_host_rank(comm);

  // set constants for config file
  lua_pushnumber(L, (lua_Number)rank);
  lua_setglobal(L, "MPI_RANK");
  lua_pushnumber(L, (lua_Number)size);
  lua_setglobal(L, "MPI_SIZE");
  lua_pushnumber(L, (lua_Number)hostrank);
  lua_setglobal(L, "HOST_RANK");

  for (int k = 0; BCs[k]; k++)
  {
    lua_pushnumber(L, (lua_Number)k);
    char name[ASD_BUFSIZ];
    snprintf(name, ASD_BUFSIZ, "BC_%s", BCs[k]);
    lua_setglobal(L, name);
  }

  lua_pushnumber(L, DIM);
  lua_setglobal(L, "DIM");

  ASD_ASSERT(lua_gettop(L) == 0);

  // evaluate config file
  if (luaL_loadfile(L, filename) || lua_pcall(L, 0, 0, 0))
    ASD_LERROR("cannot run configuration file: `%s'", lua_tostring(L, -1));

  // occa
  prefs->occa_info = asd_lua_expr_string(L, "app.occa.info", "mode = Serial");

  int have_occa_flags = asd_lua_expr_boolean(L, "app.occa.flags ~= nil", 0);
  prefs->occa_flags =
      (have_occa_flags) ? asd_lua_expr_string(L, "app.occa.flags", "") : NULL;

  prefs->occa_mode = get_occa_mode(prefs->occa_info);
  prefs->occa_kmax_mem_frac = asd_lua_expr_number(
      L, "app.occa.kmax_mem_frac", (prefs->occa_mode == SERIAL) ? 0.1 : 0.9);

  // connectivity
  prefs->brick =
      asd_lua_expr_boolean(L, "app.conn == nil or app.conn.name == nil", 0);
  prefs->conn_name =
      (prefs->brick) ? NULL : asd_lua_expr_string(L, "app.conn.name", "");

  const char *default_conn_mapping = "conn_mapping_identity";
  if (prefs->conn_name && strcmp(prefs->conn_name, "shell") == 0)
    default_conn_mapping = "conn_mapping_shell";
  prefs->conn_mapping =
      asd_lua_expr_string(L, "app.conn.mapping", default_conn_mapping);
  prefs->conn_vertex_mapping =
      asd_lua_expr_boolean(L, "app.conn.vertex_mapping and true or false", 0);

  prefs->shell_R1 = (dfloat_t)asd_lua_expr_number(L, "app.conn.shell.R1", 1);
  prefs->shell_R2 = (dfloat_t)asd_lua_expr_number(L, "app.conn.shell.R2", 2);

  prefs->brick_n[0] = (int)asd_lua_expr_integer(L, "app.conn.brick.n[1]", 1);
  prefs->brick_n[1] = (int)asd_lua_expr_integer(L, "app.conn.brick.n[2]", 1);
  prefs->brick_n[2] = (int)asd_lua_expr_integer(L, "app.conn.brick.n[3]", 1);

  prefs->brick_p[0] = (int)asd_lua_expr_boolean(L, "app.conn.brick.p[1]", 0);
  prefs->brick_p[1] = (int)asd_lua_expr_boolean(L, "app.conn.brick.p[2]", 0);
  prefs->brick_p[2] = (int)asd_lua_expr_boolean(L, "app.conn.brick.p[3]", 0);

  // mesh
  prefs->mesh_initial_refinement = asd_lua_expr_boolean(
      L, "app.mesh.initial_refinement and true or false", 0);
  prefs->mesh_start_level =
      (int)asd_lua_expr_integer(L, "app.mesh.start_level", 0);
  prefs->mesh_static_refinement =
      (int)asd_lua_expr_integer(L, "app.mesh.static_refinement", 0);
  prefs->mesh_N = (int)asd_lua_expr_integer(L, "app.mesh.N", 3);
  prefs->mesh_continuous = asd_lua_expr_boolean(L, "app.mesh.continuous", 1);

  // kernel
  prefs->kernel_include = asd_lua_expr_string(L, "app.kernel.include", "");
  prefs->kernel_bcs_advection =
      asd_lua_expr_string(L, "app.kernel.bcs.advection", "bcs_zero");
  prefs->kernel_ics = asd_lua_expr_string(L, "app.kernel.ics", "ics_zero");
  prefs->kernel_exact = asd_lua_expr_string(L, "app.kernel.exact", "ics_zero");
  prefs->kernel_KblkV = (int)asd_lua_expr_integer(L, "app.kernel.KblkV", 2);
  prefs->kernel_KblkS = (int)asd_lua_expr_integer(L, "app.kernel.KblkS", 8);
  prefs->kernel_Nt = (int)asd_lua_expr_integer(L, "app.kernel.Nt", 256);
  prefs->kernel_reduce_ldim =
      (int)asd_lua_expr_integer(L, "app.kernel.reduce.ldim", 256);
  prefs->kernel_reduce_max_copy =
      (int)asd_lua_expr_integer(L, "app.kernel.reduce.max_copy", 1024);

  // output
  prefs->output_datadir = asd_lua_expr_string(L, "app.output.datadir", ".");

  {
    struct stat sb;
    if (stat(prefs->output_datadir, &sb) != 0 &&
        mkdir(prefs->output_datadir, 0755) != 0 && errno != EEXIST)
      perror("making datadir");
  }

  prefs->output_prefix = asd_lua_expr_string(L, "app.output.prefix", APP_NAME);
  prefs->output_mfem = asd_lua_expr_boolean(L, "app.output.mfem", 1);

  // vtk
  prefs->output_vtk = asd_lua_expr_boolean(L, "app.output.vtk ~= nil", 0);
  prefs->output_vtk =
      asd_lua_expr_boolean(L, "app.output.vtk", prefs->output_vtk);
  prefs->output_vtk_binary =
      asd_lua_expr_boolean(L, "app.output.vtk.binary", 1);
  prefs->output_vtk_compress =
      asd_lua_expr_boolean(L, "app.output.vtk.compress", 1);

  ASD_ASSERT(lua_gettop(L) == 0);

  prefs->comm = comm;
  prefs->size = size;
  prefs->rank = rank;
  prefs->hostrank = hostrank;

  prefs->loglevel = loglevel;

  prefs->L = L;

  // read kernels from disk
  const size_t ilen = strlen(prefs->kernel_include);

  size_t klen;
  char *kernels = asd_read_file("okl/kernels.okl", &klen);
  prefs->kernels = asd_malloc(ilen + klen + 2);
  strcpy(prefs->kernels, prefs->kernel_include);
  strcat(prefs->kernels, "\n");
  strcat(prefs->kernels, kernels);
  asd_free(kernels);

  return prefs;
}

static void prefs_free(prefs_t *prefs)
{
  lua_close(prefs->L);
  asd_free(prefs->occa_info);
  asd_free(prefs->occa_flags);

  if (!prefs->brick)
    asd_free(prefs->conn_name);
  asd_free(prefs->conn_mapping);

  asd_free_aligned(prefs->brick_TToC);

  asd_free(prefs->kernel_include);
  asd_free(prefs->kernel_bcs_advection);
  asd_free(prefs->kernel_ics);
  asd_free(prefs->kernel_exact);
  asd_free(prefs->kernels);

  asd_free(prefs->output_datadir);
  asd_free(prefs->output_prefix);
}

static void prefs_print(prefs_t *prefs)
{
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("----- Preferences Read -----------------------------------");
  ASD_ROOT_INFO("  occa_info         = \"%s\"", prefs->occa_info);
  if (prefs->occa_flags)
    ASD_ROOT_INFO("  occa_flags        = \"%s\"", prefs->occa_flags);
  ASD_ROOT_INFO("  occa_mode         = \"%s\"", occa_modes[prefs->occa_mode]);
  ASD_ROOT_INFO("  occa_kmax_mem_frac = %g", prefs->occa_kmax_mem_frac);

  ASD_ROOT_INFO("");

  if (prefs->conn_name)
    ASD_ROOT_INFO("  conn_name    = \"%s\"", prefs->conn_name);
  else
  {
    ASD_ROOT_INFO("  conn_name    = brick");
#if DIM == 2
    ASD_ROOT_INFO("  brick_n      = {%d %d}", prefs->brick_n[0],
                  prefs->brick_n[1]);
    ASD_ROOT_INFO("  brick_p      = {%d %d}", prefs->brick_p[0],
                  prefs->brick_p[1]);
#else
    ASD_ROOT_INFO("  brick_n      = {%d %d %d}", prefs->brick_n[0],
                  prefs->brick_n[1], prefs->brick_n[2]);
    ASD_ROOT_INFO("  brick_p      = {%d %d %d}", prefs->brick_p[0],
                  prefs->brick_p[1], prefs->brick_p[2]);
#endif
  }
  if (prefs->conn_name && strcmp(prefs->conn_name, "shell") == 0)
  {
    ASD_ROOT_INFO("  shell_R1 = %f", prefs->shell_R1);
    ASD_ROOT_INFO("  shell_R2 = %f", prefs->shell_R2);
  }
  ASD_ROOT_INFO("  conn_mapping = \"%s\"", prefs->conn_mapping);
  ASD_ROOT_INFO("  conn_vertex_mapping = %d", prefs->conn_vertex_mapping);

  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  mesh_initial_refinement = %d",
                prefs->mesh_initial_refinement);
  ASD_ROOT_INFO("  mesh_static_refinement  = %d",
                prefs->mesh_static_refinement);
  ASD_ROOT_INFO("  mesh_start_level        = %d", prefs->mesh_start_level);
  ASD_ROOT_INFO("  mesh_N                  = %d", prefs->mesh_N);
  ASD_ROOT_INFO("  mesh_continuous         = %d", prefs->mesh_continuous);
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  kernel_include       = \"%s\"", prefs->kernel_include);
  ASD_ROOT_INFO("  kernel_bcs_advection = \"%s\"", prefs->kernel_bcs_advection);
  ASD_ROOT_INFO("  kernel_ics           = \"%s\"", prefs->kernel_ics);
  ASD_ROOT_INFO("  kernel_exact         = \"%s\"", prefs->kernel_exact);
  ASD_ROOT_INFO("  kernel_KblkV         = %d", prefs->kernel_KblkV);
  ASD_ROOT_INFO("  kernel_KblkS         = %d", prefs->kernel_KblkS);
  ASD_ROOT_INFO("  kernel_Nt            = %d", prefs->kernel_Nt);
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  kernel_reduce_ldim     = %d", prefs->kernel_reduce_ldim);
  ASD_ROOT_INFO("  kernel_reduce_max_copy = %d", prefs->kernel_reduce_max_copy);
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  output_datadir      = %s", prefs->output_datadir);
  ASD_ROOT_INFO("  output_prefix       = %s", prefs->output_prefix);
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  output_mfem = %d", prefs->output_mfem);
  ASD_ROOT_INFO("  output_vtk  = %d", prefs->output_vtk);
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("  output_vtk_binary   = %d", prefs->output_vtk_binary);
  ASD_ROOT_INFO("  output_vtk_compress = %d", prefs->output_vtk_compress);
  ASD_ROOT_INFO("----------------------------------------------------------");
}
// }}}

// {{{ Connectivity
/** Get a new connectivity
 *
 * \param [in]  name    p4est name or file name of connectivity;
 *                      or NULL for brick connectivity
 * \param [in]  n[DIM]  number (per dimension) of quads
 * \param [in]  p[DIM]  Boolean (per dimension) indicating periodicity of brick
 *
 * \return Initialized app structure
 *
 */
static p4est_connectivity_t *get_connectivity(char *name, int *n, int *p)
{
  p4est_connectivity_t *conn = NULL;

  if (name)
  {
    FILE *stream = fopen(name, "r");
    if (stream)
    {
      if (fclose(stream))
        ASD_ABORT("Failed fclose on: %s", name);

      conn = p4est_connectivity_read_inp(name);

      if (conn == NULL)
        ASD_ABORT("Failed to read a valid connectivity from \"%s\"", name);
    }
    else
    {
      char elename[ASD_BUFSIZ];
      snprintf(elename, ASD_BUFSIZ, "%s.ele", name);
      FILE *elestream = (DIM == 3) ? fopen(elename, "r") : NULL;
      if (elestream)
      {
        if (fclose(elestream))
          ASD_ABORT("Failed fclose on: %s", elename);

#if DIM == 3
        p8est_tets_t *tets;

        tets = p8est_tets_read(name);
        if (tets == NULL)
          ASD_ABORT("Failed to read tets: %s", name);

        p8est_tets_make_righthanded(tets);
        conn = p8est_connectivity_new_tets(tets);
        p8est_tets_destroy(tets);

        if (conn == NULL)
          ASD_ABORT("Failed get valid tet based connectivity \"%s\"", name);
#endif
      }
      else
      {

        conn = p4est_connectivity_new_byname(name);

        if (conn == NULL)
          ASD_ABORT("Failed get valid pxest connectivity \"%s\"", name);
      }
    }
  }
  else
  {
    conn = p4est_connectivity_new_brick(n[0], n[1],
#if DIM == 3
                                        n[2],
#endif
                                        p[0], p[1]
#if DIM == 3
                                        ,
                                        p[2]
#endif
                                        );
  }

  return conn;
}
// }}}

// {{{ Kernels
static occaKernelInfo common_kernelinfo_new(prefs_t *prefs, occaDevice device)
{
  const int N = prefs->mesh_N;

  const int Nq = N + 1;
#if DIM == 3
  const int Nqk = Nq;
  const int Np = Nq * Nq * Nq;
  const int Nfp = Nq * Nq;
#else
  const int Nqk = 1;
  const int Np = Nq * Nq;
  const int Nfp = Nq;
#endif
  const int Nfaces = 2 * DIM;

  occaKernelInfo info = occaCreateKernelInfo();

  occaKernelInfoAddParserFlag(info, "automate-add-barriers", "disabled");

  occaKernelInfoAddDefine(info, "iint", occaString(occa_iint_name));
  const char *const dfloat =
      (sizeof(double) == sizeof(dfloat_t)) ? "double" : "float";
  occaKernelInfoAddDefine(info, "dfloat", occaString(dfloat));
  if (sizeof(double) == sizeof(dfloat_t))
    occaKernelInfoAddDefine(info, "p_DFLOAT_DOUBLE", occaInt(1));
  else
    occaKernelInfoAddDefine(info, "p_DFLOAT_FLOAT", occaInt(1));

  occaKernelInfoAddDefine(info, "p_DFLOAT_MAX", occaDfloat(DFLOAT_MAX));

  occaKernelInfoAddDefine(info, "p_KblkV", occaIint(prefs->kernel_KblkV));
  occaKernelInfoAddDefine(info, "p_KblkS", occaIint(prefs->kernel_KblkS));
  occaKernelInfoAddDefine(info, "p_Nt", occaIint(prefs->kernel_Nt));

  occaKernelInfoAddDefine(info, "p_REDUCE_LDIM",
                          occaIint(prefs->kernel_reduce_ldim));

  occaKernelInfoAddDefine(info, "p_DIM", occaIint(DIM));

  ASD_ASSERT(sizeof(p4est_qcoord_t) <= sizeof(iint_t));
  occaKernelInfoAddDefine(info, "p_P4EST_ROOT_LEN", occaIint(P4EST_ROOT_LEN));
  occaKernelInfoAddDefine(info, "p_P4EST_HALF", occaIint(P4EST_HALF));
  occaKernelInfoAddDefine(info, "p_P4EST_FACES", occaIint(P4EST_FACES));
  occaKernelInfoAddDefine(info, "p_P4EST_EDGES", occaIint(P4EST_EDGES));

  occaKernelInfoAddDefine(info, "p_NX", occaIint(NX));

  if (prefs->brick &&
      (prefs->brick_p[0] || prefs->brick_p[1] || prefs->brick_p[2]))
    occaKernelInfoAddDefine(info, "p_PERIODIC_BRICK", occaIint(1));

  occaKernelInfoAddDefine(info, "p_FIELD_UX", occaIint(FIELD_UX));
  occaKernelInfoAddDefine(info, "p_FIELD_UY", occaIint(FIELD_UY));
  occaKernelInfoAddDefine(info, "p_FIELD_UZ", occaIint(FIELD_UZ));
  occaKernelInfoAddDefine(info, "p_FIELD_P", occaIint(FIELD_P));
  occaKernelInfoAddDefine(info, "p_FIELD_THETA", occaIint(FIELD_THETA));
  occaKernelInfoAddDefine(info, "p_NFIELDS", occaIint(NFIELDS));

  occaKernelInfoAddDefine(info, "p_VGEO_X", occaIint(VGEO_X));
  occaKernelInfoAddDefine(info, "p_VGEO_Y", occaIint(VGEO_Y));
  occaKernelInfoAddDefine(info, "p_VGEO_RX", occaIint(VGEO_RX));
  occaKernelInfoAddDefine(info, "p_VGEO_SX", occaIint(VGEO_SX));
  occaKernelInfoAddDefine(info, "p_VGEO_RY", occaIint(VGEO_RY));
  occaKernelInfoAddDefine(info, "p_VGEO_SY", occaIint(VGEO_SY));
  occaKernelInfoAddDefine(info, "p_VGEO_J", occaIint(VGEO_J));
  occaKernelInfoAddDefine(info, "p_NVGEO", occaIint(NVGEO));

#if DIM == 3
  occaKernelInfoAddDefine(info, "p_VGEO_Z", occaIint(VGEO_Z));
  occaKernelInfoAddDefine(info, "p_VGEO_TX", occaIint(VGEO_TX));
  occaKernelInfoAddDefine(info, "p_VGEO_TY", occaIint(VGEO_TY));
  occaKernelInfoAddDefine(info, "p_VGEO_RZ", occaIint(VGEO_RZ));
  occaKernelInfoAddDefine(info, "p_VGEO_SZ", occaIint(VGEO_SZ));
  occaKernelInfoAddDefine(info, "p_VGEO_TZ", occaIint(VGEO_TZ));
#endif

  occaKernelInfoAddDefine(info, "p_SGEO_NX", occaIint(SGEO_NX));
  occaKernelInfoAddDefine(info, "p_SGEO_NY", occaIint(SGEO_NY));
#if DIM == 3
  occaKernelInfoAddDefine(info, "p_SGEO_NZ", occaIint(SGEO_NZ));
#endif
  occaKernelInfoAddDefine(info, "p_SGEO_SJ", occaIint(SGEO_SJ));
  occaKernelInfoAddDefine(info, "p_NSGEO", occaIint(NSGEO));

  occaKernelInfoAddDefine(info, "p_shell_R1", occaDfloat(prefs->shell_R1));
  occaKernelInfoAddDefine(info, "p_shell_R2", occaDfloat(prefs->shell_R2));

  occaKernelInfoAddDefine(info, "p_BC_SKIP", occaIint(BC_SKIP));
  occaKernelInfoAddDefine(info, "p_BC_NONE", occaIint(BC_NONE));
  occaKernelInfoAddDefine(info, "p_BC_DEFAULT", occaIint(BC_DEFAULT));

  occaKernelInfoAddDefine(info, "p_M_PI", occaDouble(M_PI));
  occaKernelInfoAddDefine(info, "p_M_PI_4", occaDouble(M_PI_4));

  occaKernelInfoAddDefine(info, "conn_mapping",
                          occaString(prefs->conn_mapping));

  occaKernelInfoAddDefine(info, "app_bcs_advection",
                          occaString(prefs->kernel_bcs_advection));
  occaKernelInfoAddDefine(info, "app_ics", occaString(prefs->kernel_ics));
  occaKernelInfoAddDefine(info, "app_exact", occaString(prefs->kernel_exact));

  // Add rank to the kernels as a workaround for occa cache issues of
  // having multiple processes trying to use the same kernel source.
  occaKernelInfoAddDefine(info, "p_RANK", occaInt(prefs->rank));

  occaKernelInfoAddDefine(info, "p_NCORNERS", occaIint(P4EST_CHILDREN));
  occaKernelInfoAddDefine(info, "p_NHANG", occaIint(P4EST_HALF));

  occaKernelInfoAddDefine(info, "p_N", occaIint(N));
  occaKernelInfoAddDefine(info, "p_Nq", occaIint(Nq));
  occaKernelInfoAddDefine(info, "p_Nqk", occaIint(Nqk));
  occaKernelInfoAddDefine(info, "p_Nq2", occaIint(Nq * Nq));
  occaKernelInfoAddDefine(info, "p_Nq3", occaIint(Nq * Nq * Nq));
  occaKernelInfoAddDefine(info, "p_Np", occaIint(Np));
  occaKernelInfoAddDefine(info, "p_Nfaces", occaIint(Nfaces));
  occaKernelInfoAddDefine(info, "p_Nfp", occaIint(Nfp));

  return info;
}
// }}}

// {{{ Mesh
/* conversions between p4est and standard orientation */
#if DIM == 2
#define PXEST_ORIENTATION(f, nf, po) (2 * po)
#define PXEST_OTOH(o, h) (((o / 2) + (h)) % (2))
#elif DIM == 3

static const int8_t pxest_FToF_code[6][6] = {
    {0, 1, 1, 0, 0, 1}, {2, 0, 0, 1, 1, 0}, {2, 0, 0, 1, 1, 0},
    {0, 2, 2, 0, 0, 1}, {0, 2, 2, 0, 0, 1}, {2, 0, 0, 2, 2, 0}};

static const int8_t pxest_code_to_perm[3][4] = {
    {1, 2, 5, 6}, {0, 3, 4, 7}, {0, 4, 3, 7}};

static const int8_t pxest_perm_to_order[8][4] = {
    {0, 1, 2, 3}, {0, 2, 1, 3}, {1, 0, 3, 2}, {1, 3, 0, 2},
    {2, 0, 3, 1}, {2, 3, 0, 1}, {3, 1, 2, 0}, {3, 2, 1, 0}};

/** PXEST_ORIENTATION
 * f    my face
 * nf   neighbor face
 * o    p8est orientation code
 *
 * orientation 0:
 *   2---3     2---3
 *   |   | --> |   |
 *   0---1     0---1
 *   same:
 *   (a,b) --> (a,b)
 *
 * orientation 1:
 *   2---3     1---3
 *   |   | --> |   |
 *   0---1     0---2
 *   switch indices:
 *   (a,b) --> (b,a)
 *
 * orientation 2:
 *   2---3     3---2
 *   |   | --> |   |
 *   0---1     1---0
 *   reverse first index:
 *   (a,b) --> (N-a,b)
 *
 * orientation 3:
 *   2---3     0---2
 *   |   | --> |   |
 *   0---1     1---3
 *   reverse first index and switch:
 *   (a,b) --> (b,N-a)
 *
 * orientation 4:
 *   2---3     3---1
 *   |   | --> |   |
 *   0---1     2---0
 *   reverse second index and switch:
 *   (a,b) --> (N-b,a)
 *
 * orientation 5:
 *   2---3     0---1
 *   |   | --> |   |
 *   0---1     2---3
 *   reverse second index:
 *   (a,b) --> (a,N-b)
 *
 * orientation 6:
 *   2---3     2---0
 *   |   | --> |   |
 *   0---1     3---1
 *   reverse both and switch:
 *   (a,b) --> (N-b,N-a)
 *
 * orientation 7:
 *   2---3     1---0
 *   |   | --> |   |
 *   0---1     3---2
 *   reverse both:
 *   (a,b) --> (N-a,N-b)
 */
#define PXEST_ORIENTATION(f, nf, po)                                           \
  (pxest_code_to_perm[pxest_FToF_code[f][nf]][po])
#define PXEST_OTOH(o, h) (pxest_perm_to_order[o][h])
#else
#error "Bad Dimension"
#endif

typedef struct mesh
{
  // {{{ Filled before p4est_iterate
  int brick, *brick_n, *brick_p;
  p4est_topidx_t *brick_TToC;
  p4est_ghost_t *ghost;

  int N;      // order of the elements
  int Nq;     // N + 1
  int Nqk;    // dofs in the k direction (=1 for 2D, =Nq for 3D)
  int Np;     // dofs per element
  int Nfp;    // dofs per face
  int Nfaces; // faces per element

  iint_t Nmortar;     // number of mortar faces
  iint_t Ncontinuous; // number of continuous nodes
  iint_t Ncindices;   // number of continuous indices

  iint_t Ktotal;      // number of total elements (these may not all be valid)
  iint_t Klocal;      // number of local elements
  iint_t Kintra;      // number of local interior (aka not mirror) elements
  iint_t Kmirror;     // number of local mirror elements
  iint_t Kuniqmirror; // number of unique local mirror elements
  iint_t Kghost;      // number of remote elements
  // }}}

  // {{{ Filled from lnodes
  p4est_gloidx_t *DToC; // discontinuous to continuous node map

  iint_t *EToC;          // element to p4est lnodes face_code
  asd_dictionary_t CToD; // continuous to discontinuous node map

  iint_t ns;            // iteration variable for continuous nodes
  iint_t ni;            // iteration variable for continuous nodes
  p4est_gloidx_t c;     // previous continuous index
  iint_t *CToD_starts;  // start of indices for each continuous node
  iint_t *CToD_indices; // indices for continuous to discontinuous node map
  // }}}

  // {{{ Filled with p4est_iterate
  iint_t *EToL; // element to p4est level
  iint_t *EToT; // element to p4est treeid
  iint_t *EToX; // element to p4est x-qcoord
  iint_t *EToY; // element to p4est y-qcoord
  iint_t *EToZ; // element to p4est z-qcoord
  iint_t *EToB; // element to boundary condition
  iint_t *EToE; // element to neighbor element
  iint_t *EToF; // element to neighbor face
  iint_t *EToO; // element to neighbor orientation
  iint_t *EToP; // element to periodicity mask (filled only for brick)

  iint_t m;       // iteration variable for adaptive mortar faces
  iint_t *MFToEM; // mortar face to minus element
  iint_t *MFToFM; // mortar face to minus element face
  iint_t *MFToEP; // mortar face to plus elements
  iint_t *MFToFP; // mortar face to plus elements face
  iint_t *MFToOP; // mortar face to plus elements orientation

  iint_t i;     // iteration variable for IToE
  iint_t *IToE; // interior elements
  // }}}

  // {{{ Filled by looping through the mirrors and ghosts
  iint_t *MToE;  // mirror elements
  iint_t *UMToE; // unique mirror elements
  iint_t *GToE;  // ghost elements
  // }}}
} mesh_t;

static int continuous_to_discontinuous_extraction(const char *key, iint_t val,
                                                  void *arg)
{
  mesh_t *mesh = (mesh_t *)arg;

  p4est_gloidx_t c;
  int h;

  int rval = sscanf(key, "%" P4EST_GLOIDX_SCN ":%d:%*d:%*" IINT_SCN, &c, &h);
  ASD_ABORT_IF_NOT(rval == 2, "CToD key corruption");
  ASD_TRACE("key = \"%s\"; %" P4EST_GLOIDX_PRI "  %" IINT_PRI, key, c, val);

  if (c != mesh->c)
  {
    mesh->CToD_starts[mesh->ns] = mesh->ni;
    ASD_ABORT_IF_NOT(h == 0, "First CToD node is hanging");

    ++mesh->ns;
    mesh->c = c;
  }

  mesh->CToD_indices[mesh->ni] = val;
  ++mesh->ni;

  return 1;
}

static int has_prefix(const char *key, iint_t val, void *arg)
{
  int *has = (int *)arg;

  *has = 1;

  return 1;
}

static inline int fmask(int N, int n, int m, int f)
{
  int a, b, c;
  switch (f)
  {
  case 0:
    a = 0, b = n, c = m;
    break;
  case 1:
    a = N, b = n, c = m;
    break;
  case 2:
    a = n, b = 0, c = m;
    break;
  case 3:
    a = n, b = N, c = m;
    break;
  case 4:
    a = n, b = m, c = 0;
    break;
  case 5:
    a = n, b = m, c = N;
    break;
  default:
    a = 0, b = 0, c = 0;
  }

#if DIM == 2
  c = 0;
#endif

  return a + b * (N + 1) + c * (N + 1) * (N + 1);
}

#if DIM == 3
static inline int gmask(int N, int n, int g)
{

  int a, b, c;

  switch (g)
  {
  case 0:
    a = n, b = 0, c = 0;
    break;
  case 1:
    a = n, b = N, c = 0;
    break;
  case 2:
    a = n, b = 0, c = N;
    break;
  case 3:
    a = n, b = N, c = N;
    break;
  case 4:
    a = 0, b = n, c = 0;
    break;
  case 5:
    a = N, b = n, c = 0;
    break;
  case 6:
    a = 0, b = n, c = N;
    break;
  case 7:
    a = N, b = n, c = N;
    break;
  case 8:
    a = 0, b = 0, c = n;
    break;
  case 9:
    a = N, b = 0, c = n;
    break;
  case 10:
    a = 0, b = N, c = n;
    break;
  case 11:
    a = N, b = N, c = n;
    break;
  default:
    a = 0, b = 0, c = 0;
  }

  return a + b * (N + 1) + c * (N + 1) * (N + 1);
}
#endif

void get_hanging(const int N, const iint_t fc, int *hanging)
{
  const int Nq = N + 1;
  const int Np = Nq * Nq * ((DIM == 3) ? Nq : 1);

  for (int n = 0; n < Np; ++n)
    hanging[n] = 0;

  if (fc)
  {
    int hanging_face[P4EST_FACES];
#if DIM == 2
    p4est_lnodes_decode((p4est_lnodes_code_t)fc, hanging_face);
    for (int f = 0; f < P4EST_FACES; ++f)
      ASD_TRACE("hanging_face[%d] = %d", f, hanging_face[f]);

    for (int f = 0; f < P4EST_FACES; ++f)
      if (hanging_face[f] >= 0)
        for (int i = 0; i < Nq; ++i)
          hanging[fmask(N, i, 0, f)] = 1;
#else
    int hanging_edge[P4EST_EDGES];
    p8est_lnodes_decode((p8est_lnodes_code_t)fc, hanging_face, hanging_edge);
    for (int f = 0; f < P4EST_FACES; ++f)
      ASD_TRACE("hanging_face[%d] = %d", f, hanging_face[f]);
    for (int g = 0; g < P4EST_EDGES; ++g)
      ASD_TRACE("hanging_edge[%d] = %d", g, hanging_edge[g]);

    for (int f = 0; f < P4EST_FACES; ++f)
      if (hanging_face[f] >= 0)
        for (int i = 0; i < Nq; ++i)
          for (int j = 0; j < Nq; ++j)
            hanging[fmask(N, i, j, f)] = 1;

    for (int g = 0; g < P4EST_EDGES; ++g)
      if (hanging_edge[g] == 0 || hanging_edge[g] == 1)
        for (int i = 0; i < Nq; ++i)
          hanging[gmask(N, i, g)] = 1;
#endif
  }
}

static void mesh_iter_volume(p4est_iter_volume_info_t *info, void *user_data)
{
  mesh_t *mesh = user_data;

  p4est_tree_t *tree = p4est_tree_array_index(info->p4est->trees, info->treeid);
  const size_t e = tree->quadrants_offset + info->quadid;

  ASD_ASSERT(sizeof(p4est_qcoord_t) <= sizeof(iint_t));
  ASD_ASSERT(sizeof(p4est_topidx_t) <= sizeof(iint_t));

  mesh->EToL[e] = (iint_t)info->quad->level;
  mesh->EToT[e] = (iint_t)info->treeid;
  mesh->EToX[e] = (iint_t)info->quad->x;
  mesh->EToY[e] = (iint_t)info->quad->y;
#if DIM == 3
  mesh->EToZ[e] = (iint_t)info->quad->z;
#else
  mesh->EToZ[e] = IINT(0);
#endif

  p4est_quadrant_t qtemp = *info->quad;
  qtemp.p.piggy3.which_tree = info->treeid;
  qtemp.p.piggy3.local_num = (p4est_locidx_t)e;

  ssize_t midx = sc_array_bsearch(&mesh->ghost->mirrors, &qtemp,
                                  p4est_quadrant_compare_piggy);
  if (midx < 0)
    mesh->IToE[mesh->i++] = (iint_t)e;

  // {{{ Fill EToP for the periodic brick
  if (mesh->brick && (mesh->brick_p[0] || mesh->brick_p[1] || mesh->brick_p[2]))
  {
    mesh->EToP[e] = 0;

    for (int d = 0; d < DIM; ++d)
    {
      const p4est_topidx_t tc = mesh->brick_TToC[3 * info->treeid + d];
      if (mesh->brick_p[d] && (tc + 1 == mesh->brick_n[d]))
      {
        p4est_quadrant_t r;

        p4est_quadrant_face_neighbor(info->quad, 2 * d + 1, &r);
        if (!p4est_quadrant_is_inside_root(&r))
          mesh->EToP[e] |= 1 << d; // set the d'th bit
      }
    }
  }
  // }}}
}

/* This function is a modified version of the one with the same name in
 * p4est/src/p4est_mesh.c.
 */
static void mesh_iter_face(p4est_iter_face_info_t *info, void *user_data)
{
  mesh_t *mesh = user_data;

  iint_t *EToB = mesh->EToB;
  iint_t *EToE = mesh->EToE;
  iint_t *EToF = mesh->EToF;
  iint_t *EToO = mesh->EToO;

  iint_t *MFToEM = mesh->MFToEM;
  iint_t *MFToFM = mesh->MFToFM;
  iint_t *MFToEP = mesh->MFToEP;
  iint_t *MFToFP = mesh->MFToFP;
  iint_t *MFToOP = mesh->MFToOP;

  if (info->sides.elem_count == 1)
  {
    /* this face is on an outside boundary of the forest */
    p4est_iter_face_side_t *side = sc_array_index(&info->sides, 0);
    p4est_tree_t *tree =
        p4est_tree_array_index(info->p4est->trees, side->treeid);

    const p4est_locidx_t e = side->is.full.quadid + tree->quadrants_offset;
    const int8_t f = side->face;

    ASD_ASSERT(info->orientation == 0 && info->tree_boundary);
    ASD_ASSERT(0 <= side->treeid &&
               side->treeid < info->p4est->connectivity->num_trees);
    ASD_ASSERT(0 <= side->face && side->face < P4EST_FACES);
    ASD_ASSERT(!side->is_hanging && !side->is.full.is_ghost);
    ASD_ASSERT(0 <= e && e < info->p4est->local_num_quadrants);

    EToB[P4EST_FACES * e + f] = BC_DEFAULT;
    EToE[P4EST_FACES * e + f] = e;
    EToF[P4EST_FACES * e + f] = f;
    EToO[P4EST_FACES * e + f] = 0;
  }
  else
  {
    /* this face is between two quadrants */
    p4est_iter_face_side_t *side0 = sc_array_index(&info->sides, 0);
    p4est_iter_face_side_t *side1 = sc_array_index(&info->sides, 1);
    p4est_tree_t *tree0 =
        p4est_tree_array_index(info->p4est->trees, side0->treeid);
    p4est_tree_t *tree1 =
        p4est_tree_array_index(info->p4est->trees, side1->treeid);

    ASD_ASSERT(info->orientation == 0 || info->tree_boundary);
    ASD_ASSERT(info->sides.elem_count == 2);
    ASD_ASSERT(info->tree_boundary || side0->treeid == side1->treeid);
    ASD_ASSERT(!side0->is_hanging || !side1->is_hanging);

    if (!side0->is_hanging && !side1->is_hanging)
    {
      p4est_locidx_t e0, e1;
      const int8_t f0 = side0->face;
      const int8_t f1 = side1->face;

      const int8_t o0 = (int8_t)PXEST_ORIENTATION(f1, f0, info->orientation);
      const int8_t o1 = (int8_t)PXEST_ORIENTATION(f0, f1, info->orientation);

      /* same-size face neighbors */
      if (!side0->is.full.is_ghost)
        e0 = side0->is.full.quadid + tree0->quadrants_offset;
      else
        e0 = side0->is.full.quadid + info->p4est->local_num_quadrants;

      if (!side1->is.full.is_ghost)
        e1 = side1->is.full.quadid + tree1->quadrants_offset;
      else
        e1 = side1->is.full.quadid + info->p4est->local_num_quadrants;

      EToE[P4EST_FACES * e0 + f0] = e1;
      EToF[P4EST_FACES * e0 + f0] = f1;
      EToO[P4EST_FACES * e0 + f0] = o1;

      EToE[P4EST_FACES * e1 + f1] = e0;
      EToF[P4EST_FACES * e1 + f1] = f0;
      EToO[P4EST_FACES * e1 + f1] = o0;
    }
    else
    {
      /* one of the faces is hanging, rename so it's always side1 */
      if (side0->is_hanging)
      {
        ASD_SWAP_PTR(side0, side1);
        ASD_SWAP_PTR(tree0, tree1);
      }
      ASD_ASSERT(!side0->is_hanging && side1->is_hanging);

      p4est_tree_t *tree0 =
          p4est_tree_array_index(info->p4est->trees, side0->treeid);
      p4est_tree_t *tree1 =
          p4est_tree_array_index(info->p4est->trees, side1->treeid);

      // {{{ fill e0, e1, f0, f1, o1
      p4est_locidx_t e0, e1[P4EST_HALF];
      const int8_t f0 = side0->face;
      const int8_t f1 = side1->face;
      const int8_t o1 = (int8_t)PXEST_ORIENTATION(f0, f1, info->orientation);

      if (!side0->is.full.is_ghost)
        e0 = side0->is.full.quadid + tree0->quadrants_offset;
      else
        e0 = side0->is.full.quadid + info->p4est->local_num_quadrants;

      for (int h = 0; h < P4EST_HALF; ++h)
      {
        if (!side1->is.hanging.is_ghost[h])
          e1[h] = side1->is.hanging.quadid[h] + tree1->quadrants_offset;
        else
          e1[h] =
              side1->is.hanging.quadid[h] + info->p4est->local_num_quadrants;
      }
      // }}}

      // {{{ disconnect hanging faces
      EToB[P4EST_FACES * e0 + f0] = BC_SKIP;
      EToE[P4EST_FACES * e0 + f0] = e0;
      EToF[P4EST_FACES * e0 + f0] = f0;
      EToO[P4EST_FACES * e0 + f0] = 0;

      for (int h = 0; h < P4EST_HALF; ++h)
      {
        EToB[P4EST_FACES * e1[h] + f1] = BC_SKIP;
        EToE[P4EST_FACES * e1[h] + f1] = e1[h];
        EToF[P4EST_FACES * e1[h] + f1] = f1;
        EToO[P4EST_FACES * e1[h] + f1] = 0;
      }
      // }}}

      // {{{ connect hanging faces
      const iint_t m = mesh->m;

      MFToEM[m] = e0;
      MFToFM[m] = f0;
      MFToFP[m] = f1;
      MFToOP[m] = o1;

      for (int h = 0; h < P4EST_HALF; ++h)
        MFToEP[P4EST_HALF * m + h] = e1[PXEST_OTOH(o1, h)];

      ++mesh->m;
      // }}}
    }
  }
}

static mesh_t *mesh_new(prefs_t *prefs, p4est_t *pxest, p4est_ghost_t *ghost)
{
  mesh_t *mesh = asd_malloc(sizeof(mesh_t));

  const int N = mesh->N = prefs->mesh_N;
  const int Nq = mesh->Nq = N + 1;
  const int Nfaces = mesh->Nfaces = 2 * DIM;

#if DIM == 3
  mesh->Nqk = Nq;
  mesh->Nfp = Nq * Nq;
  const int Np = mesh->Np = Nq * Nq * Nq;
#else
  mesh->Nqk = 1;
  mesh->Nfp = Nq;
  const int Np = mesh->Np = Nq * Nq;
#endif

  ASD_ABORT_IF(sizeof(p4est_locidx_t) > sizeof(iint_t),
               "p4est_locidx_t not compatible with iint_t");

  // {{{ get lnodes and ghost layer to match
  p4est_lnodes_t *lnodes = p4est_lnodes_new(pxest, ghost, N);
  p4est_ghost_support_lnodes(pxest, lnodes, ghost);
  p4est_ghost_expand_by_lnodes(pxest, lnodes, ghost);
  // }}}

  const iint_t Klocal = mesh->Klocal = (iint_t)pxest->local_num_quadrants;
  const iint_t Kghost = mesh->Kghost = (iint_t)ghost->ghosts.elem_count;
  const iint_t Kuniqmirror = mesh->Kuniqmirror =
      (iint_t)ghost->mirrors.elem_count;
  const iint_t Kmirror = mesh->Kmirror =
      (iint_t)ghost->mirror_proc_offsets[pxest->mpisize];
  const iint_t Ktotal = mesh->Ktotal = Klocal + Kghost;
  const iint_t Kintra = mesh->Kintra = Klocal - Kuniqmirror;

  ASD_ASSERT(ghost->proc_offsets[pxest->mpisize] ==
             (int)ghost->ghosts.elem_count);

  mesh->ghost = ghost;
  mesh->brick = prefs->brick;
  mesh->brick_n = prefs->brick_n;
  mesh->brick_p = prefs->brick_p;
  mesh->brick_TToC = prefs->brick_TToC;

  // {{{ continuous to discontinuous node map
  asd_dictionary_init(&mesh->CToD);
  mesh->EToC = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToP = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->DToC = asd_malloc_aligned(sizeof(p4est_gloidx_t) * Ktotal * Np);
  mesh->CToD_starts = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Np);
  mesh->CToD_indices = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Np);

  // {{{ get EToC (from p4est_plex.c)
  {
    iint_t **mirror_EToC;

    for (iint_t e = 0; e < Klocal; ++e)
      mesh->EToC[e] = lnodes->face_code[e];

    mirror_EToC =
        asd_malloc_aligned(sizeof(iint_t *) * ghost->mirrors.elem_count);
    for (size_t m = 0; m < ghost->mirrors.elem_count; ++m)
    {
      p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->mirrors, m);
      mirror_EToC[m] = &mesh->EToC[q->p.piggy3.local_num];
    }
    p4est_ghost_exchange_custom(pxest, ghost, sizeof(iint_t),
                                (void **)mirror_EToC, &mesh->EToC[Klocal]);
    asd_free_aligned(mirror_EToC);
  }
  // }}}

  // {{{ get DToC (from p4est_plex.c)
  {
    p4est_gloidx_t **mirror_DToC;

    for (iint_t n = 0; n < Klocal * Np; ++n)
      mesh->DToC[n] =
          p4est_lnodes_global_index(lnodes, lnodes->element_nodes[n]);

    mirror_DToC = asd_malloc_aligned(sizeof(p4est_gloidx_t *) *
                                     ghost->mirrors.elem_count);
    for (size_t m = 0; m < ghost->mirrors.elem_count; ++m)
    {
      p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->mirrors, m);
      mirror_DToC[m] = &mesh->DToC[q->p.piggy3.local_num * Np];
    }
    p4est_ghost_exchange_custom(pxest, ghost, sizeof(p4est_gloidx_t) * Np,
                                (void **)mirror_DToC, &mesh->DToC[Klocal * Np]);
    asd_free_aligned(mirror_DToC);
  }
  // }}}

  int *hanging = asd_malloc_aligned(sizeof(int) * Np);
  // {{{  Fill CToD with local dofs
  for (iint_t e = 0; e < Klocal; ++e)
  {
    get_hanging(N, mesh->EToC[e], hanging);

    for (int n = 0; n < Np; ++n)
    {
      const iint_t id = Np * e + n;
      char key[ASD_BUFSIZ];
      snprintf(key, ASD_BUFSIZ,
               "%0" P4EST_GLOIDX_MAX_DIGITS P4EST_GLOIDX_PRI
               ":%d:%0" INT_MAX_DIGITS "d:%0" IINT_MAX_DIGITS IINT_PRI
               ":%0" INT_MAX_DIGITS "d",
               mesh->DToC[id], hanging[n], pxest->mpirank, e, n);
      asd_dictionary_insert_iint(&mesh->CToD, key, id);
    }
  }
  // }}}

  // {{{ Fill CToD with ghost dofs
  int grank = 0;
  for (iint_t g = 0; g < (iint_t)ghost->ghosts.elem_count; ++g)
  {
    while (ghost->proc_offsets[grank + 1] <= g)
    {
      ++grank;
      ASD_ASSERT(grank < pxest->mpisize);
    }

    get_hanging(N, mesh->EToC[Klocal + g], hanging);

    for (int n = 0; n < Np; ++n)
    {
      const iint_t id = Np * Klocal + Np * g + n;
      char key[ASD_BUFSIZ];
      snprintf(key, ASD_BUFSIZ, "%0" P4EST_GLOIDX_MAX_DIGITS P4EST_GLOIDX_PRI,
               mesh->DToC[id]);

      int has_local = 0;
      asd_dictionary_allprefixed_iint(&mesh->CToD, key, &has_prefix,
                                      &has_local);
      // add the ghost node if it is associated with a local node
      if (has_local)
      {
        p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->ghosts, g);
        snprintf(key, ASD_BUFSIZ,
                 "%0" P4EST_GLOIDX_MAX_DIGITS P4EST_GLOIDX_PRI
                 ":%d:%0" INT_MAX_DIGITS "d:%0" IINT_MAX_DIGITS IINT_PRI
                 ":%0" INT_MAX_DIGITS "d",
                 mesh->DToC[id], hanging[n], grank,
                 (iint_t)q->p.piggy3.local_num, n);
        asd_dictionary_insert_iint(&mesh->CToD, key, id);
      }
    }
  }
  // }}}
  asd_free_aligned(hanging);

  // {{{ Extract starts and indices from CToD
  iint_t Ncindices = mesh->Ncindices = (iint_t)mesh->CToD.num_entries;

  mesh->ns = mesh->ni = 0;
  mesh->c = -1;
  asd_dictionary_allprefixed_iint(
      &mesh->CToD, "", &continuous_to_discontinuous_extraction, mesh);

  ASD_ASSERT(Ncindices == mesh->ni);
  iint_t Ncontinuous = mesh->Ncontinuous = mesh->ns;
  mesh->CToD_starts[Ncontinuous] = Ncindices;
  // }}}

  // }}}

  mesh->i = 0;
  mesh->IToE = asd_malloc_aligned(sizeof(iint_t) * Klocal);
  mesh->MToE = asd_malloc_aligned(sizeof(iint_t) * Kmirror);
  mesh->UMToE = asd_malloc_aligned(sizeof(iint_t) * Kuniqmirror);
  mesh->GToE = asd_malloc_aligned(sizeof(iint_t) * Kghost);

  mesh->EToL = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToT = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToX = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToY = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToZ = asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToB = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  for (iint_t n = 0; n < Ktotal * Nfaces; ++n)
    mesh->EToB[n] = BC_NONE;

  mesh->EToE = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->EToF = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->EToO = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);

  mesh->m = 0;
  mesh->MFToEM = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToFM = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToEP =
      asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces * P4EST_HALF);
  mesh->MFToFP = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToOP = asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);

  p4est_iterate(pxest, ghost, mesh, mesh_iter_volume, mesh_iter_face, NULL
#if DIM == 3
                ,
                NULL
#endif
                );

  ASD_ASSERT(Kintra == mesh->i);
  mesh->Nmortar = mesh->m;

  // {{{ update mesh information with ghosts
  for (size_t g = 0; g < ghost->ghosts.elem_count; ++g)
  {
    iint_t e = (iint_t)g + pxest->local_num_quadrants;
    p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->ghosts, g);

    mesh->GToE[g] = e;

    mesh->EToL[e] = (iint_t)q->level;
    mesh->EToT[e] = (iint_t)q->p.piggy3.which_tree;
    mesh->EToX[e] = (iint_t)q->x;
    mesh->EToY[e] = (iint_t)q->y;
#if DIM == 3
    mesh->EToZ[e] = (iint_t)q->z;
#else
    mesh->EToZ[e] = IINT(0);
#endif

    for (int f = 0; f < P4EST_FACES; ++f)
    {
      mesh->EToB[P4EST_FACES * e + f] = BC_SKIP;
      mesh->EToE[P4EST_FACES * e + f] = e;
      mesh->EToF[P4EST_FACES * e + f] = f;
      mesh->EToO[P4EST_FACES * e + f] = 0;
    }

    // {{{ Fill EToP for the periodic brick
    if (mesh->brick &&
        (mesh->brick_p[0] || mesh->brick_p[1] || mesh->brick_p[2]))
    {
      mesh->EToP[e] = 0;

      for (int d = 0; d < DIM; ++d)
      {
        const p4est_topidx_t tc =
            mesh->brick_TToC[3 * q->p.piggy3.which_tree + d];
        if (mesh->brick_p[d] && (tc + 1 == mesh->brick_n[d]))
        {
          p4est_quadrant_t r;

          p4est_quadrant_face_neighbor(q, 2 * d + 1, &r);
          if (!p4est_quadrant_is_inside_root(&r))
            mesh->EToP[e] |= 1 << d; // set the d'th bit
        }
      }
    }
    // }}}
  }
  // }}}

  // {{{ update mesh information with mirrors
  for (size_t m = 0; m < ghost->mirrors.elem_count; ++m)
  {
    p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->mirrors, m);
    mesh->UMToE[m] = q->p.piggy3.local_num;
  }

  for (int m = 0; m < ghost->mirror_proc_offsets[pxest->mpisize]; ++m)
  {
    p4est_quadrant_t *q = p4est_quadrant_array_index(
        &ghost->mirrors, ghost->mirror_proc_mirrors[m]);
    mesh->MToE[m] = q->p.piggy3.local_num;
  }
// }}}

#ifdef ASD_DEBUG
  // {{{ Print communication information
  ASD_LDEBUG("Ncontinuous = %jd", (intmax_t)Ncontinuous);
  ASD_LDEBUG("");
  ASD_LDEBUG("");
  ASD_LDEBUG("Ghost neighbors");
  ASD_LDEBUG("rank  remote element number  local element number");
  ASD_LDEBUG("---------------------------");
  for (p4est_locidx_t g = 0, r = 0;
       g < (p4est_locidx_t)ghost->ghosts.elem_count; ++g)
  {
    while (ghost->proc_offsets[r + 1] <= g)
    {
      ++r;
      ASD_ASSERT(r < pxest->mpisize);
    }

    p4est_quadrant_t *q = p4est_quadrant_array_index(&ghost->ghosts, g);
    ASD_LDEBUG("%4zd   %12jd   %12jd", r, (intmax_t)q->p.piggy3.local_num,
               (intmax_t)mesh->GToE[g]);
  }

  ASD_LDEBUG("");
  ASD_LDEBUG("");
  ASD_LDEBUG("Mirrors ");
  ASD_LDEBUG("rank  local element number");
  ASD_LDEBUG("---------------------------");
  for (iint_t m = 0, r = 0; m < Kmirror; ++m)
  {
    while (ghost->mirror_proc_offsets[r + 1] <= m)
    {
      ++r;
      ASD_ASSERT(r < pxest->mpisize);
    }
    ASD_LDEBUG("%4zd   %jd", r, (intmax_t)mesh->MToE[m]);
  }
  ASD_LDEBUG("");
  ASD_LDEBUG("");
// }}}
#endif

  p4est_lnodes_destroy(lnodes);

  return mesh;
}

static void mesh_free(mesh_t *mesh)
{
  asd_free_aligned(mesh->IToE);
  asd_free_aligned(mesh->MToE);
  asd_free_aligned(mesh->UMToE);
  asd_free_aligned(mesh->GToE);

  asd_free_aligned(mesh->EToL);
  asd_free_aligned(mesh->EToT);
  asd_free_aligned(mesh->EToX);
  asd_free_aligned(mesh->EToY);
  asd_free_aligned(mesh->EToZ);
  asd_free_aligned(mesh->EToB);
  asd_free_aligned(mesh->EToE);
  asd_free_aligned(mesh->EToF);
  asd_free_aligned(mesh->EToO);

  asd_free_aligned(mesh->MFToEM);
  asd_free_aligned(mesh->MFToFM);
  asd_free_aligned(mesh->MFToEP);
  asd_free_aligned(mesh->MFToFP);
  asd_free_aligned(mesh->MFToOP);

  asd_dictionary_clear(&mesh->CToD);

  asd_free_aligned(mesh->DToC);
  asd_free_aligned(mesh->EToC);
  asd_free_aligned(mesh->EToP);
  asd_free_aligned(mesh->CToD_starts);
  asd_free_aligned(mesh->CToD_indices);
}
// }}}

// {{{ Level
typedef struct level
{
  int N;      // order of the elements
  int Nq;     // N + 1
  int Nqk;    // dofs in the k direction (=1 for 2D, =Nq for 3D)
  int Np;     // dofs per element
  int Nfp;    // dofs per face
  int Nfaces; // faces per element

  iint_t Nmortar;     // number of mortar faces
  iint_t Ncontinuous; // number of continuous nodes
  iint_t Ncindices;   // number of continuous indices

  iint_t Kmax;        // number of elements used for the allocations
  iint_t Ktotal;      // number of total elements (these may not all be valid)
  iint_t Klocal;      // number of local elements
  iint_t Kintra;      // number of local interior (aka not mirror) elements
  iint_t Kmirror;     // number of local mirror elements
  iint_t Kuniqmirror; // number of unique local mirror elements
  iint_t Kghost;      // number of remote elements

  // p4est conn information
  occaMemory o_tree_to_vertex;
  occaMemory o_tree_vertices;

  // p4est information
  int Nn;    // number of mirror and ghost neighbors
  int *NToR; // neighbor to rank

  int8_t *EToA; // refine, none, and coarsen flag for each element

  occaMemory o_IToE;  // intra elements
  occaMemory o_MToE;  // mirror elements
  occaMemory o_UMToE; // unique mirror elements
  occaMemory o_GToE;  // ghost elements

  occaMemory o_EToL;   // element to p4est level
  occaMemory o_EToT;   // element to p4est treeid
  occaMemory o_EToX;   // element to p4est x-coord
  occaMemory o_EToY;   // element to p4est y-coord
  occaMemory o_EToZ;   // element to p4est z-coord
  occaMemory o_EToB;   // element to boundary condition
  occaMemory o_EToE;   // element to neighbor element
  occaMemory o_EToF;   // element to neighbor face
  occaMemory o_EToO;   // element to neighbor orientation
  occaMemory o_EToC;   // element to p4est lnodes face_code
  occaMemory o_EToP;   // element to periodicity mask
  occaMemory o_EToOff; // element to offset used in adaptivity

  occaMemory o_CToD_starts;  // start of indices for each continuous node
  occaMemory o_CToD_indices; // indices for continuous to discontinuous node map

  occaMemory o_MFToEM; // mortar face to minus element
  occaMemory o_MFToFM; // mortar face to minus element face
  occaMemory o_MFToEP; // mortar face to plus elements
  occaMemory o_MFToFP; // mortar face to plus elements face
  occaMemory o_MFToOP; // mortar face to plus elements orientation

  // element operators
  occaMemory o_r;
  occaMemory o_w;
  occaMemory o_D;

  // element interpolation operators
  occaMemory o_Ib;
  occaMemory o_It;

  // element L^2 projection operators
  occaMemory o_Pb;
  occaMemory o_Pt;

  // fields
  occaMemory o_q;
  occaMemory o_rhsq;

  // arrays for sending & receiving q
  occaMemory o_q_buf;

  occaMemory pin_q_recv;
  occaMemory pin_q_send;

  dfloat_t *q_recv;
  dfloat_t *q_send;

  dfloat_t *q_recv_buf;
  dfloat_t *q_send_buf;

  MPI_Request *q_recv_requests;
  MPI_Request *q_send_requests;

  MPI_Status *q_recv_statuses;
  MPI_Status *q_send_statuses;

  // geometry information
  occaMemory o_vgeo;
  occaMemory o_sgeo;

  // reduction buffers
  occaMemory o_red_buf[2];

  // kernels
  occaKernel compute_X;
  occaKernel interp_X;
  occaKernel coarse_X;

  occaKernel compute_geo;
  occaKernel compute_ics;
  occaKernel compute_dt;
  occaKernel compute_energy;
  occaKernel compute_error;

  occaKernel coarsen_fields;
  occaKernel refine_and_fill_fields;

  occaKernel volume_advection;
  occaKernel mortar_advection;
  occaKernel update_advection;
  occaKernel zero_fields;

  occaKernel get_mirror_fields;
  occaKernel set_ghost_fields;

  occaKernel reduce_min;
  occaKernel reduce_sum;
} level_t;

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

static level_t *level_new(prefs_t *prefs, p4est_t *pxest,
                            p4est_ghost_t *ghost, occaDevice device)
{
  level_t *lvl = asd_malloc(sizeof(level_t));

  ASD_ABORT_IF(sizeof(p4est_locidx_t) > sizeof(iint_t),
               "p4est_locidx_t not compatible with iint_t");

  ASD_ABORT_IF(sizeof(p4est_topidx_t) > sizeof(iint_t),
               "p4est_topidx_t not compatible with iint_t");

  // {{{ send connectivity to the device
  p4est_connectivity_t *conn = pxest->connectivity;
  lvl->o_tree_vertices =
      occa_double_to_dfloat(device, NX * conn->num_vertices, conn->vertices);
  lvl->o_tree_to_vertex = occa_p4est_topidx_to_iint(
      device, P4EST_CHILDREN * conn->num_trees, conn->tree_to_vertex);
  // }}}

  // {{{ Build Operators
  get_operators(prefs->mesh_N, device, &lvl->o_r, &lvl->o_w, &lvl->o_D,
                &lvl->o_Ib, &lvl->o_It, &lvl->o_Pb, &lvl->o_Pt);
  // }}}

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

  return lvl;
}

static void level_free(level_t *lvl)
{
  occaMemoryFree(lvl->o_tree_to_vertex);
  occaMemoryFree(lvl->o_tree_vertices);

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
}
// }}}

// {{{ App
typedef struct app
{
  prefs_t *prefs;
  occaDevice device;
  occaStream copy;
  occaStream cmdx;

  p4est_connectivity_t *conn;
  p4est_t *pxest;
  p4est_ghost_t *ghost;

  level_t *lvl;
} app_t;

static int static_refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                            p4est_quadrant_t *q)
{
  return 1;
}

static int app_refine_fn(p4est_t *p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t *q)
{

  int val = 0;

  prefs_t *prefs = (prefs_t *)p4est->user_pointer;

  ASD_ASSERT(sizeof(p4est_locidx_t) <= sizeof(int));
  ASD_ASSERT(sizeof(p4est_topidx_t) <= sizeof(int));

  const int tree = which_tree, level = q->level, x = q->x, y = q->y;
#if DIM == 3
  const int z = q->z;
#else
  const int z = 0;
#endif

  int result = asd_lua_global_function_call(
      prefs->L, "app.mesh.initial_refinement", "iiiii"
                                               ">"
                                               "i",
      tree, level, x, y, z, &val);
  ASD_ABORT_IF_NOT(result == 0, "Expected 0 got %d", result);

  return val;
}

/** create app structure
 *
 * \param [in]  filename   the preference file
 * \param [in]  comm       MPI communicator
 * \param [in]  loglevel   level of logging
 *
 * \return Initialized app structure
 */
static app_t *app_new(const char *filename, MPI_Comm comm, int loglevel)
{
  app_t *app = asd_malloc(sizeof(app_t));

  //
  // Preferences
  //
  app->prefs = prefs_new(filename, comm, loglevel);
  prefs_print(app->prefs);

  //
  // OCCA
  //
  app->device = occaCreateDevice(app->prefs->occa_info);
  if (app->prefs->occa_flags)
    occaDeviceSetCompilerFlags(app->device, app->prefs->occa_flags);

  app->copy = occaDeviceCreateStream(app->device);
  app->cmdx = occaDeviceCreateStream(app->device);
  occaDeviceSetStream(app->device, app->cmdx);

  app->conn = get_connectivity(app->prefs->conn_name, app->prefs->brick_n,
                               app->prefs->brick_p);

  if (app->prefs->brick)
  {
    app->prefs->brick_TToC =
        asd_malloc_aligned(sizeof(p4est_topidx_t) * app->conn->num_trees * 3);
    const double *vertices = app->conn->vertices;
    const p4est_topidx_t *tree_to_vertex = app->conn->tree_to_vertex;

    for (p4est_topidx_t t = 0; t < app->conn->num_trees; ++t)
      for (int d = 0; d < 3; ++d)
        app->prefs->brick_TToC[3 * t + d] = (p4est_topidx_t)
            vertices[3 * tree_to_vertex[t * P4EST_CHILDREN] + d];
  }
  else
  {
    app->prefs->brick_TToC = asd_malloc_aligned(0);
  }

  /* TODO Should this be moved to the device? */
  if (app->prefs->conn_vertex_mapping)
  {
    /* let the user modify the connectivity vertices */
    ASD_ASSERT(sizeof(int) >= sizeof(p4est_topidx_t));
    for (p4est_topidx_t i = 0; i < app->conn->num_vertices; i++)
    {
      double x = app->conn->vertices[i * 3 + 0];
      double y = app->conn->vertices[i * 3 + 1];
      double z = app->conn->vertices[i * 3 + 2];
      int result =
          asd_lua_global_function_call(app->prefs->L, "app.conn.vertex_mapping",
                                       "iddd>ddd", (int)i, x, y, z, &x, &y, &z);
      if (result != 0)
        ASD_ABORT("Failed to move node %jd (%d)", (intmax_t)i);
      app->conn->vertices[i * 3 + 0] = x;
      app->conn->vertices[i * 3 + 1] = y;
      app->conn->vertices[i * 3 + 2] = z;
    }
  }

  app->pxest = p4est_new_ext(app->prefs->comm, app->conn, 0,
                             app->prefs->mesh_start_level, 1,
                             sizeof(quad_data_t), NULL, NULL);

  /* initial refinement */
  if (app->prefs->mesh_initial_refinement)
  {
    void *current_user_pointer = app->pxest->user_pointer;
    app->pxest->user_pointer = app->prefs;
    p4est_refine_ext(app->pxest, 1, -1, app_refine_fn, NULL, NULL);
    p4est_refine(app->pxest, 1, app_refine_fn, NULL);
    app->pxest->user_pointer = current_user_pointer;
  }
  p4est_balance_ext(app->pxest, P4EST_CONNECT_FULL, NULL, NULL);
  p4est_partition(app->pxest, 1, NULL);

  if (app->prefs->mesh_static_refinement > 0)
  {
    for (int r = 0; r < app->prefs->mesh_static_refinement; ++r)
      p4est_refine_ext(app->pxest, 0, -1, static_refine_fn, NULL, NULL);
    p4est_balance_ext(app->pxest, P4EST_CONNECT_FULL, NULL, NULL);
    p4est_partition(app->pxest, 1, NULL);
  }

  app->ghost = p4est_ghost_new(app->pxest, P4EST_CONNECT_FULL);

  app->lvl = level_new(app->prefs, app->pxest, app->ghost, app->device);

  return app;
}

static void app_free(app_t *app)
{
  prefs_free(app->prefs);
  asd_free(app->prefs);

  level_free(app->lvl);
  asd_free(app->lvl);

  p4est_ghost_destroy(app->ghost);
  p4est_destroy(app->pxest);
  p4est_connectivity_destroy(app->conn);

  occaStreamFree(app->copy);
  occaStreamFree(app->cmdx);

  uintmax_t bytes = occaDeviceBytesAllocated(app->device);
  uintmax_t total_bytes = occaDeviceMemorySize(app->device);
  ASD_INFO("");
  ASD_INFO("Device bytes allocated %ju (%.2f GiB) out of %ju (%.2f GiB)", bytes,
           ((double)bytes) / GiB, total_bytes, ((double)total_bytes) / GiB);

  occaDeviceFree(app->device);
}
// }}}

// {{{ mfem

static void mfem_write_mesh(const int rank, const char *outdir, const int N,
                            const iint_t K, dfloat_t *vgeo)
{
  const int Ncorners = (DIM == 2) ? 4 : 8;
  const int geo_type = (DIM == 2) ? 3 : 5;
  const int Np = (DIM == 2) ? (N + 1) * (N + 1) : (N + 1) * (N + 1) * (N + 1);

  char filename[ASD_BUFSIZ];
  snprintf(filename, ASD_BUFSIZ, "%s/mesh.%06" IINT_PRI, outdir, rank);
  ASD_VERBOSE("Writing file: '%s'", filename);

  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    ASD_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "MFEM mesh v1.0\n\n");
  fprintf(file, "dimension\n%d\n\n", DIM);

  fprintf(file, "elements\n%" IINT_PRI "\n", K);
  for (iint_t e = 0; e < K; ++e)
  {
    fprintf(file, "%" IINT_PRI " %d", e + 1, geo_type);
    for (int c = 0; c < Ncorners; ++c)
      fprintf(file, " %" IINT_PRI, Ncorners * e + c);
    fprintf(file, "\n");
  }
  fprintf(file, "\n");

  fprintf(file, "boundary\n0\n\n");
  fprintf(file, "vertices\n%" IINT_PRI "\n\n", K * Ncorners);
  fprintf(file, "nodes\n");
  fprintf(file, "FiniteElementSpace\n");
  fprintf(file, "FiniteElementCollection: L2_T1_%dD_P%d\n", DIM, N);
  fprintf(file, "VDim: %d\n", DIM);
  fprintf(file, "Ordering: 1\n");

  for (iint_t e = 0; e < K; ++e)
  {
    for (int n = 0; n < Np; ++n)
    {
      const size_t idd = n + e * Np * NVGEO;
      dfloat_t x = vgeo[idd + Np * VGEO_X];
      dfloat_t y = vgeo[idd + Np * VGEO_Y];

#if DIM == 2
      fprintf(file, "         %" DFLOAT_FMTe " %" DFLOAT_FMTe "\n", x, y);
#else
      dfloat_t z = vgeo[idd + Np * VGEO_Z];
      fprintf(file,
              "         %" DFLOAT_FMTe " %" DFLOAT_FMTe " %" DFLOAT_FMTe "\n",
              x, y, z);
#endif
    }
  }

  if (ferror(file))
  {
    ASD_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    ASD_LERROR("Error closing %s\n", filename);
  }
}

static void mfem_write_scalar(const int rank, const char *outdir,
                              const char *name, int N, iint_t K,
                              const dfloat_t *d, const int o, const int ot)
{
  const int Np = (DIM == 2) ? (N + 1) * (N + 1) : (N + 1) * (N + 1) * (N + 1);

  char filename[ASD_BUFSIZ];
  snprintf(filename, ASD_BUFSIZ, "%s/%s.%06" IINT_PRI, outdir, name, rank);
  ASD_VERBOSE("Writing file: '%s'", filename);

  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    ASD_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "FiniteElementSpace\n");
  fprintf(file, "FiniteElementCollection: L2_T1_%dD_P%d\n", DIM, N);
  fprintf(file, "VDim: 1\n");
  fprintf(file, "Ordering: 1\n");

  for (iint_t e = 0; e < K; ++e)
  {
    for (int n = 0; n < Np; ++n)
    {
      const size_t idd = n + e * Np * ot;
      dfloat_t a = d[idd + Np * o];
      fprintf(file, "         %" DFLOAT_FMTe "\n", a);
    }
  }

  if (ferror(file))
  {
    ASD_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    ASD_LERROR("Error closing %s\n", filename);
  }
}

static void mfem_write_vector(const int rank, const char *outdir,
                              const char *name, int N, iint_t K,
                              const dfloat_t *d, const int o1, const int o2,
                              const int o3, const int ot)
{
  const int Np = (DIM == 2) ? (N + 1) * (N + 1) : (N + 1) * (N + 1) * (N + 1);

  char filename[ASD_BUFSIZ];
  snprintf(filename, ASD_BUFSIZ, "%s/%s.%06" IINT_PRI, outdir, name, rank);
  ASD_VERBOSE("Writing file: '%s'", filename);

  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    ASD_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "FiniteElementSpace\n");
  fprintf(file, "FiniteElementCollection: L2_T1_%dD_P%d\n", DIM, N);
  fprintf(file, "VDim: %d\n", DIM);
  fprintf(file, "Ordering: 1\n");

  for (iint_t e = 0; e < K; ++e)
  {
    for (int n = 0; n < Np; ++n)
    {
      const size_t idd = n + e * Np * ot;
      dfloat_t a = d[idd + Np * o1];
      dfloat_t b = d[idd + Np * o2];

#if DIM == 2
      fprintf(file, "         %" DFLOAT_FMTe " %" DFLOAT_FMTe "\n", a, b);
#else
      dfloat_t c = d[idd + Np * o3];
      fprintf(file,
              "         %" DFLOAT_FMTe " %" DFLOAT_FMTe " %" DFLOAT_FMTe "\n",
              a, b, c);
#endif
    }
  }

  if (ferror(file))
  {
    ASD_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    ASD_LERROR("Error closing %s\n", filename);
  }
}

static void mfem_write_file(const int rank, const char *directory,
                            const char *prefix, const int N, const iint_t K,
                            dfloat_t *vgeo, dfloat_t *q)
{
  char outdir[ASD_BUFSIZ];
  struct stat sb;

  // create directory for the timestep data
  snprintf(outdir, ASD_BUFSIZ, "%s/%s", directory, prefix);
  if (stat(outdir, &sb) != 0 && mkdir(outdir, 0755) != 0 && errno != EEXIST)
    perror("making mfem directory");

  mfem_write_mesh(rank, outdir, N, K, vgeo);

  const char *const *scalars = FIELD_OUT_SCALARS;
  const char *const *vectors = FIELD_OUT_VECTORS;

  if (scalars)
    for (size_t s = 0; scalars[s]; ++s)
    {
      if (FIELD_OUT_SCALARS_OFF[s] >= 0)
        mfem_write_scalar(rank, outdir, scalars[s], N, K, q,
                          FIELD_OUT_SCALARS_OFF[s], NFIELDS);
      else
        mfem_write_scalar(rank, outdir, scalars[s], N, K, vgeo,
                          -FIELD_OUT_SCALARS_OFF[s] - 1, NVGEO);
    }

  if (vectors)
    for (size_t v = 0; vectors[v]; ++v)
      mfem_write_vector(rank, outdir, vectors[v], N, K, q,
                        FIELD_OUT_COMPONENTS_OFF[3 * v + 0],
                        FIELD_OUT_COMPONENTS_OFF[3 * v + 1],
                        FIELD_OUT_COMPONENTS_OFF[3 * v + 2], NFIELDS);
}

static void app_output_mfem(app_t *app, int s, const char *prefix)
{
  ASD_ROOT_INFO("");
  char fname[ASD_BUFSIZ];
  snprintf(fname, ASD_BUFSIZ, "%s%s_%06d", prefix, app->prefs->output_prefix,
           s);

  ASD_ROOT_INFO("Writing mfem files for \"%s\"", fname);

  const iint_t K = app->lvl->Klocal;
  const int N = app->lvl->N;
  const int Np = app->lvl->Np;

  /* Make sure the copy commands are placed in the command queue.  This is
   * to ensure that Q will be finished computing before we copy the data
   * to the host.
   */
  occaDeviceSetStream(app->device, app->cmdx);

  size_t local_vgeo_sz = NVGEO * K * Np * sizeof(dfloat_t);
  dfloat_t *vgeo = asd_malloc_aligned(local_vgeo_sz);
  occaCopyMemToPtr(vgeo, app->lvl->o_vgeo, local_vgeo_sz, occaNoOffset);

  size_t local_q_sz = NFIELDS * K * Np * sizeof(dfloat_t);
  dfloat_t *q = asd_malloc_aligned(local_q_sz);
  occaCopyMemToPtr(q, app->lvl->o_q, local_q_sz, occaNoOffset);

  mfem_write_file(app->prefs->rank, app->prefs->output_datadir, fname, N, K,
                  vgeo, q);

  asd_free_aligned(vgeo);
  asd_free_aligned(q);
}
// }}}

// {{{ VTK
static void strcsl(const size_t n, char *str, const char *const *list)
{
  if (n > 0)
  {
    str[0] = '\0';

    if (list != NULL)
    {
      for (size_t l = 0; list[l]; ++l)
      {
        if (strlen(list[l] + 1) <= (n - strlen(str) - 1))
          strncat(str, list[l], (n - strlen(str) - 1));
        if (list[l + 1] && 2 <= (n - strlen(str) - 1))
          strncat(str, ",", (n - strlen(str) - 1));
      }
    }
  }
}

/*
 * Integer power routine from:
 *   http://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
 */
static int ipow(int base, int exp)
{
  ASD_ASSERT(exp >= 0);

  int result = 1;
  while (exp)
  {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

#define VTU_FORMAT "%s_%05d.vtu"

/** Utility function to write binary data in VTK format.
  *
  * Currently this is just a wrapper to call a similar function in libsc.
  *
  * \param [in]  compress boolean specifying if the binary data should be
  *                       compressed.
  * \param [out] file     stream to write the data to.
  * \param [in]  data     data to write out.
  * \param [in]  size     size of the data in bytes.
  *
  * \returns 0 on success and -1 on file error.
  */
static int vtk_write_binary_data(int compress, FILE *file, char *data,
                                 size_t size)
{
  return (compress) ? sc_vtk_write_compressed(file, data, size)
                    : sc_vtk_write_binary(file, data, size);
}

static void vtk_write_file_pvtu(int size, const char *directory,
                                const char *const prefix,
                                const char *const *scalars,
                                const char *const *vectors, int binary,
                                int compress)
{
  const int endian = asd_endian();
  const char *format = (binary) ? "binary" : "ascii";

  char filename[ASD_BUFSIZ];
  snprintf(filename, ASD_BUFSIZ, "%s.pvtu", prefix);
  ASD_VERBOSE("Writing file: '%s'", filename);

  FILE *file = fopen(filename, "w");
  if (file == NULL)
  {
    ASD_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == ASD_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
  fprintf(file, "    <PPoints>\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"Position\""
                " NumberOfComponents=\"3\" format=\"%s\"/>\n",
          DFLOAT_VTK, format);
  fprintf(file, "    </PPoints>\n");

  fprintf(file, "    <PCellData Scalars=\"time,mpirank\">\n");
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"time\""
                " format=\"%s\"/>\n",
          DFLOAT_VTK, format);

  const char *int_vtk_type;
  if (sizeof(int) == sizeof(int32_t))
    int_vtk_type = "Int32";
  else if (sizeof(int) == sizeof(int64_t))
    int_vtk_type = "Int64";
  else
    int_vtk_type = "UNKNOWN";
  fprintf(file, "      <PDataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\"/>\n",
          int_vtk_type, format);
  fprintf(file, "    </PCellData>\n");

  char pointscalars[ASD_BUFSIZ];
  strcsl(ASD_BUFSIZ, pointscalars, scalars);

  char pointvectors[ASD_BUFSIZ];
  strcsl(ASD_BUFSIZ, pointvectors, vectors);

  fprintf(file, "    <PPointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
    for (size_t s = 0; scalars[s]; ++s)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " format=\"%s\"/>\n",
              DFLOAT_VTK, scalars[s], format);

  if (vectors)
    for (size_t v = 0; vectors[v]; ++v)
      fprintf(file, "      <PDataArray type=\"%s\" Name=\"%s\""
                    " NumberOfComponents=\"3\" format=\"%s\"/>\n",
              DFLOAT_VTK, vectors[v], format);

  fprintf(file, "    </PPointData>\n");

  for (int s = 0; s < size; ++s)
  {

    if (directory)
      fprintf(file, "    <Piece Source=\"%s/" VTU_FORMAT "\"/>\n", directory,
              prefix, s);
    else
      fprintf(file, "    <Piece Source=\"" VTU_FORMAT "\"/>\n", prefix, s);
  }
  fprintf(file, "  </PUnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    ASD_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    ASD_LERROR("Error closing %s\n", filename);
  }
}

/** Utility function for writing a vector data array.
 *
 * \param [out] file     stream to write the vector to.
 * \param [in]  name     name of the vector.
 * \param [in]  binary   boolean indicating if the data should be written
 *                       in binary.
 * \param [in]  compress boolean indicating if the data should be
 *                       compressed.
 * \param [in]  K        Number of elements.
 * \param [in]  Np       Number of dofs in each element.
 * \param [in]  d        data
 * \param [in]  o1       1st component offset.
 * \param [in]  o2       2nd component offset.
 * \param [in]  o3       3rd component offset. If -1 then zero will be
 *                       used for the third component.
 * \param [in]  ot       total number of components.
 */
static void vtk_write_vector(FILE *file, const char *name, int binary,
                             int compress, iint_t K, int Np, const dfloat_t *d,
                             const int o1, const int o2, const int o3,
                             const int ot)
{
  const char *format = (binary) ? "binary" : "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " NumberOfComponents=\"3\" format=\"%s\">\n",
          DFLOAT_VTK, name, format);

  if (binary)
  {
    size_t v_sz = 3 * K * Np * sizeof(dfloat_t);
    dfloat_t *v = asd_malloc_aligned(v_sz);

    for (iint_t e = 0; e < K; ++e)
    {
      for (int n = 0; n < Np; ++n)
      {
        const size_t idv = n + e * Np;
        const size_t idd = n + e * Np * ot;
        v[3 * idv + 0] = (o1 >= 0) ? d[idd + Np * o1] : 0;
        v[3 * idv + 1] = (o2 >= 0) ? d[idd + Np * o2] : 0;
        v[3 * idv + 2] = (o3 >= 0) ? d[idd + Np * o3] : 0;
      }
    }

    fprintf(file, "          ");
    int rval = vtk_write_binary_data(compress, file, (char *)v, v_sz);
    fprintf(file, "\n");
    if (rval)
      ASD_WARNING("Error encoding %s", name);

    asd_free_aligned(v);
  }
  else
  {
    for (iint_t e = 0; e < K; ++e)
    {
      for (int n = 0; n < Np; ++n)
      {
        const size_t idd = n + e * Np * ot;
        dfloat_t a = (o1 >= 0) ? d[idd + Np * o1] : 0;
        dfloat_t b = (o2 >= 0) ? d[idd + Np * o2] : 0;
        dfloat_t c = (o3 >= 0) ? d[idd + Np * o3] : 0;

        fprintf(file,
                "         %" DFLOAT_FMTe " %" DFLOAT_FMTe " %" DFLOAT_FMTe "\n",
                a, b, c);
      }
    }
  }

  fprintf(file, "        </DataArray>\n");
}

/** Utility function for writing a vector data array.
 *
 * \param [out] file     stream to write the vector to.
 * \param [in]  name     name of the vector.
 * \param [in]  binary   boolean indicating if the data should be written
 *                       in binary.
 * \param [in]  compress boolean indicating if the data should be
 *                       compressed.
 * \param [in]  K        Number of elements.
 * \param [in]  Np       Number of dofs in each element.
 * \param [in]  d        data
 * \param [in]  o        offset.
 * \param [in]  ot       total number of components.
 */
static void vtk_write_scalar(FILE *file, const char *name, int binary,
                             int compress, iint_t K, int Np, const dfloat_t *d,
                             const int o, const int ot)
{
  const char *format = (binary) ? "binary" : "ascii";

  fprintf(file, "        <DataArray type=\"%s\" Name=\"%s\""
                " format=\"%s\">\n",
          DFLOAT_VTK, name, format);

  if (binary)
  {
    size_t v_sz = K * Np * sizeof(dfloat_t);
    dfloat_t *v = asd_malloc_aligned(v_sz);

    for (iint_t e = 0; e < K; ++e)
    {
      for (int n = 0; n < Np; ++n)
      {
        const size_t idv = n + e * Np;
        const size_t idd = n + e * Np * ot;
        v[idv] = d[idd + Np * o];
      }
    }

    fprintf(file, "          ");
    int rval = vtk_write_binary_data(compress, file, (char *)v, v_sz);
    fprintf(file, "\n");
    if (rval)
      ASD_WARNING("Error encoding %s", name);

    asd_free_aligned(v);
  }
  else
  {
    for (iint_t e = 0; e < K; ++e)
    {
      for (int n = 0; n < Np; ++n)
      {
        const int idd = n + e * Np * ot;
        fprintf(file, "         %" DFLOAT_FMTe "\n", d[idd + Np * o]);
      }
    }
  }

  fprintf(file, "        </DataArray>\n");
}

static void vtk_write_file(const int rank, const int size,
                           const char *directory, const char *prefix,
                           const int binary, const int compress,
                           const dfloat_t time, const int N, const iint_t K,
                           dfloat_t *vgeo, dfloat_t *q)
{
  const int Nq = N + 1;
  const int Np = (DIM == 2) ? Nq * Nq : Nq * Nq * Nq;
  const int Ncorners = (DIM == 2) ? 4 : 8;
  const iint_t Ncells = K * ipow(N, DIM);
  const iint_t Ntotal = K * Np;
  const char *const *scalars = FIELD_OUT_SCALARS;
  const char *const *vectors = FIELD_OUT_VECTORS;

  const int endian = asd_endian();

  const char *format = (binary) ? "binary" : "ascii";

  if (rank == 0)
    vtk_write_file_pvtu(size, directory, prefix, FIELD_OUT_SCALARS,
                        FIELD_OUT_VECTORS, binary, compress);

  char filename[ASD_BUFSIZ];
  if (directory)
    snprintf(filename, ASD_BUFSIZ, "%s/" VTU_FORMAT, directory, prefix, rank);
  else
    snprintf(filename, ASD_BUFSIZ, VTU_FORMAT, prefix, rank);

  ASD_VERBOSE("Writing file: '%s'", filename);
  FILE *file = fopen(filename, "w");

  if (file == NULL)
  {
    ASD_LERROR("Could not open %s for output!\n", filename);
    return;
  }

  fprintf(file, "<?xml version=\"1.0\"?>\n");
  fprintf(file, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");

  if (binary && compress)
    fprintf(file, " compressor=\"vtkZLibDataCompressor\"");

  if (endian == ASD_BIG_ENDIAN)
    fprintf(file, " byte_order=\"BigEndian\">\n");
  else
    fprintf(file, " byte_order=\"LittleEndian\">\n");

  fprintf(file, "  <UnstructuredGrid>\n");

  fprintf(file, "    <Piece NumberOfPoints=\"%jd\" NumberOfCells=\"%jd\">\n",
          (intmax_t)Ntotal, (intmax_t)Ncells);

  // {{{ Points
  fprintf(file, "      <Points>\n");

#if DIM == 3
  vtk_write_vector(file, "Position", binary, compress, K, Np, vgeo, VGEO_X,
                   VGEO_Y, VGEO_Z, NVGEO);
#else
  vtk_write_vector(file, "Position", binary, compress, K, Np, vgeo, VGEO_X,
                   VGEO_Y, -1, NVGEO);
#endif

  fprintf(file, "      </Points>\n");
  // }}}

  // {{{ Cells
  fprintf(file, "      <Cells>\n");

  // {{{ Connectivity
  fprintf(file, "        <DataArray type=\"%s\" Name=\"connectivity\""
                " format=\"%s\">\n",
          IINT_VTK, format);
  if (binary)
  {
    size_t cells_sz = Ncells * Ncorners * sizeof(iint_t);
    iint_t *cells = asd_malloc_aligned(cells_sz);

#if DIM == 2
    for (iint_t k = 0, i = 0; k < K; ++k)
    {
      for (int m = 0; m < N; ++m)
      {
        for (int n = 0; n < N; ++n)
        {
          cells[i++] = Np * k + Nq * (m + 0) + (n + 0);
          cells[i++] = Np * k + Nq * (m + 0) + (n + 1);
          cells[i++] = Np * k + Nq * (m + 1) + (n + 0);
          cells[i++] = Np * k + Nq * (m + 1) + (n + 1);
        }
      }
    }
#else
    for (iint_t k = 0, i = 0; k < K; ++k)
    {
      for (int l = 0; l < N; ++l)
      {
        for (int m = 0; m < N; ++m)
        {
          for (int n = 0; n < N; ++n)
          {
            cells[i++] = Np * k + Nq * Nq * (l + 0) + Nq * (m + 0) + (n + 0);
            cells[i++] = Np * k + Nq * Nq * (l + 0) + Nq * (m + 0) + (n + 1);
            cells[i++] = Np * k + Nq * Nq * (l + 0) + Nq * (m + 1) + (n + 0);
            cells[i++] = Np * k + Nq * Nq * (l + 0) + Nq * (m + 1) + (n + 1);
            cells[i++] = Np * k + Nq * Nq * (l + 1) + Nq * (m + 0) + (n + 0);
            cells[i++] = Np * k + Nq * Nq * (l + 1) + Nq * (m + 0) + (n + 1);
            cells[i++] = Np * k + Nq * Nq * (l + 1) + Nq * (m + 1) + (n + 0);
            cells[i++] = Np * k + Nq * Nq * (l + 1) + Nq * (m + 1) + (n + 1);
          }
        }
      }
    }
#endif

    fprintf(file, "          ");
    int rval = vtk_write_binary_data(compress, file, (char *)cells, cells_sz);
    fprintf(file, "\n");
    if (rval)
      ASD_WARNING("Error encoding cells");

    asd_free_aligned(cells);
  }
  else
  {
#if DIM == 2
    for (iint_t k = 0; k < K; ++k)
      for (int m = 0; m < N; ++m)
        for (int n = 0; n < N; ++n)
          fprintf(file, "          %8jd %8jd %8jd %8jd\n",
                  (intmax_t)Np * k + Nq * (m + 0) + (n + 0),
                  (intmax_t)Np * k + Nq * (m + 0) + (n + 1),
                  (intmax_t)Np * k + Nq * (m + 1) + (n + 0),
                  (intmax_t)Np * k + Nq * (m + 1) + (n + 1));
#else
    for (iint_t k = 0; k < K; ++k)
      for (int l = 0; l < N; ++l)
        for (int m = 0; m < N; ++m)
          for (int n = 0; n < N; ++n)
            fprintf(
                file, "          %8jd %8jd %8jd %8jd %8jd %8jd %8jd %8jd\n",
                (intmax_t)Np * k + Nq * Nq * (l + 0) + Nq * (m + 0) + (n + 0),
                (intmax_t)Np * k + Nq * Nq * (l + 0) + Nq * (m + 0) + (n + 1),
                (intmax_t)Np * k + Nq * Nq * (l + 0) + Nq * (m + 1) + (n + 0),
                (intmax_t)Np * k + Nq * Nq * (l + 0) + Nq * (m + 1) + (n + 1),
                (intmax_t)Np * k + Nq * Nq * (l + 1) + Nq * (m + 0) + (n + 0),
                (intmax_t)Np * k + Nq * Nq * (l + 1) + Nq * (m + 0) + (n + 1),
                (intmax_t)Np * k + Nq * Nq * (l + 1) + Nq * (m + 1) + (n + 0),
                (intmax_t)Np * k + Nq * Nq * (l + 1) + Nq * (m + 1) + (n + 1));
#endif
  }
  fprintf(file, "        </DataArray>\n");
  // }}}

  // {{{ Offsets
  fprintf(file, "        <DataArray type=\"%s\" Name=\"offsets\""
                " format=\"%s\">\n",
          IINT_VTK, format);
  fprintf(file, "          ");
  if (binary)
  {
    size_t offsets_sz = Ncells * sizeof(iint_t);
    iint_t *offsets = asd_malloc_aligned(offsets_sz);

    for (iint_t i = 1; i <= Ncells; ++i)
      offsets[i - 1] = Ncorners * i;

    int rval =
        vtk_write_binary_data(compress, file, (char *)offsets, offsets_sz);
    if (rval)
      ASD_WARNING("Error encoding offsets");

    asd_free_aligned(offsets);
  }
  else
  {
    for (iint_t i = 1, sk = 1; i <= Ncells; ++i, ++sk)
    {
      fprintf(file, " %8jd", (intmax_t)(Ncorners * i));
      if (!(sk % 20) && i != Ncells)
        fprintf(file, "\n          ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  // }}}

  // {{{ Types
  fprintf(file, "        <DataArray type=\"UInt8\" Name=\"types\""
                " format=\"%s\">\n",
          format);
  fprintf(file, "          ");
  if (binary)
  {
    size_t types_sz = Ncells * sizeof(uint8_t);
    uint8_t *types = asd_malloc_aligned(types_sz);

#if DIM == 2
    for (iint_t i = 0; i < Ncells; ++i)
      types[i] = 8; /* VTK_PIXEL */
#else
    for (iint_t i = 0; i < Ncells; ++i)
      types[i] = 11;        /* VTK_VOXEL */
#endif

    int rval = vtk_write_binary_data(compress, file, (char *)types, types_sz);
    if (rval)
      ASD_WARNING("Error encoding types");

    asd_free_aligned(types);
  }
  else
  {
    for (iint_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
#if DIM == 2
      fprintf(file, " 8"); /* VTK_PIXEL */
#else
      fprintf(file, " 11"); /* VTK_VOXEL */
#endif
      if (!(sk % 20) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  // }}}

  fprintf(file, "      </Cells>\n");
  // }}}

  // {{{ Cell Data
  fprintf(file, "      <CellData Scalars=\"time,mpirank\">\n");

  // {{{ time
  fprintf(file, "        <DataArray type=\"%s\" Name=\"time\""
                " format=\"%s\">\n",
          DFLOAT_VTK, format);
  fprintf(file, "          ");
  if (binary)
  {
    size_t times_sz = Ncells * sizeof(dfloat_t);
    dfloat_t *times = asd_malloc_aligned(times_sz);

    for (iint_t i = 0; i < Ncells; ++i)
      times[i] = time;

    int rval = vtk_write_binary_data(compress, file, (char *)times, times_sz);
    if (rval)
      ASD_WARNING("Error encoding times");

    asd_free_aligned(times);
  }
  else
  {
    for (iint_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %" DFLOAT_FMTe, time);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  // }}}

  // {{{ MPI rank
  const char *int_vtk_type;
  if (sizeof(int) == sizeof(int32_t))
    int_vtk_type = "Int32";
  else if (sizeof(int) == sizeof(int64_t))
    int_vtk_type = "Int64";
  else
    int_vtk_type = "UNKNOWN";
  fprintf(file, "        <DataArray type=\"%s\" Name=\"mpirank\""
                " format=\"%s\">\n",
          int_vtk_type, format);
  fprintf(file, "          ");
  if (binary)
  {
    size_t ranks_sz = Ncells * sizeof(int);
    int *ranks = asd_malloc_aligned(ranks_sz);

    for (iint_t i = 0; i < Ncells; ++i)
      ranks[i] = rank;

    int rval = vtk_write_binary_data(compress, file, (char *)ranks, ranks_sz);
    if (rval)
      ASD_WARNING("Error encoding ranks");

    asd_free_aligned(ranks);
  }
  else
  {
    for (iint_t i = 0, sk = 1; i < Ncells; ++i, ++sk)
    {
      fprintf(file, " %6jd", (intmax_t)rank);
      if (!(sk % 8) && i != (Ncells - 1))
        fprintf(file, "\n         ");
    }
  }
  fprintf(file, "\n");
  fprintf(file, "        </DataArray>\n");
  // }}}

  fprintf(file, "      </CellData>\n");
  // }}}

  // {{{ Point Data
  char pointscalars[ASD_BUFSIZ];
  strcsl(ASD_BUFSIZ, pointscalars, scalars);

  char pointvectors[ASD_BUFSIZ];
  strcsl(ASD_BUFSIZ, pointvectors, vectors);

  fprintf(file, "      <PointData Scalars=\"%s\" Vectors=\"%s\">\n",
          pointscalars, pointvectors);

  if (scalars)
    for (size_t s = 0; scalars[s]; ++s)
    {
      if (FIELD_OUT_SCALARS_OFF[s] >= 0)
        vtk_write_scalar(file, scalars[s], binary, compress, K, Np, q,
                         FIELD_OUT_SCALARS_OFF[s], NFIELDS);
      else
        vtk_write_scalar(file, scalars[s], binary, compress, K, Np, vgeo,
                         -FIELD_OUT_SCALARS_OFF[s] - 1, NVGEO);
    }

  if (vectors)
    for (size_t v = 0; vectors[v]; ++v)
      vtk_write_vector(file, vectors[v], binary, compress, K, Np, q,
                       FIELD_OUT_COMPONENTS_OFF[3 * v + 0],
                       FIELD_OUT_COMPONENTS_OFF[3 * v + 1],
                       FIELD_OUT_COMPONENTS_OFF[3 * v + 2], NFIELDS);

  fprintf(file, "      </PointData>\n");
  // }}}

  fprintf(file, "    </Piece>\n");

  fprintf(file, "  </UnstructuredGrid>\n");
  fprintf(file, "</VTKFile>\n");

  if (ferror(file))
  {
    ASD_LERROR("Error writing to %s\n", filename);
  }

  if (fclose(file))
  {
    ASD_LERROR("Error closing %s\n", filename);
  }
}

static void app_output_vtk(app_t *app, int s, dfloat_t time, const char *prefix)
{
  ASD_ROOT_INFO("");
  char fname[ASD_BUFSIZ];
  snprintf(fname, ASD_BUFSIZ, "%s%s_fields_%06d", prefix,
           app->prefs->output_prefix, s);
  ASD_ROOT_INFO("Writing VTK file \"%s\"", fname);

  const iint_t K = app->lvl->Klocal;
  const int N = app->lvl->N;
  const int Np = app->lvl->Np;

  /* Make sure the copy commands are placed in the command queue.  This is
   * to ensure that Q will be finished computing before we copy the data
   * to the host.
   */
  occaDeviceSetStream(app->device, app->cmdx);

  size_t local_vgeo_sz = NVGEO * K * Np * sizeof(dfloat_t);
  dfloat_t *vgeo = asd_malloc_aligned(local_vgeo_sz);
  occaCopyMemToPtr(vgeo, app->lvl->o_vgeo, local_vgeo_sz, occaNoOffset);

  size_t local_q_sz = NFIELDS * K * Np * sizeof(dfloat_t);
  dfloat_t *q = asd_malloc_aligned(local_q_sz);
  occaCopyMemToPtr(q, app->lvl->o_q, local_q_sz, occaNoOffset);

  vtk_write_file(app->prefs->rank, app->prefs->size, app->prefs->output_datadir,
                 fname, app->prefs->output_vtk_binary,
                 app->prefs->output_vtk_compress, time, N, K, vgeo, q);

  asd_free_aligned(vgeo);
  asd_free_aligned(q);
}
// }}}

// {{{ Time Step Initialization
typedef struct time_step_parameters
{
  dfloat_t dt;
  int nsteps;
  int nadp;
  int neng;
  int nerr;
  int nout;
} time_step_parameters_t;

static time_step_parameters_t *app_new_time_step_parameters(app_t *app,
                                                            dfloat_t dt)
{
  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("---- Init Time Step Parameters ---------------------------");

  time_step_parameters_t *tsp = asd_malloc(sizeof(time_step_parameters_t));

  double ddt;
  int result = asd_lua_global_function_call(
      app->prefs->L, "app.time_step_parameters", "d>diiiii", (double)dt, &ddt,
      &tsp->nsteps, &tsp->nadp, &tsp->nout, &tsp->neng, &tsp->nerr);
  tsp->dt = (dfloat_t)ddt;

  ASD_ABORT_IF_NOT(result == 0,
                   "problem with lua call to 'app.time_step_parameters': "
                   "should be a function that takes dt "
                   "and returns dt, nsteps, nadp, nout, neng, nerr");

  ASD_ROOT_INFO("final time  = %" DFLOAT_FMTe, (dfloat_t)tsp->nsteps * tsp->dt);
  ASD_ROOT_INFO("dt          = %" DFLOAT_FMTe, tsp->dt);
  ASD_ROOT_INFO("nsteps      = %d", tsp->nsteps);
  ASD_ROOT_INFO("nadp        = %d", tsp->nadp);
  ASD_ROOT_INFO("nout        = %d", tsp->nout);
  ASD_ROOT_INFO("neng        = %d", tsp->neng);
  ASD_ROOT_INFO("nerr        = %d", tsp->nerr);
  ASD_ROOT_INFO("----------------------------------------------------------");

  return tsp;
}

static dfloat_t app_dt(app_t *app)
{
  const int LDIM = app->prefs->kernel_reduce_ldim;
  level_t *lvl = app->lvl;

  dfloat_t *red_dt =
      asd_malloc(sizeof(dfloat_t) * app->prefs->kernel_reduce_max_copy);

  dfloat_t dt = DFLOAT_MAX;
  dfloat_t rank_dt = DFLOAT_MAX;

  iint_t n_reduce = lvl->Klocal * lvl->Np;

  if (n_reduce > 0)
  {
    iint_t n_groups = (n_reduce + LDIM - 1) / LDIM;
    n_groups = (n_groups + 8 - 1) / 8;

    occaKernelRun(lvl->compute_dt, occaIint(lvl->Klocal), lvl->o_vgeo, lvl->o_q,
                  lvl->o_red_buf[0]);

    int r = 0;
    do
    {
      n_groups = (n_reduce + LDIM - 1) / LDIM;
      n_groups = (n_groups + 8 - 1) / 8;

      occaDim global = {n_groups, 1, 1}, local = {LDIM, 1, 1};
      occaKernelSetWorkingDims(lvl->reduce_min, 1, local, global);

      occaKernelRun(lvl->reduce_min, occaIint(n_reduce),
                    lvl->o_red_buf[(r + 0) % 2], lvl->o_red_buf[(r + 1) % 2]);

      r = (r + 1) % 2;
      n_reduce = n_groups;
    } while (n_reduce > app->prefs->kernel_reduce_max_copy);

    occaCopyMemToPtr(red_dt, lvl->o_red_buf[r], n_reduce * sizeof(dfloat_t),
                     occaNoOffset);

    for (int n = 0; n < n_reduce; ++n)
      rank_dt = ASD_MIN(rank_dt, red_dt[n]);
  }

  ASD_MPI_CHECK(
      MPI_Allreduce(&rank_dt, &dt, 1, DFLOAT_MPI, MPI_MIN, app->prefs->comm));

  asd_free(red_dt);

  return dt;
}
// }}}

// {{{ Output
static dfloat_t app_energy(app_t *app)
{
  const int LDIM = app->prefs->kernel_reduce_ldim;
  level_t *lvl = app->lvl;

  dfloat_t *red_energy =
      asd_malloc(sizeof(dfloat_t) * app->prefs->kernel_reduce_max_copy);

  dfloat_t energy = 0;
  dfloat_t rank_energy = 0;

  iint_t n_reduce = lvl->Klocal * lvl->Np;
  if (n_reduce > 0)
  {

    iint_t n_groups = (n_reduce + LDIM - 1) / LDIM;
    n_groups = (n_groups + 8 - 1) / 8;

    occaKernelRun(lvl->compute_energy, occaIint(lvl->Klocal), app->lvl->o_w,
                  lvl->o_vgeo, lvl->o_q, lvl->o_red_buf[0]);

    int r = 0;
    do
    {
      n_groups = (n_reduce + LDIM - 1) / LDIM;
      n_groups = (n_groups + 8 - 1) / 8;

      occaDim global = {n_groups, 1, 1}, local = {LDIM, 1, 1};
      occaKernelSetWorkingDims(lvl->reduce_sum, 1, local, global);

      occaKernelRun(lvl->reduce_sum, occaIint(n_reduce),
                    lvl->o_red_buf[(r + 0) % 2], lvl->o_red_buf[(r + 1) % 2]);

      r = (r + 1) % 2;
      n_reduce = n_groups;
    } while (n_reduce > app->prefs->kernel_reduce_max_copy);

    occaCopyMemToPtr(red_energy, lvl->o_red_buf[r], n_reduce * sizeof(dfloat_t),
                     occaNoOffset);

    for (int n = 0; n < n_reduce; ++n)
      rank_energy += red_energy[n];
  }

  ASD_MPI_CHECK(MPI_Allreduce(&rank_energy, &energy, 1, DFLOAT_MPI, MPI_SUM,
                              app->prefs->comm));

  asd_free(red_energy);

  return DFLOAT_SQRT(energy);
}

static dfloat_t app_error(app_t *app, dfloat_t t)
{
  const int LDIM = app->prefs->kernel_reduce_ldim;
  level_t *lvl = app->lvl;

  dfloat_t *red_error =
      asd_malloc(sizeof(dfloat_t) * app->prefs->kernel_reduce_max_copy);

  dfloat_t error = 0;
  dfloat_t rank_error = 0;

  iint_t n_reduce = lvl->Klocal * lvl->Np;

  if (n_reduce > 0)
  {

    iint_t n_groups = (n_reduce + LDIM - 1) / LDIM;
    n_groups = (n_groups + 8 - 1) / 8;

    occaKernelRun(lvl->compute_error, occaIint(lvl->Klocal), lvl->o_EToT,
                  lvl->o_w, lvl->o_vgeo, occaDfloat(t), lvl->o_q,
                  lvl->o_red_buf[0]);

    int r = 0;
    do
    {
      n_groups = (n_reduce + LDIM - 1) / LDIM;
      n_groups = (n_groups + 8 - 1) / 8;

      occaDim global = {n_groups, 1, 1}, local = {LDIM, 1, 1};
      occaKernelSetWorkingDims(lvl->reduce_sum, 1, local, global);

      occaKernelRun(lvl->reduce_sum, occaIint(n_reduce),
                    lvl->o_red_buf[(r + 0) % 2], lvl->o_red_buf[(r + 1) % 2]);

      r = (r + 1) % 2;
      n_reduce = n_groups;
    } while (n_reduce > app->prefs->kernel_reduce_max_copy);

    occaCopyMemToPtr(red_error, lvl->o_red_buf[r], n_reduce * sizeof(dfloat_t),
                     occaNoOffset);

    for (int n = 0; n < n_reduce; ++n)
      rank_error += red_error[n];
  }

  ASD_MPI_CHECK(MPI_Allreduce(&rank_error, &error, 1, DFLOAT_MPI, MPI_SUM,
                              app->prefs->comm));

  asd_free(red_error);

  return DFLOAT_SQRT(error);
}

static void app_output(app_t *app, time_step_parameters_t *tsp, int s,
                       const char *prefix)
{
  static dfloat_t initial_energy = -1;
  static dfloat_t prev_energy = -1;

  occaDeviceSetStream(app->device, app->cmdx);
  occaDeviceFinish(app->device);

  int eng = (tsp->neng > 0) && ((s % tsp->neng == 0) || (s == tsp->nsteps));
  int err = (tsp->nerr > 0) && ((s % tsp->nerr == 0) || (s == tsp->nsteps));
  int out = (tsp->nout > 0) && ((s % tsp->nout == 0) || (s == tsp->nsteps));

  if (!(eng || err || out))
    return;

  const dfloat_t time = s * tsp->dt;

  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("---- Step ----");
  ASD_ROOT_INFO("time   = %" DFLOAT_FMTe, time);
  ASD_ROOT_INFO("step   =   %d", s);

  if (eng)
  {

    const dfloat_t energy = app_energy(app);

    if (initial_energy < 0)
      initial_energy = energy;

    if (prev_energy < 0)
      prev_energy = energy;

    ASD_ROOT_INFO("           energy        = %" DFLOAT_FMTe, energy);
    ASD_ROOT_INFO("normalized energy        = %" DFLOAT_FMTe,
                  energy / initial_energy);
    ASD_ROOT_INFO("normalized energy change = %" DFLOAT_FMTe,
                  (energy - prev_energy) / initial_energy);

    prev_energy = energy;
  }

  if (err)
  {
    const dfloat_t error = app_error(app, time);
    ASD_ROOT_INFO("           error  = %" DFLOAT_FMTe, error);
    if (eng)
      ASD_ROOT_INFO("normalized error  = %" DFLOAT_FMTe, error / prev_energy);
  }

  if (out)
  {
    if (app->prefs->output_mfem)
      app_output_mfem(app, s, prefix);

    if (app->prefs->output_vtk)
      app_output_vtk(app, s, time, prefix);
  }

  ASD_ROOT_INFO("--------------");
}

// }}}

// {{{ Adapt
static void app_mark_elements_rand(app_t *app)
{
  static asd_pcg32_random_t rng;
  static int init = 0;

  if (!init)
  {
#ifdef NONDETERMINISTIC_REFINEMENT
    asd_pcg32_srandom_r(&rng, time(NULL) ^ (intptr_t)&printf, (intptr_t)&exit);
#else
    asd_pcg32_srandom_r(&rng, 42u, 64u);
#endif
    init = 1;
  }

  int a = 1;
  int h = 0;
  p4est_gloidx_t fq0 = app->pxest->global_first_quadrant[app->pxest->mpirank];
  p4est_gloidx_t fq1 =
      app->pxest->global_first_quadrant[app->pxest->mpirank + 1];
  p4est_gloidx_t fqs = app->pxest->global_first_quadrant[app->pxest->mpisize];

  for (p4est_gloidx_t e = 0; e < fqs; ++e)
  {
    if ((e) % a == 0)
    {
      a = asd_pcg32_boundedrand_r(&rng, 12) + 6;
      h = asd_pcg32_boundedrand_r(&rng, 3);
    }

    if (e >= fq0 && e < fq1)
      app->lvl->EToA[e - fq0] = (int8_t)h;
  }
}

static void replace_quads(p4est_t *p4est, p4est_topidx_t which_tree,
                          int num_outgoing, p4est_quadrant_t *outgoing[],
                          int num_incoming, p4est_quadrant_t *incoming[])
{
  const quad_data_t *outd = (quad_data_t *)outgoing[0]->p.user_data;

  for (int i = 0; i < num_incoming; ++i)
  {
    quad_data_t *ind = (quad_data_t *)incoming[i]->p.user_data;
    ind->old_level = outd->old_level;
    ind->adapt_flag = ADAPT_TOUCHED;
  }
}

static int refine_quads(p4est_t *p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t *q)
{
  const quad_data_t *d = (quad_data_t *)q->p.user_data;

  if (d->adapt_flag == ADAPT_REFINE)
    return 1;
  else
    return 0;
}

static int coarsen_quads(p4est_t *p4est, p4est_topidx_t which_tree,
                         p4est_quadrant_t *children[])
{
  int retval = 1;

  for (int i = 0; i < P4EST_CHILDREN; ++i)
  {
    const quad_data_t *d = (quad_data_t *)children[i]->p.user_data;

    if (d->adapt_flag != ADAPT_COARSEN)
      retval = 0;
  }

  return retval;
}

static void get_partition(p4est_locidx_t *local_num_quadrants, int *KrecvNn,
                          int *KrecvNToR, p4est_locidx_t *Krecv, int *KsendNn,
                          int *KsendNToR, p4est_locidx_t *Ksend, MPI_Comm comm)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  p4est_gloidx_t *global_last_quad_index =
      asd_malloc_aligned(sizeof(p4est_gloidx_t) * 2 * size);

  p4est_locidx_t *num_quadrants_in_proc =
      asd_malloc_aligned(sizeof(p4est_locidx_t) * 2 * size);

  MPI_Allgather(local_num_quadrants, 2, P4EST_MPI_LOCIDX, num_quadrants_in_proc,
                2, P4EST_MPI_LOCIDX, comm);

  global_last_quad_index[2 * 0 + 0] = num_quadrants_in_proc[2 * 0 + 0] - 1;
  global_last_quad_index[2 * 0 + 1] = num_quadrants_in_proc[2 * 0 + 1] - 1;
  for (int r = 1; r < size; ++r)
  {
    global_last_quad_index[2 * r + 0] = num_quadrants_in_proc[2 * r + 0] +
                                        global_last_quad_index[2 * (r - 1) + 0];
    global_last_quad_index[2 * r + 1] = num_quadrants_in_proc[2 * r + 1] +
                                        global_last_quad_index[2 * (r - 1) + 1];
  }

  // {{{ Code from p4est_algorithms.c to determine partition
  p4est_gloidx_t my_begin =
      (rank == 0) ? 0 : (global_last_quad_index[2 * (rank - 1) + 1] + 1);
  p4est_gloidx_t my_end = global_last_quad_index[2 * rank + 1];

  *KrecvNn = 0;
  for (int from_proc = 0; from_proc < size; ++from_proc)
  {
    p4est_gloidx_t from_begin =
        (from_proc == 0)
            ? 0
            : (global_last_quad_index[2 * (from_proc - 1) + 0] + 1);
    p4est_gloidx_t from_end = global_last_quad_index[2 * from_proc + 0];

    if (from_begin <= my_end && from_end >= my_begin)
    {
      /* from_proc sends to me but may be empty */
      Krecv[from_proc] = (p4est_locidx_t)(ASD_MIN(my_end, from_end) -
                                          ASD_MAX(my_begin, from_begin) + 1);
      ASD_ASSERT(Krecv[from_proc] >= 0);
      if (from_proc != rank)
      {
        KrecvNToR[*KrecvNn] = from_proc;
        ++(*KrecvNn);
      }
    }
    else
    {
      Krecv[from_proc] = 0;
    }
  }

  my_begin = (rank == 0) ? 0 : (global_last_quad_index[2 * (rank - 1) + 0] + 1);
  my_end = global_last_quad_index[2 * rank + 0];

  *KsendNn = 0;
  for (int to_proc = 0; to_proc < size; ++to_proc)
  {
    p4est_gloidx_t to_begin =
        (to_proc == 0) ? 0
                       : (global_last_quad_index[2 * (to_proc - 1) + 1] + 1);
    p4est_gloidx_t to_end = global_last_quad_index[2 * to_proc + 1];

    if (to_begin <= my_end && to_end >= my_begin)
    {
      /* I send to to_proc which may be empty */
      Ksend[to_proc] = (p4est_locidx_t)(ASD_MIN(my_end, to_end) -
                                        ASD_MAX(my_begin, to_begin) + 1);
      ASD_ASSERT(Ksend[to_proc] >= 0);
      if (to_proc != rank)
      {
        KsendNToR[*KsendNn] = to_proc;
        ++(*KsendNn);
      }
    }
    else
    {
      /* I don't send to to_proc */
      Ksend[to_proc] = 0;
    }
  }
  // }}}

  // Turn Krecv and Ksend into offsets
  for (int r = size; r > 0; --r)
  {
    Krecv[r] = Krecv[r - 1];
    Ksend[r] = Ksend[r - 1];
  }
  Krecv[0] = Ksend[0] = 0;

  for (int r = 0; r < size; ++r)
  {
    Krecv[r + 1] += Krecv[r];
    Ksend[r + 1] += Ksend[r];
  }

  asd_free_aligned(global_last_quad_index);
  asd_free_aligned(num_quadrants_in_proc);

#if ASD_DEBUG
  ASD_LDEBUG("");
  ASD_LDEBUG("------ get_partition dump ---------");
  ASD_LDEBUG("");

  ASD_LDEBUG("KrecvNn = %d", *KrecvNn);
  for (int n = 0; n < *KrecvNn; ++n)
    ASD_LDEBUG("KrecvNToR[%d] = %d", n, KrecvNToR[n]);
  for (int r = 0; r < (size + 1); ++r)
    ASD_LDEBUG("Krecv[%d] = %jd", r, (intmax_t)Krecv[r]);

  ASD_LDEBUG("");
  ASD_LDEBUG("");

  ASD_LDEBUG("KsendNn = %d", *KsendNn);
  for (int n = 0; n < *KsendNn; ++n)
    ASD_LDEBUG("KsendNToR[%d] = %d", n, KsendNToR[n]);
  for (int r = 0; r < (size + 1); ++r)
    ASD_LDEBUG("Ksend[%d] = %jd", r, (intmax_t)Ksend[r]);

  ASD_LDEBUG("");
  ASD_LDEBUG("-----------------------------------");
#endif
}

static void app_adapt(app_t *app, time_step_parameters_t *tsp, int s)
{
  level_t *lvl = app->lvl;
  p4est_t *pxest = app->pxest;
  prefs_t *prefs = app->prefs;
  const int Nq = lvl->Nq;
  const int Nqk = lvl->Nqk;
  const int KblkV = app->prefs->kernel_KblkV;
  const int dim = 3;
  occaDim global = {1, 1, 1}, local = (occaDim){Nq, Nq, Nqk * KblkV};

  occaDeviceSetStream(app->device, app->cmdx);
  occaDeviceFinish(app->device);

  int adp = (s > 0) && (tsp->nadp > 0) && (s % tsp->nadp == 0);

  if (!adp)
    return;

  const dfloat_t time = s * tsp->dt;

  ASD_ROOT_INFO("");
  ASD_ROOT_INFO("---- Adapt ----");

  app_mark_elements_rand(app);

  // {{{ Fill level and adapt flags in quadrant data
  // Here we put the level and adapt flag into the quadrant user_data.
  //
  // The adapt flag is used to determine which mesh modifications should be
  // performed.
  //
  // The original quadrant level is added so that after we make all of the mesh
  // modifications we can determine which elements were adapted.
  {
    p4est_locidx_t e = 0;
    for (p4est_topidx_t t = pxest->first_local_tree;
         t <= pxest->last_local_tree; ++t)
    {
      p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
      sc_array_t *tquadrants = &tree->quadrants;

      const p4est_locidx_t Q = (p4est_locidx_t)tquadrants->elem_count;
      for (p4est_locidx_t q = 0; q < Q; ++q, ++e)
      {
        p4est_quadrant_t *quad = p4est_quadrant_array_index(tquadrants, q);
        quad_data_t *d = (quad_data_t *)quad->p.user_data;

        d->old_level = quad->level;
        d->adapt_flag = lvl->EToA[e];
      }
    }
  }
  // }}}

  p4est_refine_ext(pxest, 0, -1, refine_quads, NULL, replace_quads);
  p4est_coarsen_ext(pxest, 0, 0, coarsen_quads, NULL, replace_quads);
  p4est_balance_ext(pxest, P4EST_CONNECT_FULL, NULL, replace_quads);

  iint_t *EToOff = asd_malloc_aligned(sizeof(iint_t) * (lvl->Kmax + 1));

  // {{{ Fill EToOff with information about how to coarsen the mesh
  p4est_locidx_t num_refine[2] = {0, 0}; // #refined elems {pre,post}-partition
  {

    p4est_locidx_t old_e = 0, e = 0;
    for (p4est_topidx_t t = pxest->first_local_tree;
         t <= pxest->last_local_tree; ++t)
    {
      p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
      sc_array_t *tquadrants = &tree->quadrants;

      const p4est_locidx_t Q = (p4est_locidx_t)tquadrants->elem_count;
      for (p4est_locidx_t q = 0; q < Q;)
      {
        p4est_quadrant_t *quad = p4est_quadrant_array_index(tquadrants, q);
        const quad_data_t *d = (quad_data_t *)quad->p.user_data;

        EToOff[e] = old_e;

        if (quad->level > d->old_level)
        {
          // refined (which we will do later so this is a copy right now)
          ++num_refine[0];
          q += P4EST_CHILDREN;
          ++e;
          ++old_e;
        }
        else if (quad->level < d->old_level)
        {
          // coarsened
          ++q;
          ++e;
          old_e += P4EST_CHILDREN;
        }
        else
        {
          // nothing
          ++q;
          ++e;
          ++old_e;
        }
      }
    }
    EToOff[e] = old_e;

    ASD_ASSERT(
        pxest->local_num_quadrants - (P4EST_CHILDREN - 1) * num_refine[0] == e);

    device_async_ptr_to_mem(lvl->o_EToOff, EToOff, sizeof(iint_t) * (e + 1),
                            occaNoOffset);
  }
  // }}}

  p4est_locidx_t coarsened_local_num_quadrants[2];

  coarsened_local_num_quadrants[0] =
      pxest->local_num_quadrants - (P4EST_CHILDREN - 1) * num_refine[0];

  p4est_partition(pxest, 1, NULL);

  ASD_ABORT_IF(pxest->local_num_quadrants > lvl->Kmax,
               "Kmax exceeded in adaptive step");

  // {{{ fill num_refine post
  {
    for (p4est_topidx_t t = pxest->first_local_tree;
         t <= pxest->last_local_tree; ++t)
    {
      p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
      sc_array_t *tquadrants = &tree->quadrants;

      const p4est_locidx_t Q = (p4est_locidx_t)tquadrants->elem_count;
      for (p4est_locidx_t q = 0; q < Q;)
      {
        p4est_quadrant_t *quad = p4est_quadrant_array_index(tquadrants, q);
        const quad_data_t *d = (quad_data_t *)quad->p.user_data;

        if (quad->level > d->old_level)
        {
          q += P4EST_CHILDREN;
          ++num_refine[1];
        }
        else
        {
          ++q;
        }
      }
    }
  }
  // }}}

  coarsened_local_num_quadrants[1] =
      pxest->local_num_quadrants - (P4EST_CHILDREN - 1) * num_refine[1];

  // {{{ get partition information
  int KrecvNn, KsendNn;
  int *KrecvNToR = asd_malloc_aligned(sizeof(int) * pxest->mpisize);
  int *KsendNToR = asd_malloc_aligned(sizeof(int) * pxest->mpisize);
  p4est_locidx_t *Krecv =
      asd_malloc_aligned(sizeof(p4est_locidx_t) * (pxest->mpisize + 1));
  p4est_locidx_t *Ksend =
      asd_malloc_aligned(sizeof(p4est_locidx_t) * (pxest->mpisize + 1));
  get_partition(coarsened_local_num_quadrants, &KrecvNn, KrecvNToR, Krecv,
                &KsendNn, KsendNToR, Ksend, pxest->mpicomm);

  const p4est_locidx_t Ksendremote =
      Ksend[pxest->mpisize] -
      (Ksend[pxest->mpirank + 1] - Ksend[pxest->mpirank]);
  const p4est_locidx_t Krecvremote =
      Krecv[pxest->mpisize] -
      (Krecv[pxest->mpirank + 1] - Krecv[pxest->mpirank]);
  // }}}

  for (int n = 0; n < KrecvNn; ++n)
  {
    const int r = KrecvNToR[n];

    // Adjust the offset to remove the local elements
    const p4est_locidx_t Koffset =
        Krecv[r] - ((r > pxest->mpirank)
                        ? (Krecv[pxest->mpirank + 1] - Krecv[pxest->mpirank])
                        : 0);

    MPI_Irecv(lvl->q_recv_buf + NFIELDS * lvl->Np * Koffset,
              NFIELDS * lvl->Np * (Krecv[r + 1] - Krecv[r]), DFLOAT_MPI, r, 333,
              pxest->mpicomm, lvl->q_recv_requests + n);
  }

  global =
      (occaDim){(coarsened_local_num_quadrants[0] + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->coarsen_fields, dim, local, global);
  occaKernelRun(lvl->coarsen_fields, occaIint(coarsened_local_num_quadrants[0]),
                occaIint(Ksend[pxest->mpirank]),
                occaIint(Ksend[pxest->mpirank + 1] - 1), lvl->o_EToOff,
                lvl->o_Pb, lvl->o_Pt, lvl->o_q, lvl->o_q_buf, lvl->o_rhsq);

  device_async_mem_to_ptr(lvl->q_send, lvl->o_q_buf,
                          sizeof(dfloat_t) * NFIELDS * lvl->Np * Ksendremote,
                          occaNoOffset);

  p4est_ghost_destroy(app->ghost);
  app->ghost = p4est_ghost_new(pxest, P4EST_CONNECT_FULL);

  occaDeviceFinish(app->device);

  memcpy(lvl->q_send_buf, lvl->q_send,
         sizeof(dfloat_t) * NFIELDS * lvl->Np * Ksendremote);

  for (int n = 0; n < KsendNn; ++n)
  {
    const int r = KsendNToR[n];

    // Adjust the offset to remove the local elements
    const p4est_locidx_t Koffset =
        Ksend[r] - ((r > pxest->mpirank)
                        ? (Ksend[pxest->mpirank + 1] - Ksend[pxest->mpirank])
                        : 0);

    MPI_Isend(lvl->q_send_buf + NFIELDS * lvl->Np * Koffset,
              NFIELDS * lvl->Np * (Ksend[r + 1] - Ksend[r]), DFLOAT_MPI, r, 333,
              app->pxest->mpicomm, lvl->q_send_requests + n);
  }

  // {{{ Fill EToOff with information about how to refine the mesh
  {

    p4est_locidx_t old_e = 0, e = 0;
    for (p4est_topidx_t t = pxest->first_local_tree;
         t <= pxest->last_local_tree; ++t)
    {
      p4est_tree_t *tree = p4est_tree_array_index(pxest->trees, t);
      sc_array_t *tquadrants = &tree->quadrants;

      const p4est_locidx_t Q = (p4est_locidx_t)tquadrants->elem_count;
      for (p4est_locidx_t q = 0; q < Q;)
      {
        p4est_quadrant_t *quad = p4est_quadrant_array_index(tquadrants, q);
        const quad_data_t *d = (quad_data_t *)quad->p.user_data;

        EToOff[old_e] = e;

        if (quad->level > d->old_level)
        {
          // refined
          q += P4EST_CHILDREN;
          e += P4EST_CHILDREN;
          ++old_e;
        }
        else if (quad->level < d->old_level)
        {
          // coarsened (which we already did so this becomes a copy)
          ++q;
          ++e;
          ++old_e;
        }
        else
        {
          // nothing
          ++q;
          ++e;
          ++old_e;
        }
      }
    }
    EToOff[old_e] = e;

    ASD_ASSERT(pxest->local_num_quadrants -
                   (P4EST_CHILDREN - 1) * num_refine[1] ==
               old_e);

    device_async_ptr_to_mem(lvl->o_EToOff, EToOff, sizeof(iint_t) * (e + 1),
                            occaNoOffset);
  }
  // }}}

  mesh_t *mesh = mesh_new(prefs, pxest, app->ghost);
  level_get_mesh_constants(lvl, mesh);
  level_set_working_dims(lvl, prefs);

  MPI_Waitall(KrecvNn, lvl->q_recv_requests, lvl->q_recv_statuses);
  memcpy(lvl->q_recv, lvl->q_recv_buf,
         sizeof(dfloat_t) * NFIELDS * Krecvremote * lvl->Np);

  device_async_ptr_to_mem(lvl->o_q_buf, lvl->q_recv,
                          sizeof(dfloat_t) * NFIELDS * lvl->Np * Krecvremote,
                          occaNoOffset);

  global =
      (occaDim){(coarsened_local_num_quadrants[1] + KblkV - 1) / KblkV, 1, 1};
  occaKernelSetWorkingDims(lvl->refine_and_fill_fields, dim, local, global);
  occaKernelRun(
      lvl->refine_and_fill_fields, occaIint(coarsened_local_num_quadrants[1]),
      occaIint(Krecv[pxest->mpirank]), occaIint(Krecv[pxest->mpirank + 1] - 1),
      lvl->o_EToOff, lvl->o_Ib, lvl->o_It, lvl->o_q_buf, lvl->o_rhsq, lvl->o_q);

  // Since we used rhsq as a buffer we need to zero it for the time stepping
  occaKernelRun(lvl->zero_fields, occaIint(lvl->Ktotal), lvl->o_rhsq);

  occaDeviceSetStream(app->device, app->copy);
  level_get_mesh(lvl, mesh, prefs, pxest, app->ghost, app->device);

  asd_free_aligned(KrecvNToR);
  asd_free_aligned(KsendNToR);
  asd_free_aligned(Krecv);
  asd_free_aligned(Ksend);

  MPI_Waitall(KsendNn, lvl->q_send_requests, lvl->q_send_statuses);

  occaDeviceSetStream(app->device, app->copy);
  occaDeviceFinish(app->device);
  occaDeviceSetStream(app->device, app->cmdx);
  occaDeviceFinish(app->device);

  asd_free_aligned(EToOff);
  mesh_free(mesh);
  asd_free(mesh);

  ASD_ROOT_INFO("time   = %" DFLOAT_FMTe, time);
  ASD_ROOT_INFO("step   =   %d", s);
  ASD_ROOT_INFO("--------------");
}
// }}}

// {{{ Run
static void app_advection_step(app_t *app, int step, dfloat_t dt)
{
  level_t *lvl = app->lvl;
  const int Nq = lvl->Nq;
  const int Nqk = lvl->Nqk;
  const int KblkV = app->prefs->kernel_KblkV;

  const int dim = 3;
  occaDim global = {1, 1, 1}, local = (occaDim){Nq, Nq, Nqk * KblkV};

  for (int stage = 0; stage < NRKSTAGES; ++stage)
  {
    const dfloat_t rkt = step * dt + RKC[stage] * dt;

    /* [host] post the MPI receives */
    for (int n = 0; n < lvl->Nn; ++n)
    {
      const int r = lvl->NToR[n];
      const int Kr =
          app->ghost->proc_offsets[r + 1] - app->ghost->proc_offsets[r];
      MPI_Irecv(lvl->q_recv_buf +
                    NFIELDS * lvl->Np * app->ghost->proc_offsets[r],
                NFIELDS * lvl->Np * Kr, DFLOAT_MPI, r, 666, app->pxest->mpicomm,
                lvl->q_recv_requests + n);
    }

    /* [cmdx] wait to ensure that the update is finished */
    occaDeviceSetStream(app->device, app->cmdx);
    occaDeviceFinish(app->device);

    /* [copy] get mirror fields and copy to host */
    occaDeviceSetStream(app->device, app->copy);
    occaKernelRun(lvl->get_mirror_fields, occaIint(lvl->Kmirror), lvl->o_MToE,
                  lvl->o_q, lvl->o_q_buf);
    device_async_mem_to_ptr(lvl->q_send, lvl->o_q_buf,
                            sizeof(dfloat_t) * NFIELDS * lvl->Kmirror * lvl->Np,
                            occaNoOffset);

    /* [cmdx] launch the volume kernel on internal elements minus the mirrors */
    occaDeviceSetStream(app->device, app->cmdx);
    global = (occaDim){(lvl->Kintra + KblkV - 1) / KblkV, 1, 1};
    occaKernelSetWorkingDims(lvl->volume_advection, dim, local, global);
    occaKernelRun(lvl->volume_advection, occaIint(lvl->Kintra), lvl->o_IToE,
                  occaDfloat(rkt), lvl->o_D, lvl->o_w, lvl->o_EToT, lvl->o_vgeo,
                  lvl->o_sgeo, lvl->o_EToB, lvl->o_EToE, lvl->o_EToF,
                  lvl->o_EToO, lvl->o_q, lvl->o_rhsq);

    /* [copy] wait on the mirror copy */
    occaDeviceSetStream(app->device, app->copy);
    occaDeviceFinish(app->device);

    /* [host] post the MPI send, wait on the receives, then copy to local mem
     */
    MPI_Waitall(lvl->Nn, lvl->q_send_requests, lvl->q_send_statuses);
    memcpy(lvl->q_send_buf, lvl->q_send,
           sizeof(dfloat_t) * NFIELDS * lvl->Kmirror * lvl->Np);

    for (int n = 0; n < lvl->Nn; ++n)
    {
      const int r = lvl->NToR[n];
      const int Kr = app->ghost->mirror_proc_offsets[r + 1] -
                     app->ghost->mirror_proc_offsets[r];
      MPI_Isend(lvl->q_send_buf +
                    NFIELDS * lvl->Np * app->ghost->mirror_proc_offsets[r],
                NFIELDS * lvl->Np * Kr, DFLOAT_MPI, r, 666, app->pxest->mpicomm,
                lvl->q_send_requests + n);
    }

    MPI_Waitall(lvl->Nn, lvl->q_recv_requests, lvl->q_recv_statuses);
    memcpy(lvl->q_recv, lvl->q_recv_buf,
           sizeof(dfloat_t) * NFIELDS * lvl->Kghost * lvl->Np);

    /* [copy] send ghosts back to the device and wait for it to be there */
    device_async_ptr_to_mem(lvl->o_q_buf, lvl->q_recv,
                            NFIELDS * sizeof(dfloat_t) * lvl->Kghost * lvl->Np,
                            occaNoOffset);
    occaKernelRun(lvl->set_ghost_fields, occaIint(lvl->Kghost), lvl->o_GToE,
                  lvl->o_q_buf, lvl->o_q);
    occaDeviceFinish(app->device);

    occaDeviceSetStream(app->device, app->cmdx);
    /* [cmdx] launch the volume on mirrors and ghosts */
    global = (occaDim){(lvl->Kuniqmirror + KblkV - 1) / KblkV, 1, 1};
    occaKernelSetWorkingDims(lvl->volume_advection, dim, local, global);
    occaKernelRun(lvl->volume_advection, occaIint(lvl->Kuniqmirror),
                  lvl->o_UMToE, occaDfloat(rkt), lvl->o_D, lvl->o_w,
                  lvl->o_EToT, lvl->o_vgeo, lvl->o_sgeo, lvl->o_EToB,
                  lvl->o_EToE, lvl->o_EToF, lvl->o_EToO, lvl->o_q, lvl->o_rhsq);

    /* [cmdx] launch mortar kernel */
    for (iint_t facegroup = 0; facegroup < P4EST_FACES / 2; ++facegroup)
      occaKernelRun(lvl->mortar_advection, occaIint(facegroup),
                    occaIint(lvl->Nmortar), lvl->o_w, lvl->o_Pb, lvl->o_Pt,
                    lvl->o_MFToEM, lvl->o_MFToFM, lvl->o_MFToEP, lvl->o_MFToFP,
                    lvl->o_MFToOP, lvl->o_vgeo, lvl->o_sgeo, lvl->o_q,
                    lvl->o_rhsq);

    /* [cmdx] launch the update kernels */
    occaKernelRun(lvl->update_advection, occaIint(lvl->Klocal), occaDfloat(dt),
                  occaDfloat(RKA[(stage + 1) % NRKSTAGES]),
                  occaDfloat(RKB[stage]), lvl->o_q, lvl->o_rhsq);
  }
}

static void run(app_t *app)
{
  level_t *lvl = app->lvl;

  // ics
  occaKernelRun(lvl->compute_ics, occaIint(lvl->Ktotal), lvl->o_EToT,
                lvl->o_vgeo, occaDfloat(0), lvl->o_q);

  // zero rhs
  occaKernelRun(lvl->zero_fields, occaIint(lvl->Ktotal), lvl->o_rhsq);

  dfloat_t dt = app_dt(app);
  time_step_parameters_t *tsp = app_new_time_step_parameters(app, dt);

  for (int step = 0; step < tsp->nsteps; ++step)
  {
    app_output(app, tsp, step, "");

    app_adapt(app, tsp, step);

    app_output(app, tsp, step, "postadapt_");

    ASD_ROOT_LDEBUG("starting step %d", step);

    app_advection_step(app, step, tsp->dt);

    ASD_ROOT_LDEBUG("done with step %d", step);
  }

  app_output(app, tsp, tsp->nsteps, "");
  asd_free(tsp);
  ASD_ROOT_INFO("----------------------------------------------------------");
}
// }}}

// {{{ Main
static void usage()
{
  const char *help_text =
      "  " APP_NAME " [options] prefs_file\n"
      "\n"
      "  there are four possible options to this program, some of which \n"
      "  have multiple names:\n"
      "\n"
      "    -h -? --help --HELP\n"
      "    -d --debug\n"
      "    -D --devices\n"
      "    -V --version\n"
      "    -v --verbose  (which may be repeated for more verbosity)\n"
      "\n";
  ASD_ROOT_INFO(help_text);
}

int main(int argc, char *argv[])
{
  int status = EXIT_SUCCESS;
  int rank;
  MPI_Comm comm = MPI_COMM_WORLD;

  ASD_MPI_CHECK(MPI_Init(&argc, &argv));
  ASD_MPI_CHECK(MPI_Comm_rank(comm, &rank));

  //
  // parse command line
  //
  void *options = asd_gopt_sort(
      &argc, (const char **)argv,
      asd_gopt_start(asd_gopt_option('h', 0, asd_gopt_shorts('h', '?'),
                                     asd_gopt_longs("help", "HELP")),
                     asd_gopt_option('d', 0, asd_gopt_shorts('d'),
                                     asd_gopt_longs("debug")),
                     asd_gopt_option('D', 0, asd_gopt_shorts('D'),
                                     asd_gopt_longs("devices")),
                     asd_gopt_option('V', 0, asd_gopt_shorts('V'),
                                     asd_gopt_longs("version")),
                     asd_gopt_option('v', ASD_GOPT_REPEAT, asd_gopt_shorts('v'),
                                     asd_gopt_longs("verbose"))));

  if (asd_gopt(options, 'd'))
    debug(comm);

  if (asd_gopt(options, 'h'))
  {
    usage();
    goto finalize;
  }

  if (asd_gopt(options, 'V'))
    ASD_ROOT_INFO("app Version: %s", "unknown");

  if (asd_gopt(options, 'D'))
    occaPrintAvailableDevices();

  int verbosity = (int)asd_gopt(options, 'v');

  if (argc != 2)
  {
    ASD_LERROR("Unexpected number of arguments.");
    usage();
    status = EXIT_FAILURE;
    goto finalize;
  }

  //
  // initialize
  //

  int loglevel = init_libs(comm, verbosity);
  app_t *app = app_new(argv[1], comm, loglevel);

  print_precision();

  //
  // run
  //
  run(app);

  //
  // cleanup
  //

  app_free(app);
  asd_free(app);

finalize:
  sc_finalize();
  ASD_MPI_CHECK(MPI_Finalize());

  asd_gopt_free(options);

  return status;
}
// }}}
