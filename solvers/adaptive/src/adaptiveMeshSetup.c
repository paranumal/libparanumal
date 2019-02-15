#include "adaptive.h"

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
  mesh_t *mesh = (mesh_t *)user_data;

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
  mesh_t *mesh = (mesh_t *)user_data;

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
    p4est_iter_face_side_t *side = (p4est_iter_face_side_t *)sc_array_index(&info->sides, 0);
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
    p4est_iter_face_side_t *side0 = (p4est_iter_face_side_t *)sc_array_index(&info->sides, 0);
    p4est_iter_face_side_t *side1 = (p4est_iter_face_side_t *)sc_array_index(&info->sides, 1);
    p4est_tree_t *tree0 = (p4est_tree_t *)
        p4est_tree_array_index(info->p4est->trees, side0->treeid);
    p4est_tree_t *tree1 = (p4est_tree_t *)
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
        p4est_iter_face_side_t * sidetmp;
        sidetmp = side0;
        side0 = side1;
        side1 = sidetmp;

        p4est_tree_t * treetmp;
        treetmp = tree0;
        tree0 = tree1;
        tree1 = treetmp;
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

mesh_t *adaptiveMeshSetup(p4est_t *pxest, p4est_ghost_t *ghost, int *brick_n, int *brick_p, int *brick_TToC, int N)
{
  mesh_t *mesh = (mesh_t*)asd_malloc(sizeof(mesh_t));

  mesh->N = N;
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
  mesh->brick = 1;
  mesh->brick_n = brick_n;
  mesh->brick_p = brick_p;
  mesh->brick_TToC = brick_TToC;

  // {{{ continuous to discontinuous node map
  asd_dictionary_init(&mesh->CToD);
  mesh->EToC = (iint_t*)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToP = (iint_t*)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->DToC = (p4est_gloidx_t*)asd_malloc_aligned(sizeof(p4est_gloidx_t) * Ktotal * Np);
  mesh->CToD_starts = (iint_t*)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Np);
  mesh->CToD_indices = (iint_t*)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Np);

  // {{{ get EToC (from p4est_plex.c)
  {
    iint_t **mirror_EToC;

    for (iint_t e = 0; e < Klocal; ++e)
      mesh->EToC[e] = lnodes->face_code[e];

    mirror_EToC =
        (iint_t**)asd_malloc_aligned(sizeof(iint_t *) * ghost->mirrors.elem_count);
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

    mirror_DToC = (p4est_gloidx_t **)asd_malloc_aligned(sizeof(p4est_gloidx_t *) *
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

  int *hanging = (int*) asd_malloc_aligned(sizeof(int) * Np);
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
  mesh->IToE = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Klocal);
  mesh->MToE = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Kmirror);
  mesh->UMToE = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Kuniqmirror);
  mesh->GToE = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Kghost);

  mesh->EToL = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToT = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToX = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToY = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToZ = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal);
  mesh->EToB = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  for (iint_t n = 0; n < Ktotal * Nfaces; ++n)
    mesh->EToB[n] = BC_NONE;

  mesh->EToE = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->EToF = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->EToO = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);

  mesh->m = 0;
  mesh->MFToEM = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToFM = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToEP =
      (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces * P4EST_HALF);
  mesh->MFToFP = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);
  mesh->MFToOP = (iint_t *)asd_malloc_aligned(sizeof(iint_t) * Ktotal * Nfaces);

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

void mesh_free(mesh_t *mesh)
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

  asd_free(mesh);
}
