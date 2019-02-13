#ifndef MESH_H
#define MESH_H 1


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

mesh_t *mesh_new(p4est_t *pxest, p4est_ghost_t *ghost, int *brick_n, int *brick_p, int *brick_TToC, int N);

void mesh_free(mesh_t *mesh);

// }}}

#endif
