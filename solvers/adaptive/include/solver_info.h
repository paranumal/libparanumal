#ifndef SOLVER_INFO_H
#define SOLVER_INFO_H 1

// {{{ Solver Info
// Number of dimensions in p4est connectivity
#define NX 3

typedef struct quad_data
{
  int8_t old_level;
  int8_t adapt_flag;
} quad_data_t;

#define FIELD_UX 0
#define FIELD_UY 1
#define FIELD_UZ 2
#define NFIELDS 3
const char *const FIELD_OUT_SCALARS[] = {NULL};
const int FIELD_OUT_SCALARS_OFF[] = {};
const char *const FIELD_OUT_VECTORS[] = {"u", NULL};
const char *const FIELD_OUT_COMPONENTS[] = {"ux", "uy", "uz", NULL};
const int FIELD_OUT_COMPONENTS_OFF[] = {FIELD_UX, FIELD_UY, FIELD_UZ};

#if DIM == 3
#define VGEO_RX 0
#define VGEO_SX 1
#define VGEO_TX 2
#define VGEO_RY 3
#define VGEO_SY 4
#define VGEO_TY 5
#define VGEO_RZ 6
#define VGEO_SZ 7
#define VGEO_TZ 8
#define VGEO_J 9
#define VGEO_X 10
#define VGEO_Y 11
#define VGEO_Z 12
#define NVGEO 13
const char *const VGEO_NAMES[] = {"rx", "sx", "tx", "ry", "sy", "ty", "rz",
                                  "sz", "tz", "J",  "x",  "y",  "z",  NULL};


#define GGEO_RR 0
#define GGEO_RS 1
#define GGEO_RT 2
#define GGEO_SS 3
#define GGEO_ST 4
#define GGEO_TT 5
#define GGEO_JW 6
#define NGGEO 7

#else
#define VGEO_RX 0
#define VGEO_SX 1
#define VGEO_RY 2
#define VGEO_SY 3
#define VGEO_J 4
#define VGEO_X 5
#define VGEO_Y 6
#define NVGEO 7
const char *const VGEO_NAMES[] = {"rx", "sx", "ry", "sy", "J", "x", "y", NULL};
#endif

#if DIM == 3
#define SGEO_NX 0
#define SGEO_NY 1
#define SGEO_NZ 2
#define SGEO_SJ 3
#define NSGEO 4
const char *const SGEO_NAMES[] = {"nx", "ny", "nz", "sJ", NULL};
#else
#define SGEO_NX 0
#define SGEO_NY 1
#define SGEO_SJ 2
#define NSGEO 3
const char *const SGEO_NAMES[] = {"nx", "ny", "sJ", NULL};
#endif

/* BC types */
#define BC_SKIP -1
#define BC_NONE 0
#define BC_DEFAULT 1
const char *const BCs[] = {"NONE", "DIRICHLET", NULL};

/* BC types */
#define ADAPT_COARSEN 0
#define ADAPT_NONE 1
#define ADAPT_REFINE 2
#define ADAPT_TOUCHED 3
#define NADAPT 4
// }}}

#endif
