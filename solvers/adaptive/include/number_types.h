#ifndef NUMBER_TYPES_H
#define NUMBER_TYPES_H 1

#include "types.h"

// {{{ Number Types

#define INT8_MAX_DIGITS "3"
#define INT_MAX_DIGITS "10"
#define INT32_MAX_DIGITS "10"
#define INT64_MAX_DIGITS "19"

#define USE_IINT_SIGNED
#ifdef USE_IINT_SIGNED
typedef int iint_t;
#define occaIint(x) occaInt(x)
#define IINT(x) (x)
#define occa_iint_name "int"
#define IINT_VTK "Int32"
#define IINT_SCN SCNd32
#define IINT_PRI PRId32
#define IINT_MAX INT32_MAX
#define IINT_MAX_DIGITS INT32_MAX_DIGITS
#define asd_dictionary_insert_iint asd_dictionary_insert_int
#define asd_dictionary_allprefixed_iint asd_dictionary_allprefixed_int
#else
typedef unsigned int iint_t;
#define occaIint(x) occaUInt((iint_t)x)
#define IINT(x) ASD_APPEND(x, u)
#define occa_iint_name "unsigned int"
#define IINT_VTK "UInt32"
#define IINT_SCN SCNu32
#define IINT_PRI PRIu32
#define IINT_MAX UINT32_MAX
#define IINT_MAX_DIGITS INT32_MAX_DIGITS
#define asd_dictionary_insert_iint asd_dictionary_insert_uint
#define asd_dictionary_allprefixed_iint asd_dictionary_allprefixed_uint
#endif

#define P4EST_LOCIDX_MAX_DIGITS INT32_MAX_DIGITS
#define P4EST_GLOIDX_MAX_DIGITS INT64_MAX_DIGITS
#define P4EST_GLOIDX_SCN SCNd64
#define P4EST_GLOIDX_PRI PRId64

#ifdef USE_DFLOAT_DOUBLE

typedef double dfloat_t;
#define occaDfloat occaDouble
#define DFLOAT_MAX DBL_MAX
#define DFLOAT_FMTe "24.16e"
#define DFLOAT_MPI MPI_DOUBLE
#define DFLOAT_VTK "Float64"
#define DFLOAT_SQRT sqrt


#else
typedef float dfloat_t;
#define occaDfloat occaFloat
#define DFLOAT_MAX FLT_MAX
#define DFLOAT_FMTe "24.16e"
#define DFLOAT_MPI MPI_FLOAT
#define DFLOAT_VTK "Float32"
#define DFLOAT_SQRT sqrtf
#endif

#define GB (1e9)
#define GiB (1024 * 1024 * 1024)
// }}}

#endif
