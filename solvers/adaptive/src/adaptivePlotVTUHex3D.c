/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Ali Karakus, Lucas Wilcox

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "adaptive.h"

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
    dfloat_t *v = (dfloat_t*)asd_malloc_aligned(v_sz);

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
    dfloat_t *v = (dfloat_t*)asd_malloc_aligned(v_sz);

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
    iint_t *cells = (iint_t*)asd_malloc_aligned(cells_sz);

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
    iint_t *offsets = (iint_t*)asd_malloc_aligned(offsets_sz);

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
    uint8_t *types = (uint8_t*)asd_malloc_aligned(types_sz);

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
    dfloat_t *times = (dfloat_t*)asd_malloc_aligned(times_sz);

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
    int *ranks = (int*)asd_malloc_aligned(ranks_sz);

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

void adaptivePlotVTUHex3D(adaptive_t *adaptive, level_t *level,
                          int s, dfloat_t time, const char *prefix,
                          occa::memory &o_fields)
{
  ASD_ROOT_INFO("");
  char fname[ASD_BUFSIZ];
  snprintf(fname, ASD_BUFSIZ, "%s%s_fields_%06d", prefix, "adaptive_", s);
  ASD_ROOT_INFO("Writing VTK file \"%s\"", fname);

  const iint_t K = level->Klocal;
  const int N = level->N;
  const int Np = level->Np;

  /* Make sure the copy commands are placed in the command queue.  This is
   * to ensure that Q will be finished computing before we copy the data
   * to the host.
   */
  //  occaDeviceSetStream(app->device, app->cmdx);

  size_t local_vgeo_sz = NVGEO * K * Np * sizeof(dfloat_t);
  dfloat_t *vgeo = (dfloat_t*)asd_malloc_aligned(local_vgeo_sz);
  level->o_vgeo.copyTo(vgeo, local_vgeo_sz, 0); // 0 no offset

  size_t local_q_sz = NFIELDS * K * Np * sizeof(dfloat_t);
  dfloat_t *fields = (dfloat_t*)asd_malloc_aligned(local_q_sz);
  o_fields.copyTo(fields, local_q_sz, 0); // 0 no offset

  vtk_write_file(adaptive->rank, adaptive->size,
                 ".", // output dir
                 fname,
                 1, // output vtk binary
                 1, // output vtk compress
                 time, N, K, vgeo, fields);

  asd_free_aligned(vgeo);
  asd_free_aligned(fields);
}
// }}}
