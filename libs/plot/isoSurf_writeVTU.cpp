/*

The MIT License (MIT)

Copyright (c) 2017-2023 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "isoSurf.hpp"

namespace libp {

void isoSurf_t::WriteVTU(std::ofstream& out, double plot_time) const
{
  LIBP_ABORT("WriteVTU: output stream is not ready.", out.fail());

  // references to trimesh data
  const memory<dfloat>& refN = this->iso_node_data;
  const memory<int>&    refT = this->iso_tri_ids;

  int Nnodes = refN.length() / 4;
  int Ntris  = refT.length() / 3;

  //---------------------------------------------
  // 1. write header
  //---------------------------------------------
  bool BigEndian = false;
  out << "<?xml version=\"1.0\"?>\n";
  if (BigEndian) {
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
  } else {
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  }
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << Nnodes << "\" NumberOfCells=\"" << Ntris << "\">\n";

  //#############################################
  if (0==Ntris && 0==m_rank) {
    // Check Paraview: does rank 0 still need to write an empty file?
    out << "      <Cells>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"></DataArray>\n"
        << "      </Cells>\n"
        << "      <PointData Scalars=\"Color\">\n"
        << "        <DataArray type=\"Float32\" Name=\"Color\"></DataArray>\n"
        << "      </PointData>\n"
        << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";
    out.flush();
    return;
  }
  //#############################################

  int noff = (3 + 1);     // {x,y,z, color} 4 data per node:

  //---------------------------------------------
  // 2. write scalar data
  //---------------------------------------------
  out << "      <PointData Scalars=\"Color\">\n";
  out << "        <DataArray type=\"Float32\" Name=\"Color\" format=\"ascii\">\n";
  for (dlong n = 0; n < Nnodes; ++n) {
    out << "          " << refN[n*noff + 3] << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </PointData>\n";

  //---------------------------------------------
  // 3. write node coordinates
  //---------------------------------------------
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  double pxn, pyn, pzn;
  for (dlong n = 0; n < Nnodes; ++n) {
    int idn = n * noff;         // load coords for node n
    pxn = refN[idn + 0];
    pyn = refN[idn + 1];
    pzn = refN[idn + 2];
    out << "          " << pxn << ' ' << pyn << ' ' << pzn << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";

  //---------------------------------------------
  // 4. write element connectivity/offsets/types
  //---------------------------------------------
  out << "      <Cells>\n";
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  int n1, n2, n3;
  for (dlong e = 0; e < Ntris; ++e) {
    n1 = refT[e * 3 + 0];   // 3 node ids for this tri
    n2 = refT[e * 3 + 1];
    n3 = refT[e * 3 + 2];
    out << "        " << n1 << ' ' << n2 << ' ' << n3 << "\n";
  }
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  dlong cnt = 0;
  for (dlong e = 0; e < Ntris; ++e) { cnt += 3; out << "        " << cnt << "\n"; }
  out << "        </DataArray>\n";
  out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (dlong e = 0; e < Ntris; ++e) { out << "          5\n"; }
  out << "        </DataArray>\n";
  out << "      </Cells>\n";

  //---------------------------------------------
  // 5. write footer
  //---------------------------------------------
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";
  out.flush();
}


void isoSurf_t::WriteVTU_parallel(const std::string& filename, const MPI_Comm comm) const
{
  //-------------------------------------------------------
  // TODO: enable merging of output from multiple ranks 
  // into n_groups of .vtu files
  //-------------------------------------------------------

  MPI_Info info;
  int ierr = MPI_Info_create(&info);
  LIBP_ABORT("MPI_Info_create error.", ierr != MPI_SUCCESS);

  MPI_File fh;
  ierr = MPI_File_open(comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
  LIBP_ABORT("MPI_File_open error.", ierr != MPI_SUCCESS);

  ierr = MPI_File_set_size(fh, 0); // delete the file contents
  LIBP_ABORT("MPI_File_set_size error.", ierr != MPI_SUCCESS);

  // Barrier required: ensure filesize is zero'd 
  // before any cores start to write to the file
  ierr = MPI_Barrier(comm);
  LIBP_ABORT("MPI_Barrier error.", ierr != MPI_SUCCESS);

  ierr = MPI_Info_free(&info);
  LIBP_ABORT("MPI_Info_free error.", ierr != MPI_SUCCESS);

  nnMSG(1, "[proc:%02d] exporting isosurf data to file: %s\n", m_rank, filename.c_str());

  // Define header size for broadcasting to comm group.
  // uint header_size;
  // uint64_t footer_offset;

  // write header
  if (0 == m_rank) {

    nnMSG(1, "[proc:%02d] writing header for group\n", m_rank);

    // std::stringstream ss;
    // WriteVTU_header(ss);
    // header_size = ss.str().size();
  }

  // ierr = MPI_Bcast(&header_size, 1, MPI_UNSIGNED, 0, comm);
  // LIBP_ABORT("MPI_Bcast error.", ierr != MPI_SUCCESS);

  // std::stringstream ss;
  // WriteVTU_main(patches, get_dataset_names(), get_nonscalar_data_ranges(), vtk_flags, ss);
  //
  // Find write offsets.
  // const std::uint64_t size_on_proc = ss.str().size();
  //       std::uint64_t prefix_sum   = 0;
  // ierr = MPI_Exscan(&size_on_proc, &prefix_sum, 1, MPI_UINT64_T, MPI_SUM, comm);
  // LIBP_ABORT("MPI_Exscan error.", ierr != MPI_SUCCESS);
  //
  // Locate specific offset for each rank.
  // const MPI_Offset offset = static_cast<MPI_Offset>(header_size) + prefix_sum;
  //
  //  if (ss.str().size() > mpi_max_int_count) {
  //    LIBP_ABORT("MPI_File_write_at_all(...): exceeded mpi_max_int_count.", true);
  //  }
  //  ierr = MPI_File_write_at_all(fh, offset, ss.str().c_str(), ss.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
  //  LIBP_ABORT("MPI_File_write_at(...) error.", ierr != MPI_SUCCESS);

  if (m_rank == m_nranks - 1) {

    nnMSG(1, "[proc:%02d] writing footer for group\n", m_rank);

    // footer_offset = size_on_proc + offset;
    //
    // std::stringstream ss;
    // WriteVTU_footer(ss);
    // const unsigned int footer_size = ss.str().size();
    // ...

  }

  ierr = MPI_File_sync(fh);
  LIBP_ABORT("MPI_File_sync error.", ierr != MPI_SUCCESS);

  ierr = MPI_File_close(&fh);
  LIBP_ABORT("MPI_File_close error.", ierr != MPI_SUCCESS);
}


void isoSurf_t::WritePVTU(std::ostream& out, const vecS& piece_names, double plot_time)
{
  LIBP_ABORT("Output stream is not ready.", out.fail());

  // Match order of sections in WriteVTU()

  // (1. header)
  out << "<?xml version=\"1.0\"?>\n";
//out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
  out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  
  // simulation time for this timestep
  out << "    <FieldData>\n"
      << "      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"
      << plot_time << "</DataArray>\n"
      << "    </FieldData>\n";

  // (2. scalar data)
  out << "    <PPointData Scalars=\"Color\">\n";
  out << "      <PDataArray type=\"Float32\" Name=\"Color\" format=\"ascii\"/>\n";
  out << "    </PPointData>\n";

  // (3. node coordinates)
  out << "    <PPoints>\n";
  out << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
  out << "    </PPoints>\n";

  // filename for each piece
  for (const auto& piece_name : piece_names) {
    out << "    <Piece Source=\"" << piece_name << "\"/>\n";
  }
  out << "  </PUnstructuredGrid>\n";
  out << "</VTKFile>\n";
  out.flush();

  // check the stream is still ok
  LIBP_ABORT("Output stream is not Ok.", out.fail());
}


// helpers
uint needed_digits(const uint max_number) {
  if (max_number > 0) {
    return static_cast<int>(std::ceil(std::log10(std::fabs(max_number + 0.1))));
  }
  return 1;
}

stdS padded_int_string(const uint value, const uint ndigits) {
  std::string sz = std::to_string(value);       // 123 -> "123"
  if (sz.size() < ndigits) {                    //     -> "00123"
    // pad with leading zeroes
    const std::string padz(ndigits - sz.size(), '0');
    sz.insert(0, padz);
  }
  return sz;
}


//---------------------------------------------------------
stdS isoSurf_t::ExportVTU(
  const stdS& plot_dir, const stdS& base_name,
  double plot_time, uint plot_num,
  uint plot_digits, uint n_groups
)
//---------------------------------------------------------
{
  const uint n_files_written = (n_groups == 0 || n_groups >= m_nranks) 
    ? m_nranks    // each rank writes one file
    : n_groups;   // merge output into n_groups files

  // assert(n_files_written >= 1);
  int digit_range = n_files_written - 1;  // e.g. 32 ranks => {0:31}
  const uint n_digits = needed_digits(std::max(0, digit_range));
  const uint color = m_rank % n_files_written;
  const stdS filename = plot_dir + base_name + "_" + padded_int_string(plot_num, plot_digits)
                                             + "." + padded_int_string(color, n_digits) + ".vtu";

  if (n_groups == 0 || n_groups >= m_nranks) {
    // every rank writes one file
    std::ofstream output(filename);
    LIBP_ABORT("Could not open file: " << filename, !output);
    WriteVTU(output, plot_time);
  } 
  else if (n_groups == 1) {
    // merge data from all ranks into a single file
    MPI_Comm mpi_communicator = comm.comm();
    //comm_t mpi_communicator = comm;
    this->WriteVTU_parallel(filename, mpi_communicator);
  }
  else {

    // TODO: merge output into n_groups files

#if (1)
    MPI_Comm mpi_communicator = comm.comm();
    MPI_Comm comm_group;
    int ierr = MPI_Comm_split(mpi_communicator, color, m_rank, &comm_group);
    LIBP_ABORT("error calling MPI_Comm_split: " << ierr, ierr != MPI_SUCCESS);
    this->WriteVTU_parallel(filename, comm_group);
    ierr = MPI_Comm_free(&comm_group);
    LIBP_ABORT("error calling MPI_Comm_free: " << ierr, ierr != MPI_SUCCESS);
#else
    comm_t comm_group;
    comm_group = comm.Split(color, rank);
    this->WriteVTU_parallel(filename, comm_group);
    comm.Free();
#endif

  }

  // write pvtu record
  // Note: not using plot_dir. Adding a directory name
  // seems to confuse Paraview's file verification routine?
  const stdS pvtu_filename = base_name + "_" + padded_int_string(plot_num, plot_digits) + ".pvtu";

  // Gather number of tris on each rank, then add
  // only non-empty <Pieces> to piece_names:

  // rank 0 gathers number of triangles on each rank
  int localNtris = this->iso_tri_ids.length() / 3;
  memory<hlong> globalNtris(m_nranks, 0);
  comm.Gather(localNtris, globalNtris, 0);

  if (m_rank == 0) {
    std::vector<stdS> piece_names;
    for (uint i = 0; i < n_files_written; ++i) {

      LIBP_ABORT("error: too many VTU pieces!", i >= m_nranks);

      if (globalNtris[i] > 0) {

        nnMSG(1, "[proc:%02d] VTU <Piece> %d has %d triangles\n", m_rank, i, globalNtris[i]);

        const stdS filename = base_name + "_" + padded_int_string(plot_num, plot_digits) 
                                        + "." + padded_int_string(i, n_digits) + ".vtu";
        piece_names.emplace_back(filename);
      }
      else {
        // empty <Piece>
        nnMSG(1, "[proc:%02d] VTU <Piece> %d has %d triangles (EMPTY)\n", m_rank, i, globalNtris[i]);
      }
    }

    std::ofstream pvtu_output(plot_dir + pvtu_filename);
    WritePVTU(pvtu_output, piece_names, plot_time);
  }

  return pvtu_filename;
}

} //namespace libp
