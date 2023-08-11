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

// use to locate occa kernel (isoSurface3D.okl)
#define ISO_DIR LIBP_DIR"/libs/plot/"

namespace libp {

isoSurfSettings_t::isoSurfSettings_t(comm_t _comm) {
  newSetting("ISOSURFACE MODE", "ISO_MODE_SURFACE", "Extract 3D surface or planar slices");
  newSetting("ISOSURFACE FILE NAME",   "iso", "Output name for isosurface plots");

  newSetting("ISOSURFACE FIELD ID",    "3", "Field for isosurfing");
  newSetting("ISOSURFACE COLOR ID",    "4", "Field for coloring isosurface");
  newSetting("ISOSURFACE NUM LEVELS",  "1", "Number of isosurface levels");
  newSetting("ISOSURFACE MAX NUMTRIS", "5e6", "Max num tris from isosurface kernel");

  newSetting("ISOSURFACE CONTOUR MIN", "1.0", "Isosurface min value");
  newSetting("ISOSURFACE CONTOUR MAX", "1.0", "Isosurface max value");
  newSetting("ISOSURFACE VOXEL SIZE",  "0.001", "grid voxel size for merging small triangles");
  newSetting("ISOSURFACE VERTEX TOL",  "0.001", "tolerance for merging close vertices");
}

void isoSurfSettings_t::report() {
  if (comm.rank() == 0) {
    std::cout << "isoSurf Settings:\n\n";
    reportSetting("ISOSURFACE MODE");
    reportSetting("ISOSURFACE FILE NAME");

    reportSetting("ISOSURFACE FIELD ID");
    reportSetting("ISOSURFACE COLOR ID");
    reportSetting("ISOSURFACE NUM LEVELS");
    reportSetting("ISOSURFACE MAX NUMTRIS");

    reportSetting("ISOSURFACE CONTOUR MIN");
    reportSetting("ISOSURFACE CONTOUR MAX");
    reportSetting("ISOSURFACE VOXEL SIZE");
    reportSetting("ISOSURFACE VERTEX TOL");
  }
}


void isoSurf_t::Setup(platform_t& _platform, 
                      isoSurfSettings_t& _settings,
                      const mesh_t& mesh, 
                      properties_t& kernelInfo,
                      comm_t _comm)
{
  platform = _platform;
  settings = _settings;
  props = kernelInfo;
  comm = _comm.Dup();

  // initialize the member MyMesh2 pointer
  this->m_isoMesh = nullptr;

  m_rank = mesh.rank;
  m_nranks = mesh.size;

  m_is3D = (3 == mesh.dim) ? true : false;
  if (!m_is3D) {
    nnMSG(1, "[proc:%02d] isoSurface export requires 3D mesh\n", m_rank);
    return;
  }

  // copy user-defined clipBox/plotRegion
  for (int r = 0; r < 2; ++r) {
    for (int c = 0; c < 3; ++c) {
      m_clipBox[r][c] = mesh.m_plotRegion[r][c];
    }
  }

  if (settings.compareSetting("ISOSURFACE MODE", "SURFACE")) {
    this->isoMode = ISO_MODE_SURF;    // general 3D isosurface
  } else {
    this->isoMode = ISO_MODE_PLANE;   // axis-aligned planar slices
  }

  settings.getSetting("ISOSURFACE FIELD ID", this->isoField);
  settings.getSetting("ISOSURFACE COLOR ID", this->isoColorField);
  settings.getSetting("ISOSURFACE NUM LEVELS", this->isoNlevels);

  // use dfloat to handle exponential format ("5e6") for integer value
  dfloat tmpVal=0;
  settings.getSetting("ISOSURFACE MAX NUMTRIS", tmpVal);
  this->isoMaxNumTris = (int) tmpVal;

  settings.getSetting("ISOSURFACE CONTOUR MAX", this->isoMaxVal);
  settings.getSetting("ISOSURFACE CONTOUR MIN", this->isoMinVal);
  settings.getSetting("ISOSURFACE VOXEL SIZE", this->isoVoxelSize);
  settings.getSetting("ISOSURFACE VERTEX TOL", this->isoVertexTol);
  if (isoMinVal > isoMaxVal) std::swap(isoMinVal, isoMaxVal);

  // number of fields exported by vorticity/valib kernel
  plotNvort = 2;  // {qw, |v|} ... qw: "valib criterion"
  plotNfields = mesh.dim + plotNvort;

  props["defines/" "p_plotNvort"] = plotNvort;     //         {qw,|v|}
  props["defines/" "p_plotNfields"] = plotNfields; // {x,y,z}+{qw,|v|}

  int plotNthreads = std::max(mesh.Np, std::max(mesh.plotNp, mesh.plotNelements));
  props["defines/" "p_plotNthreads"] = plotNthreads;
  props["defines/" "p_plotNp"] = mesh.plotNp;
  props["defines/" "p_plotNelements"] = mesh.plotNelements;

  // set isosurface levels
  isoLevels.malloc(isoNlevels);
  if (isoNlevels > 1) {
    dfloat delta = (isoMaxVal - isoMinVal) / (dfloat)(isoNlevels - 1);
    for (int lev = 0; lev < isoNlevels; ++lev) {
      isoLevels[lev] = isoMinVal + (lev * delta);
    }
  }
  else {
    isoLevels[0] = 0.5 * (isoMinVal + isoMaxVal);
  }

  // transfer isosurface levels to device
  o_isoLevels = platform.malloc<dfloat>(isoLevels);

  // transfer num tris in trimesh from device to host
  isoNtris.malloc(1);
  isoNtris[0] = 0;   // must be zero!
  o_isoNtris = platform.malloc<int>(1);

  // allocate buffer to receive trimesh data from device
  //            ( {x,y,z}   {qw,|v|} ) * (3 per tri) * (Ntris)
  isoMaxQData = (mesh.dim + plotNvort) * 3 * isoMaxNumTris;
  isoQ.malloc(isoMaxQData);
  o_isoQ = platform.malloc<dfloat>(isoQ);

  o_plotEToV = platform.malloc<int>(mesh.plotEToV);
  o_plotInterp = platform.malloc<dfloat>(mesh.plotInterp);

  stdS oklFilePrefix = ISO_DIR "okl/";
  stdS oklFileSuffix = ".okl";

  stdS fileName = oklFilePrefix + "isoSurface3D" + oklFileSuffix;
  stdS kernelName = "isoSurface3D";

  isoSurfaceKernel = platform.buildKernel(fileName, kernelName, props);
  nnMSG(1, "[proc:%02d] isoSurfaceKernel() ready. Nlevels = %d\n", m_rank, isoNlevels);
}


void isoSurf_t::doPlot(const mesh_t& mesh,
                       const deviceMemory<dfloat>& o_q,
                       const deviceMemory<dfloat>& o_Vort,
                       int plot_num)
{
  nnMSG(1, "[proc:%02d] calling isoSurf_t::doPlot()\n", m_rank);

  // switch (isoField) -------------------------
  // case 0: break;  // slice plane normal to x
  // case 1: break;  // slice plane normal to y
  // case 2: break;  // slice plane normal to z
  // case 3: break;  // qw / rho / grad(rho)
  // case 4: break;  // magnitude |curl|
  // -------------------------------------------

  // Note: MUST zero isoNtris[0] before calling the kernel
  o_isoNtris.free();
  isoNtris[0] = 0;
  o_isoNtris = platform.malloc<int>(isoNtris);

  // number of elements in clipped plot region on this rank
  int Nclip = mesh.clipK_local;
  int Ntris1 = 0, Ntris2 = 0;
  hlong fld_Nno = 0, fld_Tri = 0;

  // If using a clipped plot region, some ranks 
  // may have no elements in that region.

  if (Nclip > 1) {

    isoSurfaceKernel(
      Nclip,              // #elements in user-defined clip region
      mesh.o_clipElems,   // ids of elements in clip region
      isoField,           // which field to use for isosurfing
      isoColorField,      // which field to use for color
      isoMaxNumTris,      // limit max number of triangles generated
      isoNlevels,         // number of isosurface levels
      o_isoLevels,        // array of isosurface levels
      mesh.o_x,
      mesh.o_y,
      mesh.o_z,
      o_q,                // {rho,u,v,w,p} [nel][5][nno]
      o_Vort,             // {qw,|v|}      [nel][2][nno]
      o_plotInterp,
      o_plotEToV,
      o_isoNtris,   // output: number of generated triangles
      o_isoQ);      // output: {x,y,z,q0,q1}*3*Ntris

    // find number of generated triangles
    o_isoNtris.copyTo(isoNtris);
    isoNtris[0] = std::min(isoNtris[0], isoMaxNumTris);
    Ntris1 = isoNtris[0];
    nnMSG(1, "[proc:%02d] generated %d iso-triangles\n", m_rank, Ntris1);

    // o_isoQ allocation is sized to store "isoMaxNumTris"
    // Only copy data for "Ntris" extracted by isosurfer
    //     len = Ntris1 * (mesh.dim + isoNfields) * 3;
    size_t len = Ntris1 * (plotNfields) * 3;
    o_isoQ.copyTo(isoQ, len);

    //-----------------------------------------
    // compact and smooth trimesh data
    //-----------------------------------------
    cleanSurface(Ntris1, plotNvort);

    // final isosurface trimesh
    fld_Nno = iso_node_data.length()/4;   // 4 data per node {x,y,z,color}
    fld_Tri = iso_tri_ids.length()  /3;   // 3 nodes per triangle
  }
  else {
    nnMSG(1, "[proc:%02d] no iso-triangles from this rank.\n", m_rank);
  }

  int nranks = mesh.size;

  //---------------------------------------------
  // [gmsh] set ELEMENT offset for this rank
  //---------------------------------------------
  memory<hlong> triOffset(nranks + 1, 0);
  comm.Allgather(fld_Tri, triOffset + 1);
  for (int rr = 0; rr < nranks; ++rr) {
    triOffset[rr + 1] = triOffset[rr] + triOffset[rr + 1];
  }
  int E_offs = triOffset[m_rank];

  //---------------------------------------------
  // [gmsh] set NODE offset for this rank
  //---------------------------------------------
  memory<hlong> nodeOffset(nranks + 1, 0);
  if (nranks > 1) {
    comm.Allgather(fld_Nno, nodeOffset + 1);
    for (int rr = 0; rr < nranks; ++rr) {
      nodeOffset[rr + 1] = nodeOffset[rr] + nodeOffset[rr + 1];
    }
  }
  int N_offs = nodeOffset[m_rank];

  int nno_1 = Ntris1 * 3;
  int nno_2 = fld_Nno;
  nnMSG(1, "\n[proc:%02d] isoSurf smoother ------------ "
           "\n before: %7d tris (%7d verts)"
           "\n  after: %7d tris (%7d verts)\n\n", 
    m_rank, Ntris1, nno_1, fld_Tri, nno_2);

  // output field files
  stdS base_name, fname;
  settings.getSetting("ISOSURFACE FILE NAME", base_name);

#ifdef _MSC_VER
  stdS plot_dir("D:/TW/bin/plot/");
#else
  stdS plot_dir("./");    // FIXME: output directory
#endif

  bool bBinary = false;
  if (0) {                // GMSH
    fname = plot_dir + nnSTR("%s_f%06d_%03d.msh", base_name.c_str(), plot_num, m_rank);
    writeGmsh(mesh, fname, g_tstep, N_offs, E_offs, plot_num, g_time, bBinary);
  }
  else {                  // PARAVIEW
    int plot_digits = 5;  // zero-padded frame numbers
    int n_groups = 0;     // each rank writes its piece of the isosurface
  //int n_groups = 1;     // TODO: merge output from all ranks into 1 group
  //int n_groups = n;     // TODO: merge output from #ranks into n groups
    ExportVTU(plot_dir, base_name, g_time, plot_num, plot_digits, n_groups);
  }
}

} //namespace libp
