using MPI
using Pxest.p8est

if !MPI.Initialized()
  MPI.Init()
  const mpicomm = MPI.COMM_WORLD
  const mpirank = MPI.Comm_rank(mpicomm)
end

let
  N = 3 # polynomial degree

  conn = p8est.Connectivity(1,1,2)
  pxest = p8est.PXEST(conn; min_lvl=0)
  p8est.refine!(pxest; maxlevel=3) do which_tree, quadrant
    if rand() > 0.9
      return Cint(1)
    else
      return Cint(0)
    end
  end
  p8est.balance!(pxest)
  p8est.partition!(pxest)
  p8est.ghost!(pxest)
  p8est.lnodes!(pxest, N)
  mesh = p8est.Mesh(pxest)

  # dump VTK
  vtk_dir = "vtk_files"
  vtk_base = "mesh_p8est"
  mpirank == 0 ? mkpath(vtk_dir) : nothing
  MPI.Barrier(mpicomm)
  p8est.vtk_write_file(pxest, string(vtk_dir, "/", vtk_base))
  if mpirank == 0
    mv(string(vtk_dir, "/", vtk_base, ".pvtu"), string(vtk_base, ".pvtu"),
      force=true)
    mv(string(vtk_dir, "/", vtk_base, ".visit"), string(vtk_base, ".visit"),
      force=true)
  end
end

if !isinteractive()
  # Run gc to make sure cleanup happens before MPI is finalized
  GC.gc()
  MPI.Finalize()
end
