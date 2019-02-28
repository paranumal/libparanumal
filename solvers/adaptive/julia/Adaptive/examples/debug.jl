using MPI
using Pxest.p8est

!MPI.Initialized() && MPI.Init()

let
  mpicomm = MPI.COMM_WORLD
  mpirank = MPI.Comm_rank(mpicomm)

  N = 3 # polynomial degree

  conn = p8est.Connectivity(1,1,2)
  pxest = p8est.PXEST(conn; min_lvl=0)


  refine_level = 1
  p8est.refine!(pxest) do which_tree, quadrant
    qid = ccall(:p8est_quadrant_child_id, Cint,
                (Ref{Pxest.p8est.pxest_quadrant_t},),
                quadrant)
    add = (qid == 0 || qid == 3 || qid == 5 || qid == 6) ? 1 : 0

    refine =  quadrant.level < refine_level + add

    refine ? Cint(1) : Cint(0)
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
