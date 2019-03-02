function writemesh(base_name, x, y, z; fields=(), realelems=1:size(x)[end])
  (Nqr, Nqs, Nqt, ~) = size(x)
  (Nr, Ns, Nt) = (Nqr-1, Nqs-1, Nqt-1)
  Nsubcells = Nr * Ns * Nt
  cells = Array{MeshCell{Array{Int,1}}, 1}(undef, Nsubcells * length(realelems))
  ind = LinearIndices((1:Nqr, 1:Nqs, 1:Nqt))
  for e ∈ realelems
    offset = (e-1) * Nqr * Nqs * Nqt
    for k = 1:Nt
      for j = 1:Ns
        for i = 1:Nr
          cells[i + (j-1) * Nr + (k-1) * Nr * Ns + (e-1) * Nsubcells] =
            MeshCell(VTKCellTypes.VTK_VOXEL,
                     offset .+ ind[i:i+1, j:j+1, k:k+1][:])
        end
      end
    end
  end

  vtkfile = vtk_grid("$(base_name)", @view(x[:]), @view(y[:]), @view(z[:]),
                     cells; compress=false)
  for (name, v) ∈ fields
    vtk_point_data(vtkfile, v, name)
  end
  outfiles = vtk_save(vtkfile)
end
