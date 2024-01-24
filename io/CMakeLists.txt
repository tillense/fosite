add_library(io 
  fileio_generic.f90
  #fileio_npy.f90
  fileio_binary.f90
  fileio_gnuplot.f90
  #fileio_hdf5.f90
  fileio_vtk.f90
  #fileio_netcdf.f90
  fileio_base.f90
  #fileio_xdmf.f90
)
target_link_libraries(io timedisc)

