add_library(boundary 
#  boundary_axis.f90
  boundary_nogradients.f90
#  boundary_farfield.f90
#  boundary_custom.f90
#  boundary_noh.f90
#  boundary_fixed.f90
#  boundary_noslip.f90
#  boundary_extrapolation.f90
#  boundary_reflecting.f90
  boundary_generic.f90
#  boundary_absorbing.f90
#  boundary_folded.f90
  boundary_periodic.f90
  boundary_base.f90
#  boundary_dmr.f90
#  boundary_shearing.f90
)
target_link_libraries(boundary physics)
