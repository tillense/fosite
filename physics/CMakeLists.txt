add_library(physics
  constants_base.f90
  constants_cgs.f90
  constants_generic.f90
  constants_geometrical.f90
  constants_SI.f90
  physics_base.f90
  physics_generic.f90
  physics_eulerisotherm.f90
  physics_euler.f90
)
target_link_libraries(physics mesh)
