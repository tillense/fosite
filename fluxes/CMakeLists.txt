add_library(fluxes
  fluxes_generic.f90
  fluxes_base.f90
  reconstruction_generic.f90
  reconstruction_base.f90
  reconstruction_constant.f90
  reconstruction_linear.f90
  fluxes_kt.f90
  #  fluxes_hll.f90
  #  fluxes_hllc.f90
  #  fluxes_exact.f90
)
target_link_libraries(fluxes boundary)
