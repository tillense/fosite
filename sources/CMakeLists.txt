add_library(sources
#  multipole_cylindrical.f90
  sources_base.f90
#  sources_cooling.f90
  sources_generic.f90
  sources_gravity.f90
#  sources_diskthomson.f90
#  sources_rotframe.f90
  sources_c_accel.f90
#  sources_sgs.f90
#  sources_diskcooling.f90
#  sources_wave_damping.f90
#  sources_forcing.f90
#  sources_stellarheating.f90
#  sources_planetheating.f90
#  sources_planetcooling.f90
  sources_shearbox.f90
  sources_viscosity.f90
#  multipole_spherical.f90
#  multipole_generic.f90
  gravity_base.f90
  gravity_pointmass.f90
  gravity_generic.f90
  gravity_sboxspectral.f90
#  gravity_monopol.f90
#  gravity_binary.f90
#  gravity_multigrid.f90
#  gravity_spectral.f90
#  gravity_potential.f90
)
target_link_libraries(sources fluxes)
