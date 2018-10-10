set (EXES
  gauss2d
  gauss3d
  sedov2d
  sedov3d
  KHI3
  KHI2D
  vortex3d
  bondi2d
  bondi3d
#  collapse
#  noh2d
#  noh3d
  ode
#  pringle
  shear
  riemann1d-x
  riemann1d-2d-x
  riemann1d-y
  riemann1d-z
  riemann2d
  riemann2d-xz
  poiseuille
  sblintheo
#  riemann3d
#  rngtest
  RTI
#  specialfunctions
#  vortex2d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
  add_test(NAME ${p} COMMAND ${p})
endforeach (p)
