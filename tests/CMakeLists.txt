set (EXES
  bondi2d
  bondi3d
  gauss2d
  gauss3d
  KHI3D
  KHI2D
  ode
  poiseuille
  pringle
  riemann1d
  riemann2d
  rngtest
  rootstest
  RTI
  sblintheo
  sblintheo3D
  sedov2d
  sedov3d
  shear
  specialfunctions
  vortex3d
#  collapse
#  noh2d
#  noh3d
  riemann3d
#  vortex2d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
  add_test(NAME ${p} COMMAND ${p})
endforeach (p)
