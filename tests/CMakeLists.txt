set (EXES
  gauss2d
  gauss3d
  sedov2d
  sedov3d
  KHI3
  vortex3d
#  bondi2d
#  bondi3d
#  collapse
#  noh2d
#  noh3d
#  ode
#  pringle
#  riemann1d
#  riemann2d
#  riemann3d
#  rngtest
#  RTI
#  specialfunctions
#  vortex2d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
  add_test(NAME ${p} COMMAND ${p})
endforeach (p)
