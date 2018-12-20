set (EXES
#  bondi2d
#  bondi3d
#  collapse
#  gauss2d
    gauss3d
#  KHI3
#  noh2d
#  noh3d
#  ode
#  pringle
#  riemann1d
#  riemann2d
#  riemann3d
#  rngtest
#  RTI
#  sedov2d
#  sedov3d
#  specialfunctions
#  vortex2d
   vortex2d3d
#  vortex3d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
  add_test(NAME ${p} COMMAND ${p})
endforeach (p)
