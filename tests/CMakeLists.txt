set (EXES
  bondi2d
  bondi3d
  gauss2d
  gauss3d
  KHI3D
  KHI2D
  KHI2D_sphere
  ode
  planet2d
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
  collapse
#  noh2d
#  noh3d
  riemann3d
#  vortex2d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
  if(PARALLEL MATCHES ON)
    add_test(NAME ${p} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${p}> ${MPIEXEC_POSTFLAGS})
    set_property(TEST ${p} PROPERTY PROCESSORS 2)
  else()
    add_test(NAME ${p} COMMAND $<TARGET_FILE:${p}>)
  endif()
endforeach (p)
