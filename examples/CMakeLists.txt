set (EXES
  #  KHI
  #  acousticwave
    agndisk
  #  betadisk1d
  #  bindisk
  #  bowshock
  #  dmr
  #  forcing-turb
  #  gresho
  #  implosion
  #  kepdisk
  #  KHI2
  #  linearwave
  keplerianvortex
  mestel
  mmsn
  orbitingcylinders
  pringle3d
  restart
  sgdisk
  shearingsheet
  planet2d
  #  planeteu
  #  poiseuille
  #  restart
  #  RMI
  #  selfsimdisk1d
  #  shearingbox
  #  statdisk
  #  TCI
    ttauridisk
  #  vortex2d_sphere
  #  wind3d
)
foreach (p ${EXES})
  add_executable (${p} "${p}.f90")
  target_link_libraries(${p} fosite)
endforeach (p)
