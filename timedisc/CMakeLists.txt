add_library(timedisc 
  timedisc_generic.f90
  timedisc_modeuler.f90
  timedisc_rkfehlberg.f90
  timedisc_base.f90
  timedisc_cashkarp.f90
  timedisc_dormand_prince.f90
  timedisc_ssprk.f90
)
target_link_libraries(timedisc sources)
