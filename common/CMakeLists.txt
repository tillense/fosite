add_library(common
  logging_base.f90
  array.f90
  common_dict.f90
  field_base.f90
  tap.f90
  tap.h
)
target_link_libraries(common numtools)
