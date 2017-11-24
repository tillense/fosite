# CMake project file for FOSITE

cmake_minimum_required (VERSION 2.8.10)
project (FOSITE NONE)
enable_language (Fortran)
enable_testing()
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the build directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

MESSAGE( STATUS "Compiler: " ${CMAKE_Fortran_COMPILER})

if (Fortran_COMPILER_NAME MATCHES "(.*)gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS "-fdefault-real-8 -cpp -ffree-line-length-none -std=f2003 -DVERSION='\"v2\"'")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O2 -g -fcheck=all")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "Insufficient gfortran version")
  endif()
elseif (Fortran_COMPILER_NAME MATCHES "(.*)ifort")
  # ifort
  set (CMAKE_Fortran_FLAGS "-cpp -r8 -stand f03 -DVERSION='\"v2\"'")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "(.*)sxf03")
  set (CMAKE_Fortran_FLAGS "-Ep -default_real_size=8 -pi auto -DVERSION='\"v2\"'")
  set (CMAKE_Fortran_FLAGS_RELEASE "")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  message(FATAL_ERROR "Fortran 77 is oldschool but still not supported.")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "(.*)gfortran")


# build executables
set (LIBDIRS numtools common boundary fluxes io mesh physics sources timedisc)

#set (EXECUTABLES "NORMA.exe" ${NMPROGRAMS})
#set (SCRIPTS "gen_pert.sh" "pert_multi_mode.sh")

foreach (p ${LIBDIRS})
  include_directories(/.$(CMAKE_BINARY_DIR)/${p})
endforeach (p)
include_directories(/.$(CMAKE_BINARY_DIR))
include_directories(/.$(CMAKE_SOURCE_DIR)/common) #because of tap.h

foreach (p ${LIBDIRS})
  add_subdirectory (${p})
endforeach (p)

add_library(fosite fosite.f90)
target_link_libraries(fosite ${LIBDIRS})

#foreach (p ${EXAMPLES})
#  add_executable (${p} "${p}.f90")
#  target_link_libraries(${p} fosite)
#endforeach (p)

#add_subdirectory(examples)
add_subdirectory(tests)

# install executables and scripts
#install (TARGETS ${EXECUTABLES} 
#         RUNTIME DESTINATION "bin")
#install (PROGRAMS ${SCRIPTS}
#         DESTINATION "bin") 
