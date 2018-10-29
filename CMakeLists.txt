# CMake project file for FOSITE
cmake_minimum_required (VERSION 2.8.10)
project (FOSITE NONE)

# get the latest git version tag and commit hash
# adds "-dirty" if working tree is not clean
execute_process(
  COMMAND git describe --tag --dirty
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_options("-DVERSION=\"${GIT_VERSION}\"")


# list of known fortran compilers
if( NOT DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER_INIT
          nfort sxf03 gfortran ifort
          ifc af95 af90 efc f95 pathf2003 pathf95 pgf95 pgfortran lf95 xlf95
          fort gfortran-4 g95 f90 pathf90 pgf90 xlf90 epcf90 fort77
          frt pgf77 xlf fl32 af77 g77 f77
          )
elseif(DEFINED ENV{FC})
  set(CMAKE_Fortran_COMPILER $ENV{FC})
endif()
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

#------------------------ comand line options -------------------------#
option(PARALLEL "PARALLEL" OFF)
MESSAGE(STATUS "Running in parallel: " ${PARALLEL})

#-------------------- choose compiler & set flags ---------------------#
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
MESSAGE( STATUS "Compiler: " ${CMAKE_Fortran_COMPILER})

if (Fortran_COMPILER_NAME MATCHES "(.*)gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS "${FCFLAGS} -fdefault-real-8 -cpp -ffree-line-length-none")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_PROFILING "-funroll-all-loops -fno-f2c -O3 -pg -no-pie")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O2 -g -fcheck=all")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "Insufficient gfortran version")
  endif()
elseif (Fortran_COMPILER_NAME MATCHES "(.*)mpifort")
  # mpifort
  message("Building Fosite in parallel mode")
  set (CMAKE_Fortran_FLAGS "${FCFLAGS} -fdefault-real-8 -cpp -ffree-line-length-none -DPARALLEL -DHAVE_MPI_MOD")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O2 -g -fcheck=all -Wunused")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "Insufficient gfortran version")
  endif()
elseif (Fortran_COMPILER_NAME MATCHES "(.*)ifort")
  # ifort
  set (CMAKE_Fortran_FLAGS "${FCFLAGS} -cpp -r8 -stand f03")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  set (CMAKE_Fortran_FLAGS_PROFILING "-O3 -pg")
elseif (Fortran_COMPILER_NAME MATCHES "(.*)nfort")
  set (CMAKE_Fortran_FLAGS "${FCFLAGS} -report-all -proginf -fdefault-real=8 -x f95-cpp-input")
  set (CMAKE_Fortran_FLAGS_RELEASE "-finline-max-depth=5 -finline-functions -floop-collapse")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  set (CMAKE_Fortran_FLAGS_PROFILING "-finline-max-depth=5 -finline-functions -floop-collapse -ftrace")
  set (CMAKE_EXE_LINKER_FLAGS "-proginf")
  set (CMAKE_EXE_LINKER_FLAGS_PROFILING "-ftrace")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  message(FATAL_ERROR "Fortran 77 is oldschool but still not supported.")
#  set (LINK_FLAGS "")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "(.*)gfortran")

# set preprocessor macro to the vector length of the target system
# DEFAULT: 1
# TODO: add other target systems
if (CMAKE_SYSTEM_NAME MATCHES "(.*)SX-AURORA")
  add_compile_options(-DVECTOR_LENGTH="256")
else ()#(CMAKE_SYSTEM_NAME MATCHES "(.*)SX-AURORA")
  add_compile_options(-DVECTOR_LENGTH="1")
endif ()#(CMAKE_SYSTEM_NAME MATCHES "(.*)SX-AURORA")

MARK_AS_ADVANCED(CMAKE_Fortran_FLAGS_PROFILING
    CMAKE_EXE_LINKER_FLAGS_PROFILING)

#message(STATUS "${CMAKE_Fortran_COMPILER_ID}")

#----------------------- search FFTW -----------------------------------------#
set(CMAKE_FIND_ROOT_PATH ${FFTW_DIR})
# search for library

find_library(
  FFTW_LIB
	NAMES fftw3
	DOC "Fastest Fourier Transform of the West."
)
if(PARALLEL)
  find_library(
    FFTWMPI_LIB
	  NAMES fftw3_mpi
    DOC "Distributed parallel implementation of FFTW."
  )
endif()

# search for headers
if(PARALLEL)
  set (fftw_incfile fftw3-mpi.f03)
else()
  set (fftw_incfile fftw3.f03)
endif()

find_path(
  FFTW_INCLUDE
  NAMES ${fftw_incfile}
)

if(NOT FFTW_LIB)
  message(WARNING "No FFTW Library found.")
else()
  if(NOT FFTW_INCLUDE)
    message(WARNING "No FFTW header file found.")
  endif()
  message(STATUS "FFTW library: ${FFTW_LIB}")
  if(PARALLEL)
    message(STATUS "FFTW MPI library: ${FFTWMPI_LIB}")
  endif()
  message(STATUS "FFTW headers: ${FFTW_INCLUDE}")
  add_definitions(-DHAVE_FFTW)
endif()

#---------------------- set include & link -----------------------------------#
# build executables
set (LIBDIRS numtools common boundary fluxes io mesh physics sources timedisc )

foreach (p ${LIBDIRS})
  include_directories(/.$(CMAKE_BINARY_DIR)/${p})
endforeach (p)
include_directories(/.$(CMAKE_BINARY_DIR))
include_directories(/.$(CMAKE_SOURCE_DIR)/common) #because of tap.h
if(FFTW_INCLUDE)
  include_directories(${FFTW_INCLUDE})
endif()

foreach (p ${LIBDIRS})
  add_subdirectory (${p})
endforeach (p)

add_library(fosite fosite.f90)
target_link_libraries(fosite ${LIBDIRS})
if(FFTW_LIB)
  target_link_libraries(fosite ${FFTW_LIB})
  if(PARALLEL)
    target_link_libraries(fosite ${FFTWMPI_LIB} m)
  endif()
endif()

add_subdirectory(examples)
add_subdirectory(tests)

#--------------------- doxygen documentation ----------------------------------#
option(BUILD_DOC "Build documentation" ON)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)

  # a quite clumpsy way to replace the placeholders in the doxygen document in
  # the right way but it was the fastest I found
  # 1. use glob (generates a list with ; as deliminter which doxygen cannot read)
  # 2. generate new custom string with every element in new line with \ at end
  FILE(GLOB htmllist ${CMAKE_CURRENT_SOURCE_DIR}/doc/html/*.html)
  foreach(line IN LISTS htmllist)
    string(APPEND htmlplaceholder "${line} \\ \n")
  endforeach()
  FILE(GLOB jslist ${CMAKE_CURRENT_SOURCE_DIR}/doc/html/js/*.js)
  foreach(line IN LISTS jslist)
    string(APPEND jsplaceholder "${line} \\ \n")
  endforeach()
  FILE(GLOB fontlist ${CMAKE_CURRENT_SOURCE_DIR}/doc/html/fonts/*)
  foreach(line IN LISTS fontlist)
    string(APPEND fontplaceholder "${line} \\ \n")
  endforeach()

  # copy doxygen file and replace all @strings@ within Doxygen.in
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

  # add target "make doc"
  add_custom_target(doc
    # main command "doxygen Doxyfile" to run doxygen
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMENT "Generating API documentation with Doxygen"

    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/doc/adjust_names.sh ${CMAKE_CURRENT_BINARY_DIR}

    VERBATIM)
endif(DOXYGEN_FOUND)

# install executables and scripts
#install (TARGETS ${EXECUTABLES}
#         RUNTIME DESTINATION "bin")
#install (PROGRAMS ${SCRIPTS}
#         DESTINATION "bin")
