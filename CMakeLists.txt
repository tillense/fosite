# CMake project file for FOSITE
cmake_minimum_required (VERSION 3.14.1)
project (FOSITE NONE)

# detect git
find_package(Git)
if(Git_FOUND)
  # get the latest git version tag and commit hash
  # adds "-dirty" if working tree is not clean
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tag --dirty
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  add_definitions("-DVERSION=\"${GIT_VERSION}\"")
else()
  add_definitions("-DVERSION=\"0.8.2\"") 
endif()

# DIRTY HACK!!! For some reason cmake does not find libz although it is installed.
# Seems to be a known bug, see https://askubuntu.com/questions/1244299/cmake-can-not-find-zlib-on-ubuntu-20-04
set(ZLIB_LIBRARY "/usr/lib/x86_64-linux-gnu/libz.so")
# check whether LLVM is available
find_package(LLVM CONFIG)
if(LLVM_FOUND)
  message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
  message(STATUS "Using LLVMConfig.cmake in: ${LLVM_DIR}")
  set(CMAKE_Fortran_COMPILER "${LLVM_INSTALL_PREFIX}/bin/flang")
endif()

# check and set the build type
set (CMAKE_BUILD_TYPE_NAMES RELEASE DEBUG PROFILING)
set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Should be one of RELEASE, DEBUG, PROFILING; default: RELEASE")
if (NOT CMAKE_BUILD_TYPE IN_LIST CMAKE_BUILD_TYPE_NAMES)
  message(FATAL_ERROR "Unknown build type, should be one of RELEASE, DEBUG, PROFILING")
endif()

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the build directory."
      FORCE)

#------------------------ command line options -------------------------#
option(PARALLEL "PARALLEL" OFF)
option(MPI_USE_SENDRECV "MPI_USE_SENDRECV" ON)
if (PARALLEL MATCHES "ON")
  message(STATUS "Configuring parallel build using MPI ...")
  if(NOT MPI_Fortran_COMPILER)
    unset(MPI_Fortran_COMPILER CACHE)
    find_program(MPI_Fortran_COMPILER NAMES $ENV{FC} mpiifort mpinfort mpifort mpif90 PATHS ${MPI_DIR} ENV MPI_HOME PATH_SUFFIXES bin DOC "MPI Fortran compiler wrapper")
  endif()
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  enable_language (Fortran)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  if (MPI_Fortran_HAVE_F90_MODULE)
    add_definitions(-DHAVE_MPI_MOD)
    message(STATUS "Using MPI Fortran 90 module")
  elseif (MPI_Fortran_HAVE_F77_HEADER)
    add_definitions(-DHAVE_MPIF_H)
    message(STATUS "Using MPI legacy Fortran 77 header")
  else()
    message(FATAL_ERROR "Neither Fortran 90 module nor Fortran 77 header found")
  endif()
  message(WARNING
      "There is a known problem causing segfaults with OpenMPI \
      (experienced with version 3.1.3), depending on the amount of processes \
      used. If you experience any problems with OpenMPI \
      please consider other versions or MPI implementations like, e.g., MPICH.")
  add_definitions(-DPARALLEL)
  set(MPI_LIBRARY_SUFFIX "_mpi") # used, e.g., for ASL and FFTW libraries

  if(MPI_USE_SENDRECV MATCHES "ON")
    add_definitions(-DMPI_USE_SENDRECV)
  endif()
else (PARALLEL MATCHES "ON")
  enable_language (Fortran)
  set(MPI_LIBRARY_SUFFIX "") # no MPI suffix for serial builds
endif()

# workaround for unknown compiler ID NEC
if (CMAKE_Fortran_COMPILER MATCHES "(.*)nfort")
  message(STATUS "Fortran compiler ID set to NEC")
  set(CMAKE_Fortran_COMPILER_ID "NEC")
endif()
enable_testing()

set (X86_VECTOR_EXTENSION NONE CACHE STRING "Enable x86 vector extensions, possible values are: NONE AVX AVX2 AVX-512")

#-------------------- choose compiler & set flags ---------------------#
#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
option(IPO "enable/disable interprocedural / linktime optimizations" OFF)
option(CHECK_LEAKS "enable/disable runtime check for memory leaks (gfortran only)" OFF)
set (CMAKE_Fortran_FLAGS "${FCFLAGS}")
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # GNU compiler, i.e. gfortran
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.4)
    message(FATAL_ERROR "gfortran must be at least version 8.4")
  endif()
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 12.4)
    message(FATAL_ERROR "gfortran must be at most version 12.4")
  endif()
  add_compile_options(-fdefault-real-8 -cpp -ffree-line-length-none)
  if (X86_VECTOR_EXTENSION MATCHES "NONE")
    add_compile_options(-mno-avx)
    message(STATUS "compiling with x86 advanced vector instructions (avx) disabled")
  else ()
    if (X86_VECTOR_EXTENSION STREQUAL "AVX")
      add_compile_options(-mavx)
    elseif (X86_VECTOR_EXTENSION STREQUAL "AVX2")
      add_compile_options(-mavx2)
    elseif (X86_VECTOR_EXTENSION STREQUAL "AVX-512")
      add_compile_options(-mavx512f)
    else ()
      message(FATAL_ERROR "X86_VECTOR_EXTENSION must be one of {NONE,AVX,AVX2,AVX-512}")
    endif()
    message(STATUS "compiling with x86 advanced vector instructions (avx): " ${X86_VECTOR_EXTENSION})
  endif()
  add_compile_options("$<$<CONFIG:RELEASE>:-funroll-all-loops;-fno-f2c;-O3>")
  add_compile_options("$<$<CONFIG:PROFILING>:-funroll-all-loops;-fno-f2c;-O3;-pg;-no-pie;-fopt-info-vec>")
  add_compile_options("$<$<CONFIG:DEBUG>:-fno-f2c;-O2;-g;-fsanitize=address;-Wall;-ffpe-trap=invalid,zero,overflow>")
  add_link_options("$<$<CONFIG:PROFILING>:-pg>")
  add_link_options("$<$<CONFIG:DEBUG>:-fsanitize=address>")
  if (CHECK_LEAKS MATCHES "ON")
     add_link_options(-fsanitize=leak)
  endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # intel Fortran compiler, i.e. ifort
  add_compile_options(-cpp -r8 -stand f08)
  if (X86_VECTOR_EXTENSION MATCHES "NONE")
    add_compile_options(-mno-avx)
    message(STATUS "compiling with x86 advanced vector instructions (avx) disabled")
  else ()
    message(STATUS "using x86 advanced vector instructions: " ${X86_VECTOR_EXTENSION})
    if (X86_VECTOR_EXTENSION STREQUAL "AVX")
      add_compile_options(-xavx)
    elseif (X86_VECTOR_EXTENSION STREQUAL "AVX2")
      add_compile_options(-xcore-avx2)
    elseif (X86_VECTOR_EXTENSION STREQUAL "AVX-512")
      add_compile_options(-xcommon-avx512)
    else ()
      message(FATAL_ERROR "X86_VECTOR_EXTENSION must be one of {NONE,AVX,AVX2,AVX-512}")
    endif()
  endif()
  add_compile_options("$<$<CONFIG:RELEASE>:-O3>")
  add_compile_options("$<$<CONFIG:PROFILING>:-O3;-pg;-qopt-report=1 -qopt-report-phase=vec>")
  add_compile_options("$<$<CONFIG:DEBUG>:-O0;-W1;-g;-CB;-CS;-traceback>")
  add_link_options("$<$<CONFIG:DEBUG>:-g;-CB;-CS;-traceback>")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  add_compile_options(-fdefault-real-8 -Mpreprocess -Mextend -Mfreeform)
  add_compile_options("$<$<CONFIG:RELEASE>:>")
  add_compile_options("$<$<CONFIG:PROFILING>:-fprofile-instr-generate>")
  add_compile_options("$<$<CONFIG:DEBUG>:-g;-Mstandard;-fsanitize=address>")
  if (X86_VECTOR_EXTENSION MATCHES "AVX")
    add_compile_options(-fvectorize)
  endif()
  add_link_options("$<$<CONFIG:DEBUG>:-fsanitize=address>")
  if (CHECK_LEAKS MATCHES "ON")
     add_link_options(-fsanitize=leak)
  endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NEC")
  add_compile_options(-report-all -proginf -fdefault-real=8 -x f95-cpp-input)
  add_compile_options("$<$<CONFIG:RELEASE>:-finline-max-depth=5;-finline-max-function-size=150;-finline-functions;-floop-collapse>")
  add_compile_options("$<$<CONFIG:DEBUG>:-O1;-g;-fcheck=all;-finline-max-depth=5;-finline-functions>")
  add_compile_options("$<$<CONFIG:PROFILING>:-finline-max-depth=5;-finline-max-function-size=150;-finline-functions;-floop-collapse;-ftrace>")
  add_link_options("$<$<CONFIG:RELEASE>:-proginf>")
  add_link_options("$<$<CONFIG:PROFILING>:-ftrace>")
else (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # unknown compiler
  message (STATUS "Found unsupported Fortran compiler: " ${CMAKE_Fortran_COMPILER})
  add_compile_options("$<$<CONFIG:RELEASE>:>")
  add_compile_options("$<$<CONFIG:PROFILING>:>")
  add_compile_options("$<$<CONFIG:DEBUG>:>")
endif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

include(CheckIPOSupported)
if(IPO MATCHES ON)
  # enable link time optimizations
  check_ipo_supported(RESULT ipo_supported OUTPUT output)
  message(STATUS "Checking whether the compiler supports IPO: ${ipo_supported}")
  if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO requested but not supported by the compiler: ${output}")
  endif()
endif()

# set preprocessor macro to the vector length of the target system
# DEFAULT: 1
if (CMAKE_SYSTEM_NAME MATCHES "(.*)SX-Aurora")
  add_definitions(-DVECTOR_LENGTH="256" -DNECSXAURORA)
else ()#(CMAKE_SYSTEM_NAME MATCHES "(.*)SX-AURORA")
  if (X86_VECTOR_EXTENSION STREQUAL "AVX")
    add_definitions(-DVECTOR_LENGTH="16")
  elseif (X86_VECTOR_EXTENSION STREQUAL "AVX2")
    add_definitions(-DVECTOR_LENGTH="16")
  elseif (X86_VECTOR_EXTENSION STREQUAL "AVX-512")
    add_definitions(-DVECTOR_LENGTH="32")
  else()
    add_definitions(-DVECTOR_LENGTH="1")
  endif()
endif ()#(CMAKE_SYSTEM_NAME MATCHES "(.*)SX-AURORA")

set(DEBUGLEVEL ERROR CACHE STRING "Set verbosity for debugging, possible values are: ERROR WARN INFO")
if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  if (DEBUGLEVEL MATCHES "ERROR")
    add_compile_definitions(DEBUG=1)
  elseif (DEBUGLEVEL MATCHES "WARN")
    add_compile_definitions(DEBUG=2)
  elseif (DEBUGLEVEL MATCHES "INFO")
    add_compile_definitions(DEBUG=3)
  else()
    message(FATAL_ERROR "Unknown verbosity level for debugging")
  endif()
endif()

#----------------------- search ASL ------------------------------------------#
find_library(ASL_LIB NAMES "asl${MPI_LIBRARY_SUFFIX}_sequential" DOC "NEC ASL library"
  HINTS ${ASL_DIR} ${FFTW_DIR} ${VE_LIBRARY_PATH} ${NVE_NLC_PATH} ENV NLC_HOME PATH_SUFFIXES lib lib64)
if (NOT ASL_LIB MATCHES ASL_LIB-NOTFOUND)
  message(STATUS "Found NEC ASL library: ${ASL_LIB}")
  find_path(ASL_MODULE_PATH NAMES asl_unified.mod DOC "ASL Fortran module path"
    HINTS ${ASL_DIR} ${FFTW_DIR} ${VE_LIBRARY_PATH} ${NVE_NLC_PATH} ENV NLC_HOME
    PATH_SUFFIXES include/mod${MPI_LIBRARY_SUFFIX})
  if (NOT ASL_MODULE_PATH MATCHES ASL_MODULE_PATH-NOTFOUND)
    message(STATUS "Found NEC ASL Fortran module in: ${ASL_MODULE_PATH}")
    include_directories(${ASL_MODULE_PATH})
  endif()
  set(FFTW_NAME_PREFIX "asl")
  set(ASL_FOUND TRUE)
else()
  set(FFTW_NAME_PREFIX "")
  set(ASL_FOUND FALSE)
endif()

  
#----------------------- search GSL ------------------------------------------#
option(USE_GSL "Add GSL/FGSL support" OFF)
if(USE_GSL MATCHES ON)
  find_package(GSL REQUIRED)

  if((GSL_VERSION VERSION_EQUAL "2.1") OR (GSL_VERSION VERSION_EQUAL "2.2.1"))
	set(FGSL_CHECKOUT_TAG "v1.1.0")
  elseif(GSL_VERSION VERSION_EQUAL "2.3")
	set(FGSL_CHECKOUT_TAG "v1.2.0")
  elseif(GSL_VERSION VERSION_EQUAL "2.4")
	set(FGSL_CHECKOUT_TAG "v1.3.0")
  elseif(GSL_VERSION VERSION_EQUAL "2.5")
	set(FGSL_CHECKOUT_TAG "v1.4.0")
  elseif(GSL_VERSION VERSION_EQUAL "2.6")
	set(FGSL_CHECKOUT_TAG "1.5.0") # missing 'v' is correct here!
  elseif(GSL_VERSION VERSION_GREATER_EQUAL "2.7")
	set(FGSL_CHECKOUT_TAG "master")
  else()
	message(FATAL_ERROR "GSL version not supported")
  endif()
  message(STATUS "Checking out FGSL version: ${FGSL_CHECKOUT_TAG}")

  set(FGSL_SOURCE_DIR ${CMAKE_SOURCE_DIR}/contrib/fgsl)
  set(FGSL_INSTALL_DIR ${CMAKE_BINARY_DIR}/contrib/fgsl)
  set(FGSL_LIBRARY_DIR ${FGSL_INSTALL_DIR}/lib)
  get_filename_component(gsl_library_suffix ${GSL_LIBRARY} LAST_EXT)
  # check whether we use static linking for GSL
  if (${gsl_library_suffix} MATCHES ".a")
	  set(FGSL_LIBRARY ${FGSL_LIBRARY_DIR}/libfgsl.a)
  else()
	  set(FGSL_LIBRARY ${FGSL_LIBRARY_DIR}/libfgsl.so)
  endif()
  include_directories(${GSL_INCLUDE_DIRS} ${FGSL_INSTALL_DIR}/include/fgsl)
  # add gsl fortran interface as external project
  if(NOT EXISTS "${FGSL_SOURCE_DIR}/configure.ac")
    message(STATUS "fgsl: found empty source dir -> fetching code")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --force contrib/fgsl
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  endif()
  execute_process(
        COMMAND ${GIT_EXECUTABLE} show-ref --quiet --verify -- "refs/heads/fosite_build"
	RESULTS_VARIABLE CREATE_FOSITE_BRANCH
	WORKING_DIRECTORY ${FGSL_SOURCE_DIR}
  )
  if(${CREATE_FOSITE_BRANCH})
	execute_process(
		COMMAND ${GIT_EXECUTABLE} checkout tags/${FGSL_CHECKOUT_TAG} -b fosite_build
		WORKING_DIRECTORY ${FGSL_SOURCE_DIR}
		ERROR_VARIABLE GIT_ERROR
		OUTPUT_QUIET)
  else()
	execute_process(
		COMMAND ${GIT_EXECUTABLE} checkout fosite_build
		WORKING_DIRECTORY ${FGSL_SOURCE_DIR}
		ERROR_VARIABLE GIT_ERROR
		OUTPUT_QUIET)
	if(NOT ${GIT_ERROR})
		execute_process(
			COMMAND ${GIT_EXECUTABLE} reset --hard ${FGSL_CHECKOUT_TAG}
			WORKING_DIRECTORY ${FGSL_SOURCE_DIR}
			OUTPUT_QUIET)
	endif()
  endif()
  if(${GIT_ERROR})
	message(FATAL_ERROR "fgsl: git fails to checkout branch fosite_build")
  endif()

  if(NOT EXISTS ${FGSL_SOURCE_DIR}/configure)
    if(NOT EXISTS ${FGSL_SOURCE_DIR}/m4)
      file(MAKE_DIRECTORY ${FGSL_SOURCE_DIR}/m4)
    endif()
    execute_process(
      COMMAND autoreconf -i
      WORKING_DIRECTORY ${FGSL_SOURCE_DIR}
    )
  endif()
  include(ExternalProject)
  ExternalProject_Add(fgsl
    PREFIX ${FGSL_INSTALL_DIR}
    SOURCE_DIR ${FGSL_INSTALL_DIR}/src
    BUILD_IN_SOURCE 0
    BINARY_DIR ${FGSL_INSTALL_DIR}/build
    INSTALL_DIR ${FGSL_INSTALL_DIR}
    DOWNLOAD_DIR ${FGSL_INSTALL_DIR}
    STAMP_DIR ${FGSL_INSTALL_DIR}/stamp
    TMP_DIR ${FGSL_INSTALL_DIR}/tmp
    DOWNLOAD_COMMAND  ${CMAKE_COMMAND} -E remove_directory "${FGSL_INSTALL_DIR}/src" && ${CMAKE_COMMAND} -E create_symlink "${FGSL_SOURCE_DIR}" "${FGSL_INSTALL_DIR}/src"
    CONFIGURE_COMMAND ../src/configure FC=${CMAKE_Fortran_COMPILER} -q --prefix=${FGSL_INSTALL_DIR}
    BUILD_COMMAND make
    INSTALL_COMMAND make install)

   add_compile_definitions(HAVE_FGSL)
else()
  set(GSL_FOUND FALSE)
endif()

#----------------------- search FFTW -----------------------------------------#
option(USE_FFTW "Add FFTW support" OFF)
if(DEFINED FFTW_DIR)
  if(IS_DIRECTORY ${FFTW_DIR})
    set(USE_FFTW ON)
  else()
    message(FATAL_ERROR "${FFTW_DIR} should be a directory!")
  endif()
endif()
if(USE_FFTW MATCHES ON)
  # on NEC SX-Aurora fftw is linked with the asl library
  if (CMAKE_SYSTEM_NAME MATCHES "(.*)SX-Aurora" AND NOT ASL_FOUND)
    message(FATAL_ERROR "FFTW support requested but no ASL library found.")
  endif()

  # find FFTW library
  find_library(FFTW_LIB NAMES ${FFTW_NAME_PREFIX}fftw3 DOC "FFTW library"
    HINTS ${FFTW_DIR} ${VE_LIBRARY_PATH} ${NVE_NLC_PATH} ENV FFTW_DIR
    PATH_SUFFIXES lib lib64
    REQUIRED)
  message(STATUS "Found FFTW library: ${FFTW_LIB}")

  # find FFTW MPI library for parallel build
  if(PARALLEL MATCHES ON)
    set(FFTW_FORTRAN_MODULE_NAME "fftw3-mpi.f03")
    find_library(FFTW_MPI_LIB NAMES ${FFTW_NAME_PREFIX}fftw3${MPI_LIBRARY_SUFFIX} DOC "FFTW MPI library"
      HINTS ${FFTW_DIR} ${VE_LIBRARY_PATH} ${NVE_NLC_PATH} ENV FFTW_DIR
      PATH_SUFFIXES lib lib64
      REQUIRED)
    message(STATUS "Found FFTW MPI library: ${FFTW_MPI_LIB}")
  else()
    set(FFTW_FORTRAN_MODULE_NAME "fftw3.f03")
  endif()

  # find FFTW Fortran module
  find_path(FFTW_MODULE_PATH NAMES ${FFTW_NAME_PREFIX}${FFTW_FORTRAN_MODULE_NAME} DOC "FFTW Fortran module path"
    HINTS ${FFTW_DIR} ${VE_LIBRARY_PATH} ${NVE_NLC_PATH} ENV FFTW_DIR
    PATH_SUFFIXES include
    REQUIRED)
  message(STATUS "Found FFTW Fortran module: ${FFTW_MODULE_PATH}/${FFTW_NAME_PREFIX}${FFTW_FORTRAN_MODULE_NAME}")
  include_directories(${FFTW_MODULE_PATH})
  add_compile_definitions(HAVE_FFTW)

# if (CMAKE_SYSTEM_NAME MATCHES "(.*)SX-Aurora")
#   set (FFTW_INCLUDE ${NVE_NLC_PATH}/include/)
#   if (NOT PARALLEL)
#     set (ASL_INCLUDE ${NVE_NLC_PATH}/include/mod)
#   elseif (PARALLEL)
#     set (ASL_INCLUDE ${NVE_NLC_PATH}/include/mod_mpi)
#   endif()
#   message(STATUS "ASL headers: ${ASL_INCLUDE}")
# endif()

  set(FFTW_FOUND TRUE)
else()
  set(FFTW_FOUND FALSE)
endif()

#---------------------- add fosite subdirs  -----------------------------------#
# build executables
set (LIBDIRS numtools common boundary fluxes io mesh physics sources timedisc )

foreach (p ${LIBDIRS})
  include_directories(/.$(CMAKE_BINARY_DIR)/${p})
endforeach (p)
include_directories(/.$(CMAKE_BINARY_DIR))
include_directories(/.$(CMAKE_SOURCE_DIR)/common) #because of tap.h

foreach (p ${LIBDIRS})
  add_subdirectory (${p})
endforeach (p)

# create target fosite
add_library(fosite fosite.f90)
target_link_libraries(fosite ${LIBDIRS})

# add asl
if(ASL_FOUND)
  target_link_libraries(fosite ${ASL_LIB})
endif()

# add gsl
if(GSL_FOUND)
  target_link_libraries(fosite ${GSL_LIBRARIES} ${FGSL_LIBRARY})
  add_dependencies(numtools fgsl)
endif()

# add fftw
if(FFTW_FOUND)
  if(PARALLEL MATCHES ON)
    target_link_libraries(fosite ${FFTW_MPI_LIB})
  endif()
  target_link_libraries(fosite ${FFTW_LIB})
else()
  message(WARNING "FFTW support DISABLED! \
    The gravity spectral solvers will not be available. \
    Add -DUSE_FFTW=ON to the cmake command line to enable FFTW support.")
endif()

#--------------------- add paths for standard simulation setups----------------#
add_subdirectory(examples)
add_subdirectory(tests)

# add directory with user defined simulation setups
set(SIMULATION_PATH "" CACHE PATH "Add directories with user defined simulation setups")
foreach (dir ${SIMULATION_PATH})
  if (EXISTS "${dir}/CMakeLists.txt")
    get_filename_component(fncomp "${dir}" NAME)
    message("-- found simulation source path ${dir}")
    message("-- adding simulation binary path ${CMAKE_BINARY_DIR}/simulations/${fncomp}")
    add_subdirectory(${dir} ${CMAKE_BINARY_DIR}/simulations/${fncomp})
  endif()
endforeach(dir)

#--------------------- doxygen documentation ----------------------------------#
option(BUILD_DOC "Build documentation" ON)
if (BUILD_DOC MATCHES ON)
  find_package(Doxygen 1.8.12)
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
endif()
