# Toolchain for SX Aurora Tsubasa cards.
# Use with: cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain-SX-Aurora.cmake ..
# this one is important
SET(CMAKE_SYSTEM_NAME SX-Aurora)
SET(CMAKE_SYSTEM_PROCESSOR tsubasa)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

# standard path to nec vector engine tools and libraries
SET(NVE_ROOT /opt/nec/ve CACHE PATH "NEC VE tools and library path") 
# standard path to nec compilers and linkers
SET(NVE_TOOL_PATH ${NVE_ROOT}/bin CACHE PATH "NEC VE compiler and linker path")
# default version and path of the NEC numeric library collection (nlc)
IF(DEFINED ENV{NLC_VERSION})
  SET(NVE_NLC_VERSION $ENV{NLC_VERSION})
ELSE()
  SET(NVE_NLC_VERSION 2.1.0)
ENDIF()
IF(DEFINED ENV{NLC_HOME})
  SET(NVE_NLC_PATH $ENV{NLC_HOME} CACHE PATH "NEC NLC library path")
ELSE()
  SET(NVE_NLC_PATH ${NVE_ROOT}/nlc/${NVE_NLC_VERSION} CACHE PATH "NEC NLC library path")
ENDIF()

# specify the cross compilers
#SET(CMAKE_C_COMPILER ${NVE_TOOL_PATH}/ncc CACHE FILEPATH "NEC VE C compiler")
#SET(CMAKE_CXX_COMPILER ${NVE_TOOL_PATH}/nc++ CACHE FILEPATH "NEC VE C++ compiler")
#SET(CMAKE_Fortran_COMPILER ${NVE_TOOL_PATH}/nfort CACHE FILEPATH "NEC VE Fortran compiler")
SET(CMAKE_C_COMPILER ${NVE_TOOL_PATH}/ncc)
SET(CMAKE_CXX_COMPILER ${NVE_TOOL_PATH}/nc++)
SET(CMAKE_Fortran_COMPILER ${NVE_TOOL_PATH}/nfort)
SET(MPI_C_COMPILER $ENV{NMPI_ROOT}/bin/mpincc CACHE FILEPATH "NEC VE MPI C compiler")
SET(MPI_CXX_COMPILER $ENV{NMPI_ROOT}/bin/mpinc++ CACHE FILEPATH "NEC VE MPI C++ compiler")
SET(MPI_Fortran_COMPILER $ENV{NMPI_ROOT}/bin/mpinfort CACHE FILEPATH "NEC VE MPI Fortran compiler")

# Sollten eigentlich automatisch gefunden/gesetzt werden wenn man die Compiler definiert hat.
SET(CMAKE_AR ${NVE_TOOL_PATH}/nar CACHE FILEPATH "NEC VE archiver")
SET(CMAKE_NM ${NVE_TOOL_PATH}/nnm CACHE FILEPATH "NEC VE list archive symbols")
SET(CMAKE_LD ${NVE_TOOL_PATH}/nld CACHE FILEPATH "NEC VE linker")
SET(CMAKE_RANLIB ${NVE_TOOL_PATH}/nranlib CACHE FILEPATH "NEC VE archive index generator")

# Modifiziere das LINK commmando, da <FLAGS> auch die compile flags enthält und dann der c-preprocessor auf die obj files angewendet wird.
# Original:
#SET(CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <FLAGS> <CMAKE_Fortran_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
SET(CMAKE_Fortran_LINK_EXECUTABLE "<CMAKE_Fortran_COMPILER> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")


# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

