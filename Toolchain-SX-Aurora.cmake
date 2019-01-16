# Toolchain for SX Aurora Tsubasa cards.
# Use with: cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain-SX-Aurora.cmake ..
# this one is important
SET(CMAKE_SYSTEM_NAME SX-Aurora)
SET(CMAKE_SYSTEM_PROCESSOR tsubasa)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

set(DEFAULT_NEC_TOOL_PATH /opt/nec/ve)

# specify the cross compilers
SET(CMAKE_C_COMPILER ${DEFAULT_NEC_TOOL_PATH}/bin/ncc)
SET(CMAKE_CXX_COMPILER ${DEFAULT_NEC_TOOL_PATH}/bin/nc++)
SET(CMAKE_Fortran_COMPILER ${DEFAULT_NEC_TOOL_PATH}/bin/nfort)
SET(MPI_C_COMPILER $ENV{NMPI_ROOT}/bin/mpincc)
SET(MPI_CXX_COMPILER $ENV{NMPI_ROOT}/bin/mpinc++)
SET(MPI_Fortran_COMPILER $ENV{NMPI_ROOT}/bin/mpinfort)

# Sollten eigentlich automatisch gefunden/gesetzt werden wenn man die Compiler definiert hat.
SET(CMAKE_AR ${DEFAULT_NEC_TOOL_PATH}/bin/nar)
SET(CMAKE_NM ${DEFAULT_NEC_TOOL_PATH}/bin/nnm)
SET(CMAKE_LD ${DEFAULT_NEC_TOOL_PATH}/bin/nld)
SET(CMAKE_RANLIB ${DEFAULT_NEC_TOOL_PATH}/bin/nranlib)

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

