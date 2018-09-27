# Toolchain for SX Aurora Tsubasa cards.
# Use with: cmake -DCMAKE_TOOLCHAIN_FILE=../Toolchain-SX-Aurora.cmake ..
# this one is important
SET(CMAKE_SYSTEM_NAME SX-Aurora)
SET(CMAKE_SYSTEM_PROCESSOR tsubasa)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

set(ve /opt/nec/ve)

# specify the cross compiler
SET(CMAKE_C_COMPILER ${ve}/bin/ncc)
SET(CMAKE_CXX_COMPILER ${ve}/bin/nc++)
SET(CMAKE_Fortran_COMPILER ${ve}/bin/nfort)

# Sollten eigentlich automatisch gefunden/gesetzt werden wenn man die Compiler definiert hat.
SET(CMAKE_AR ${ve}/bin/nar)
SET(CMAKE_NM ${ve}/bin/nnm)
SET(CMAKE_LD ${ve}/bin/nld)
SET(CMAKE_RANLIB ${ve}/bin/nranlib)

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

