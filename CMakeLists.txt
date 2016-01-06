cmake_minimum_required(VERSION 3.1)
#enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)
list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

if (NOT DEFINED USE_EVERYTRACE)
	set(USE_EVERYTRACE NO)
endif()

if (USE_EVERYTRACE)
	find_package(Everytrace REQUIRED)
	add_definitions(-DUSE_EVERYTRACE)
	include_directories(${EVERYTRACE_INCLUDE_DIR})
endif()

# -----------------------------------------------------
#https://cmake.org/pipermail/cmake/2007-February/012796.html
FIND_PACKAGE(Doxygen)

# ----------------------------------------------------------

include (setup_rpath)

# --------------------- Define our toolchain

find_package(gtest REQUIRED)
find_package(ibmisc REQUIRED)
find_package(Blitz++ REQUIRED)
find_package(netcdf4_c REQUIRED)
find_package(netcdf4_cxx REQUIRED)

include_directories(${BLITZ++_INCLUDE_DIR})

include_directories(${NETCDF4_C_INCLUDE_DIR})
include_directories(${NETCDF4_CXX_INCLUDE_DIR})
include_directories(${IBMISC_INCLUDE_DIR})

set(ALL_LIBS ${EVERYTRACE_C_LIBRARY} ${IBMISC_LIBRARY} ${NETCDF4_CXX_LIBRARY})

# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(bin)
add_subdirectory(test)

