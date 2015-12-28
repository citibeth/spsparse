cmake_minimum_required(VERSION 3.1)
#enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include (setup_rpath)

# --------------------- Define our toolchain

find_package(gtest REQUIRED)

# -------- Process subdirectories of the build
add_subdirectory(slib)
#add_subdirectory(bin)
add_subdirectory(test)
