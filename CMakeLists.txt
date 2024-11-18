cmake_minimum_required (VERSION 3.10)

# use vcpkg to build
# set(VCPKG_ROOT "$ENV{VCPKG_ROOT}")
# message(STATUS "VCPKG_ROOT: ${VCPKG_ROOT}")
# set(CMAKE_TOOLCHAIN_FILE "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
# message(STATUS "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")

project(MeshReconstructor)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE "Release")
add_definitions(-DDEBUG_VERBOSE) # for test verbose
add_definitions(-DDEBUG_OUTPUT) # for test output
add_compile_definitions(OUTPUT_PATH="${PROJECT_SOURCE_DIR}/data/") # set output file path
# add_compile_options(/openmp:llvm)

# Eigen
find_package(Eigen3 REQUIRED)
IF(Eigen3_FOUND)
    message(STATUS "Eigen3 FOUNDED ${EIGEN3_VERSION_STRING} ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
ELSE()
    message(FATAL_ERROR "Eigen3 NOT FOUNDED")
ENDIF()

# Boost
find_package(Boost REQUIRED)
IF(Boost_FOUND)
    message(STATUS "Boost FOUNDED ${Boost_VERSION_STRING} ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "Boost NOT FOUNDED")
ENDIF()

# CGAL
find_package(CGAL 5.6 REQUIRED)
IF(CGAL_FOUND)
    message(STATUS "CGAL FOUNDED ${CGAL_VERSION_STRING} ${CGAL_INCLUDE_DIRS}")
    include_directories(${CGAL_INCLUDE_DIRS}) 
ELSE()
    message(FATAL_ERROR "CGAL NOT FOUNDED")
ENDIF()

# others
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/external)

set(SRC_DIRS ${PROJECT_SOURCE_DIR}/src)
aux_source_directory(${SRC_DIRS} SRC_SOURCE)
add_executable(MeshReconstructor ${SRC_SOURCE})
