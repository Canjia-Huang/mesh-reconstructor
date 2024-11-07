cmake_minimum_required (VERSION 3.10)

set(CMAKE_TOOLCHAIN_FILE, "$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake")
message(STATUS "CMAKE_TOOLCHAIN_FILE: ${CMAKE_TOOLCHAIN_FILE}")

project(MeshReconstructor)

set(CMAKE_BUILD_TYPE "Release")
add_definitions(-DRPD_VERBOSE) # for test verbose

add_compile_definitions(OUTPUT_PATH="${PROJECT_SOURCE_DIR}/data/")
# add_compile_options(/openmp:llvm)

# Eigen
find_package(Eigen3 REQUIRED)
IF(Eigen3_FOUND)
    message(STATUS "Eigen3 FOUNDED ${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
ELSE()
    message(FATAL_ERROR "Eigen3 NOT FOUNDED")
ENDIF()

# Boost
find_package(Boost 1.77.0 REQUIRED)
IF(Boost_FOUND)
    message(STATUS "Boost FOUNDED ${Boost_VERSION_STRING}")
    include_directories(${Boost_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "Boost NOT FOUNDED")
ENDIF()

# CGAL
find_package(CGAL 5.6 REQUIRED)
IF(CGAL_FOUND)
    message(STATUS "CGAL FOUNDED ${CGAL_VERSION_STRING}")
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