cmake_minimum_required (VERSION 3.8)

project("MeshReconstructor")

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

aux_source_directory(src/RPD/. RPD_SRCS)

add_executable(MeshReconstructor
	${RPD_SRCS})

