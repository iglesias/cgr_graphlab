project(cgr_graphlab)

cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
set(ROS_BUILD_TYPE Release)

message(STATUS "Compiling CGR sources...")
message(STATUS "Compilers found: ${CMAKE_CXX_COMPILER_LIST}" )
message(STATUS "Using compiler: ${CMAKE_CXX_COMPILER}" )
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Arch: ${CMAKE_SYSTEM_PROCESSOR}")

set(CMAKE_CXX_FLAGS "-mmmx -msse -msse2 -msse4.1 -save-temps -DHAS_MPI -fopenmp")

rosbuild_init()
#uncomment if you have defined messages
rosbuild_genmsg()
#uncomment if you have defined services
rosbuild_gensrv()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

find_package(Lua51 REQUIRED)

include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/shared/math)
include_directories(${PROJECT_SOURCE_DIR}/src/shared/util)
include_directories(${PROJECT_SOURCE_DIR}/src/map)
include_directories(${PROJECT_SOURCE_DIR}/src/cgr)
include_directories(${PROJECT_SOURCE_DIR}/src/fspf)

#FIXME make machine independent
include_directories("/home/iglesias/workspace/boost_1_53_0")
include_directories("/usr/include/mpi")

include_directories(${LUA_INCLUDE_DIR})
include_directories(${USB_INCLUDE_DIR})

# find libraries
find_library(GRAPHLAB_LIBRARY NAMES graphlab libgraphlab)
find_library(MPI_LIBRARY NAMES mpi libmpi)
find_library(MPIXX_LIBRARY NAMES mpi++ libmpi++)
find_library(TCMALLOC_LIBRARY NAMES tcmalloc libtcmalloc)
find_library(ZOOKEEPER_MT_LIBRARY NAMES zookeeper_mt libzookeeper_mt)
find_library(ZOOKEEPER_ST_LIBRRARY NAMES zookeeper_st libzookeeper_st)
find_library(EVENT_LIBRARY NAMES event libevent)
find_library(EVENT_PTHREADS_LIBRARY NAMES event_pthreads libevent_pthreads)
#FIXME make machine independent
find_library(BOOST_CONTEXT_LIBRARY NAMES boost_context libboost_context
  HINTS "/home/iglesias/workspace/boost_1_53_0/stage/lib")

set(libs popt pthread ${LUA_LIBRARIES} rt
  ${MPI_LIBRARY} ${MPIXX_LIBRARY} ${TCMALLOC_LIBRARY}
  ${EVENT_LIBRARY} ${EVENT_PTHREADS_LIBRARY}
  ${ZOOKEEPER_MT_LIBRARY} ${ZOOKEEPER_ST_LIBRRARY}
  ${BOOST_CONTEXT_LIBRARY}
  ${GRAPHLAB_LIBRARY})

set(SHARED_SRCS
  src/shared/math/gvector.h
  src/shared/math/gvector.cpp
  src/shared/math/geomalgo.h
  src/shared/math/geometry.h
  src/shared/util/proghelp.cc
  src/shared/util/terminal_utils.cpp
  src/shared/util/configreader.cpp
  src/shared/util/watch_files.cpp)

set(MAP_SRCS src/map/vector_map.cpp)

set(CGR_SRCS src/cgr/vectorparticlefilter.cpp)

set(FSPF_SRCS
  src/fspf/plane_filtering.cpp
  src/fspf/plane_polygon.cpp
  src/fspf/grahams_scan.cpp)

#create libraries to eliminate redundant compile steps

add_library(shared_libraries ${SHARED_SRCS})
target_link_libraries(shared_libraries shared_libraries ${libs})

add_library(map_libraries ${MAP_SRCS})
add_dependencies(map_libraries shared_libraries )
target_link_libraries(map_libraries shared_libraries ${libs})

add_library(fspf_libraries ${FSPF_SRCS} )
add_dependencies(fspf_libraries shared_libraries )
target_link_libraries(fspf_libraries shared_libraries ${libs})

add_library(cgr_libraries ${CGR_SRCS})
add_dependencies(cgr_libraries shared_libraries map_libraries)
target_link_libraries(cgr_libraries shared_libraries map_libraries ${libs})

#executable target
set(target cgr_localization)
rosbuild_add_executable(${target} src/cgr/localization_main.cpp )
add_dependencies(${target} cgr_libraries fspf_libraries shared_libraries map_libraries rosbuild_precompile)
target_link_libraries(${target} shared_libraries map_libraries cgr_libraries fspf_libraries ${libs} roslib)
