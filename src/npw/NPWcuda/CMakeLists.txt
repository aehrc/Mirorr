#-------------------------------------------------------
# This is the root npwcuda CMakeLists file.
# Version 2.6 is the minimum.
#-------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#-------------------------------------------------------
# Project Name
#-------------------------------------------------------
PROJECT(NPW_CUDA)


#-------------------------------------------------------
# Version string should not include patch level.  The major.minor is
# enough to distinguish available features of the toolkit.
# There will be a #define made for each of these static variables
#-------------------------------------------------------
SET(NPW_CUDA_VERSION_MAJOR "1")
SET(NPW_CUDA_VERSION_MINOR "5")
SET(NPW_CUDA_VERSION_PATCH "0")
SET(NPW_CUDA_VERSION_STRING "${NPW_CUDA_VERSION_MAJOR}.${NPW_CUDA_VERSION_MINOR}")


#-------------------------------------------------------
# Set basic build locations.
# Pre-defined variables are NPW_CUDA_SOURCE_DIR and NPW_CUDA_BINARY_DIR
#-------------------------------------------------------
#SET (LIBRARY_OUTPUT_PATH ${NPW_CUDA_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
#SET (PLUGIN_OUTPUT_PATH ${NPW_CUDA_BINARY_DIR}/plugins CACHE PATH "Single output directory for building all plugins.")
#SET (DEPENDENCY_PATH ${NPW_CUDA_BINARY_DIR}/depend CACHE PATH "Temporary directory to hold downloaded dependencies.")
#SET (EXECUTABLE_OUTPUT_PATH ${NPW_CUDA_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
#SET (TEST_EXECUTABLE_OUTPUT_PATH ${NPW_CUDA_BINARY_DIR}/bin/test CACHE PATH "Single output directory for building all test programs.")
#
#SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${NPW_CUDA_BINARY_DIR}/bin/test/*)
#SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${NPW_CUDA_BINARY_DIR}/DartTestfile.txt)


#-------------------------------------------------------
# CMake root location to that FIND_PACKAGE uses to find packages
#
# Tell CMake's FIND_PACKAGE command where to look for our Find<pkg>
# modules.
#-------------------------------------------------------
IF(UNIX)
    SET(CMAKE_MODULE_PATH "${NPW_CUDA_SOURCE_DIR}/utilities/CMake;/usr/local/share/CMake/Modules;${CMAKE_MODULE_PATH}")
    SET(CMAKE_ROOT $ENV{CMAKE_ROOT})
    IF(NOT CMAKE_ROOT)
       SET(CMAKE_ROOT "/usr/share/CMake")
    ENDIF(NOT CMAKE_ROOT)
ENDIF(UNIX)
IF(WIN32)
  SET(CMAKE_MODULE_PATH "${NPW_CUDA_SOURCE_DIR}/utilities/CMake;${CMAKE_MODULE_PATH}")
ENDIF(WIN32)

#-------------------------------------------------------
# Enable/Disable Shared Libraries
#-------------------------------------------------------
OPTION(BUILD_SHARED_LIBS "Build NPWCudaFramework with shared libraries." ON)
SET(NPW_CUDA_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

#-------------------------------------------------------
# Enable/Disable Build Examples
#-------------------------------------------------------

#-------------------------------------------------------
# Limit the applications/libraries that are built for each OS here
#-------------------------------------------------------
#OPTION(DOWNLOAD_DEP "Download dependent libraries if not present" ON)

#OPTION(USE_CUDA "Use NVidia CUDA" ON)
#OPTION(USE_DEVIL "Use Devil imaging library" ON)
#OPTION(USE_GZSTREAM "Use GZStream" ON)
#OPTION(USE_ZLIB "Use Zlib" ON)
#OPTION(USE_RPATH "Use an embedded rpath for binaries" OFF)

#IF (NOT USE_RPATH)
#    SET(CMAKE_SKIP_BUILD_RPATH TRUE)
#ELSE(NOT USE_RPATH)
#    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
#ENDIF(NOT USE_RPATH)


#-------------------------------------------------------
# Finding external packages
#-------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)  # P-Threads
find_package(CUDA REQUIRED)
find_package(SDL REQUIRED)
find_package(cudpp REQUIRED)

#-------------------------------------------------------
# Compiler flags for relevent platform
IF(UNIX)
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG")
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os -DNDEBUG")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -W -Wall -Wno-unused-function -Wno-unused-parameter -Wconversion -Wno-long-long -fno-strict-aliasing" CACHE STRING "Debug CXX flags" FORCE)
    SET(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} -g --coverage")
    SET(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_DEBUG} -g -pg")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2 -DNDEBUG")
    SET(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -Os -DNDEBUG")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O2 -g")
    SET(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_DEBUG} -g -pg")
ENDIF(UNIX)

#-------------------------------------------------------
# Configure NPWCuda
LINK_DIRECTORIES (${LIBRARY_OUTPUT_PATH})

#-------------------------------------------------------
# Configure ZLIB.
#IF(USE_ZLIB)
#     FIND_PACKAGE(ZLIB)
#ENDIF(USE_ZLIB)

#-------------------------------------------------------
# Dispatch the build into the proper subdirectories.
#-------------------------------------------------------
#SUBDIRS(npw)

#-------------------------------------------------------
# Create the list of include directories needed for MILX header files.
#-------------------------------------------------------
SET(NPW_CUDA_INCLUDE_DIRS_BUILD_TREE ${NPW_CUDA_BINARY_DIR}/src/ ${NPW_CUDA_SOURCE_DIR}/include)

#-------------------------------------------------------
# Setup the Include directories for the cmake configured headers
#-------------------------------------------------------
#INCLUDE_DIRECTORIES(
#  ${NPW_CUDA_INCLUDE_DIRS_BUILD_TREE}
#)

#install( 
#  CODE "execute_process(
#	  COMMAND ant 
#          WORKING_DIRECTORY ${NPW_CUDA_SOURCE_DIR}/ant)"
#)




# TODO
#
#   - find cuda using appropriate make file
#
#   - include cudpp properly
#/home/bro86j/Dev/milx-view/lib/src/milxRegistration/NPWcuda

INCLUDE_DIRECTORIES(
    ${NPW_CUDA_CORE_INCLUDES}
    ${CUDA_TOOLKIT_INCLUDE}
    ${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWopengl
    ${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWcuda/ext-cudpp_1.1.1/cudpp/include
    ${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWcuda/ext-cudpp_1.1.1/common/inc
)

SET( SRCS_CUDA
     NPWCudaExtractGeometryKernel.cu
     NPWCudaRenderGeometryDirectlyFromCUDAKernel.cu
)

SET(NPWCudaCoreLibraries
    ${NPWCudaCoreLibraries}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CUDA_LIBRARIES}
    cudpp
)

SET(CUDA_NVCC_FLAGS	--compiler-options;-fpermissive;-arch=sm_20)

CUDA_ADD_LIBRARY(NPWCuda_CUDA ${SRCS_CUDA})
SET_TARGET_PROPERTIES(NPWCuda_CUDA PROPERTIES SOVERSION 1.0)
TARGET_LINK_LIBRARIES(
        NPWCuda_CUDA 
        ${NPWCudaCoreLibraries}
        GLEW
        SDL 
)

ADD_LIBRARY( dbNPWindows
  NPWCudaHistogramCreator.cxx
  NPWCudaOpenGLRenderer.cxx
  NPWCudaUtilities.cxx
  NPWCudaWindow.cxx
  )

TARGET_LINK_LIBRARIES( dbNPWindows
                       NPWCuda_CUDA )



