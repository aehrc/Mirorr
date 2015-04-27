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
# Finding external packages
#-------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)  # P-Threads
find_package(CUDA REQUIRED)
find_package(SDL REQUIRED)

INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})    # TODO: None of the include seem to satisfy CLion...
INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE}) # TODO:

find_package(cudpp CONFIG)
if (NOT ${cudpp_FOUND})
  message(SEND_ERROR "CUDPP could NOT be found on the system")
  message("In date of May 2015, it is strongly recommanded to install the 2.2 version")
  message("directly from the https://github.com/drhenault/cudpp repository. The install ")
  message("process in the main repository (http://cudpp.github.io/) is currently broken.")
  message("")
  message("Quick install cues:")
  message("  cd ~/myDevDirectory")
  message("  git clone git@github.com:drhenault/cudpp.git cudpp")
  message("  cd cudpp")
  message("  git checkout 2.2-drh-install-fixed")
  message("  git submodule init")
  message("  git submodule update")
  message("  mkdir build; cd build")
  message("  [do the regular CMake thing]")
  message("  make -j 24")
  message("  sudo make install  # IMPORTANT: need to be installed")
  message("")
  message("See also: https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-Linux")
endif ()

#-------------------------------------------------------
# Compiler flags for relevent platform
#-------------------------------------------------------
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
#-------------------------------------------------------
LINK_DIRECTORIES (${LIBRARY_OUTPUT_PATH})

#-------------------------------------------------------
# Define & link targets
#-------------------------------------------------------
INCLUDE_DIRECTORIES(
    ${NPW_CUDA_CORE_INCLUDES}
    ${CUDA_TOOLKIT_INCLUDE}
    ${CMAKE_CURRENT_SOURCE_DIR}/../NPWopengl
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

ADD_LIBRARY( NPWindowsCUDA
  NPWCudaHistogramCreator.cxx
  NPWCudaOpenGLRenderer.cxx
  NPWCudaUtilities.cxx
  NPWCudaWindow.cxx
  )

TARGET_LINK_LIBRARIES( NPWindowsCUDA
                       NPWCuda_CUDA )



