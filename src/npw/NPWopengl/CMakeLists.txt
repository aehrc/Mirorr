#-------------------------------------------------------
# This is the root npwcuda CMakeLists file.
# Version 2.6 is the minimum.
#-------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#-------------------------------------------------------
# Project Name
#-------------------------------------------------------
PROJECT(NPW_OPENGL)

#-------------------------------------------------------
# Enable/Disable Shared Libraries
#-------------------------------------------------------
OPTION(BUILD_SHARED_LIBS "Build NPWindowOpenGL with shared libraries." ON)

#-------------------------------------------------------
# Finding external packages
#-------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(SDL REQUIRED)

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
# Define & link targets
#-------------------------------------------------------
# INCLUDE_DIRECTORIES(${milxview_SOURCE_DIR}/lib/include/milxRegistration/NPWmontecarlo)

ADD_LIBRARY( NPWindowsOpenGL
      NPWGLImage.cxx
      NPWGeometry.cxx
      NPWJointHistogramCreator.cxx
      NPWOpenGLRenderer.cxx
      NPWWindow.cxx
      NPWVertices.cxx
      )

TARGET_LINK_LIBRARIES( NPWindowsOpenGL
                       GLEW
                       SDL
                       )



