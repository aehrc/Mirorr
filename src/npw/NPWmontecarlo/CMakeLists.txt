#-------------------------------------------------------
# This is the root npwcuda CMakeLists file.
# Version 2.6 is the minimum.
#-------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#-------------------------------------------------------
# Project Name
#-------------------------------------------------------
PROJECT(NPW_MONTECARLO)

#-------------------------------------------------------
# Finding external packages
#-------------------------------------------------------
#find_package(OpenGL REQUIRED)
#find_package(SDL REQUIRED)
find_package(CImg REQUIRED)

#-------------------------------------------------------
# Define & link targets
#-------------------------------------------------------
# INCLUDE_DIRECTORIES(${milxview_SOURCE_DIR}/lib/include/milxRegistration/NPWmontecarlo)

ADD_LIBRARY( NPWindowsMontecarlo
      db3DNPWmontecarlo.cxx
      utils.cxx
      )

TARGET_LINK_LIBRARIES( NPWindowsMontecarlo
                       )



