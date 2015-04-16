cmake_minimum_required(VERSION 2.8)
project( Mirorr )

option( USE_NPW "Compile in np windows stuff developed by Daan Broekhuizen" OFF )
option( USE_OPENCL "Compile in GPU-based block matching developed by J. Coatelen, M. Golebiewski and D. R. Henault" OFF )
option( VTK_MIRORR "Debug Mirorr using VTK visualisation" OFF )
option( USE_MIRORR_NON_SYSTEM_BOOST "Use non-system path boost installation. " OFF )

add_definitions( "-Wall -ansi" )

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
include_directories( ${ITK_INCLUDE_DIRS} )

if(${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR} VERSION_GREATER 3.20)
  set(ITK_IO_LIBRARIES ITKIOBioRad ITKIOBMP ITKIOCSV ITKIOGDCM ITKIOGE ITKIOGIPL ITKIOHDF5 ITKIOImageBase ITKIOIPL ITKIOJPEG ITKIOLSM ITKIOMesh ITKIOMeta ITKIONIFTI ITKIONRRD ITKIOPNG ITKIOSiemens ITKIOSpatialObjects ITKIOStimulate ITKIOTIFF ITKIOTransformBase ITKIOTransformHDF5 ITKIOTransformInsightLegacy ITKIOTransformMatlab ITKIOVTK ITKIOXML ITKReview)
else(${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR} VERSION_GREATER 3.20)
  set(ITK_IO_LIBRARIES ITKIO)
endif(${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR} VERSION_GREATER 3.20)

if( VTK_MIRORR )
  add_definitions("-DVTK_MIRORR")
  find_package( VTK REQUIRED)
  include(${VTK_USE_FILE})
  set(VTK_MIRORR_SRC vtkMirorrUtils.cxx)
endif( VTK_MIRORR )

if( USE_NPW )
  set(MILX_VIEW_LIB_DIR "/home/riv019/dev/milx-view/build/lib/" CACHE PATH "Path milx-view compiled library directory).")
  link_directories(${MILX_VIEW_LIB_DIR})
endif( USE_NPW )

if( USE_OPENCL )
  include_directories(/usr/local/cuda/include/)
endif( USE_OPENCL )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") #Turn off those anoying warning... (most of them won't solve: they're part of ITK V4)
add_executable( mirorr mirorr.cxx ${VTK_MIRORR_SRC}) 
target_link_libraries( mirorr
  ${ITK_LIBRARIES}
  boost_program_options boost_filesystem boost_system boost_timer
)

if (USE_MIRORR_NON_SYSTEM_BOOST)
  set(Boost_USE_STATIC_RUNTIME    OFF)
  set(Boost_NO_SYSTEM_PATHS  TRUE)
  set(Boost_NO_BOOST_CMAKE   TRUE)
  set(Boost_REALPATH TRUE)
  find_package(Boost 1.49.0 REQUIRED COMPONENTS program_options filesystem system)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries( mirorr ${Boost_LIBRARIES} )
    set(CMAKE_EXE_LINKER_FLAGS " -Wl,-rpath-link,${Boost_LIBRARY_DIRS} -Wl,-rpath,${Boost_LIBRARY_DIRS}")
  else(Boost_FOUND)
    message("Boost was not found. Expect linking errors.")
  endif(Boost_FOUND)
endif(USE_MIRORR_NON_SYSTEM_BOOST)

if( USE_OPENCL )
  find_package(OpenCL REQUIRED)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENCL") #USE_NPW_CUDA
  target_link_libraries(mirorr ${OPENCL_LIBRARIES})
endif( USE_OPENCL )

#Compile in np windows stuff developed by Daan Broekhuizen
if( USE_NPW )
  message("NPW is not yet available. Coming soon.")
  #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_NPW")
  #TARGET_LINK_LIBRARIES(mirorr dbNPWindows GLU)
endif( USE_NPW )

if( VTK_MIRORR )
  target_link_libraries(mirorr ${VTK_LIBRARIES})
endif( VTK_MIRORR )

