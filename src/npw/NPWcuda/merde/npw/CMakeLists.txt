


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
    #${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWcuda/ext-cudpp_1.1.1/lib/libcutil.a
    #${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWcuda/ext-cudpp_1.1.1/lib/libcudpp_i686.a
    ${milxview_SOURCE_DIR}/lib/common/src/milxRegistration/NPWcuda/ext-cudpp_1.1.1/lib/libcudpp_x86_64.a
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

