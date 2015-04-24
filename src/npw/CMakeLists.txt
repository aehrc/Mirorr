project(NPW)

option(NPW_USE_CUDA "Use NVidia CUDA in NPW" ON)
option(NPW_USE_OPENGL "Use OpenGL in NPW" OFF)

if(NPW_USE_CUDA)
  add_subdirectory(NPWcuda)
endif(NPW_USE_CUDA)

if(NPW_USE_OPENGL)
  add_subdirectory(NPWopengl)
endif(NPW_USE_OPENGL)


