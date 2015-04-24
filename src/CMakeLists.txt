set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/CMakeConfig)

option(USE_NPW "Activate support for Non Parametric Window (NPW) matching" ON)

add_subdirectory(mirorr)

if (USE_NPW)
  add_subdirectory(npw)
endif(USE_NPW)

