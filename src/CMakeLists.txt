set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/CMakeConfig)

# Not yet released:
# option(USE_NPW "Activate support for Non Parametric Window (NPW) matching" ON)
set(USE_NPW OFF)

if (USE_NPW)
  add_subdirectory(npw)
endif(USE_NPW)

add_subdirectory(mirorr)


