IF(WIN32)
  SET(CIMG_ROOT $ENV{CIMG_ROOT})
  IF(NOT DEFINED CIMG_ROOT)
    MESSAGE(FATAL_ERROR "System environment variable, CIMG_ROOT not set. Please set it to continue...")
  ENDIF(NOT DEFINED CIMG_ROOT)

  FIND_PATH( CIMG_INCLUDE_DIR CImg.h $ENV{CIMG_ROOT}
  DOC "Include directory for CImg.h (Nyala Denoiser)"
  )
ELSE(WIN32)
  FIND_PATH( CIMG_INCLUDE_DIR CImg.h
  PATHS $ENV{CIMG_ROOT}/include /usr/include /usr/local/include /usr/local/CImg /usr/local/include/CImg
  DOC "Include directory for CImg.h (Nyala Denoiser)"
  )
ENDIF(WIN32)

IF (CIMG_INCLUDE_DIR)
  SET(CIMG_FOUND 1)
  SET(CIMG_INCLUDES "${CIMG_INCLUDE_DIR}")
ELSE(CIMG_INCLUDE_DIR)
  SET(CIMG_FOUND 0)
  SET(CIMG_INCLUDE_DIR "")
ENDIF (CIMG_INCLUDE_DIR)

MARK_AS_ADVANCED(CIMG_INCLUDE_DIR)
MARK_AS_ADVANCED(CIMG_INCLUDES)
