

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  Thyra_Amesos2Types.hpp
  Thyra_Amesos2LinearOpWithSolve.hpp
  Thyra_Amesos2LinearOpWithSolveFactory.hpp
  )

SET(SOURCES
  Thyra_Amesos2Types.cpp
  Thyra_Amesos2LinearOpWithSolve.cpp
  Thyra_Amesos2LinearOpWithSolveFactory.cpp
  )

TRIBITS_ADD_LIBRARY(
  stratimikosamesos2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
