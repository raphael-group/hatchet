cmake_minimum_required( VERSION 2.8 )

project( ALLELIC-IMF )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

# first we can indicate the documentation build as an option and set it to ON by default
#option(BUILD_DOC "Build documentation" OFF)

##if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
#  # check if Doxygen is installed
#  find_package(Doxygen)
#  if (DOXYGEN_FOUND)
#      # set input and output files
#      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
#      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#      # request to configure the file
#      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
#      message("Doxygen build started")

#      # note the option ALL which allows to build the docs together with the application
#      add_custom_target( doc_doxygen ALL
#          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          COMMENT "Generating API documentation with Doxygen"
#          VERBATIM )
#  else (DOXYGEN_FOUND) 
#    message("Doxygen need to be installed to generate the doxygen documentation")
#  endif (DOXYGEN_FOUND)
##endif()

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

include_directories( "src" )

find_package( GUROBI REQUIRED )

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable ( solve
    src/solve.cpp
    src/ilp-min.h
    src/ilp-min.cpp
    src/gurobi-utils.h
    src/utils.h
    src/utils.cpp
    src/argparse.h
    src/argparse.cpp
    src/input_instance.h
    src/input_instance.cpp
    src/coordinate_descent.h
    src/coordinate_descent.cpp
    src/worker.h
    src/worker.cpp
    src/bbc_instance.h
    src/bbc_instance.cpp
)
target_include_directories( solve PUBLIC "src" ${GUROBI_INCLUDE_DIR} )
target_link_libraries( solve ${GUROBI_LIBRARIES} )
