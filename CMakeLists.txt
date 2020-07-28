cmake_minimum_required( VERSION 2.8 )

project( HATCHet )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread" )

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
