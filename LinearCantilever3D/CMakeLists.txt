project( LinearCantilever3D )


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


# --- libMesh commands
# Check if there is a env. variable for the libMesh
set(LIBMESH_DIR $ENV{LIBMESH_DIR})

# If not, use the default path
if ( "${LIBMESH_DIR}" STREQUAL "" )
   set (LIBMESH_DIR "/usr/local" )
endif()

# Set the path to the libmesh-config binary
set(LIBMESH_CONFIG "${LIBMESH_DIR}/bin/libmesh-config")

# Get the libMesh's flags, includes, libs ...
execute_process(COMMAND ${LIBMESH_CONFIG} --cppflags OUTPUT_VARIABLE LIBMESH_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --cxxflags OUTPUT_VARIABLE LIBMESH_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --include OUTPUT_VARIABLE LIBMESH_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --libs OUTPUT_VARIABLE LIBMESH_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${LIBMESH_CONFIG} --version OUTPUT_VARIABLE LIBMESH_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

# Append them to the CMake flags and linkers
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LIBMESH_LIBS}" ) 

message(STATUS "libMesh version: ${LIBMESH_VERSION}")
message(STATUS "libMesh dir: ${LIBMESH_DIR}")


# Creating entries for target: LinearCantilever3D
# ############################

add_executable( LinearCantilever3D LinearCantilever3D.cxx )

add_to_cached_list( CGAL_LIBRARY_TARGETS )

# Link the executable to CGAL and third-party libraries
target_link_libraries(LinearCantilever3D ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

