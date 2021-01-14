cmake_minimum_required(VERSION 3.14)
project(panglos LANGUAGES CXX)

########################################################################
# Dependencies
########################################################################
set( REPOSITORIES "release"
    CACHE STRING
    "Options for where to fetch repositories: develop, release, local"
    )

if( REPOSITORIES STREQUAL "develop" )
    include( cmake/develop_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "release" )
    include( cmake/release_dependencies.cmake )

elseif( REPOSITORIES STREQUAL "local" )
    include( cmake/local_dependencies.cmake )

endif()

########################################################################
# Compile Options
########################################################################
# Create panglos library target
add_library(panglos STATIC src/panglos.cpp)

# Add public interface includes directory
target_include_directories(panglos PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Require C++17 standard
target_compile_features(panglos PRIVATE cxx_std_17)

# Set compile options for library
target_compile_options(panglos PRIVATE -W -Wall -Wextra -Wpedantic -Weffc++)
target_compile_options(panglos PRIVATE $<$<CONFIG:DEBUG>:-g>)
target_compile_options(panglos PRIVATE $<$<CONFIG:RELEASE>:-O2>)
target_compile_options(panglos PRIVATE $<$<BOOL:PNDL_GO_FASTER>:-Ofast -ffast-math>)

# Link dependencies
target_link_libraries(panglos INTERFACE interpolation
                              INTERFACE ENDFtk
                              INTERFACE GSL
)

########################################################################
# Test Application
########################################################################
# TODO