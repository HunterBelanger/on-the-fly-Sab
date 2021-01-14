cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

#######################################################################
# Declare project dependencies
#######################################################################

FetchContent_Declare(interpolation
    GIT_REPOSITORY https://github.com/njoy/interpolation
    GIT_TAG        2a76934a148bf379ab594f6cdd2cdf4c8c28e447
)

FetchContent_Declare(ENDFtk
    GIT_REPOSITORY https://github.com/njoy/ENDFtk
    GIT_TAG        v0.1.0
)

FetchContent_Declare(GSL
    GIT_REPOSITORY https://github.com/microsoft/GSL
    GIT_TAG        v3.1.0
)

#######################################################################
# Load dependencies
#######################################################################

FetchContent_MakeAvailable(
    interpolation
    ENDFtk
    GSL
)