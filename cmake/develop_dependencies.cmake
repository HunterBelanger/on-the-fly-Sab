cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

#######################################################################
# Declare project dependencies
#######################################################################

FetchContent_Declare(interpolation
    GIT_REPOSITORY https://github.com/HunterBelanger/interpolation
    GIT_TAG        9710c077477315617cceceda6f3995bc90def3df
)

FetchContent_Declare(ENDFtk
    GIT_REPOSITORY https://github.com/njoy/ENDFtk
    GIT_TAG        origin/develop
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