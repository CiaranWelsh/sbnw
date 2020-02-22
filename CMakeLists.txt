cmake_minimum_required(VERSION 3.15)
project(SBNW)


# version information
set(SBNW_MAJOR_VER 0)
set(SBNW_MINOR_VER 0)
set(SBNW_PATCHLEVEL 1)

set(SBNW_VER_STRING "${SBNW_MAJOR_VER}.${SBNW_MINOR_VER}.${SBNW_PATCHLEVEL}")

# setup googletest

# build options
set(SBNW_LINK_TO_STATIC_LIBSBML ON)
set(BUILD_DOCS OFF)
set(SBNW_SPYDER_DIST OFF CACHE BOOL "Install Spyder plugin files")
set(BUILD_STATIC_LIB ON CACHE BOOL "Build the static version of the library")
set(ENABLE_PYTHON3_BINDINGS OFF CACHE BOOL "Enable Python 3 bindings")
set(LINK_WITH_LIBSBML 1 CACHE BOOL "Link with libSBML? (effectively a requirement)")


# compiler position independent code (PIC)
#   (necessary for building shared libraries https://stackoverflow.com/questions/5311515/gcc-fpic-option)
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(FPIC_FLAG "-fPIC")
else ()
    set(FPIC_FLAG)
endif ()

#Important options & paths
set(SAGITTARIUS_BUILD_ROOT ${CMAKE_CURRENT_BINARY_DIR})

# set debug level dependent on build type
if (${CMAKE_BUILD_TYPE} STREQUAL [Rr]elease)
    set(SBNW_DEBUG_LEVEL 1 CACHE STRING "Verbosity of debug checks")
else ()
    set(SBNW_DEBUG_LEVEL 2 CACHE STRING "Verbosity of debug checks")
endif ()


#User Configuration
if (${WIN32})
    set(SBNW_WINDOWS_BUILD ON)
    if (MINGW)
        set(SBNW_MINGW ON)
        set(SBNW_NTDDI_VERSION "NTDDI_WIN7" CACHE STRING "Latest Windows version to include support for")
        add_definitions(-DNTDDI_VERSION=${SBNW_NTDDI_VERSION})
    else (MINGW)
        set(SBNW_MINGW OFF)
    endif (MINGW)

    #Internal Options
    set(SBNW_DEBUG_LEVEL 2 CACHE STRING "Verbosity of debug checks, must be greater than zero")
    set(SBNW_LINK_TO_STATIC_LIBSBML ON CACHE BOOL "Link to a static version of libSBML where available")

else ()
    #Internal Options
    set(ENABLE_FULL_RPATH 1 CACHE BOOL "Enable full RPath in generated executable code")
    set(SBNW_LINK_TO_STATIC_LIBSBML OFF CACHE BOOL "Link to a static version of libSBML where available")

    #C/C++ compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -Wall -Wno-inline") # -pedantic -Wextra
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif ()

#MagickWand
set(LINK_WITH_MAGICK OFF CACHE BOOL "Link with Image Magick?")
if (LINK_WITH_MAGICK)
    set(MAGICK_INCLUDE "/usr/include/ImageMagick")
    set(MAGICK_LIBDIR "/usr/lib")
    set(MAGICK_LIBS MagickWand)
endif ()


#libSBML
set(LIBSBML_PREFIX "C:/Program Files (x86)/libsbml" CACHE PATH "libSBML root directory (contains lib, include, etc.)")
set(LIBSBML_INCLUDE_DIR ${LIBSBML_PREFIX}/include)
set(LIBSBML_BIN_DIR ${LIBSBML_PREFIX}/bin)
set(LIBSBML_LIB_DIR ${LIBSBML_PREFIX}/lib)
set(LIBSBML_LIBRARY ${LIBSBML_LIB_DIR}/libsbml.lib)
set(LIBSBML_STATIC_LIBRARY ${LIBSBML_LIB_DIR}/libsbml-static.lib)

set(LIBXML2_PREFIX "C:/Program Files (x86)/libxml2" CACHE PATH "Root directory of libxml2 deployment")
set(LIBXML2_LIB_DIR ${LIBXML2_PREFIX}/lib)

add_subdirectory(lib/googletest)
add_subdirectory(sbnw)
add_subdirectory(test)


message(STATUS "
LIBSBML_PREFIX             ${LIBSBML_PREFIX}
LIBSBML_INCLUDE_DIR        ${LIBSBML_INCLUDE_DIR}
LIBSBML_BIN_DIR            ${LIBSBML_BIN_DIR}
LIBSBML_LIB_DIR            ${LIBSBML_LIB_DIR}
LIBSBML_STATIC_LIBRARY     ${LIBSBML_STATIC_LIBRARY}
LIBXML2_PREFIX             ${LIBXML2_PREFIX}
LIBXML2_LIB_DIR            ${LIBXML2_LIB_DIR}
        ")







