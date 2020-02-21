cmake_minimum_required(VERSION 3.15)
project(SBNW)


# version information
set(SBNW_MAJOR_VER  0)
set(SBNW_MINOR_VER  0)
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
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(FPIC_FLAG "-fPIC")
else()
  set(FPIC_FLAG)
endif()

#Important options & paths
set(SAGITTARIUS_BUILD_ROOT ${CMAKE_CURRENT_BINARY_DIR})

# set debug level dependent on build type
if (${CMAKE_BUILD_TYPE} STREQUAL [Rr]elease)
    set(SBNW_DEBUG_LEVEL 1 CACHE STRING "Verbosity of debug checks")
else()
    set(SBNW_DEBUG_LEVEL 2 CACHE STRING "Verbosity of debug checks")
endif()




#User Configuration
if(${WIN32})
    set(SBNW_WINDOWS_BUILD ON)
    if(MINGW)
        set(SBNW_MINGW ON)
        set(SBNW_NTDDI_VERSION "NTDDI_WIN7" CACHE STRING "Latest Windows version to include support for")
        add_definitions( -DNTDDI_VERSION=${SBNW_NTDDI_VERSION} )
    else(MINGW)
        set(SBNW_MINGW OFF)
    endif(MINGW)

    #Internal Options
    set(SBNW_DEBUG_LEVEL 2 CACHE STRING "Verbosity of debug checks, must be greater than zero")
    set(SBNW_LINK_TO_STATIC_LIBSBML ON CACHE BOOL "Link to a static version of libSBML where available")

else()
    #Internal Options
    set(ENABLE_FULL_RPATH 1 CACHE BOOL "Enable full RPath in generated executable code")
    set(SBNW_LINK_TO_STATIC_LIBSBML OFF CACHE BOOL "Link to a static version of libSBML where available")

    #C/C++ compiler flags
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -Wall -Wno-inline") # -pedantic -Wextra
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

#MagickWand
set(LINK_WITH_MAGICK OFF CACHE BOOL "Link with Image Magick?")
if(LINK_WITH_MAGICK)
    set(MAGICK_INCLUDE "/usr/include/ImageMagick")
    set(MAGICK_LIBDIR "/usr/lib")
    set(MAGICK_LIBS MagickWand)
endif()


# GTest
#find_package(gtest)
#if(GTEST_FOUND)
#  set(WITH_GTEST ON CACHE BOOL "Use gtest framework")
#else()
#  set(WITH_GTEST OFF CACHE BOOL "Use gtest framework")
#endif()

#set(ENABLE_PYTHON FALSE CACHE BOOL "Enable Python bindings")
#if(ENABLE_PYTHON)
#  find_package(PythonInterp)
#  if(PYTHONINTERP_FOUND)
#    set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_STRING})
#  else()
#    message(ERROR FATAL "Unable to find Python")
#  endif()
#endif()

# Find Python libs
#find_package(PythonLibs)
#message(STATUS "Using Python include directory ${PYTHON_INCLUDE_DIR}")

#Dependencies

#SBML
#find_package(libsbml 5.18.0 REQUIRED
#        NAMES libsbml libsbml-static
#        PATHS
#            "C:/Program Files (x86)/libsbml"
#            "C:/Program Files (x86)/libsbml/lib"
#            "C:/Program Files (x86)/libsbml/lib/cmake"
#        )


#libSBML
set(LIBSBML_PREFIX "C:/Program Files (x86)/libsbml" CACHE PATH "libSBML root directory (contains lib, include, etc.)")
set(LIBSBML_INCLUDE_DIR ${LIBSBML_PREFIX}/include)
set(LIBSBML_BIN_DIR ${LIBSBML_PREFIX}/BIN)
set(LIBSBML_LIB_DIR ${LIBSBML_PREFIX}/lib)
set(LIBSBML_STATIC_LIBRARY ${LIBSBML_LIB_DIR}/libsbml-static.lib)

set(LIBXML2_PREFIX "C:/Program Files (x86)/libxml2" CACHE PATH "Root directory of libxml2 deployment")
set(LIBXML2_LIB_DIR ${LIBXML2_PREFIX}/lib)

#link_directories(${LIBSBML_LIB_DIR})
#link_directories(${LIBXML2_LIBDIR})
#
#set(LIBSBML_LIBS libsbml libxml2)



#find_library(LIBSBML_LIBRARY sbml libsbml HINTS ${LIBSBML_PREFIX}/lib)
#find_library(LIBSBML_STATIC_LIBRARY sbml-static REQUIRED libsbml-static HINTS ${LIBSBML_PREFIX}/lib)
#find_path(LIBSBML_INCLUDE_DIR sbml/SBMLTypes.h HINTS ${LIBSBML_PREFIX}/include)
#message(STATUS "LIBSBML_INCLUDE_DIR ${LIBSBML_INCLUDE_DIR}")
#include_directories(${LIBSBML_INCLUDE_DIR})
#
##Image Magick
#if(LINK_WITH_MAGICK)
#    set(SBNW_USE_MAGICK 1)
#    include_directories(${MAGICK_INCLUDE})
#else()
#    set(SBNW_USE_MAGICK 0)
#endif()
#
#if(ENABLE_FULL_RPATH) #http://www.cmake.org/Wiki/CMake_RPATH_handling
#    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
#
#    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#
#    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#
#    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#
#    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#    IF("${isSystemDir}" STREQUAL "-1")
#        SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#    ENDIF("${isSystemDir}" STREQUAL "-1")
#endif()
#
##add include directories
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${CMAKE_CURRENT_BINARY_DIR}) #we have to do this to get the generated header
#

add_subdirectory(lib/googletest)
#add_subdirectory(graphfab)
add_subdirectory(test)
#
##add_subdirectory(sandbox)
##
##add_subdirectory(icons)
##
##add_subdirectory(testcases)
##
##if(WITH_GTEST)
##  add_subdirectory(test)
##endif()
