project(SagittariusSrc)

set(SOURCES
        core/SagittariusAssert.c
        core/SagittariusCommon.cpp
        core/SagittariusException.cpp
        diag/error.cpp
        draw/tikz.cpp
        io/io.cpp
        interface/layout.cpp
        layout/arrowhead.cpp
        layout/box.cpp
        layout/canvas.cpp
        layout/fr.cpp
        layout/point.cpp
        math/cubic.cpp
        math/geom.cpp
        math/transform.cpp
        network/network.cpp
        sbml/autolayoutSBML.cpp
        util/string.c
        )

set(HEADERS
        core/BaseWORDs.h
        core/config.h
        core/Debug.hpp
        core/DefaultBool.hpp
        core/SagittariusAssert.h
        core/SagittariusCommon.h
        core/SagittariusConfigure.h
        core/SagittariusCore.h
        core/SagittariusException.hpp
        core/SagittariusPlatform.h
        core/SagittariusPlatformWin.h
        core/SagittariusPrefetch.h
        diag/error.h
#        draw/magick.h
        io/io.h
        interface/layout.h
        layout/arrowhead.h
        layout/box.h
        layout/canvas.h
        layout/curve.h
        layout/fr.h
        layout/layoutall.h
        layout/point.h
        math/allen.h
        math/dist.h
        math/geom.h
        math/min_max.h
        math/rand_unif.h
        math/round.h
        math/sig.h
        math/sign_mag.h
        math/transform.h
        network/network.h
        sbml/autolayoutSBML.h
        util/string.h
        )

configure_file(core/config.h core/config.h)

##Image Magick-related sources
if (LINK_WITH_MAGICK)
    set(MAGICK_SOURCES
            draw/magick.cpp)
    set(MAGICK_HEADERS
            draw/magick.h)
else ()
    set(MAGICK_SOURCES) # Clear magick sources
    set(MAGICK_HEADERS)
endif ()

set(SBNW_SOURCES ${SOURCES} ${MAGICK_SOURCES})
set(SBNW_HEADERS ${HEADERS} ${MAGICK_HEADERS})


# Sphinx
configure_file(sphinx/source/conf.py sphinx/source/conf.py)
configure_file(sphinx/source/index.rst sphinx/source/index.rst COPYONLY)
configure_file(sphinx/Makefile sphinx/Makefile COPYONLY)

# ** Shared Library **
add_library(sbnw SHARED ${SBNW_SOURCES} ${HEADERS})
target_include_directories(sbnw PUBLIC ${CMAKE_SOURCE_DIR} ${LIBSBML_INCLUDE_DIR})

set_target_properties(sbnw PROPERTIES COMPILE_DEFINITIONS SBNW_CLIENT_BUILD=0)

#Link in dependencies

##SBML
#     if(SBNW_LINK_TO_STATIC_LIBSBML)
target_link_libraries(sbnw ${LIBSBML_STATIC_LIBRARY}) # libxml2 libbz2 zdll libiconv Ws2_32


# ** Static Library **
if (BUILD_STATIC_LIB)
    add_library(sbnw_static STATIC ${SBNW_SOURCES})

    # Link in dependencies
    # SBML
    target_link_libraries(sbnw_static ${LIBSBML_STATIC_LIBRARY})
endif ()

#Library dist
install(TARGETS sbnw RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

#Doxygen
#Other files besides headers that contain documentation
if (${BUILD_DOCS})
    set(DOXYGEN_OTHER_FILES)
    set(SBNW_DOXY_INPUT "")
    foreach (f ${HEADERS} ${DOXYGEN_OTHER_FILES})
        set(SBNW_DOXY_INPUT "${SBNW_DOXY_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${f}")
    endforeach ()
    configure_file(../doxy/doxy.cfg ../doxy/doxy.cfg)
    add_custom_target(docs COMMAND doxygen ../doxy/doxy.cfg)
ENDIF ()