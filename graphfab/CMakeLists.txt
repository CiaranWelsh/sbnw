project(SagittariusSrc)

set(SOURCES
        core/SagittariusAssert.c
        core/SagittariusCommon.cpp
        core/SagittariusException.cpp
        error.cpp
        draw/tikz.cpp
        io.cpp
        layout.cpp
        arrowhead.cpp
        box.cpp
        canvas.cpp
        fr.cpp
        ../sbnw/point.cpp
        math/cubic.cpp
        math/geom.cpp
        math/transform.cpp
        network.cpp
        autolayoutSBML.cpp
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
        error.h
        draw/magick.h
        io.h
        layout.h
        arrowhead.h
        box.h
        canvas.h
        curve.h
        fr.h
        layoutall.h
        ../sbnw/point.h
        math/allen.h
        math/dist.h
        math/geom.h
        math/min_max.h
        math/rand_unif.h
        math/round.h
        math/sig.h
        math/sign_mag.h
        math/transform.h
        network.h
        autolayoutSBML.h
        util/string.h
        )

set(SOURCES
        #        core/BaseWORDs.h
        #        core/config.h
        #        core/Debug.hpp
        #        core/DefaultBool.hpp
        #        core/SagittariusAssert.h
        #        core/SagittariusCommon.h
        #        core/SagittariusConfigure.h
        #        core/SagittariusCore.h
        #        core/SagittariusException.hpp
        #        core/SagittariusPlatform.h
        #        core/SagittariusPlatformWin.h
        #        core/SagittariusPrefetch.h
        #        draw/magick.h
        ../sbnw/NotARealClass.cpp
#        error.h
#        io.h
#        layout.h
#        arrowhead.h
#        box.h
#        canvas.h
#        curve.h
#        fr.h
#        layoutall.h
#        point.h
#        math/allen.h
#        math/dist.h
#        math/geom.h
#        math/min_max.h
#        math/rand_unif.h
#        math/round.h
#        math/sig.h
#        math/sign_mag.h
#        math/transform.h
#        network.h
#        autolayoutSBML.h
        )

configure_file(core/config.h core/config.h)


# Sphinx
configure_file(sphinx/source/conf.py sphinx/source/conf.py)
configure_file(sphinx/source/index.rst sphinx/source/index.rst COPYONLY)
configure_file(sphinx/Makefile sphinx/Makefile COPYONLY)


# ** Shared Library **
add_library(sbnw SHARED ${SOURCES} ../test/NotARealTest.cpp)
#Link in dependencies
message(STATUS "LIBSBML_INCLUDE_DIR ${LIBSBML_INCLUDE_DIR}")
message(STATUS "LIBSBML_STATIC_LIBRARY ${LIBSBML_STATIC_LIBRARY}")

target_include_directories(sbnw PUBLIC ${CMAKE_SOURCE_DIR} ${LIBSBML_INCLUDE_DIR})

target_link_libraries(sbnw ${LIBSBML_STATIC_LIBRARY})


# ** Static Library **
if (BUILD_STATIC_LIB)
    add_library(sbnw_static STATIC ${SOURCES})
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