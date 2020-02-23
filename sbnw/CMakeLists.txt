set(SOURCES
        network.h network.cpp
        arrowhead.h     arrowhead.cpp
        curve.h
        point.h         point.cpp
        transform.h     transform.cpp
        box.h           box.cpp
        sign_mag.h
        allen.h
        constants.h
        autolayoutSBML.h autolayoutSBML.cpp
        error.h         error.cpp
        layout.h        layout.cpp
        geom.h          geom.cpp
        export.h
        canvas.h        canvas.cpp
#        cubic.h         cubic.cpp
        rand_unif.h
        io.h            io.cpp
        dist.h

        )

# ** Shared Library **
add_library(sbnw SHARED ${SOURCES} )
target_include_directories(sbnw PUBLIC ${CMAKE_SOURCE_DIR} ${LIBSBML_INCLUDE_DIR})
target_link_libraries(sbnw ${LIBSBML_LIBRARY})

# ** Static Library **
add_library(sbnw-static STATIC ${SOURCES})
target_include_directories(sbnw-static PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${LIBSBML_INCLUDE_DIR})
target_link_libraries(sbnw-static ${LIBSBML_STATIC_LIBRARY})
