set(SOURCES
        Network.h Network.cpp
        )

# ** Shared Library **
#add_library(sbnw SHARED ${SOURCES} )
#target_include_directories(sbnw PUBLIC ${CMAKE_SOURCE_DIR} ${LIBSBML_INCLUDE_DIR})
#target_link_libraries(sbnw ${LIBSBML_LIBRARY})

# ** Static Library **
add_library(sbnw-static STATIC ${SOURCES})
target_include_directories(sbnw-static PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${LIBSBML_INCLUDE_DIR})
target_link_libraries(sbnw-static ${LIBSBML_STATIC_LIBRARY})
