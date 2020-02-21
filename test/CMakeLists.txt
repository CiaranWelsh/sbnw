add_executable(Test NotARealTest.cpp)

target_include_directories(Test PUBLIC
        # graphfab sources
#        ${CMAKE_SOURCE_DIR}/graphfab
        # include googletest sources
        ${CMAKE_SOURCE_DIR}/lib/googletest/googletest/include
#        ${CMAKE_CURRENT_BINARY_DIR}
        )

message(STATUS "CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(Test PUBLIC
        gtest gtest_main
#        sbnw
        )