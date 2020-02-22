







add_executable(Test
#        NotARealTest.cpp
        NetworkTests.cpp
        )

target_include_directories(Test PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/lib/googletest/googletest/include
        )

message(STATUS "CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(Test PUBLIC
        gtest gtest_main
        sbnw-static
        )















