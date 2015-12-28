# https://cmake.org/cmake/help/v3.0/module/FindGTest.html

# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/slib)

SET(ALL_LIBS ${GTEST_BOTH_LIBRARIES} sparsesparse)

add_executable(test_xiter test_xiter.cpp)
target_link_libraries(test_xiter ${ALL_LIBS})
add_test(AllTests test_xiter)

add_executable(test_sparsesparse test_sparsesparse.cpp)
target_link_libraries(test_sparsesparse ${ALL_LIBS})
add_test(AllTests test_sparsesparse)
