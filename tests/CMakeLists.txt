# https://cmake.org/cmake/help/v3.0/module/FindGTest.html

# https://github.com/snikulov/google-test-examples/blob/master/CMakeLists.txt

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/slib)

SET(ALL_LIBS ${GTEST_BOTH_LIBRARIES} ${EXTERNAL_LIBS} spsparse)

foreach(TEST xiter array netcdf multiply_sparse)
	add_executable(test_${TEST} test_${TEST}.cpp)
	target_link_libraries(test_${TEST} ${ALL_LIBS})
	add_test(AllTests test_${TEST})
endforeach()

