# ------------ Header Files

install (FILES
	${CMAKE_CURRENT_SOURCE_DIR}/sparsesparse.hpp
	${CMAKE_CURRENT_BINARY_DIR}/xiter.hpp
         DESTINATION include)


# ----------- sparsesparse_c
add_library(sparsesparse SHARED
	sparsesparse.cpp)
install(TARGETS sparsesparse DESTINATION lib)
