include_directories(${PROJECT_SOURCE_DIR}/slib)


# ----------- spsparse_c
add_library(spsparse SHARED
	${EVERYTRACE_REFADDR}	# Blank if not USE_EVERYTRACE
	spsparse/spsparse.cpp)
install(TARGETS spsparse DESTINATION lib)

# ------------ Header Files
# install SpSparse headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")
