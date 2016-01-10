include_directories(${PROJECT_SOURCE_DIR}/slib)


# ----------- spsparse_c
set(SPSPARSE_SOURCE spsparse/spsparse.cpp)

add_library(spsparse SHARED
	${EVERYTRACE_REFADDR}	# Blank if not USE_EVERYTRACE
	${SPSPARSE_SOURCE})

target_link_libraries(spsparse ${EXTERNAL_LIBS})
install(TARGETS spsparse DESTINATION lib)

# ------------ Header Files
# install SpSparse headers
file(GLOB SPSPARSE_HEADERS "spsparse/*.hpp")

# http://www.vtk.org/Wiki/CMake_FAQ#How_to_convert_a_semicolon_separated_list_to_a_whitespace_separated_string.3F
foreach(arg ${SPSPARSE_HEADERS})
  set(SPSPARSE_HEADERS_SPACES "${SPSPARSE_HEADERS_SPACES} ${arg}")
endforeach(arg ${SPSPARSE_HEADERS})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/slib/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")


# ------------- API Docs
IF (DOXYGEN_FOUND)

	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
		${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

	add_custom_target(apidoc
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)


#ADD_CUSTOM_TARGET(apidoc_forced
#  COMMAND ${CMAKE_COMMAND} -E echo_append "Building API Documentation..."
#  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
#  COMMAND ${CMAKE_COMMAND} -E echo "Done."
#  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

ENDIF (DOXYGEN_FOUND)
