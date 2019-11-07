add_executable( beams ${CMAKE_CURRENT_SOURCE_DIR}/beams.cpp )
target_link_libraries( beams admmelastic ${RENDER_LIBS} )
