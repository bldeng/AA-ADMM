add_executable( beams ${CMAKE_CURRENT_SOURCE_DIR}/beams.cpp )
target_link_libraries( beams admmelastic ${RENDER_LIBS} )

add_executable( plinkopony ${CMAKE_CURRENT_SOURCE_DIR}/plinkopony.cpp )
target_link_libraries( plinkopony admmelastic ${RENDER_LIBS} )

add_executable( plinkohit ${CMAKE_CURRENT_SOURCE_DIR}/plinkohit.cpp )
target_link_libraries( plinkohit admmelastic ${RENDER_LIBS} )

add_executable( windyflag ${CMAKE_CURRENT_SOURCE_DIR}/windyflag.cpp )
target_link_libraries( windyflag admmelastic ${RENDER_LIBS} )
