add_executable(ghost_test test.c)
target_link_libraries(ghost_test PUBLIC ghost)
target_include_directories(ghost_test PUBLIC ghost)