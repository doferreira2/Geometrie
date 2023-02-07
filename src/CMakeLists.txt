add_executable(genre "genre.cpp")
add_executable(measures "measures.cpp")

target_link_libraries(genre PRIVATE CGAL)
target_link_libraries(measures PRIVATE CGAL)