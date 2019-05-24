cmake_minimum_required(VERSION 3.14)
project(fdm_bar C)
set(CMAKE_C_STANDARD 11)
add_executable(fdm_bar fdm_bar.c)
target_link_libraries(fdm_bar gsl)


