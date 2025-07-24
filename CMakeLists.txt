cmake_minimum_required(VERSION 3.24)

add_library(signalsmith-dsp INTERFACE)
target_include_directories(signalsmith-dsp INTERFACE ${CMAKE_CURRENT_LIST_DIR}/include)