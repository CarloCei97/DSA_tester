cmake_minimum_required(VERSION 3.27)
project(DS_A_Tester)

set(CMAKE_CXX_STANDARD 23)

find_package(Matplot++)
if (Matplot++_FOUND)
    add_executable(DS_A_Tester main.cpp
            src/DownsamplingAlgorithms.cpp
            include/DownsamplingAlgorithms.h
            src/CSVHandler.cpp
            include/CSVHandler.h)
    target_link_libraries(DS_A_Tester PRIVATE Matplot++::matplot)
endif ()

#target_include_directories(downsampling PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})