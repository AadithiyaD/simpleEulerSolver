cmake_minimum_required(VERSION 3.10)
project(mySolver CXX)

# Ensure c++17 is used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the include folder
include_directories(include)

# List all source files
set(SOURCES
    main.cpp
    src/postProcessor.cpp
    src/preProcessor.cpp
    src/solver.cpp)

add_executable(mySolver ${SOURCES})
