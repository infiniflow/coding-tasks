cmake_minimum_required(VERSION 3.2)
project(cpp17)

# Settings
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(CMAKE_CXX_FLAGS "-std=c++17")

# System Checking
find_package(Threads)

add_subdirectory(src bin)

