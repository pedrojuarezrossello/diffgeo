project(DIFFGEO)
cmake_minimum_required(VERSION 3.12)

# Install
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src DESTINATION include)

# C++17
set(CMAKE_CXX_STANDARD 17)

# DIFFGEO
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

# Set appropriate optimization flags
# For debugging, use `cmake -DCMAKE_BUILD_TYPE=Debug .`
SET(CMAKE_BUILD_TYPE Release)



