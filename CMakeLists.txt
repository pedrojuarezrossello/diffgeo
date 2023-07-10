cmake_minimum_required(VERSION 3.12)

project("diffgeo" VERSION 1.0.1
             DESCRIPTION "A header-only library that provides a computational framework to engage with differential geometric objects in C++."
             HOMEPAGE_URL "https://github.com/pedrojuarezrossello/diffgeo")
			 
include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
			
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)