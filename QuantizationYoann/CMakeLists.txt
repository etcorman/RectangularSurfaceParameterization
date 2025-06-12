project(Quantization)

cmake_minimum_required(VERSION 3.5.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU OR CMAKE_CXX_COMPILER_ID MATCHES Clang)
	set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
	add_compile_options(
		-Wall
		-Wextra
		-fdiagnostics-color
	)
	set_source_files_properties(src/refine/cgls.cpp PROPERTIES COMPILE_FLAGS -march=native)
endif()

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

file(GLOB SOURCES
	src/*.cpp
	src/refine/*.cpp
)
include_directories(src)

add_executable(${PROJECT_NAME} ${SOURCES})