cmake_minimum_required(VERSION 3.6)
project(dsc3_seg)
set(CMAKE_CXX_STANDARD 14)

set(EXEC dsc3_seg)

include_directories(CGLA)
include_directories(CImg)
include_directories(is_mesh)
include_directories(SOIL)
include_directories(src)
include_directories(DEMO)
include_directories(./)
include_directories(SCGenerator)
include_directories(TetGen)

aux_source_directory(src SRC_DSC)
aux_source_directory(SOIL SRC_SOIL)
aux_source_directory(is_mesh SRC_MESH)
aux_source_directory(DEMO SRC_DEMO)
aux_source_directory(CGLA SRC_CGLA)
aux_source_directory(SCGenerator SRC_MG)
aux_source_directory(TetGen SRC_GEN)


# Open GL
find_package (OpenGL)
find_package (GLEW REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS})

# Lib png
#find_package(PNG REQUIRED libpng16)
#include_directories(${PNG_INCLUDE_DIR})
# search for pkg-config
include (FindPkgConfig)
if (NOT PKG_CONFIG_FOUND)
    message (FATAL_ERROR "pkg-config not found")
endif ()

# check for libpng
pkg_check_modules (LIBPNG libpng16 REQUIRED)
if (NOT LIBPNG_FOUND)
    message(FATAL_ERROR "You don't seem to have libpng16 development libraries installed")
else ()
    include_directories (${LIBPNG_INCLUDE_DIRS})
    link_directories (${LIBPNG_LIBRARY_DIRS})
    link_libraries (${LIBPNG_LIBRARIES})
endif ()

# Link library

add_executable(${EXEC} ${SRC_DSC} ${SRC_SOIL} ${SRC_MESH} ${SRC_DEMO} ${SRC_CGLA} ${SRC_MG} ${SRC_GEN} ${LIBPNG_LINK_FLAGS})
target_link_libraries(${EXEC} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
