CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
PROJECT(VIR_BodyTracking)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_BUILD_TYPE Release)
set(PYTHON_INCLUDE_DIRS ~/anaconda3/envs/face/include/python3.8)
set(PYTHON_LIBRARIES ~/anaconda3/envs/face/lib/libpython3.8.so)


option(LINK_SHARED_ZED "Link with the ZED SDK shared executable" ON)
option(ENABLE_GUI "Enable OpenGL/ OpenCV display" ON)

if (NOT LINK_SHARED_ZED AND MSVC)
    message(FATAL_ERROR "LINK_SHARED_ZED OFF : ZED SDK static libraries not available on Windows")
endif()

find_package(ZED 4 REQUIRED)
find_package(CUDA REQUIRED)
find_package(Python3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(k4a QUIET)
if (ENABLE_GUI)
    add_definitions(-DENABLE_GUI)
    find_package(GLUT REQUIRED)
    find_package(GLEW REQUIRED)
    SET(OpenGL_GL_PREFERENCE GLVND)
    find_package(OpenGL REQUIRED)
	find_package(OpenCV REQUIRED)
endif()

# libigl
option(LIBIGL_WITH_TETGEN            "Use TetGen"         ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"          ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
include(libigl)
set(LIBIGL_DIR ${LIBIGL_INCLUDE_DIR}/..)

# openMP
FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${ZED_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${LIBIGL_DIR})
include_directories(${LIBIGL_DIR}/external)
include_directories(${LIBIGL_DIR}/external/tetgen)

link_directories(${ZED_LIBRARY_DIR})
link_directories(${CUDA_LIBRARY_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS})

IF(NOT WIN32)
    SET(SPECIAL_OS_LIBS "pthread")
ENDIF()

if (ENABLE_GUI)
    include_directories(${GLEW_INCLUDE_DIRS})
    include_directories(${GLUT_INCLUDE_DIR})
	include_directories(${OpenCV_INCLUDE_DIRS})
    
    link_directories(${GLEW_LIBRARY_DIRS})
    link_directories(${GLUT_LIBRARY_DIRS})
    link_directories(${OpenGL_LIBRARY_DIRS})
	link_directories(${OpenCV_LIBRARY_DIRS})

    SET(GUI_LIB ${OPENGL_LIBRARIES} ${OpenCV_LIBRARIES} ${GLUT_glut_LIBRARY} ${GLEW_LIBRARIES})
endif()

FILE(GLOB_RECURSE SRC_FILES src/*.c*)
FILE(GLOB_RECURSE HDR_FILES include/*.h*)

add_executable(${PROJECT_NAME} main.cpp ${HDR_FILES} ${SRC_FILES})

if (LINK_SHARED_ZED)
    SET(ZED_LIBS ${ZED_LIBRARIES} ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY})
else()
    SET(ZED_LIBS ${ZED_STATIC_LIBRARIES} ${CUDA_CUDA_LIBRARY} ${CUDA_LIBRARY})
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
# Link
if (k4a_FOUND)
    message("k4a found")
    target_compile_definitions(${PROJECT_NAME} PRIVATE K4A_FOUND)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} 
        ${ZED_LIBS} 
        ${UTILS_LIB} 
        ${SPECIAL_OS_LIBS} 
        ${OpenCV_LIBRARIES} 
        ${PYTHON_LIBRARIES} 
        igl::core
        igl::tetgen
        igl::opengl_glfw
        ${GUI_LIB}
        k4a::k4a)
else()
    message("k4a not found")
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} 
        ${ZED_LIBS} 
        ${UTILS_LIB} 
        ${SPECIAL_OS_LIBS} 
        ${OpenCV_LIBRARIES} 
        ${PYTHON_LIBRARIES} 
        igl::core
        igl::tetgen
        igl::opengl_glfw
        ${GUI_LIB}
        )
endif()