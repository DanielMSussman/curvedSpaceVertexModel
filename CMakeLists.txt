cmake_minimum_required(VERSION 3.0.0)
project(sphericalFlocks)

add_definitions(-DDIMENSION=3)

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} -g -frounding-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -frounding-math -O3 -std=c++11")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --expt-relaxed-constexpr
                                        -arch=sm_30
                                        -gencode=arch=compute_30,code=sm_30
                                        -gencode=arch=compute_50,code=sm_50
                                        -gencode=arch=compute_52,code=sm_52
                                        ")
                                        #")#-gencode arch=compute_30,code=sm_30 -Wno-deprecated-gpu-targets ")#-Xptxas -O3 -Xcompiler -O3")

if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -Wno-deprecated-declarations")
    set(Qt5_DIR "/Users/dmsussma/Qt/5.11.2/clang_64/lib/cmake/Qt5")
    set(Qt5OpenGL_DIR "/Users/dmsussma/Qt/5.11.2/clang_64/lib/cmake/Qt5OpenGL")
    find_package(Qt5OpenGL REQUIRED)
    #find_package(Qt5 COMPONENTS Core Gui Widgets 3DCore 3DExtras 3DRender 3DInput OpenGL)
     find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL)
     find_package(OpenGL)
 elseif(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    #set(Qt5_DIR "/home/user/Qt/5.11.2/gcc_64/lib/cmake/Qt5")
    find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL)
    find_package(OpenGL)
endif()

add_definitions( ${QT_DEFINITIONS} )

find_package(CUDA REQUIRED)
find_package(CGAL REQUIRED)


if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -lineinfo -Xptxas --generate-line-info")
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

add_subdirectory(src)

include_directories(inc ${CUDA_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS}/GL src src/utilities src/simulation src/model src/forces src/updaters /usr/include/GL)
link_directories(${CUDA_LIBRARY_DIRS})

set(myLibs
    force
    forceGPU
    updater
    updaterGPU
    simulation
    model
    modelGPU
    utility
    utilityGPU
    )
message("libraries = " "${myLibs}")

qt5_wrap_ui(UI_HEADERS mainwindow.ui)
set(SOURCES mainwindow.cpp oglwidget.cpp)
set(HEADERS
       mainwindow.h
       oglwidget.h
    )

# list the names of cpp files corresponding to linked executables you'd like...
foreach(ARG guiVersion sphericalVertex)
add_executable("${ARG}.out" "${ARG}.cpp" ${SOURCES} ${HEADERS} ${UI_HEADERS})
target_link_libraries("${ARG}.out" ${myLibs} ${CUDA_LIBRARIES} ${OPENGL_LIBRARIES}
    Qt5::Widgets
    Qt5::Core
    Qt5::Gui
    CGAL::CGAL
    )
qt5_use_modules("${ARG}.out" Widgets)
endforeach()
