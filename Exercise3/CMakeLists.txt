cmake_minimum_required(VERSION 3.13)

project(exercise_4)

# Options
set(CMAKE_CXX_STANDARD 14)
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../Libs CACHE PATH "这是lib的路径")
set(Eigen3_DIR ${LIBRARY_DIR}/eigen-master/build CACHE PATH "Path to installed Eigen")
set(glog_DIR ${LIBRARY_DIR}/glog-master/build CACHE PATH "Path to installed glog")
set(Ceres_DIR ${LIBRARY_DIR}/ceres_install/lib/cmake/Ceres CACHE PATH "Path to installed Ceres")
#生成 CeresConfig.cmake 文件并不在构建阶段完成，而是在“安装”阶段生成。即便你执行了 cmake --build . 或 make，
#构建的文件不会自动包括 CeresConfig.cmake。为了生成并安装配置文件，必须执行以下命令之一：make install
# 可以用find 来寻找文件路径 当你 make install 以后 如果找不到 xxxconfig。cmake文件时候
set(CMAKE_PREFIX_PATH ${LIBRARY_DIR}/ceres_install/lib/cmake/Ceres )

find_package(Eigen3 REQUIRED)
find_package(glog REQUIRED)
find_package(Ceres REQUIRED)

add_definitions("-D_DISABLE_EXTENDED_ALIGNED_STORAGE")

set(UTILS
    utils/io.h
    utils/points.h
)
#回忆录 这里的注释是有条件的， 想运行那个就把其他的注释掉 下面的link也一样
#此外 exercise3的build文件夹可以删了重新建，此外数据data的路径就别计较了， 运行的时候就放在build下面就行 别找麻烦

#add_executable(gaussian ${UTILS} gaussian.cpp)
add_executable(dragon ${UTILS} dragon.cpp)
#add_executable(surface ${UTILS} surface.cpp)

#target_link_libraries(gaussian ceres glog::glog)
target_link_libraries(dragon ceres)
#target_link_libraries(surface ceres)

if(WIN32)
    target_compile_definitions(gaussian PUBLIC NOMINMAX _USE_MATH_DEFINES)
    target_compile_definitions(dragon PUBLIC NOMINMAX _USE_MATH_DEFINES)
    target_compile_definitions(surface PUBLIC NOMINMAX _USE_MATH_DEFINES)

    set_property(TARGET gaussian PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)
    set_property(TARGET dragon PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)
    set_property(TARGET surface PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)
       
    get_target_property(GLOG_DLL_PATH_DEBUG glog::glog IMPORTED_LOCATION_DEBUG)
    get_target_property(GLOG_DLL_PATH_RELEASE glog::glog IMPORTED_LOCATION_RELEASE)
    
    add_custom_command(TARGET gaussian POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GLOG_DLL_PATH_DEBUG}"
            ${CMAKE_BINARY_DIR})
    
    add_custom_command(TARGET gaussian POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${GLOG_DLL_PATH_RELEASE}"
            ${CMAKE_BINARY_DIR})
endif(WIN32)
