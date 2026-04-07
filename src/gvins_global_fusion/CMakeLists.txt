cmake_minimum_required(VERSION 2.8.3)
project(gvins_global_fusion)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  sensor_msgs
  nav_msgs
  geometry_msgs
  visualization_msgs
  gnss_comm
)

find_package(Ceres REQUIRED)

add_subdirectory(./ThirdParty/GeographicLib/)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ./ThirdParty/GeographicLib/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

catkin_package(
  CATKIN_DEPENDS roscpp std_msgs sensor_msgs nav_msgs geometry_msgs visualization_msgs gnss_comm
)

add_executable(gvins_global_fusion_node
  src/globalOptNode.cpp
  src/globalOpt.cpp)

target_link_libraries(gvins_global_fusion_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} libGeographiccc)

install(TARGETS gvins_global_fusion_node RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})
install(DIRECTORY launch DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
install(DIRECTORY models DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
