cmake_minimum_required(VERSION 2.8.3)
project(hybrid_trajectory_logger)
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(catkin REQUIRED COMPONENTS roscpp nav_msgs geometry_msgs)

include_directories(${catkin_INCLUDE_DIRS})

catkin_package(CATKIN_DEPENDS roscpp nav_msgs geometry_msgs)

add_executable(hybrid_trajectory_logger_node src/hybrid_trajectory_logger_node.cpp)
target_link_libraries(hybrid_trajectory_logger_node ${catkin_LIBRARIES})

install(TARGETS hybrid_trajectory_logger_node RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})
install(DIRECTORY launch DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})
