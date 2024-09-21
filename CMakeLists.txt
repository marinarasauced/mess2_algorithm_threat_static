cmake_minimum_required(VERSION 3.8)
project(mess2_algorithm_threat_static)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)

include_directories(include)

add_library(threat SHARED src/threat.cpp)

install(
  TARGETS threat
  EXPORT export_mess2_algorithm_threat_static
  LIBRARY DESTINATION lib
)

install(
  DIRECTORY include/
  DESTINATION include
)

ament_export_include_directories(include)
ament_export_targets(export_mess2_algorithm_threat_static)
ament_package()
