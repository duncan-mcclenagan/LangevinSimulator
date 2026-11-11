cmake_minimum_required(VERSION 3.20)

project(NoneqControl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(NoneqControl
  BumpDistribution.cpp
  EntropyCalculator.cpp
  NoneqControl.cpp
  PhaseCloud.cpp
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/data/bumps-1x1-0.csv"
  "${CMAKE_CURRENT_BINARY_DIR}/data/bumps-1x1-0.csv"
  COPYONLY
)
