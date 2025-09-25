cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(G4_PS_converter LANGUAGES CXX)

FIND_PACKAGE(ROOT REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(INFO "Setting build to release")
  set(CMAKE_CXX_FLAGS_RELEASE "-O2")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(INFO "Setting build to debug")
  set(CMAKE_CXX_FLAGS_DEBUG "-Wall -DDEBUG -DDEBUG2 -Wextra -Wpedantic -g -Og -DDEBUG_REL_ABUNDANCE")
endif()

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
FIND_PACKAGE(Geant4 REQUIRED)


include(${Geant4_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/inc)
link_directories(${CMAKE_SOURCE_DIR}/lib)

#add_library(G4_PS_hist STATIC inc/G4_PS_hist.hpp lib/G4_PS_hist.cpp)
#add_library(G4_PS_hist STATIC inc/G4_hist.hpp lib/G4_hist.cpp)

file(GLOB SOURCES "lib/*.cpp")



#target_link_libraries(G4_PS_hist PUBLIC ROOT::TMVA)

add_executable(G4_PS_creator G4_PS_converter.cpp ${SOURCES})
target_link_libraries(G4_PS_creator PUBLIC ROOT::TMVA ${Geant4_LIBRARIES})

