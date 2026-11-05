cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

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
  set(CMAKE_CXX_FLAGS_DEBUG
    "-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion \
     -Wformat=2 -Wnull-dereference -Wdouble-promotion \
     -Wcast-align -Wcast-qual -Woverloaded-virtual \
     -Wnon-virtual-dtor -Wold-style-cast -Wzero-as-null-pointer-constant \
     -Wimplicit-fallthrough -Wduplicated-cond -Wduplicated-branches \
     -Wlogical-op -Wuseless-cast \
     -DDEBUG -DDEBUG2 -DDEBUG_REGEX \
     -g3 -Og \
     -fno-omit-frame-pointer \
     -fsanitize=address,undefined"
)
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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SOURCES "lib/*.cpp")

# Static library with the converter logic; linked by the G4_PS_creator
# executable and by the unit tests in test/.
add_library(g4_ps_lib STATIC ${SOURCES})
target_link_libraries(g4_ps_lib PUBLIC ROOT::TMVA ${Geant4_LIBRARIES})

add_executable(G4_PS_creator G4_PS_converter.cpp)
target_link_libraries(G4_PS_creator PRIVATE g4_ps_lib)

enable_testing()
add_subdirectory(test)

