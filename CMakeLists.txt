cmake_minimum_required(VERSION 2.8)
project( viabLabLib )
set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++11")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -O3 -fopenmp -Wno-deprecated -Wall")	 
  set(CMAKE_LINKER_FLAGS "")
    
find_package(Boost)
message("-- BOOST INCLUDE DIR : ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})

include_directories( include ${BOOST_INCLUDE_DIRS}   )

file(GLOB_RECURSE Viab_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS data/*.h ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS include/*.h )
 
file(GLOB_RECURSE Viab_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS src/*.cpp) 

ADD_EXECUTABLE(viabLabExe  ${Viab_SOURCES} ${Viab_HEADERS})  # creation de la librairie  
 