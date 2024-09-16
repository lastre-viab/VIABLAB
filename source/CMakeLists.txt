cmake_minimum_required(VERSION 2.8)
project( viabLabLib )
#To build ViabLab as a library, set the variable BUILD_LIB to ON
option(BUILD_LIB "Builds ViabLab as Library" OFF)   
	
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++11")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -O3 -fopenmp -Wall  -flat_namespace,-undefined dynamic_lookup")	# <--- O3 MODE
  set(CMAKE_LINKER_FLAGS "")
  find_package(Boost)

ELSE()
  set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++11")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wno-deprecated -Wall")	 
  set(CMAKE_LINKER_FLAGS "")

  set(Boost_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../../../libs/boost_1_75_0")
ENDIF()
#

set(SPDLOG_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/spdlog/include")
message("-- SPD LOG   DIR :  ${SPDLOG_INCLUDE_DIR}  ")
include_directories(${SPDLOG_INCLUDE_DIR})
add_subdirectory("${CMAKE_SOURCE_DIR}/spdlog")
target_link_libraries(spdlog_header_only INTERFACE spdlog_header_only)

include_directories(${Boost_INCLUDE_DIRS})
message("-- BOOST INCLUDE DIR : ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})

include_directories( include ${BOOST_INCLUDE_DIRS}   )

file(GLOB_RECURSE Viab_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS data/*.h ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS include/*.h )
 
file(GLOB_RECURSE Viab_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS src/*.cpp) 

if(BUILD_LIB)
message("Building ViabLab as library  ")
ADD_LIBRARY(viabLabLib STATIC ${Viab_SOURCES} ${Viab_HEADERS})  # creation de la librairie  
else()
message("Building ViabLab as executable ")
ADD_EXECUTABLE(viabLabExe  ${Viab_SOURCES} ${Viab_HEADERS})  # Executable creation
endif(BUILD_LIB)
unset(BUILD_LIB CACHE)


 