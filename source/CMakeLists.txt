cmake_minimum_required(VERSION 3.5)
project( viabLabLib )
#To build ViabLab as a library, set the variable BUILD_LIB to ON
option(BUILD_LIB "Builds ViabLab as Library" OFF)   
	
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++17")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wall -Wextra -O3 -flat_namespace,-undefined dynamic_lookup -Wno-deprecated -Wno-unused")	# <--- O3 MODE
  set(CMAKE_LINKER_FLAGS "-flto") # Linking time optimisations
  set(CMAKE_EXECUTABLE_ENABLE_EXPORTS 1)
ELSE()
  set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++17")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wall -Wextra -g -Wno-deprecated -Wno-unused")	 
  set(CMAKE_LINKER_FLAGS "-flto") # Linking time optimisations
  set(CMAKE_EXECUTABLE_ENABLE_EXPORTS 1)
ENDIF() 

# Uniformity for any system
set(CMAKE_SHARED_LIBRARY_PREFIX_CXX "")
set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX ".so")

set(SPDLOG_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/spdlog/include")
message("-- SPD LOG   DIR :  ${SPDLOG_INCLUDE_DIR}  ")
include_directories(${SPDLOG_INCLUDE_DIR})
add_subdirectory("${CMAKE_SOURCE_DIR}/spdlog")

message("-- BOOST INCLUDE DIR : ${Boost_INCLUDE_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

include_directories( include ${BOOST_INCLUDE_DIRS}   )

file(GLOB_RECURSE Viab_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS data/*.h ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS include/*.h )
 
file(GLOB_RECURSE Viab_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS src/*.cpp)


# Find all CPP files in the data directory
file(GLOB MODEL_FILES "data/*.cpp")

# Create a shared library for each CPP file
foreach(MODEL_FILE ${MODEL_FILES})
    get_filename_component(LIB_NAME ${MODEL_FILE} NAME_WE)  # Get the file name without extension
    message("Found model : ${MODEL_FILE}")
    add_library(${LIB_NAME} SHARED ${MODEL_FILE})           # Create a shared library
    target_compile_options(${LIB_NAME} PRIVATE -fPIC)   # Add -fPIC flag
endforeach()

if(BUILD_LIB)
    message("Building ViabLab as library  ")
    ADD_LIBRARY(viabLabLib STATIC ${Viab_SOURCES} ${Viab_HEADERS})  # creation de la librairie  
else()
    message("Building ViabLab as executable ")
    ADD_EXECUTABLE(viabLabExe  ${Viab_SOURCES} ${Viab_HEADERS})  # Executable creation

    # Link all shared libraries to the executable
    target_link_libraries(viabLabExe ${CMAKE_DL_LIBS} spdlog::spdlog_header_only)
    # Set the output directory for the shared libraries
    set_target_properties(${LIB_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

endif(BUILD_LIB)

unset(BUILD_LIB CACHE)
