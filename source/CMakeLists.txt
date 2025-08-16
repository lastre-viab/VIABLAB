cmake_minimum_required(VERSION 3.5)
project( viabLab )
#To build ViabLab as a library, set the variable BUILD_LIB to ON
option(BUILD_LIB "Builds ViabLab as static Library" OFF)   
	
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin) # The target system is MAC OS
  set(CMAKE_C_COMPILER   gcc)
  set(CMAKE_CXX_COMPILER c++)
   
  set(OPTION "-std=c++17")
  message(".. OPTIONS: ${OPTION}")
  message("-- [-O3] MODE.")				# <--- O3 MODE
  set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wall -Wextra -O3 -flat_namespace,-undefined dynamic_lookup -Wno-deprecated -Wno-unused")	# <--- O3 MODE
  set(CMAKE_EXE_LINKER_FLAGS "-flto") # Linking time optimisations
  set(CMAKE_EXECUTABLE_ENABLE_EXPORTS TRUE)
 
  
ELSE()
  if(WIN32) #The target system is Windows, including Win64
    set(CMAKE_C_COMPILER   gcc)
    set(CMAKE_CXX_COMPILER c++)
    set(OPTION "-std=c++17")
    message(".. OPTIONS: ${OPTION}")
    message("-- [-O3] MODE.")				# <--- O3 MODE
    set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wall -Wextra -g -Wno-deprecated -Wno-unused")	 
    set(CMAKE_EXE_LINKER_FLAGS "-flto -Wl,--no-as-needed -ldl") # Linking time optimisations
    set(CMAKE_SHARED_LIBRARY_ENABLE_EXPORTS TRUE)
    set(CMAKE_EXECUTABLE_ENABLE_EXPORTS TRUE)
  ELSE()
    set(CMAKE_C_COMPILER   gcc)
    set(CMAKE_CXX_COMPILER c++)
    set(OPTION "-std=c++17")
    message(".. OPTIONS: ${OPTION}")
    message("-- [-O3] MODE.")				# <--- O3 MODE
    set(CMAKE_CXX_FLAGS "${OPTION} -fopenmp -Wall -Wextra -g -Wno-deprecated -Wno-unused")	 
    set(CMAKE_EXE_LINKER_FLAGS "-flto -Wl,--no-as-needed") # Linking time optimisations
    set(CMAKE_SHARED_LIBRARY_ENABLE_EXPORTS TRUE)
    set(CMAKE_EXECUTABLE_ENABLE_EXPORTS TRUE)
  ENDIF()
ENDIF() 
# Uniformity for any system
set(CMAKE_SHARED_LIBRARY_PREFIX_CXX "")
set(CMAKE_SHARED_LIBRARY_SUFFIX_CXX ".so")

find_package(spdlog REQUIRED)
find_package(Boost REQUIRED)

file(GLOB_RECURSE Viab_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS data/*.h ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS include/*.h )
 
file(GLOB_RECURSE Viab_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} FOLLOW_SYMLINKS src/*.cpp)


# Find all CPP files in the data directory
file(GLOB MODEL_FILES "data/*.cpp")


if(BUILD_LIB)
    message("Building ViabLab as library  ")
    ADD_LIBRARY(viabLabLib STATIC ${Viab_SOURCES} ${Viab_HEADERS})  # creation de la librairie  
	
else()
    message("Building ViabLab as executable ")
	option(VIABLAB_LIBRARY_STATIC "Compile foo_library as static library" OFF)
	ADD_EXECUTABLE(${PROJECT_NAME}  ${Viab_SOURCES} ${Viab_HEADERS})  # Executable creation
    set_target_properties(${PROJECT_NAME} PROPERTIES ENABLE_EXPORTS 1)
    # Link all shared libraries to the executable
	 message("Dynamic libs value from cmake ${CMAKE_DL_LIBS}  ")
    target_link_libraries(${PROJECT_NAME} dl spdlog::spdlog_header_only)
    #target_link_libraries(${PROJECT_NAME} PRIVATE dl spdlog::spdlog $<$<BOOL:${MINGW}>:ws2_32>)
	target_compile_definitions(${PROJECT_NAME} PRIVATE VIABLAB_LIBRARY_BUILD)

    if(WIN32)
        #Invoke ar to generate a .dll.a from the .obj files, required to link the executable
        add_custom_command(TARGET ${PROJECT_NAME} PRE_LINK
            COMMAND sh -c "${CMAKE_AR} cr lib${PROJECT_NAME}.dll.a $$(find . -name '*.obj')"
            COMMENT "Generating lib${PROJECT_NAME}.dll.a for external linking")
    endif(WIN32)


    # Create a shared library for each CPP file
    foreach(MODEL_FILE ${MODEL_FILES})
        get_filename_component(LIB_NAME ${MODEL_FILE} NAME_WE)  # Get the file name without extension
        message("Found model : ${MODEL_FILE}")
    	set(CMAKE_SHARED_LIBRARY_ENABLE_EXPORTS TRUE)
        add_library(${LIB_NAME} SHARED ${MODEL_FILE})           # Create a shared library
        target_compile_options(${LIB_NAME} PRIVATE -fPIC)   # Add -fPIC flag
    	set_target_properties(${LIB_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        target_link_libraries(${LIB_NAME} ${PROJECT_NAME})
    endforeach()
    
endif(BUILD_LIB)

unset(BUILD_LIB CACHE)
