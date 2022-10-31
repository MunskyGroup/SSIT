CMAKE_MINIMUM_REQUIRED(VERSION 3.14)
PROJECT(FSP_CVODE_MEX)

option(USE_OPENMP "allow the libraries to use OpenMP (Note: Sundials must be compiled with OpenMP support!)" OFF)

set(SOLVER_LIBS FspCvode FspCvodeForwardSens
        )

set(MEX_LIBS
        FspCVodeMex FspSensCVodeMex
        )

set(TEST_PROGRAMS
        test_cvode_fsp
        test_cvodes_fsp
        )

set(EXT_LIBRARIES
        -lsundials_cvodes
        -lsundials_nvecserial
        )

if (APPLE)
    if (USE_OPENMP)
        message("OpenMP is not supported for Matlab on MacOS.")
        SET(USE_OPENMP OFF CACHE BOOL "allow the libraries to use OpenMP (Note: Sundials must be compiled with OpenMP support!)" FORCE)
    endif(USE_OPENMP)

else()
    if (USE_OPENMP)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSENVECOMP")
        message("Use OpenMP")
        set(EXT_LIBRARIES "${EXT_LIBRARIES} -lsundials_nvecopenmp")
    endif(USE_OPENMP)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

FIND_PACKAGE(Threads REQUIRED)
FIND_PACKAGE(MATLAB REQUIRED)

# Set the directory to put our compiled files
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

# Where to put .a files (archive for static libraries)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Where to put .so files (dynamic libraries)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# Where to put the executable files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Set the folders to include
include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

foreach( LIB ${SOLVER_LIBS})
    add_library(${LIB} SHARED ${CMAKE_SOURCE_DIR}/src/${LIB}.c ${CMAKE_SOURCE_DIR}/src/${LIB}.h)
    target_link_libraries(${LIB} ${EXT_LIBRARIES})
endforeach()

foreach ( PROG ${TEST_PROGRAMS})
    add_executable(${PROG} ${CMAKE_SOURCE_DIR}/test/${PROG}.c)
    target_link_libraries(${PROG} ${SOLVER_LIBS} ${EXT_LIBRARIES})
endforeach()

foreach( LIB ${MEX_LIBS})
    matlab_add_mex(NAME ${LIB} SRC ${CMAKE_SOURCE_DIR}/${LIB}.c LINK_TO ${EXT_LIBRARIES} ${SOLVER_LIBS} R2018a)
endforeach()



