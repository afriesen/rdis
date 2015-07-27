cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW) # set so link_directories() treats paths relative to the source dir

project(rdis)



# ----------- parameters for the build -----------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

# parameters for deciding which libraries and functionality to include in the build
set(USE_PATOH TRUE CACHE BOOL "Use the PaToH hypergraph partitioning library")
set(USE_HMETIS FALSE CACHE BOOL "Use the HMetis hypergraph partitioning library")
set(USE_LEVMAR FALSE CACHE BOOL "Use the levmar optimization library" )

# if your boost install is not being successfully found, then set the path to it here
# (or set BOOST_INCLUDEDIR and BOOST_LIBRARYDIR)
#set(BOOST_ROOT <your boost install directory>)

# you can set the following variables to specify paths to external libraries
#   EXT_LIB_DIR # general path(s) for external libraries
#   PATOH_DIR # path for PaToH library
#   HMETIS_DIR # path to HMetis library
#   LEVMAR_DIR # path to levmar library


# ------ set up the build (ideally, shouldn't need to edit below here) ------

set(REQ_CXX_FLAGS "-std=c++11 -fPIC -Wno-deprecated-register")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${REQ_CXX_FLAGS} ${MY_CXX_FLAGS}")

include_directories(.)
include_directories(src)
include_directories(external/include)

if ( NOT ROSETTA_EXISTS )
    add_definitions(-DBOOSTNS=boost)
    add_definitions(-DBOOSTPATH=<boost)
    add_definitions(-DNOROSETTA)

    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_MULTITHREADED OFF)

    # find Boost 1.55.0 (or later) on this system
    find_package(Boost 1.55.0 REQUIRED
        COMPONENTS chrono filesystem program_options system)

    message(STATUS "boost include dir: ${Boost_INCLUDE_DIR}")
    message(STATUS "boost libs dir: ${Boost_LIBRARY_DIR}")

    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})

else()
    include(rosetta_options.cmake)
endif()


add_definitions(-DBOOST1_55_0_NO_CXX11_VARIADIC_TEMPLATES)
add_definitions(-DBOOST1_55_0_CONTAINER_PERFECT_FORWARDING)
add_definitions(-DBOOST1_55_0_MPL_CFG_NO_ADL_BARRIER_NAMESPACE)


set( MAIN_SOURCES
    src/common.cpp
    src/common.h
    src/Component.cpp
    src/Component.h
    src/ComponentCache.cpp
    src/ComponentCache.h
    src/ConnectivityGraph.cpp
    src/ConnectivityGraph.h
    src/EulerTourTree.cpp
    src/EulerTourTree.h
    src/Factor.cpp
    src/Factor.h
    src/IntrusivePtrPool.h
    src/NonlinearProductFactor.cpp
    src/NonlinearProductFactor.h
    src/OptimizableFunction.cpp
    src/OptimizableFunction.h
    src/OptimizableFunctionGenerator.cpp
    src/OptimizableFunctionGenerator.h
    src/Optimizer.cpp
    src/Optimizer.h
    src/SimpleSumFactor.cpp
    src/SimpleSumFactor.h
    src/PolynomialFunction.cpp
    src/PolynomialFunction.h
    src/RDISOptimizer.cpp
    src/RDISOptimizer.h
    src/OptimizerVariableData.cpp
    src/OptimizerVariableData.h
    src/Semiring.h
    src/SplayTree.h
    src/State.cpp
    src/State.h
    src/SubspaceOptimizer.cpp
    src/SubspaceOptimizer.h
    src/Variable.cpp
    src/Variable.h
    src/VariableDomain.cpp
    src/VariableDomain.h
    src/optimizers/BCDOptimizer.cpp
    src/optimizers/BCDOptimizer.h
    src/optimizers/CGDSubspaceOptimizer.cpp
    src/optimizers/CGDSubspaceOptimizer.h
    src/optimizers/EMOptimizer.cpp
    src/optimizers/EMOptimizer.h
    src/optimizers/EMSubspaceOptimizer.cpp
    src/optimizers/EMSubspaceOptimizer.h
    src/optimizers/LMSubspaceOptimizer.cpp
    src/optimizers/LMSubspaceOptimizer.h
    src/util/numeric.cpp
    src/util/numeric.h
    src/util/signal_handler.cpp
    src/util/signal_handler.h
    src/util/utility.cpp
    src/util/utility.h
    )


# --------- find the external libraries to build with ---------

set( EXTERNAL_SOURCES external/include/minimize_nrc.h )
set( EXTERNAL_LIBS ${Boost_LIBRARIES} )
set( MY_EXTINCL_DIR external/include )

# (weakly) check whether this is mac os x or linux (this can fail, but should be fine for now)
if ( APPLE )
    set( PLATFORM_LIB_DIR external/macosx )
elseif ( UNIX )
    set( PLATFORM_LIB_DIR external/linux )
    list( APPEND EXTERNAL_LIBS pthread rt )
else()
    message( WARNING "Unsupported platform (only OS X and Linux are supported), so building will probably fail" )
endif()

# get whether this system is 32 or 64 bit
EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
message( STATUS "Detected architecture: ${ARCHITECTURE}" )

if ( PLATFORM_LIB_DIR )
    if( ${ARCHITECTURE} STREQUAL "x86_64" )
        set( PLATFORM_LIB_DIR "${PLATFORM_LIB_DIR}/64" )
    else()
        set( PLATFORM_LIB_DIR "${PLATFORM_LIB_DIR}/32" )
    endif()

    set( MY_EXTLIB_DIR ${PLATFORM_LIB_DIR} )
endif()

message( STATUS "using platform-specific external lib dir: ${MY_EXTLIB_DIR}" )

if ( USE_PATOH )
    find_file( PATOH_INCL patoh.h PATHS ${MY_EXTINCL_DIR} ${EXT_INCL_DIR} ${PATOH_DIR} PATH_SUFFIXES patoh include )
    find_library( PATOH_LIB patoh libpatoh PATHS ${MY_EXTLIB_DIR} ${EXT_LIB_DIR} ${PATOH_DIR} PATH_SUFFIXES patoh lib )

    if ( PATOH_INCL AND PATOH_LIB )
        message(STATUS "Found PaToH lib at ${PATOH_LIB}")
        message(STATUS "Found PaToH include at ${PATOH_INCL} and lib at ${PATOH_LIB}")
        add_definitions(-DUSE_PATOH)
        get_filename_component(PATOH_INCLUDE_DIR ${PATOH_INCL} DIRECTORY ABSOLUTE)
        include_directories( ${PATOH_INCLUDE_DIR} )
        list(APPEND EXTERNAL_SOURCES ${PATOH_INCL} )
        list(APPEND EXTERNAL_LIBS ${PATOH_LIB} )
    else()
        message(FATAL_ERROR "Did not find PaToH include or library, so cannot use it")
    endif()
endif()


if ( USE_HMETIS )
    find_file( HMETIS_INCL hmetis.h PATHS ${MY_EXTINCL_DIR} ${EXT_INCL_DIR} ${HMETIS_DIR} PATH_SUFFIXES hmetis include )
    find_library( HMETIS_LIB hmetis libhmetis PATHS ${MY_EXTLIB_DIR} ${EXT_LIB_DIR} ${HMETIS_DIR} PATH_SUFFIXES hmetis lib )

    if ( HMETIS_LIB AND HMETIS_INCL )
        message(STATUS "Found HMetis include at ${HMETIS_INCL} and lib at ${HMETIS_LIB}")
        add_definitions(-DUSE_HMETIS)
        get_filename_component(HMETIS_INCLUDE_DIR ${HMETIS_INCL} DIRECTORY ABSOLUTE)
        include_directories( ${HMETIS_INCLUDE_DIR} )
        list(APPEND EXTERNAL_SOURCES ${HMETIS_INCL} )
        list(APPEND EXTERNAL_LIBS ${HMETIS_LIB} )
    else()
        message(FATAL_ERROR "Did not find HMetis include or library, so cannot use it")
    endif()
endif()


if ( USE_LEVMAR )
    find_file( LEVMAR_INCL levmar.h PATHS ${MY_EXTINCL_DIR} ${EXT_INCL_DIR} ${LEVMAR_DIR} PATH_SUFFIXES levmar include )
    find_library( LEVMAR_LIB levmar liblevmar PATHS ${MY_EXTLIB_DIR} ${EXT_LIB_DIR} ${LEVMAR_DIR} PATH_SUFFIXES levmar lib )

    if ( LEVMAR_INCL AND LEVMAR_LIB )
        message(STATUS "Found levmar include at ${LEVMAR_INCL} and library at ${LEVMAR_LIB}")
        add_definitions(-DUSE_LEVMAR)
        find_package(LAPACK REQUIRED)
        find_package(BLAS REQUIRED)
        get_filename_component(LEVMAR_INCLUDE_DIR ${LEVMAR_INCL} DIRECTORY ABSOLUTE)
        include_directories( ${LEVMAR_INCLUDE_DIR} )
        list(APPEND EXTERNAL_SOURCES ${LEVMAR_INCL})
        list(APPEND EXTERNAL_LIBS ${LEVMAR_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
    else()
        message(FATAL_ERROR "Did not find levmar include or library, so cannot use it")
    endif()
endif()


# ------ define the build targets ------

# RDIS library
add_library( rdis STATIC ${MAIN_SOURCES} ${EXTERNAL_SOURCES} )

# debug / test executable
add_executable( testRDIS src/main.cpp )
target_link_libraries( testRDIS PUBLIC rdis ${EXTERNAL_LIBS} )

# high-dim sinusoid executable
add_executable( optSinusoid src/optimize_sinusoid.cpp )
target_link_libraries( optSinusoid PUBLIC rdis ${EXTERNAL_LIBS} )

# bundle adjustment executable
add_subdirectory( src/bundleadjust )

## protein folding executable
if(ROSETTA_EXISTS)
    add_subdirectory( src/proteinfolding )
endif()
