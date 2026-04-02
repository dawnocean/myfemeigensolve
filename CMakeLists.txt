cmake_minimum_required(VERSION 3.16)

# ── Find MFEM and use its compiler BEFORE project() ─────────────────────────
set(MFEM_DIR "$ENV{MFEM_DIR}" CACHE PATH "MFEM installation prefix")
if(NOT MFEM_DIR)
    file(GLOB _mfem_candidates "$ENV{HOME}/spack/opt/spack/linux-*/mfem-*/")
    list(SORT _mfem_candidates ORDER DESCENDING)
    if(_mfem_candidates)
        list(GET _mfem_candidates 0 MFEM_DIR)
    endif()
endif()

if(NOT EXISTS "${MFEM_DIR}/share/mfem/config.mk")
    message(FATAL_ERROR "Cannot find MFEM config.mk. Set MFEM_DIR to the MFEM install prefix.")
endif()

# Extract the MPI compiler wrapper from config.mk and set it before project()
file(READ "${MFEM_DIR}/share/mfem/config.mk" _mfem_config)
string(REGEX MATCH "MFEM_CXX[ \t]*=[ \t]*([^\n]*)" _ "${_mfem_config}")
set(CMAKE_CXX_COMPILER "${CMAKE_MATCH_1}" CACHE STRING "C++ compiler" FORCE)

project(eigensolve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "Using MFEM at: ${MFEM_DIR}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")

# Parse remaining MFEM config variables
string(REGEX MATCH "MFEM_TPLFLAGS[ \t]*=[ \t]*([^\n]*)" _ "${_mfem_config}")
set(MFEM_TPLFLAGS "${CMAKE_MATCH_1}")

string(REGEX MATCH "MFEM_EXT_LIBS[ \t]*=[ \t]*([^\n]*)" _ "${_mfem_config}")
set(MFEM_EXT_LIBS "${CMAKE_MATCH_1}")

# ── nlohmann/json (header-only) ──────────────────────────────────────────────
include(FetchContent)
FetchContent_Declare(
    json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
)
FetchContent_MakeAvailable(json)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(eigensolve
    src/main.cpp
    src/mesh_generator.cpp
    src/waveguide_solver.cpp
    src/eigensolver.cpp
)

target_include_directories(eigensolve PRIVATE
    ${MFEM_DIR}/include
)

# Parse TPL include flags
string(REGEX MATCHALL "-I[^ ]+" _tpl_incs "${MFEM_TPLFLAGS}")
foreach(_inc ${_tpl_incs})
    string(REPLACE "-I" "" _inc_dir "${_inc}")
    target_include_directories(eigensolve PRIVATE "${_inc_dir}")
endforeach()

# Link MFEM and its dependencies
target_link_libraries(eigensolve PRIVATE
    nlohmann_json::nlohmann_json
    -L${MFEM_DIR}/lib -lmfem
)

separate_arguments(_ext_libs UNIX_COMMAND "${MFEM_EXT_LIBS}")
target_link_libraries(eigensolve PRIVATE ${_ext_libs})

# ── Optional SLEPc eigensolver ────────────────────────────────────────────────
option(EIGENSOLVE_USE_SLEPC "Enable SLEPc eigensolver backend" OFF)

if(EIGENSOLVE_USE_SLEPC)
    # Auto-detect PETSc and SLEPc from spack if not specified
    set(PETSC_DIR "$ENV{PETSC_DIR}" CACHE PATH "PETSc installation prefix")
    set(SLEPC_DIR "$ENV{SLEPC_DIR}" CACHE PATH "SLEPc installation prefix")

    if(NOT PETSC_DIR)
        file(GLOB _petsc_candidates "$ENV{HOME}/spack/opt/spack/linux-*/petsc-*/")
        list(SORT _petsc_candidates ORDER DESCENDING)
        if(_petsc_candidates)
            list(GET _petsc_candidates 0 PETSC_DIR)
        endif()
    endif()

    if(NOT SLEPC_DIR)
        file(GLOB _slepc_candidates "$ENV{HOME}/spack/opt/spack/linux-*/slepc-*/")
        list(SORT _slepc_candidates ORDER DESCENDING)
        if(_slepc_candidates)
            list(GET _slepc_candidates 0 SLEPC_DIR)
        endif()
    endif()

    if(NOT EXISTS "${PETSC_DIR}/include/petsc.h")
        message(FATAL_ERROR "PETSc not found at ${PETSC_DIR}. Set PETSC_DIR.")
    endif()
    if(NOT EXISTS "${SLEPC_DIR}/include/slepc.h")
        message(FATAL_ERROR "SLEPc not found at ${SLEPC_DIR}. Set SLEPC_DIR.")
    endif()

    message(STATUS "PETSc found at: ${PETSC_DIR}")
    message(STATUS "SLEPc found at: ${SLEPC_DIR}")

    target_compile_definitions(eigensolve PRIVATE EIGENSOLVE_USE_SLEPC)
    target_include_directories(eigensolve PRIVATE
        ${PETSC_DIR}/include
        ${SLEPC_DIR}/include
    )
    target_link_libraries(eigensolve PRIVATE
        -L${SLEPC_DIR}/lib -lslepc
        -L${PETSC_DIR}/lib -lpetsc
    )
    target_sources(eigensolve PRIVATE src/slepc_eigensolver.cpp)
endif()

# ── Optional ARPACK eigensolver ─────────────────────────────────────────────
option(EIGENSOLVE_USE_ARPACK "Enable ARPACK/PARPACK eigensolver backend" OFF)

if(EIGENSOLVE_USE_ARPACK)
    set(ARPACK_DIR "$ENV{ARPACK_DIR}" CACHE PATH "ARPACK-ng installation prefix")

    if(NOT ARPACK_DIR)
        file(GLOB _arpack_candidates "$ENV{HOME}/spack/opt/spack/linux-*/arpack-ng-*/")
        list(SORT _arpack_candidates ORDER DESCENDING)
        if(_arpack_candidates)
            list(GET _arpack_candidates 0 ARPACK_DIR)
        endif()
    endif()

    if(NOT EXISTS "${ARPACK_DIR}/lib/libparpack.so" AND
       NOT EXISTS "${ARPACK_DIR}/lib/libparpack.a")
        message(FATAL_ERROR "PARPACK not found at ${ARPACK_DIR}. Set ARPACK_DIR.")
    endif()

    message(STATUS "ARPACK-ng found at: ${ARPACK_DIR}")

    target_compile_definitions(eigensolve PRIVATE EIGENSOLVE_USE_ARPACK)
    target_include_directories(eigensolve PRIVATE ${ARPACK_DIR}/include)
    target_link_libraries(eigensolve PRIVATE
        -L${ARPACK_DIR}/lib -lparpack -larpack
    )
    target_sources(eigensolve PRIVATE src/arpack_eigensolver.cpp)
endif()

# ── RPATH ────────────────────────────────────────────────────────────────────
set(_rpath "${MFEM_DIR}/lib")
if(EIGENSOLVE_USE_SLEPC)
    set(_rpath "${_rpath};${PETSC_DIR}/lib;${SLEPC_DIR}/lib")
endif()
if(EIGENSOLVE_USE_ARPACK)
    set(_rpath "${_rpath};${ARPACK_DIR}/lib")
endif()
set_target_properties(eigensolve PROPERTIES
    INSTALL_RPATH "${_rpath}"
    BUILD_RPATH   "${_rpath}"
)
