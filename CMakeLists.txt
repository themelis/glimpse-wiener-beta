cmake_minimum_required(VERSION 3.0)

# Paths to external libraries
set(NICAEA_URL http://www.cosmostat.org/wp-content/uploads/2014/12/nicaea_2.5.tgz)
set(NICAEA_MD5 27c4d6d46f234133fc69130e6460c2d4)

set(SPARSE2D_URL https://github.com/CosmoStat/Sparse2D.git)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(FindPkgConfig)
include(ExternalProject)

project(glimpse)

# Include external libraries
find_package(GSL       REQUIRED)
find_package(CCFits    REQUIRED)
find_package(FFTW      REQUIRED)
find_package(NFFT      REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Boost     COMPONENTS program_options REQUIRED)

# Include external projects
## NICAEA
ExternalProject_Add(NICAEA
        PREFIX NICAEA
        URL ${NICAEA_URL}
        URL_MD5 ${NICAEA_MD5}
        PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/nicaea.patch)
        set(NICAEA_LIBRARY_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/lib)
        set(NICAEA_INCLUDE_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/include)
        set(NICAEA_LIBRARIES -lnicaea)

## Sparse2D
ExternalProject_Add(SPARSE2D
         PREFIX SPARSE2D
         GIT_REPOSITORY ${SPARSE2D_URL}
         GIT_TAG  22e24e0
         CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern)
        set(SPARSE2D_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib)
        set(SPARSE2D_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include)
        set(SPARSE2D_LIBRARIES -lsparse2d)
        
include_directories(${SPARSE2D_INCLUDE_DIR}
		    ${GSL_INCLUDE_DIR}
		    ${CCFITS_INCLUDE_DIR}
		    ${CFITSIO_INCLUDE_DIR}
		    ${Boost_INCLUDE_DIR}
		    ${NICAEA_INCLUDE_DIR}
		    ${NFFT_INCLUDE_DIR}
		    ${FFTW_INCLUDE_DIR}
		    ${ARMADILLO_INCLUDE_DIRS})
		    
link_directories(${SPARSE2D_LIBRARY_DIR}
		    ${GSL_LIBRARY_DIR}
		    ${CCFITS_LIBRARY_DIR}
		    ${CFITSIO_LIBRARY_DIR}
		    ${Boost_LIBRARY_DIR}
		    ${NICAEA_LIBRARY_DIR}
		    ${NFFT_LIBRARY_DIR}
		    ${FFTW_LIBRARY_DIR})

set(GLIMPSE_SRC src/survey.cpp
		src/redshift_distribution.cpp
		src/field.cpp
		src/wiener_filtering.cpp
		src/surface_reconstruction.cpp
		src/density_reconstruction.cpp
		src/starlet_2d.cpp
		src/wavelet_transform.cpp
		src/gpu_utils.c)
		
# VERSIONING from git
execute_process(
  COMMAND git describe
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(
  ${CMAKE_SOURCE_DIR}/src/version.h.in
  ${CMAKE_SOURCE_DIR}/src/version.h
)

find_package(CUDA)

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}  -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fpermissive -fopenmp -std=c++11 -Wno-write-strings -Wno-literal-suffix")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG_FITS")

if(${CUDA_FOUND})
    message("Compiling CUDA accelerated reconstruction code, with 3D support")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCUDA_ACC")
    cuda_add_executable(glimpse src/glimpse.cpp ${GLIMPSE_SRC} src/spg.cu src/spg.cpp)
    cuda_add_cufft_to_target(glimpse)
else(${CUDA_FOUND})
    message("Compiling without CUDA acceleration")
    add_executable(glimpse src/glimpse.cpp ${GLIMPSE_SRC})
endif(${CUDA_FOUND})

add_dependencies(glimpse NICAEA SPARSE2D)

target_link_libraries(glimpse ${NICAEA_LIBRARIES}
                              ${SPARSE2D_LIBRARIES}
                              ${Boost_LIBRARIES} 
                              ${GSL_LIBRARIES}
                              ${NFFT_LIBRARIES}
                              ${FFTW_LIBRARIES}
                              ${CCFITS_LIBRARY}
                              ${CFITSIO_LIBRARY}
                              ${ARMADILLO_LIBRARIES})
