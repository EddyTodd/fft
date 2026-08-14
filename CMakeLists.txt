cmake_minimum_required(VERSION 3.23)
project(fftlab VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(FFTLAB_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(FFTLAB_WARNINGS_AS_ERRORS "Treat fftlab warnings as errors" ${PROJECT_IS_TOP_LEVEL})
option(FFTLAB_BUILD_TESTS "Build deterministic FFT correctness tests" ${PROJECT_IS_TOP_LEVEL})

include(cmake/ProjectOptions.cmake)

set(FFTLAB_PUBLIC_HEADERS
  include/fftlab/arbitrary_algorithms.hpp
  include/fftlab/arbitrary_plan.hpp
  include/fftlab/codelet.hpp
  include/fftlab/convolution_plan.hpp
  include/fftlab/fft.hpp
  include/fftlab/fftlab.hpp
  include/fftlab/good_thomas_plan.hpp
  include/fftlab/kernel.hpp
  include/fftlab/mixed_plan.hpp
  include/fftlab/oracle.hpp
  include/fftlab/plan.hpp
  include/fftlab/planner.hpp
  include/fftlab/power2_algorithms.hpp
  include/fftlab/types.hpp
  include/fftlab/version.hpp)

add_library(fftlab src/kernel.cpp)
add_library(fftlab::fftlab ALIAS fftlab)
target_compile_features(fftlab PUBLIC cxx_std_23)
set_target_properties(fftlab PROPERTIES
  CXX_EXTENSIONS OFF
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME fftlab)
target_include_directories(fftlab PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_sources(fftlab PUBLIC
  FILE_SET public_headers TYPE HEADERS
  BASE_DIRS ${PROJECT_SOURCE_DIR}/include
  FILES ${FFTLAB_PUBLIC_HEADERS})
fftlab_apply_warnings(fftlab)
fftlab_apply_sanitizers(fftlab)

if(BUILD_TESTING AND FFTLAB_BUILD_TESTS)
  add_subdirectory(tests)
endif()

include(cmake/Install.cmake)
