cmake_minimum_required(VERSION 3.23)
project(fftlab VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)

option(FFTLAB_BUILD_TESTS "Build correctness tests" ${PROJECT_IS_TOP_LEVEL})
option(FFTLAB_BUILD_EXAMPLES "Build usage examples" ${PROJECT_IS_TOP_LEVEL})
option(FFTLAB_INSTALL "Enable install/package rules" ${PROJECT_IS_TOP_LEVEL})

add_library(fftlab src/kernel.cpp)
add_library(fftlab::fftlab ALIAS fftlab)
target_compile_features(fftlab PUBLIC cxx_std_23)
set_target_properties(fftlab PROPERTIES
  CXX_EXTENSIONS OFF
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
target_include_directories(fftlab PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

if(FFTLAB_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(FFTLAB_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(FFTLAB_INSTALL)
  include(CMakePackageConfigHelpers)

  install(TARGETS fftlab
    EXPORT fftlabTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(EXPORT fftlabTargets
    FILE fftlabTargets.cmake
    NAMESPACE fftlab::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab)

  file(WRITE "${PROJECT_BINARY_DIR}/fftlabConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/fftlabTargets.cmake\")\n")
  write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/fftlabConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
  install(FILES
    "${PROJECT_BINARY_DIR}/fftlabConfig.cmake"
    "${PROJECT_BINARY_DIR}/fftlabConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab)
  install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/fftlab)
endif()
