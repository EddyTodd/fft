cmake_minimum_required(VERSION 3.23)
project(fftlab VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(FFTLAB_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(FFTLAB_WARNINGS_AS_ERRORS "Treat fftlab warnings as errors" OFF)
option(FFTLAB_BUILD_TESTS "Build deterministic FFT correctness tests" ${PROJECT_IS_TOP_LEVEL})

function(fftlab_apply_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
        if(FFTLAB_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
        if(FFTLAB_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()

function(fftlab_apply_sanitizers target)
    if(FFTLAB_ENABLE_SANITIZERS AND NOT MSVC)
        target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
endfunction()

add_library(fftlab src/kernel.cpp)
add_library(fftlab::fftlab ALIAS fftlab)
target_compile_features(fftlab PUBLIC cxx_std_23)
set_target_properties(fftlab PROPERTIES
    CXX_EXTENSIONS OFF
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    EXPORT_NAME fftlab
)
target_include_directories(fftlab
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_sources(fftlab PUBLIC
    FILE_SET public_headers
    TYPE HEADERS
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
    FILES
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
        include/fftlab/version.hpp
)
fftlab_apply_warnings(fftlab)
fftlab_apply_sanitizers(fftlab)

if(BUILD_TESTING AND FFTLAB_BUILD_TESTS)
    include(cmake/HeaderSelfContainment.cmake)
    add_public_header_self_containment_target(
        fftlab-header-self-containment
        fftlab
        "${CMAKE_CURRENT_SOURCE_DIR}/include")
    fftlab_apply_warnings(fftlab-header-self-containment)

    add_executable(fftlab_core_correctness tests/correctness.cpp)
    target_link_libraries(fftlab_core_correctness PRIVATE fftlab::fftlab)
    target_compile_features(fftlab_core_correctness PRIVATE cxx_std_23)
    set_target_properties(fftlab_core_correctness PROPERTIES CXX_EXTENSIONS OFF)

    add_executable(fftlab_kernel_correctness tests/kernel_test.cpp)
    target_link_libraries(fftlab_kernel_correctness PRIVATE fftlab::fftlab)
    target_compile_features(fftlab_kernel_correctness PRIVATE cxx_std_23)
    set_target_properties(fftlab_kernel_correctness PROPERTIES CXX_EXTENSIONS OFF)

    foreach(target fftlab_core_correctness fftlab_kernel_correctness)
        fftlab_apply_warnings(${target})
        fftlab_apply_sanitizers(${target})
    endforeach()

    add_test(NAME fftlab_core_correctness COMMAND fftlab_core_correctness)
    add_test(NAME fftlab_kernel_correctness COMMAND fftlab_kernel_correctness)
    if(NOT FFTLAB_ENABLE_SANITIZERS)
        include(cmake/PackageConsumerTest.cmake)
        register_package_consumer_test(
            fftlab.package-consumer
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/package-consumer")
    endif()
endif()

install(TARGETS fftlab
    EXPORT fftlabTargets
    FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

configure_package_config_file(
    cmake/fftlabConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/fftlabConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/fftlabConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(EXPORT fftlabTargets
    FILE fftlabTargets.cmake
    NAMESPACE fftlab::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/fftlabConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/fftlabConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab
)
install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/fftlab)
