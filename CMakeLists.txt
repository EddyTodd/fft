cmake_minimum_required(VERSION 3.20)
project(fft VERSION 3.0.0 LANGUAGES CXX)

option(FFT_NATIVE "Enable host-specific ISA optimizations" OFF)
option(FFT_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(FFT_LTO "Enable interprocedural optimization when supported" OFF)

add_library(fftlab STATIC
    src/common.cpp
    src/power2.cpp
    src/arbitrary.cpp
    src/research.cpp
    src/planned.cpp
)
add_library(fftlab::fftlab ALIAS fftlab)
target_include_directories(fftlab
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE src
)
target_compile_features(fftlab PUBLIC cxx_std_23)
set_target_properties(fftlab PROPERTIES CXX_EXTENSIONS OFF)

add_executable(fft src/main.cpp)
target_link_libraries(fft PRIVATE fftlab::fftlab)
target_compile_features(fft PRIVATE cxx_std_23)
set_target_properties(fft PROPERTIES CXX_EXTENSIONS OFF)

add_executable(fft-plan src/plan_main.cpp)
target_link_libraries(fft-plan PRIVATE fftlab::fftlab)
target_compile_features(fft-plan PRIVATE cxx_std_23)
set_target_properties(fft-plan PROPERTIES CXX_EXTENSIONS OFF)

foreach(target fftlab fft fft-plan)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    endif()
endforeach()

if(FFT_NATIVE)
    if(MSVC)
        target_compile_options(fftlab PRIVATE /arch:AVX2)
    else()
        target_compile_options(fftlab PRIVATE -march=native)
    endif()
endif()

if(FFT_SANITIZERS AND NOT MSVC AND CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")
    foreach(target fftlab fft fft-plan)
        target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endforeach()
endif()

if(FFT_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_ok OUTPUT ipo_error)
    if(NOT ipo_ok)
        message(FATAL_ERROR "IPO/LTO requested but unsupported: ${ipo_error}")
    endif()
    foreach(target fftlab fft fft-plan)
        set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endforeach()
endif()

enable_testing()
add_test(NAME fft_self_test COMMAND fft --self-test)
add_test(NAME fft_accuracy_smoke COMMAND fft --verify --algorithm split-radix --size 64 --signal random)
add_test(NAME fft_rader_smoke COMMAND fft --verify --algorithm rader --size 31 --signal dynamic-range)
add_test(NAME fft_planned_self_test COMMAND fft-plan --self-test)
add_test(NAME fft_planned_smoke COMMAND fft-plan --benchmark --size 256 --samples 5 --target-ms 0.2)
