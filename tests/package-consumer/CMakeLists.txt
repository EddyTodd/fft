cmake_minimum_required(VERSION 3.23)
project(fftlab_package_consumer LANGUAGES CXX)

enable_testing()
find_package(fftlab 1 CONFIG REQUIRED)

add_executable(fftlab-package-consumer main.cpp)
target_link_libraries(fftlab-package-consumer PRIVATE fftlab::fftlab)
target_compile_features(fftlab-package-consumer PRIVATE cxx_std_23)
set_target_properties(fftlab-package-consumer PROPERTIES CXX_EXTENSIONS OFF)
add_test(NAME fftlab.package-runtime COMMAND fftlab-package-consumer)
