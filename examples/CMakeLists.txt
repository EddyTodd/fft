add_executable(fftlab-example-basic basic.cpp)
target_link_libraries(fftlab-example-basic PRIVATE fftlab::fftlab)
fftlab_configure_private_target(fftlab-example-basic)

if(BUILD_TESTING)
  add_test(NAME fftlab.example.basic COMMAND fftlab-example-basic)
endif()
