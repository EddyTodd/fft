function(add_fftlab_test name source)
  add_executable(fftlab-test-${name} ${source})
  target_link_libraries(fftlab-test-${name} PRIVATE fftlab::fftlab)
  add_test(NAME fftlab.${name} COMMAND fftlab-test-${name})
endfunction()

add_fftlab_test(core test_correctness.cpp)
add_fftlab_test(kernel test_kernel.cpp)
add_fftlab_test(plan-contract test_plan_contract.cpp)
