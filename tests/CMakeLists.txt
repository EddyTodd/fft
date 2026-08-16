include(${PROJECT_SOURCE_DIR}/cmake/HeaderSelfContainment.cmake)
add_public_header_self_containment_target(
  fftlab-header-self-containment
  fftlab
  "${PROJECT_SOURCE_DIR}/include")
fftlab_apply_warnings(fftlab-header-self-containment)

function(fftlab_add_test name source)
  set(target "fftlab-test-${name}")
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE fftlab::fftlab)
  fftlab_configure_private_target(${target})
  add_test(NAME "fftlab.${name}" COMMAND ${target})
endfunction()

fftlab_add_test(core test_correctness.cpp)
fftlab_add_test(kernel test_kernel.cpp)
fftlab_add_test(plan-contract test_plan_contract.cpp)

if(FFTLAB_BUILD_PACKAGE_TESTS AND NOT FFTLAB_ENABLE_SANITIZERS)
  include(${PROJECT_SOURCE_DIR}/cmake/PackageConsumerTest.cmake)
  register_package_consumer_test(
    fftlab.package-consumer
    "${PROJECT_SOURCE_DIR}/tests/package-consumer")
endif()
