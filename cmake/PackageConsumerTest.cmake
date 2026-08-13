cmake_minimum_required(VERSION 3.23)

function(register_package_consumer_test name consumer_source_dir)
  if(NOT PROJECT_IS_TOP_LEVEL)
    return()
  endif()

  add_test(
    NAME ${name}
    COMMAND ${CMAKE_COMMAND}
      -DPACKAGE_CONSUMER_RUN=ON
      "-DPACKAGE_CONSUMER_PROJECT_BINARY_DIR=${CMAKE_BINARY_DIR}"
      "-DPACKAGE_CONSUMER_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/package-consumer-prefix"
      "-DPACKAGE_CONSUMER_SOURCE_DIR=${consumer_source_dir}"
      "-DPACKAGE_CONSUMER_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/package-consumer-build"
      "-DPACKAGE_CONSUMER_CONFIG=$<CONFIG>"
      "-DPACKAGE_CONSUMER_GENERATOR=${CMAKE_GENERATOR}"
      "-DPACKAGE_CONSUMER_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
      "-DPACKAGE_CONSUMER_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
      "-DPACKAGE_CONSUMER_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
      "-DPACKAGE_CONSUMER_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}"
      "-DPACKAGE_CONSUMER_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}"
      "-DPACKAGE_CONSUMER_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}"
      "-DPACKAGE_CONSUMER_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}"
      -P "${CMAKE_CURRENT_FUNCTION_LIST_FILE}"
  )
endfunction()

if(PACKAGE_CONSUMER_RUN)
  foreach(required
      PACKAGE_CONSUMER_PROJECT_BINARY_DIR
      PACKAGE_CONSUMER_INSTALL_PREFIX
      PACKAGE_CONSUMER_SOURCE_DIR
      PACKAGE_CONSUMER_BINARY_DIR)
    if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
      message(FATAL_ERROR "Missing package-consumer variable: ${required}")
    endif()
  endforeach()

  function(_package_consumer_execute label)
    execute_process(
      COMMAND ${ARGN}
      RESULT_VARIABLE result
      OUTPUT_VARIABLE output
      ERROR_VARIABLE error
    )
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "${label} failed (${result})\n${output}\n${error}")
    endif()
  endfunction()

  file(REMOVE_RECURSE
    "${PACKAGE_CONSUMER_INSTALL_PREFIX}"
    "${PACKAGE_CONSUMER_BINARY_DIR}"
  )

  set(install_command
    "${CMAKE_COMMAND}" --install "${PACKAGE_CONSUMER_PROJECT_BINARY_DIR}"
    --prefix "${PACKAGE_CONSUMER_INSTALL_PREFIX}"
  )
  if(PACKAGE_CONSUMER_CONFIG)
    list(APPEND install_command --config "${PACKAGE_CONSUMER_CONFIG}")
  endif()
  _package_consumer_execute("package install" ${install_command})

  set(configure_command
    "${CMAKE_COMMAND}"
    -S "${PACKAGE_CONSUMER_SOURCE_DIR}"
    -B "${PACKAGE_CONSUMER_BINARY_DIR}"
    "-DCMAKE_PREFIX_PATH=${PACKAGE_CONSUMER_INSTALL_PREFIX}"
  )
  if(PACKAGE_CONSUMER_GENERATOR)
    list(APPEND configure_command -G "${PACKAGE_CONSUMER_GENERATOR}")
  endif()
  if(PACKAGE_CONSUMER_GENERATOR_PLATFORM)
    list(APPEND configure_command -A "${PACKAGE_CONSUMER_GENERATOR_PLATFORM}")
  endif()
  if(PACKAGE_CONSUMER_GENERATOR_TOOLSET)
    list(APPEND configure_command -T "${PACKAGE_CONSUMER_GENERATOR_TOOLSET}")
  endif()
  if(PACKAGE_CONSUMER_TOOLCHAIN_FILE)
    list(APPEND configure_command "-DCMAKE_TOOLCHAIN_FILE=${PACKAGE_CONSUMER_TOOLCHAIN_FILE}")
  elseif(PACKAGE_CONSUMER_CXX_COMPILER)
    list(APPEND configure_command "-DCMAKE_CXX_COMPILER=${PACKAGE_CONSUMER_CXX_COMPILER}")
  endif()
  if(PACKAGE_CONSUMER_CONFIG)
    list(APPEND configure_command "-DCMAKE_BUILD_TYPE=${PACKAGE_CONSUMER_CONFIG}")
  endif()
  if(PACKAGE_CONSUMER_OSX_ARCHITECTURES)
    list(APPEND configure_command "-DCMAKE_OSX_ARCHITECTURES=${PACKAGE_CONSUMER_OSX_ARCHITECTURES}")
  endif()
  if(PACKAGE_CONSUMER_OSX_SYSROOT)
    list(APPEND configure_command "-DCMAKE_OSX_SYSROOT=${PACKAGE_CONSUMER_OSX_SYSROOT}")
  endif()
  if(PACKAGE_CONSUMER_OSX_DEPLOYMENT_TARGET)
    list(APPEND configure_command "-DCMAKE_OSX_DEPLOYMENT_TARGET=${PACKAGE_CONSUMER_OSX_DEPLOYMENT_TARGET}")
  endif()
  _package_consumer_execute("consumer configure" ${configure_command})

  set(build_command "${CMAKE_COMMAND}" --build "${PACKAGE_CONSUMER_BINARY_DIR}")
  if(PACKAGE_CONSUMER_CONFIG)
    list(APPEND build_command --config "${PACKAGE_CONSUMER_CONFIG}")
  endif()
  _package_consumer_execute("consumer build" ${build_command})
endif()
