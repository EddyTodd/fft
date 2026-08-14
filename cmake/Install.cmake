install(TARGETS fftlab
  EXPORT fftlabTargets
  FILE_SET public_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/cmake/fftlabConfig.cmake.in
  ${PROJECT_BINARY_DIR}/fftlabConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/fftlabConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(EXPORT fftlabTargets
  FILE fftlabTargets.cmake
  NAMESPACE fftlab::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab)

install(FILES
  ${PROJECT_BINARY_DIR}/fftlabConfig.cmake
  ${PROJECT_BINARY_DIR}/fftlabConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fftlab)

install(FILES ${PROJECT_SOURCE_DIR}/LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/fftlab)
