cmake_minimum_required(VERSION 3.23)

get_filename_component(RELEASE_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

function(_release_read path out_var)
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "Required release metadata file is missing: ${path}")
  endif()
  file(READ "${path}" content)
  set(${out_var} "${content}" PARENT_SCOPE)
endfunction()

function(_release_require_contains path needle label)
  _release_read("${path}" content)
  string(FIND "${content}" "${needle}" found)
  if(found EQUAL -1)
    message(FATAL_ERROR "${label} does not contain expected release value: ${needle}")
  endif()
endfunction()

_release_read("${RELEASE_ROOT}/CMakeLists.txt" cmake_text)
string(REGEX MATCH
  "project\\([^\\)]*VERSION[ \t\r\n]+([0-9]+\\.[0-9]+\\.[0-9]+)"
  project_declaration
  "${cmake_text}")
if(NOT project_declaration)
  message(FATAL_ERROR "Could not derive semantic version from top-level project(...) declaration")
endif()
set(release_version "${CMAKE_MATCH_1}")

file(GLOB version_headers "${RELEASE_ROOT}/include/*/version.hpp")
list(LENGTH version_headers version_header_count)
if(NOT version_header_count EQUAL 1)
  message(FATAL_ERROR "Expected exactly one public include/*/version.hpp; found ${version_header_count}")
endif()
list(GET version_headers 0 version_header)

_release_require_contains("${version_header}" "\"${release_version}\"" "public version header")
_release_require_contains("${RELEASE_ROOT}/CITATION.cff" "version: ${release_version}" "CITATION.cff")
_release_require_contains("${RELEASE_ROOT}/CHANGELOG.md" "## [${release_version}]" "CHANGELOG.md")

_release_read("${RELEASE_ROOT}/CITATION.cff" citation_text)
string(REGEX MATCH
  "date-released:[ \t]+[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]"
  citation_date
  "${citation_text}")
if(NOT citation_date)
  message(FATAL_ERROR "CITATION.cff must contain date-released: YYYY-MM-DD")
endif()

if(EXISTS "${RELEASE_ROOT}/pyproject.toml")
  _release_require_contains(
    "${RELEASE_ROOT}/pyproject.toml"
    "version = \"${release_version}\""
    "pyproject.toml")
endif()

if(EXISTS "${RELEASE_ROOT}/benchctl/__init__.py")
  _release_require_contains(
    "${RELEASE_ROOT}/benchctl/__init__.py"
    "__version__ = \"${release_version}\""
    "benchctl/__init__.py")
endif()

if(EXISTS "${RELEASE_ROOT}/SUBJECTS.lock")
  find_program(RELEASE_GIT git REQUIRED)
  file(STRINGS "${RELEASE_ROOT}/SUBJECTS.lock" subject_lines REGEX "^[A-Za-z0-9_-]+[ \t]+[0-9a-f]+")
  foreach(line IN LISTS subject_lines)
    string(REGEX MATCH "^([^ \t]+)[ \t]+([0-9a-f]+)[ \t]+" parsed "${line}")
    if(NOT parsed)
      message(FATAL_ERROR "Malformed SUBJECTS.lock entry: ${line}")
    endif()
    set(subject_name "${CMAKE_MATCH_1}")
    set(subject_sha "${CMAKE_MATCH_2}")
    string(LENGTH "${subject_sha}" subject_sha_length)
    if(NOT subject_sha_length EQUAL 40)
      message(FATAL_ERROR "SUBJECTS.lock SHA for ${subject_name} is not 40 hex characters")
    endif()
    execute_process(
      COMMAND "${RELEASE_GIT}" -C "${RELEASE_ROOT}" rev-parse "HEAD:subjects/${subject_name}"
      RESULT_VARIABLE gitlink_result
      OUTPUT_VARIABLE gitlink_sha
      ERROR_VARIABLE gitlink_error
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT gitlink_result EQUAL 0)
      message(FATAL_ERROR "Could not resolve gitlink for ${subject_name}: ${gitlink_error}")
    endif()
    if(NOT gitlink_sha STREQUAL subject_sha)
      message(FATAL_ERROR
        "SUBJECTS.lock mismatch for ${subject_name}: lock=${subject_sha}, gitlink=${gitlink_sha}")
    endif()
  endforeach()
endif()

message(STATUS "Release metadata verified for version ${release_version}")
