# Create an SBOM file for a package and copy it back to the current build
# directory. Optionally, verify conformance to NTIA standards using the
# ntia-conformance-checker python package before copying back if python3
# is available.

function (create_source_sbom archive)
  # TODO: DEMO
  if (TRUE)
    execute_process (
      COMMAND "/home/jhenderson/git/hdf5/config/install/sbommer.sh" "${archive}" "${archive}.spdx.json"
      RESULT_VARIABLE _hdf5_sbom_verify_result
      ERROR_VARIABLE _hdf5_sbom_verify_error
    )
    if (NOT _hdf5_sbom_verify_result EQUAL 0)
      if ("${_hdf5_sbom_verify_error}" STREQUAL "")
        set (_err_msg "${_hdf5_sbom_verify_result}")
      else ()
        set (_err_msg "${_hdf5_sbom_verify_error}")
      endif ()

      message (NOTICE "-- Couldn't create SBOM file: ${_err_msg}")
      return ()
    endif ()
    return ()
  endif ()

  # TODO
  message (FATAL_ERROR "Not Yet")

  

endfunction ()

function (verify_sbom sbom_file)
  # Find Python3 3.5 or greater (require support for virtual environments
  # with pip installed by default + a compatiblity change to the command
  # for creating them)
  find_package (Python3 QUIET COMPONENTS Interpreter)
  if (NOT Python3_FOUND OR (Python3_VERSION VERSION_LESS 3.5))
    message (NOTICE "-- Not verifying SBOM file ${sbom_file}; python3 not found or version is < 3.5")
    return ()
  endif ()

  set (_system_python_exec "${Python3_EXECUTABLE}")

  # Get the base directory from CPACK_PACKAGE_FILES
  list (GET CPACK_PACKAGE_FILES 0 _package_files_0)
  if ("${_package_files_0}" STREQUAL "")
    message (NOTICE "-- Couldn't get base directory for package ${_package_files_0}")
    return ()
  endif ()

  cmake_path (GET _package_files_0 PARENT_PATH _package_files_root_path)
  set (VENV_PATH "${_package_files_root_path}/sbom_verify_venv")

  # Create a virtual environment for temporary use
  message (STATUS "Setting up python virtual environment")
  execute_process (
    COMMAND "${Python3_EXECUTABLE}" -m venv "${VENV_PATH}"
    RESULT_VARIABLE _hdf5_sbom_verify_result
    ERROR_VARIABLE _hdf5_sbom_verify_error
  )
  if (NOT _hdf5_sbom_verify_result EQUAL 0)
    if ("${_hdf5_sbom_verify_error}" STREQUAL "")
      set (_err_msg "${_hdf5_sbom_verify_result}")
    else ()
      set (_err_msg "${_hdf5_sbom_verify_error}")
    endif ()

    message (NOTICE "-- Couldn't create python virtual environment: ${_err_msg}")
    return ()
  endif ()

  # Re-find the python3 package, preferring to use our virtual environment
  set (ENV{VIRTUAL_ENV} "${VENV_PATH}")
  set (Python3_FIND_VIRTUALENV ONLY)
  unset(Python3_EXECUTABLE)
  find_package (Python3 QUIET COMPONENTS Interpreter)
  if (NOT Python3_FOUND OR (Python3_VERSION VERSION_LESS 3.5) OR ("${Python3_EXECUTABLE}" STREQUAL "${_system_python_exec}"))
    message (NOTICE "-- Couldn't find python3 virtual environment")

    # Remove the virtual environment
    execute_process (
      COMMAND "${CMAKE_COMMAND}" -E rm -Rf "${VENV_PATH}"
    )

    return ()
  endif ()

  # Install the ntia-conformance-checker python package for SBOM NTIA
  # conformance checking
  message (STATUS "Installing ntia-conformance-checker python package")
  execute_process (
    COMMAND "${Python3_EXECUTABLE}" -m pip install ntia-conformance-checker -qqq
    RESULT_VARIABLE _hdf5_sbom_verify_result
    ERROR_VARIABLE _hdf5_sbom_verify_error
  )
  if (NOT _hdf5_sbom_verify_result EQUAL 0)
    if ("${_hdf5_sbom_verify_error}" STREQUAL "")
      set (_err_msg "${_hdf5_sbom_verify_result}")
    else ()
      set (_err_msg "${_hdf5_sbom_verify_error}")
    endif ()

    message (NOTICE "-- Couldn't install ntia-conformance-checker python package: ${_err_msg}")

    # Remove the virtual environment
    execute_process (
      COMMAND "${CMAKE_COMMAND}" -E rm -Rf "${VENV_PATH}"
    )

    return ()
  endif ()

  # Perform conformance checking
  message (STATUS "Verifying SBOM file ${sbom_file}")
  if (NOT WIN32)
    set (_sbom_checker "${VENV_PATH}/bin/sbomcheck")
  else ()
    set (_sbom_checker "${VENV_PATH}/Scripts/sbomcheck")
  endif ()
  execute_process (
    COMMAND "${_sbom_checker}" -s spdx2 -c ntia -r print "${sbom_file}"
    RESULT_VARIABLE _hdf5_sbom_verify_result
    ERROR_VARIABLE _hdf5_sbom_verify_error
  )
  if (NOT _hdf5_sbom_verify_result EQUAL 0)
    if ("${_hdf5_sbom_verify_error}" STREQUAL "")
      set (_err_msg "${_hdf5_sbom_verify_result}")
    else ()
      set (_err_msg "${_hdf5_sbom_verify_error}")
    endif ()

    message (FATAL_ERROR "-- Verification failed for SBOM file ${sbom_file}: ${_err_msg}")
  endif ()

  # Remove the virtual environment
  message (STATUS "Removing python virtual environment")
  execute_process (
    COMMAND "${CMAKE_COMMAND}" -E rm -Rf "${VENV_PATH}"
  )
endfunction ()

#-------------------------------------------------------------------------------
# Main logic
#-------------------------------------------------------------------------------
if (NOT CPACK_OUTPUT_FILE_PATH)
  message (FATAL_ERROR "This script relies on CPACK_OUTPUT_FILE_PATH being set")
endif ()

find_program (JQ_EXECUTABLE jq)
if (NOT JQ_EXECUTABLE)
  message (FATAL_ERROR "-- Creating SBOM files requires the jq program")
  return ()
endif ()

# Get the output directory to copy SBOM files back to
cmake_path (GET CPACK_OUTPUT_FILE_PATH PARENT_PATH _output_file_path)

if (CPACK_SOURCE_INSTALLED_DIRECTORIES)
  set (_cpack_package_type "source")
else ()
  set (_cpack_package_type "binary")
endif ()

foreach (_cpack_package "${CPACK_PACKAGE_FILES}")
  message (STATUS "Creating SBOM file for ${_cpack_package_type} package ${_cpack_package}")

  if (_cpack_package_type STREQUAL "source")
    create_source_sbom ("${_cpack_package}")
  else ()
    # TODO
    message (FATAL_ERROR "Binary packages not supported yet")
    create_binary_sbom ("${_cpack_package}")
  endif ()

  verify_sbom ("${_cpack_package}.spdx.json")

  # Copy the SBOM file back to the build directory
  execute_process (
    COMMAND "${CMAKE_COMMAND}" -E copy "${_cpack_package}.spdx.json" "${_output_file_path}"
  )
endforeach ()
