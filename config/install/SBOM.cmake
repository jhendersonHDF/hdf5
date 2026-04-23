# Create an SBOM file for a package and copy it back to the current build
# directory. Optionally, verify conformance to NTIA standards using the
# ntia-conformance-checker (https://github.com/spdx/ntia-conformance-checker)
# python package before copying back if python3 is available.

#---------------------------------------------------------------------------------------------------
# Create an SBOM document for the given source package archive
#
# - archive
#   Full path to the input source package archive
#
# - sbom_output_name
#   Full path of the output SBOM document
#
# - ret_var
#   Used for error checking by the caller (TRUE on success/FALSE if an error occurred)
#---------------------------------------------------------------------------------------------------
function (create_source_sbom archive sbom_output_name ret_var)
  set (${ret_var} FALSE)

  # Reject non-tgz and non-zip archives for now
  cmake_path (GET archive EXTENSION _base_filename_ext)
  if (NOT "${_base_filename_ext}" MATCHES ".tar.gz|.tgz|.zip")
    message (NOTICE "-- Skipping creation of SBOM for unsupported archive format ${_base_filename_ext}")
    return ()
  endif ()

  cmake_path (GET archive PARENT_PATH _base_dir)
  cmake_path (GET archive FILENAME _base_filename)
  set (_tmpdir_name "${_base_dir}/tmpdir")
  set (_tmp_archive_name "${_tmpdir_name}/${_base_filename}")

  # TODO: Remove
  message (STATUS "Creating temporary directory at ${_base_dir}")
  # TODO: Remove

  # Create temporary directory for processing archive
  execute_process (
    COMMAND "${CMAKE_COMMAND}" "-E" "make_directory" "${_tmpdir_name}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create temporary directory ${_tmpdir_name}: ${_err_msg}")
    return ()
  endif ()

  # TODO: Remove
  message (STATUS "Copying archive into temporary directory")
  # TODO: Remove

  execute_process (
    COMMAND "${CMAKE_COMMAND}" "-E" "copy" "${archive}" "${_tmpdir_name}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't copy archive into temporary directory ${_tmpdir_name}: ${_err_msg}")
    return ()
  endif ()

  # TODO: Remove
  message (STATUS "Extracting archive into temporary directory")
  # TODO: Remove

  execute_process (
    COMMAND "${CMAKE_COMMAND}" "-E" "tar" "x" "${_tmp_archive_name}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
    WORKING_DIRECTORY "${_tmpdir_name}"
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create temporary directory ${_tmpdir_name}: ${_err_msg}")
    return ()
  endif ()

  # TODO: need to programatically determine root directory name (hard-code for now)
  # TODO: check for errors from function
  set (_source_dir "${_tmpdir_name}/hdf5-2.2.0")
  create_source_sbom_helper ("${archive}" "${_source_dir}" "${sbom_output_name}" _helper_ret_var)

  # TODO: Remove
  message (STATUS "Removing temporary directory ${_tmpdir_name}")
  # TODO: Remove

  # TODO: Commented out for now for debugging until process works
#  execute_process (
#    COMMAND "${CMAKE_COMMAND}" "-E" "rm" "-R" "${_tmpdir_name}"
#    RESULT_VARIABLE _tmp_result_var
#    ERROR_VARIABLE _tmp_error_var
#  )
#  if (NOT _tmp_result_var EQUAL 0)
#    if ("${_tmp_error_var}" STREQUAL "")
#      set (_err_msg "${_tmp_result_var}")
#    else ()
#      set (_err_msg "${_tmp_error_var}")
#    endif ()
#
#    message (NOTICE "-- Couldn't remove temporary directory ${_tmpdir_name}: ${_err_msg}")
#    return ()
#  endif ()

  if (_helper_ret_var)
    set (${ret_var} TRUE)
  endif ()
endfunction ()

#---------------------------------------------------------------------------------------------------
# Performs the actual work of creating an SBOM document for a source package archive
#
# - archive
#   Full path to the input source package archive
#
# - dir
#   Full path to the temporary directory used for intermediate processing
#
# - sbom_output_name
#   Full path of the output SBOM document
#
# - ret_var
#   Used for error checking by the caller (TRUE on success/FALSE if an error occurred)
#---------------------------------------------------------------------------------------------------
function (create_source_sbom_helper archive dir sbom_output_name ret_var)
  set (${ret_var} FALSE)

  cmake_path (GET archive FILENAME _archive_name)
  cmake_path (GET sbom_output_name FILENAME _sbom_document_name)

  # Get timestamp for document's "created" field
  string (TIMESTAMP _created_at "%Y-%m-%dT%H:%M:%SZ" UTC)
  if (_created_at STREQUAL "")
    message (NOTICE "-- Couldn't obtain timestamp for final output file")
    return ()
  endif ()

  set (_spdx_version "SPDX-2.3")
  set (_spdx_data_license "CC0-1.0")
  set (_spdx_license_list_version "3.22")
  set (_package_name "HDF5")
  set (_package_license "BSD-3-Clause")
  set (_package_originator "Organization: The HDF Group")
  set (_package_supplier "Organization: The HDF Group")
  set (_spdx_package_id "SPDXRef-PACKAGE-hdf5")
  set (_document_name "HDF5 ${CPACK_PACKAGE_VERSION} SBOM")
  set (_document_namespace "https://www.hdfgroup.org/HDF5/${_sbom_document_name}")
  set (_download_location "https://support.hdfgroup.org/downloads/index.html") # TODO: Correct? Allow override?
  set (_creator_org "The HDF Group")
  set (_creator_tool "CMake ${CMAKE_VERSION}") 
  set (_creator_person "") # TODO: Use user ID or something consistent? 'creator_person=$(id -un)'
  set (_package_purl "pkg:generic/org.hdfgroup/hdf5@${CPACK_PACKAGE_VERSION}")

  # Process any overrides from user-set CMake variables
  if (DEFINED CPACK_HDF5_SBOM_CREATOR_ORGANIZATION)
    set (_creator_org "${CPACK_HDF5_SBOM_CREATOR_ORGANIZATION}")
  endif ()
  if (DEFINED CPACK_HDF5_SBOM_CREATOR_PERSON)
    set (_creator_person "${CPACK_HDF5_SBOM_CREATOR_PERSON}")
  endif ()
  if (DEFINED CPACK_HDF5_SBOM_CREATOR_TOOL)
    set (_creator_tool "${CPACK_HDF5_SBOM_CREATOR_TOOL}")
  endif ()
  if (DEFINED CPACK_HDF5_SBOM_DOCUMENT_NAMESPACE)
    set (_document_namespace "${CPACK_HDF5_SBOM_DOCUMENT_NAMESPACE}/${_sbom_document_name}")
  endif ()
  if (DEFINED CPACK_HDF5_SBOM_LICENSE_LIST_VERSION)
    set (_spdx_license_list_version "${CPACK_HDF5_SBOM_LICENSE_LIST_VERSION}")
  endif ()
  if (DEFINED CPACK_HDF5_SBOM_DATA_LICENSE)
    set (_spdx_data_license "${CPACK_HDF5_SBOM_DATA_LICENSE}")
  endif ()

  # Calculate sha256 for package archive file
  execute_process (
    COMMAND "${CMAKE_COMMAND}" "-E" "sha256sum" "${archive}"
    OUTPUT_VARIABLE _package_sha256
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create final output file: ${_err_msg}")
    return ()
  endif ()

  # Omit filename from sha256 result
  string (REPLACE " " ";" _package_sha256 "${_package_sha256}")
  list (GET _package_sha256 0 _package_sha256)

  # Create temporary files used in assembling final output file
  set (_files_jsonl "${dir}/files.jsonl")
  set (_relationships_jsonl "${dir}/relationships.jsonl")
  set (_package_json "${dir}/package.json")
  set (_sha1_list "${dir}/file-sha1s.txt")
  execute_process (
    COMMAND "${CMAKE_COMMAND}" "-E" "touch" "${_files_jsonl}" "${_relationships_jsonl}" "${_package_json}" "${_sha1_list}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create final output file: ${_err_msg}")
    return ()
  endif ()

  # TODO: Create 'files' file and 'relationships' file

  # TODO: Calculate package verification code
  # Calculate the SPDX package verification code as the SHA1 of the sorted
  # file SHA1 values 
  set (_package_verification_code "")

  # Create 'package' JSON file
  execute_process (
    COMMAND "${JQ_EXECUTABLE}" "-cn"
            "--arg" "spdxid" "${_spdx_package_id}"
            "--arg" "checksum" "${_package_sha256}"
            "--arg" "download_location" "${_download_location}"
            "--arg" "purl" "${_package_purl}"
            "--arg" "license" "${_package_license}"
            "--arg" "name" "${_package_name}"
            "--arg" "originator" "${_package_originator}"
            "--arg" "package_file_name" "${_archive_name}"
            "--arg" "verification_code" "${_package_verification_code}"
            "--arg" "supplier" "${_package_supplier}"
            "--arg" "version" "${CPACK_PACKAGE_VERSION}"
            "\
{\
 SPDXID: $spdxid,\
 checksums: [\
   {algorithm: \"SHA256\", checksumValue: $checksum}\
 ],\
 downloadLocation: $download_location,\
 externalRefs: [\
   {\
     referenceCategory: \"PACKAGE_MANAGER\",\
     referenceLocator: $purl,\
     referenceType: \"purl\"\
   }\
 ],\
 filesAnalyzed: true,\
 licenseConcluded: $license,\
 name: $name,\
 originator: $originator,\
 packageFileName: $package_file_name,\
 packageVerificationCode: {\
   packageVerificationCodeValue: $verification_code\
 },\
 primaryPackagePurpose: \"SOURCE\",\
 supplier: $supplier,\
 versionInfo: $version\
}"
    OUTPUT_FILE "${_package_json}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create 'packages' JSON file: ${_err_msg}")
    return ()
  endif ()

  # TODO: Append to 'relationships' file

  # Assemble final output file
  execute_process (
    COMMAND "${JQ_EXECUTABLE}" "-n"
            "--arg" "spdx_version" "${_spdx_version}"
            "--arg" "document_spdx_id" "SPDXRef-DOCUMENT"
            "--arg" "document_name" "${_document_name}"
            "--arg" "document_namespace" "${_document_namespace}"
            "--arg" "created" "${_created_at}"
            "--arg" "creator_org" "Organization: ${_creator_org}"
            "--arg" "creator_person" "Person: ${_creator_person}"
            "--arg" "creator_tool" "Tool: ${_creator_tool}"
            "--arg" "license_list_version" "${_spdx_license_list_version}"
            "--arg" "data_license" "${_spdx_data_license}"
            "--slurpfile" "packages" "${_package_json}"
            "--slurpfile" "files" "${_files_jsonl}"
            "--slurpfile" "relationships" "${_relationships_jsonl}"
            "\
{\
 spdxVersion: $spdx_version,\
 SPDXID: $document_spdx_id,\
 name: $document_name,\
 documentNamespace: $document_namespace,\
 creationInfo: {\
   created: $created,\
   creators: [\
     $creator_org,\
     $creator_person,\
     $creator_tool\
   ],\
   licenseListVersion: $license_list_version\
 },\
 dataLicense: $data_license,\
 packages: $packages,\
 files: $files,\
 relationships: $relationships\
}"
    OUTPUT_FILE "${sbom_output_name}"
    RESULT_VARIABLE _tmp_result_var
    ERROR_VARIABLE _tmp_error_var
    COMMAND_ECHO STDOUT # TODO: Debugging
  )
  if (NOT _tmp_result_var EQUAL 0)
    if ("${_tmp_error_var}" STREQUAL "")
      set (_err_msg "${_tmp_result_var}")
    else ()
      set (_err_msg "${_tmp_error_var}")
    endif ()

    message (NOTICE "-- Couldn't create final output file: ${_err_msg}")
    return ()
  endif ()

  set (${ret_var} TRUE)

  message (STATUS "Wrote SPDX JSON SBOM to ${sbom_output_name}")
endfunction ()

#---------------------------------------------------------------------------------------------------
# Given an existing SBOM file, use the ntia-conformance-checker python package
# (https://github.com/spdx/ntia-conformance-checker) to verify that the SBOM document conforms to
# the NTIA's "minimum elements" specification requirements.
#
# - sbom_file
#   Full path to the input SBOM document
#---------------------------------------------------------------------------------------------------
function (verify_sbom sbom_file)
  # Find Python3 3.5 or greater (require support for virtual environments
  # with pip installed by default + a compatibility change to the command
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

#---------------------------------------------------------------------------------------------------
# Main logic
#---------------------------------------------------------------------------------------------------
if (NOT CPACK_OUTPUT_FILE_PATH)
  message (FATAL_ERROR "This script relies on CPACK_OUTPUT_FILE_PATH being set")
endif ()

find_program (JQ_EXECUTABLE jq)
if (NOT JQ_EXECUTABLE)
  message (FATAL_ERROR "-- Creating SBOM files requires the 'jq' program")
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
  unset (_tmp_ret)

  message (STATUS "Creating SBOM file for ${_cpack_package_type} package ${_cpack_package}")

  if (_cpack_package_type STREQUAL "source")
    create_source_sbom ("${_cpack_package}" "${_cpack_package}.spdx.json" _tmp_ret)
  else ()
    # TODO
    message (FATAL_ERROR "Binary packages not supported yet")
    create_binary_sbom ("${_cpack_package}")
  endif ()

  if (_tmp_ret)
    verify_sbom ("${_cpack_package}.spdx.json")

    # Copy the SBOM file back to the build directory
    execute_process (
      COMMAND "${CMAKE_COMMAND}" -E copy "${_cpack_package}.spdx.json" "${_output_file_path}"
    )
  endif ()
endforeach ()
