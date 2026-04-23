# When building source packages, remove the packaging install prefix so
# the packages are left with only the top-level directory.
if (CPACK_SOURCE_INSTALLED_DIRECTORIES)
  unset (CPACK_PACKAGING_INSTALL_PREFIX)
endif ()
