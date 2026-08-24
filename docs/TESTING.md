# HDF5 Testing Information

This file provides information about testing HDF5.

---

## Table of Contents

* [Common ctest options](#common-ctest-options)
* [TestExpress](#testexpress)
* [Testing timeouts](#testing-timeouts)
* [Using CMake Regex Options for Testing](#using-cmake-regex-options-for-testing)
* [Java FFM Testing](#java-ffm-testing)
  * [FFM Test Organization](#ffm-test-organization)
  * [Running FFM Tests](#running-ffm-tests)

---

## Common ctest options

Running all of HDF5's tests can be performed with the simple command:

```Shell
ctest --test-dir <build_dir>
```

where `<build_dir>` is the path to an HDF5 build directory which contains tests to run. This command
runs testing with a single process, printing out the test name and a PASS/FAIL status for each test.
Some helpful additional options that can be provided to `ctest` include:

- `-V` - Enables verbose output from tests
- `--output-on-failure` - Enables verbose output from tests that fail; useful for debugging failures
  while generally keeping overall test output minimal
- `--parallel` - Runs testing in parallel; a specific number of processes to use for running tests
  can be given or can be omitted to use as many processes as available (with CMake 3.28 and earlier
  it is required to specify a number of processes to use)
- `-R` / `-E` - Run or exclude specific tests based on a regular expression

See [ctest](https://cmake.org/cmake/help/latest/manual/ctest.1.html) for more information on using
`ctest` for testing.

## TestExpress

HDF5 testing includes a mechanism for expediting testing, referred to as "TestExpress". Tests which
support this functionality typically check the current TestExpress level and try to adjust the
amount of testing performed in order to fit within a specific test runtime:

    0 - Tests should take as long as necessary
    1 - Tests should take no more than 30 minutes
    2 - Tests should take no more than 10 minutes
    3 - Tests should take no more than 1 minute

HDF5's CMake build system sets the TestExpress level to 3 by default; more exhaustive testing is
regularly performed in CI. The TestExpress level can be changed when configuring HDF5 (by setting
the `HDF_TEST_EXPRESS` option), or can be set to a different value at test runtime by setting the
`HDF5TestExpress` environment variable.

## Testing timeouts

HDF5 testing uses a default value of 1200 seconds before a test program is timed out. To adjust this
value, pass the `--timeout` option when running `ctest`:

```Shell
# Run testing with timeout of 2 minutes (applies to each individual test program)
ctest --timeout 120 .
```

## Using CMake Regex Options for Testing

`ctest` provides the ability to use regular expressions to control which tests are executed. HDF5's
CMake build logic adds prefixes to most test names which can be used to select specific groups of
tests to run or to exclude. These prefixes are:

SERIAL
- CPP for C++ tests.
- HL_CPP for high-level C++ tests.
- FORTRAN for Fortran tests.
- HL_FORTRAN for high-level Fortran tests.
- HL for high-level tests.
- JUnit for Java tests.
- H5WATCH for tests that use the h5watch SWMR program.
- SWMR for tests that use the SWMR feature.
- h5_api for the API tests.
- VOL for tests that use the VOL feature.
- VFD for tests that use the VFD feature.
- H5TEST for the library tests.
- EX for examples
- H5CLEAR for the h5clear tool tests.
- H5COPY for the h5copy tool tests.
- H5DIFF for the h5diff tool tests.
- H5DUMP for the h5dump tool tests.
- H5FC for the h5fc tool tests.
- H5IMPORT for the h5import tool tests.
- H5JAM for the h5jam tool tests.
- H5LS for the h5ls tool tests.
- H5MKGRP for the h5mkgrp tool tests.
- H5REPACK for the h5repack tool tests.
- H5REPART for the h5repart tool tests.
- H5STAT for the h5stat tool tests.
- PERFORM for performance tests.

PARALLEL
- MPI_TEST for parallel tests.
- MPI_TEST_FORT for just parallel Fortran tests.

To run tests with a specific prefix, use the `--tests-regex` (or `-R`) option with `ctest`. For
example, to run only the C++ tests, use:

```bash
ctest --test-dir <build_dir> --tests-regex "CPP"
```

To select multiple groups of tests with different prefixes, use the `|` operator to separate the
prefixes. For example, to run both the C++ and Fortran tests, use:

```bash
ctest --test-dir <build_dir> --tests-regex "CPP|FORTRAN"
```

To exclude tests with a specific prefix, use the `--exclude-regex` (or `-E`) option with `ctest`. For
example, to run all tests except for the C++ tests, use:

```bash
ctest --test-dir <build_dir> --exclude-regex "CPP"
```

For more information on using regular expressions with `ctest`, see the
[ctest documentation](https://cmake.org/cmake/help/latest/manual/ctest.1.html#regular-expressions).

## Java FFM Testing

HDF5 2.0.0 and newer includes testing for the experimental Java Foreign Function & Memory (FFM) API
wrappers. HDF5 must be configured with FFM support, which requires Java 25 or newer, for these tests
to be available.

### FFM Test Organization

Tests are organized by HDF5 module in `java/jtest/`:

Module Test Files:
  - TestH5ffm.java     - General library operations
                         H5open, H5close, memory management, version info
  - TestH5Affm.java    - Attribute operations
  - TestH5Dffm.java    - Dataset operations
  - TestH5Effm.java    - Error handling
  - TestH5Fffm.java    - File operations
  - TestH5FDffm.java   - File drivers
  - TestH5Gffm.java    - Group operations
  - TestH5Iffm.java    - Identifier management
  - TestH5Lffm.java    - Link operations
  - TestH5Offm.java    - Object operations
  - TestH5Pffm.java    - Property lists
  - TestH5PLffm.java   - Plugin management
  - TestH5Rffm.java    - References
  - TestH5Sffm.java    - Dataspace operations
  - TestH5Tffm.java    - Datatype operations
  - TestH5VLffm.java   - VOL connector
  - TestH5Zffm.java    - Filter operations

### Running FFM Tests

To run all FFM tests:
```bash
ctest -R "JUnitFFM" -V
```

To run specific module tests:
```bash
ctest -R "JUnit-TestH5Affm" -V    # Attributes
ctest -R "JUnit-TestH5Pffm" -V    # Properties
ctest -R "JUnit-TestH5Tffm" -V    # Datatypes
ctest -R "JUnit-TestH5Sffm" -V    # Dataspaces
```
