# Instructions for the installation of HDF5 Software

This file provides instructions for installing the HDF5 software.

For additional help with installing, questions can be posted to the HDF Forum or sent to the HDF
Helpdesk:

* **HDF Forum:** [https://forum.hdfgroup.org/](https://forum.hdfgroup.org/)
* **HDF Helpdesk:** [https://help.hdfgroup.org/](https://help.hdfgroup.org/)

Questions may also be sent to <help@hdfgroup.org>.

---

## Table of Contents

* [1. Prerequisites](#1-prerequisites)
* [2. Obtaining HDF5](#2-obtaining-hdf5)
* [3. Building HDF5 with CMake command mode](#3-building-hdf5-with-cmake-command-mode)
  * [Quick build steps](#quick-build-steps)
  * [Customizing the build](#customizing-the-build)
  * [HDF5 compression and data filters](#hdf5-compression-and-data-filters)
  * [Parallel HDF5](#parallel-hdf5)
  * [Thread-safe HDF5](#thread-safe-hdf5)
  * [Read-only S3 VFD](#read-only-s3-vfd)
* [4. Alternative methods of building HDF5 with CMake](#4-alternative-methods-of-building-hdf5-with-cmake)
* [5. Cross-compiling considerations](#5-cross-compiling-considerations)
* [6. Packaging HDF5 with CPack](#6-packaging-hdf5-with-cpack)
  * [WiX Toolset](#wix-toolset)
  * [Nullsoft Scriptable Install System](#nullsoft-scriptable-install-system)
* [7. Additional information](#7-additional-information)
  * [Backward compatibility](#backward-compatibility)
  * [User-defined CMake Options for HDF5 Libraries](#user-defined-cmake-options-for-hdf5-libraries)
  * [CDash server](#cdash-server)

---

## 1. Prerequisites

**Required software**

- A C compiler supporting the C11 standard

- [CMake](https://cmake.org/)
  - HDF5 requires CMake 3.26 or newer for building

- **Windows-only**
  - HDF5 requires [Visual Studio](https://visualstudio.microsoft.com/) for building. Visual Studio
    17 2022 and newer are officially supported; older versions may work as well.

**Optional software**

The following software may need to be available during the build process if specific HDF5 features
are enabled:

- **A Fortran compiler** -
  HDF5's Fortran language wrappers and example programs require a Fortran compiler supporting at
  least the 2003 standard.

- **A C++ compiler** -
  HDF5's C++ language wrappers and example programs require a C++ compiler supporting at least the
  C++11 standard.

- **A Java JDK** -
  HDF5's Java language wrappers and example programs require a Java JDK version 11 or newer for
  building.

- **A Python3 interpreter** -
  HDF5's Python example programs require a Python3 interpreter.

- **MPI** -
  HDF5 support for parallel I/O requires MPI libraries and development files to be available on the
  system. The MPI implementation must support at least the MPI 3.0 standard.

- **Threading support** -
  Threading support (C11 threads, Win32 threads or pthreads) is required to build HDF5 with thread
  safety or concurrency support. Threading support is also required to build HDF5's Subfiling VFD.

- **zlib-ng** OR **zlib** -
  HDF5 includes a predefined filter which uses the [Deflate](https://en.wikipedia.org/wiki/Deflate)
  algorithm for compressing data. This filter requires either zlib-ng or zlib libraries and
  development files to be available on the system.

- **libaec** -
  HDF5 includes a predefined filter which uses the extended-Rice lossless algorithm for compressing
  data. This filter requires libaec libraries and development files to be available on the system.

- **various other libraries for data filter plugins** -
  The [HDF5 plugins project](https://github.com/HDFGroup/hdf5_plugins) provides various data filter
  plugins which each have their own requirements for building.

- **OpenSSL** -
  HDF5's "digitally signed plugins" feature requires OpenSSL 1.1.0 or newer.

- **aws-c-s3** -
  HDF5's ROS3 VFD requires the [aws-c-s3](https://github.com/awslabs/aws-c-s3) library.

- **libhdfs** -
  HDF5's HDFS VFD requires libhdfs.

- **Doxygen** -
  Doxygen is required for building HDF5's doxygen documentation.

## 2. Obtaining HDF5

The HDF Group provides source code and pre-compiled binaries from the HDF5 GitHub releases page:
[https://github.com/HDFGroup/hdf5/releases](https://github.com/HDFGroup/hdf5/releases). The source
code for the latest HDF5 release may be obtained from that page or directly from

- [https://github.com/HDFGroup/hdf5/releases/latest/download/hdf5.tar.gz](hdf5.tar.gz)
- [https://github.com/HDFGroup/hdf5/releases/latest/download/hdf5.zip](hdf5.zip)

The `hdf5.tar.gz` or `hdf5.zip` archive should be extracted somewhere on disk in preparation for
building from source.

The latest HDF5 development source code may be obtained with git:

```Shell
git clone https://github.com/HDFGroup/hdf5.git
```

## 3. Building HDF5 with CMake command mode

The following steps assume that the HDF5 source code has already been obtained by one of the methods
listed above and exists in a directory referred to below as `<source_dir>`. The path where HDF5 will
be installed to is referred to as `<install_dir>`. Additional CMake options which can be specified to
customize the HDF5 configuration are listed as `<options>`.

> [!NOTE]
> HDF5 does not support direct in-source builds, so the steps below create a build directory in the
> source tree to contain the build files generated by CMake.

### Quick build steps

**Linux / MacOS / non-Windows / If using a single-configuration generator** - from a terminal (bash, Terminal, etc.):

```Shell
cd <source_dir>
# Configure HDF5 with any options specified
cmake -B build/ -S . -DCMAKE_INSTALL_PREFIX=<install_dir> <options>
# Build HDF5 in parallel
cmake --build build/ --parallel
# Optional, test HDF5 before installing
ctest --test-dir build/
# Install HDF5
cmake --install build/
```

**Windows / If using a multi-configuration generator** - from a terminal (PowerShell, Git Bash, etc.):

```Shell
cd <source_dir>
# Configure HDF5 with any options specified
cmake -B build/ -S . -DCMAKE_INSTALL_PREFIX=<install_dir> <options>
# Build HDF5 in parallel
cmake --build build/ --config Release --parallel
# Optional, test HDF5 before installing
ctest --test-dir build/ --build-config Release
# Install HDF5
cmake --install build/ --config Release
```

This creates a build directory called `build`, configures HDF5 for a Release build with an
installation directory of `<install_dir>` and any other options specified for `<options>`, builds HDF5
in parallel using a default number of concurrent processes and finally installs HDF5 after the build
is finished.

Note that the steps for building on Windows or when using a multi-configuration CMake generator are
nearly identical to the steps for other platforms or when using a single-configuration CMake
generator. The only difference is that the configuration type to build, test and install is given
with the `--config Release` / `--build-config Release` options. Since multi-configuration CMake
generators, such as the Visual Studio generators for Windows, can generate build files for multiple
configuration types simultaneously, the desired configuration type is specified explicitly. For more
information on CMake generators, see [cmake-generators](https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html).

As shown above, the HDF5 build can be tested with `ctest` before installation. The `--parallel` option
may be passed to the `ctest` command to run tests in parallel. For more information on HDF5 testing,
see [TESTING.md](./TESTING.md).

See [USING_HDF5_CMake.md](./USING_HDF5_CMake.md) for information on using the HDF5 installation
after it has been built.

### Customizing the build

TODO: reference `<options>` in some way
HDF5 has many different features and settings that can be enabled by specifying them on the command
line in the form `-D<OPTION>=<VALUE>` when configuring the library. For example, the following
command, modified slightly from above, configures HDF5 with Fortran support:

```Shell
cmake -B build/ -S . -DCMAKE_INSTALL_PREFIX=<install_dir> -DHDF5_BUILD_FORTRAN=ON
```

For a list of the different options that HDF5 can be configured with, see
[INSTALL_options](./INSTALL_options.md).

See CMake's [User Interaction Guide](https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#guide:User%20Interaction%20Guide) for general information on building with CMake.

TODO: mention presets for setting configuration options

### HDF5 compression and data filters

For specific notes on enabling support for compression and other data filters in HDF5, see
[INSTALL_filters.md](./INSTALL_filters.md).

### Parallel HDF5

HDF5 can be configured to use MPI and MPI-IO for parallelism on a distributed multi-processor
system by specifying

    -DHDF5_ENABLE_PARALLEL=ON

when configuring the library. Refer to [INSTALL_HPC.md](./INSTALL_HPC.md) for detailed information.

### Thread-safe HDF5

HDF5 can be configured to be thread-safe by specifying

    -DHDF5_ENABLE_THREADSAFE=ON

when configuring the library. For further information, see
https://support.hdfgroup.org/documentation/hdf5/latest/thread-safe-lib.html.

### Read-only S3 VFD

See [INSTALL_S3.md](./INSTALL_S3.md) for instructions on how to configure HDF5 with support for
reading files from S3 using the read-only S3 VFD.

## 4. Alternative methods of building HDF5 with CMake

The following are alternative ways of building HDF5 with CMake for different platforms and use
cases:

| File | Description |
|------|-------------|
| [INSTALL_presets.md](./INSTALL_presets.md) | Building HDF5 with CMake presets |
| [INSTALL_GUI_TUI.md](./INSTALL_GUI_TUI.md) | Building HDF5 with `cmake-gui` or `ccmake` |
| [INSTALL_VS.md](./INSTALL_VS.md) | Building HDF5 with Visual Studio (Windows) |
| [INSTALL_script.md](./INSTALL_script.md) | Building HDF5 with a `ctest` script (HDF5 maintainers) |

## 5. Cross-compiling considerations

Cross-compiling has several consequences for CMake:

* CMake cannot automatically detect the target platform.
* CMake cannot find libraries and headers in the default system directories.
* Executables built during cross compiling cannot be executed.

Cross-compiling support means that CMake separates information about the
build platform and target platform and gives the user mechanisms to solve
cross-compiling issues without additional requirements such as running
virtual machines, etc.

CMake uses a toolchain of utilities to compile, link libraries, create
archives, and other tasks to drive the build. The toolchain utilities
available are determined by the languages enabled.

CMake stores info about the current toolchain in the variables `CMAKE_C_COMPILER`, `CMAKE_CXX_COMPILER`.
They contain paths to the C and C++ compilers respectively. This is usually
enough on desktop platforms. In the case of embedded systems, a custom
linker and assembler setting may be needed. In more complex projects
you may need to additionally specify binaries to other parts of the toolchain
(size, ranlib, objcopy…). All these tools should be set in the corresponding
variables: `CMAKE_AR`, `CMAKE_ASM_COMPILER`, `CMAKE_LINKER`, `CMAKE_OBJCOPY`, and `CMAKE_RANLIB`.

As for the host and target operating systems, CMake stores their names in the
following variables:

* `CMAKE_HOST_SYSTEM_NAME` – name of the platform, on which CMake is running (host platform). On major
operating systems this is set to the `Linux`, `Windows` or `Darwin` (macOS) value.

* `CMAKE_SYSTEM_NAME` – name of the platform, for which we are building (target platform). By default,
this value is the same as `CMAKE_HOST_SYSTEM_NAME`, which means that we are building for the local
platform (no cross-compilation).

Put the toolchain variables into a separate file (e.g. `<toolchain_name>.cmake`) and set the
`CMAKE_TOOLCHAIN_FILE` variable to the path of that file. If `cmake` is invoked with the command line
parameter `--toolchain path/to/file` or `-DCMAKE_TOOLCHAIN_FILE=path/to/file` the file will be loaded
early to set values for the compilers. The `CMAKE_CROSSCOMPILING` variable is set to true when CMake
is cross-compiling.

### Structure of the toolchain file

In fact, the toolchain file doesn’t have any structure. You can put anything you
want there. But the best practice is to define at least these settings:
path to the toolchain binaries (C compiler, C++ compiler, linker, etc.),
name of the target platform (and optionally target processor architecture),
required compilation and linking flags on that particular platform, and
toolchain sysroot settings.

It is recommended that you set the `CMAKE_FIND_ROOT_PATH` variable to a path where
you have an exact copy of the root filesystem you have on your target device (with
libraries and binaries pre-compiled for the target processor).

References:

https://cmake.org/cmake/help/latest/manual/cmake-toolchains.7.html<br/>
https://gitlab.com/embeddedlinux/libs/platform<br/>
https://discourse.cmake.org/t/cross-compile-for-aarch64-on-ubuntu/2161/10<br/>
https://stackoverflow.com/questions/54539682/how-to-set-up-cmake-to-cross-compile-with-clang-for-arm-embedded-on-windows?rq=1<br/>
https://developer.android.com/ndk/guides/cmake<br/>

## 6. Packaging HDF5 with CPack

This section covers using CMake's [cpack](https://cmake.org/cmake/help/latest/manual/cpack.1.html#manual:cpack(1))
functionality to create binary and source packages for HDF5.

> [!NOTE]
> Windows developers should install NSIS or WiX to create an install image with CPack.

TODO
    * On **Windows (with WiX installed)**, execute `hdf5-2.X.Y-win-vs2022_cl.msi` or
      `hdf5-2.X.Y-win-vs2022_intel.msi`. By default this program will install the hdf5
      library into the `C:\Program Files` directory and will create the
      following directory structure:

        ```
        HDF_Group
        --HDF5
        ----2.X.Y
        ------bin
        ------include
        ------lib
        --------plugin
        ------cmake
        ```

    * On **Linux**, change to the install destination directory (create if doesn't exist) and execute: `<source_dir>/build/HDF5-2.X.Y-Linux.sh`. After accepting the license, the script will prompt:

        ```
        By default the HDF5 will be installed in:
        "<current directory>/HDF5-2.X.Y-Linux"
        Do you want to include the subdirectory HDF5-2.X.Y-Linux?
        Saying no will install in: "<current directory>" [Yn]:
        ```

        Note that the script will create the following directory structure
        relative to the install point:

        ```
        HDF_Group
        --HDF5
        ----2.X.Y
        ------bin
        ------include
        ------lib
        --------plugin
        ------share
        ```

    * On **macOS**, locate `HDF5-2.X.Y-Darwin.dmg` file in the build folder. Click on the  file to proceed with installation. After accepting the license, there will be a folder with the following structure:

        ```
        HDF_Group
        --HDF5
        ----2.X.Y
        ------bin
        ------include
        ------lib
        --------plugin
        ------share
        ```

### WiX Toolset

WiX--the Windows Installer XML toolset--lets developers create installers for Windows Installer, the
Windows installation engine. See http://wixtoolset.org.

### Nullsoft Scriptable Install System

The Nullsoft Scriptable Install System (NSIS) is an open source installation
system. It was created by the WinAmp authors to distribute that application,
but it is now a general-purpose system which anyone might use. NSIS installers
recognize `/S` for silent installation and `/D=dir` to specify the
output directory, which is where the program will be installed. These
options are case-sensitive, so be sure to type them in upper case.

## 7. Additional information

### Backward compatibility

The HDF5 library can be configured to set all versioned API functions to their respective versions
available in a particular release of the library by specifying the `HDF5_DEFAULT_API_VERSION` CMake
option when configuring HDF5. This allows existing code to be compiled with a newer version of the
library without requiring immediate changes to the application source code, as long as the application
does not use any symbols removed by a newer release of HDF5. For more information on this mechanism,
see <a href="https://support.hdfgroup.org/documentation/hdf5/latest/api-compat-macros.html">API Compatibility Macros</a>.

### User-defined CMake Options for HDF5 Libraries

HDF5's CMake logic includes support for user-defined CMake macros and options. The file
`UserMacros.cmake` has an example of the technique. Replace the template code with your macro in the
`UserMacros.cmake` file. Then enable the option in the CMake configuration, build and test process.

### CDash server

Build and test results can be submitted to our HDF5 CDash server. The CDash server for community
submissions of HDF5 is at https://my.cdash.org/index.php?project=HDF5. We ask that all submissions
include the configuration information and contact information in the CTest Notes Files upload step.
See the current reports on CDash for examples.

Please follow the convention that "NIGHTLY" submissions maintain the same configuration every time.
"EXPERIMENTAL" submissions can be expected to be different for each submission.
