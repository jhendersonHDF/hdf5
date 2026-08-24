# Building HDF5 with a ctest script

This file provides instructions for building HDF5 with a `ctest` script. This method is intended for
HDF5 maintainers and automated testing; most HDF5 users should choose a different installation
method provided in [INSTALL.md](./INSTALL.md).

These instructions assume that a directory, referred to below as `<source_dir>`, will be created to
contain the HDF5 source code and all the files needed to build it. This installation method will use
the default settings in the `config/cmake/cacheinit.cmake` file.

## Individual files needed as mentioned in this document

HDF5 source code:
[hdf5.tar.gz](https://github.com/HDFGroup/hdf5/releases/latest)

Script files from https://github.com/HDFGroup/hdf5/blob/develop/config/cmake/scripts:
* `CTestScript.cmake`  -- CMake build script
* `HDF5config.cmake`   -- CMake configuration script
* `HDF5options.cmake`  -- CMake configuration options script

HDF5 filter plugins:
* **Plugins:** [hdf5_plugins-<version>.tar.gz](https://github.com/HDFGroup/hdf5_plugins/releases/latest)

External libraries:
* **ZLIB:** [zlib-1.3.2.tar.gz](https://github.com/madler/zlib/releases/download/v1.3.2/zlib-1.3.2.tar.gz)
* **ZLIBNG:** [2.3.3.tar.gz](https://github.com/zlib-ng/zlib-ng/archive/refs/tags/2.3.3.tar.gz)
* **LIBAEC:** [libaec-1.1.6.tar.gz](https://github.com/MathisRosenhauer/libaec/releases/download/v1.1.6/libaec-1.1.6.tar.gz)

## Build Scripts for Windows or Linux

To build HDF5 with the SZIP and ZLIB external libraries you will need to:

1. Download/copy the individual files mentioned above to `<source_dir>`. Do not uncompress the
   `tar.gz` files.

2. Change to the directory `<source_dir>`. The `CTestScript.cmake` file should not be modified.

3. Edit the platform configuration file, `HDF5options.cmake`, if you want to change the default build
environment. The file is a compilation of the most used options and by commenting/uncommenting lines
the options can easily be changed.

4. From the `<source_dir>` directory execute the CTest Script with the following options:

    * **64-bit Windows with Visual Studio 2026:**

      `ctest -S HDF5config.cmake,BUILD_GENERATOR=VS2026 -C Release -VV -O hdf5.log`

    * **64-bit Windows with Visual Studio 2022:**

      `ctest -S HDF5config.cmake,BUILD_GENERATOR=VS202264 -C Release -VV -O hdf5.log`

    * **Linux and Mac:**

      `ctest -S HDF5config.cmake,BUILD_GENERATOR=Unix -C Release -VV -O hdf5.log`

    The command above will configure, build, test, and create an install package in
    the `<source_dir>` folder. It will have the format: `HDF5-2.X.Y-<platform>.<zip or tar.gz>`.

    On Unix, `<platform>` will be "Linux". A similar `.sh` file will also be created. On Windows, `<platform>`
    will be "win-vs2026_cl" or "win-vs2022_intel". If you have an installer on your system, you will also see a
    similar file that ends in either `.exe` (NSIS) or `.msi` (WiX).

    Notes on the command line options.

    * The `-S` option uses the script version of ctest.

    * The value for the `-C` option (as shown above, `-C Release`) must match the setting for
    `CTEST_CONFIGURATION_TYPE` in the platform configuration file.

    * The `-VV` option is for most verbose; use `-V` for less verbose.

    * The `-O hdf5.log` option saves the output to a log file `hdf5.log`.

5. To install, `X.Y` is the current release version.

    * **On Windows (with WiX):** Execute `hdf5-2.X.Y-win-vs2026_cl.msi` or `hdf5-2.X.Y-win-vs2022_intel.msi`.
    By default this program will install the HDF5 library into the `C:\Program Files` directory and will create
    the following directory structure:

          HDF_Group
          --HDF5
          ----2."X"."Y"
          ------bin
          ------include
          ------lib
          --------plugin
          ------cmake

    * **On Linux:** Change to the install destination directory (create it if doesn't exist) and execute
    `<path-to>/<source_dir>/HDF5-2.X.Y-Linux.sh`. After accepting the license, the script will prompt:

          By default the HDF5 will be installed in:
              <current directory>/HDF5-2.X.Y-Linux
          Do you want to include the subdirectory HDF5-2.X.Y-Linux?
          Saying no will install in: "<current directory>" [Yn]:

      Note that the script will create the following directory structure
      relative to the install point:

          HDF_Group
          --HDF5
          ----2."X"."Y"
          ------bin
          ------include
          ------lib
          --------plugin
          ------share

    * **On macOS:** Locate `HDF5-2.X.Y-Darwin.dmg` file in the `<source_dir>` folder. Click
      on the file to proceed with installation. After accepting the license,
      there will be a folder with the following structure:

          HDF_Group
          --HDF5
          ----2."X"."Y"
          ------bin
          ------include
          ------lib
          --------plugin
          ------share

    By default the installation will create the `bin`, `include`, `lib` and `cmake`
    folders in the `<install destination directory>/HDF_Group/HDF5/2.X.Y` directory.
    The `<install destination directory>` depends on the build platform:

    * Windows will set the default to `C:/Program Files/HDF_Group/HDF5/2.X.Y`

    * Linux will set the default to `<source_dir>/HDF_Group/HDF5/2.X.Y`

    The default can be changed by adding `,INSTALLDIR=<my new dir>` to the
    `ctest -S HDF5config.cmake...` command. For example on Linux:
    `ctest -S HDF5config.cmake,INSTALLDIR=/usr/local/myhdf5,BUILD_GENERATOR=Unix -C Release -VV -O hdf5.log`.
