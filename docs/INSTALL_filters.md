# Build HDF5 with support for data filters

This file provides instructions for building HDF5 with support for various data filters, including
zlib and libaec. To see what filter support was configured into a built HDF5 library, check the
`I/O filters (external):` line in the `libhdf5.settings` file. Also, the `hdf5-config.cmake` file contains
the list of configuration options used to build the library.

---

## Table of Contents

* [CMake options](#cmake-options)
* [Using system-installed libraries](#using-system-installed-libraries)
* [Using libraries built with CMake FetchContent](#using-libraries-built-with-cmake-fetchcontent)
  * [CMake options for `HDF5_ALLOW_EXTERNAL_SUPPORT=GIT`](#cmake-options-for-hdf5_allow_external_supportgit)
  * [CMake options for `HDF5_ALLOW_EXTERNAL_SUPPORT=TGZ`](#cmake-options-for-hdf5_allow_external_supporttgz)
* [Troubleshooting](#troubleshooting)

---

## CMake options
<a name="cmake_options"></a>

HDF5 currently has built-in support for zlib, zlib-ng and libaec (szip), as well as support for a collection of dynamically-loaded filter libraries from the [HDF5 filter plugins](https://github.com/HDFGroup/hdf5_plugins) repository, which are controlled by the following CMake options:

  - `HDF5_ENABLE_ZLIB_SUPPORT` (`ON`/`OFF` Default: `OFF`) - Enable zlib support
  - `HDF5_ENABLE_SZIP_SUPPORT` (`ON`/`OFF` Default: `OFF`) - Enable libaec (szip) support
  - `HDF5_USE_ZLIB_NG` (`ON`/`OFF` Default: `OFF`) - Enable zlib-ng support
  - `HDF5_ENABLE_PLUGIN_SUPPORT` (`ON`/`OFF` Default: `OFF`) - Enable support for plugins from the HDF5 filter plugins repository
  - `HDF5_ALLOW_EXTERNAL_SUPPORT` (`NO`/`GIT`/`TGZ` Default: `NO`) - Specify whether support for libraries from external sources should be enabled (see [Using libraries built with CMake FetchContent](#fetchcontent_building))
  - `TGZPATH` (`<path>` Default: top level of HDF5 source directory) - When `HDF5_ALLOW_EXTERNAL_SUPPORT` is `TGZ` and a filter's `_USE_LOCALCONTENT` variable is `ON` (see below), specifies the directory to look in for a library's source code compressed file.

By default, setting `HDF5_ENABLE_ZLIB_SUPPORT` to `ON` will enable support for zlib. Setting the additional option `HDF5_USE_ZLIB_NG` to `ON` will cause HDF5 to instead use zlib-ng in place of zlib.

Additional CMake options for controlling how specific filters are found/used include the following:

[zlib(-ng) options](./INSTALL_options.md#zlib_filter_options)  
[libaec options](./INSTALL_options.md#libaec_filter_options)  
[HDF5 filter plugins options](./INSTALL_options.md#plugins_filter_options)  

> [!NOTE]
> The `HDF5_USE_ZLIB_STATIC`, `ZLIB_USE_EXTERNAL` and `ZLIB_USE_LOCALCONTENT` options also apply to zlib-ng when `HDF5_USE_ZLIB_NG` is set to `ON`.

## Using system-installed libraries

To build HDF5 with data filter support using libraries that are installed on the system, setting the relevant `HDF5_ENABLE_XXX_SUPPORT` [CMake options](#cmake_options) to `ON` should be enough in most cases. However, CMake can be instructed to look for specific libraries if desired or if the needed libraries are installed to directories that aren't in CMake's default search paths:

  - For zlib, set the CMake variable `ZLIB_ROOT` to the top level directory of a zlib installation
  - For zlib-ng, set the CMake variable `ZLIBNG_ROOT` to the top level directory of a zlib installation
  - For libaec (szip), set the CMake variable `SZIP_ROOT` or `libaec_ROOT` to the top level directory of a libaec installation

<b>Example:</b> Configure HDF5 with zlib and libaec support from system libraries

```bash
cmake -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_ENABLE_SZIP_SUPPORT=ON <hdf5_source_dir>
```

<b>Example:</b> Configure HDF5 with zlib installed to a non-system directory

```bash
cmake -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DZLIB_ROOT=/path/to/my/zlib/installation <hdf5_source_dir>
```

<b>Example:</b> Configure HDF5 with zlib-ng static libraries installed to a non-system directory

```bash
cmake -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_USE_ZLIB_NG=ON -DHDF5_USE_ZLIB_STATIC=ON -DZLIBNG_ROOT=/path/to/my/zlib-ng/installation <hdf5_source_dir>
```

## Using libraries built with CMake FetchContent
<a name="fetchcontent_building"></a>

> [!WARNING]
> HDF5 does not currently namespace libraries that are built using CMake's FetchContent functionality. Use caution when installing HDF5 to a system directory if any libraries are built this way, as there is risk of overwriting pre-existing libraries on the system. For this reason, it is recommended to install filter libraries on the system and allow HDF5 to find them at build time instead. This method is only provided as a convenience for building HDF5 binaries and for testing of HDF5.

When building HDF5, one may choose to build the libraries for data filters alongside HDF5 by using CMake's FetchContent functionality.
This will either download the source code for libraries or use compressed source code files from the filesystem and add the source code to the HDF5 build and installation processes. To set up this method, several CMake variables need to be set when configuring HDF5. First, the CMake variable `HDF5_ALLOW_EXTERNAL_SUPPORT` must be set to the value `GIT` or `TGZ`. The value `GIT` will instruct the HDF5 build process to retrieve the source code for libraries from git URLs, whereas the value `TGZ` will instruct the build process to retrieve the source code for libraries from compressed files from either a web URL or a path on the filesystem.

Next, external support for the desired data filters needs to be configured into HDF5. This simply involves setting the relevant `HDF5_ENABLE_XXX_SUPPORT` and `XXX_USE_EXTERNAL` [CMake options](#cmake_options) to `ON`.

Finally, any CMake options for how the source code should be retrieved and built must be set. The default values set for most of these options can be found in [CacheURLs.cmake](../config/CacheURLs.cmake) and can be changed by setting values for the CMake variables listed there. Note that the `HDF5_USE_XXX_STATIC` options affect which type of library, shared or static, is built when building external libraries.

> [!TIP]
> CMake options for a library can generally be passed through when building them with FetchContent by defining those options when configuring HDF5 (e.g., add `-DZLIB_BUILD_TESTING=OFF` to the list of HDF5 configuration options), but be aware that HDF5 may override the values for some of these variables.

### CMake options for `HDF5_ALLOW_EXTERNAL_SUPPORT=GIT`

When the `HDF5_ALLOW_EXTERNAL_SUPPORT` CMake option is specified as `GIT`, all filter libraries will be downloaded from git URLs using a specified git tag, commit hash or branch name. For example, using the default values for libaec will cause CMake to download the libaec source from `https://github.com/MathisRosenhauer/libaec.git` (`LIBAEC_GIT_URL`) using the git tag `v1.1.6` (`LIBAEC_GIT_TAG` as of the time of writing). Each filter's `_GIT_URL` and `_GIT_TAG` CMake variables may have values specified for them to change where the source code is downloaded from.

<b>Example:</b> Configure HDF5 with zlib-ng from system libraries and libaec (szip) from the default git repository

```bash
cmake -DHDF5_ALLOW_EXTERNAL_SUPPORT=GIT -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_USE_ZLIB_NG=ON -DHDF5_ENABLE_SZIP_SUPPORT=ON -DSZIP_USE_EXTERNAL=ON <hdf5_source_dir>
```

<b>Example:</b> Configure HDF5 with zlib-ng, libaec and the HDF5 filter plugins from the default git repositories

```bash
cmake -DHDF5_ALLOW_EXTERNAL_SUPPORT=GIT -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_USE_ZLIB_NG=ON -DZLIB_USE_EXTERNAL=ON -DHDF5_ENABLE_SZIP_SUPPORT=ON -DSZIP_USE_EXTERNAL=ON -DHDF5_ENABLE_PLUGIN_SUPPORT=ON -DPLUGIN_USE_EXTERNAL=ON <hdf5_source_dir>
```

<b>Example:</b> Configure HDF5 with libaec from a custom git repository

```bash
cmake -DHDF5_ALLOW_EXTERNAL_SUPPORT=GIT -DHDF5_ENABLE_SZIP_SUPPORT=ON -DSZIP_USE_EXTERNAL=ON -DLIBAEC_GIT_URL="https://github.com/user/libaec.git" -DLIBAEC_GIT_TAG="my_libaec_tag" <hdf5_source_dir>
```

### CMake options for `HDF5_ALLOW_EXTERNAL_SUPPORT=TGZ`

When the `HDF5_ALLOW_EXTERNAL_SUPPORT` CMake option is specified as `TGZ`, the way in which filter libraries are retrieved depends on the value for that filter's `_USE_LOCALCONTENT` CMake option. For each filter that has its `_USE_LOCALCONTENT` variable set to `OFF`, CMake will download a compressed file with the relevant library's source code from a URL formed by combining the filter's `_TGZ_ORIGPATH` and `_TGZ_NAME` CMake variables. For example, the default values of `LIBAEC_TGZ_ORIGPATH` and `LIBAEC_TGZ_NAME` are `https://github.com/MathisRosenhauer/libaec/releases/download/v<libaec_version>` and `libaec-<libaec_version>.tar.gz`, respectively. As of the time of writing, the current default version of libaec is 1.1.6, resulting in a URL of https://github.com/MathisRosenhauer/libaec/releases/download/v1.1.6/libaec-1.1.6.tar.gz. CMake will download this file, uncompress it and use the contents as source code when building the library alongside HDF5. Each filter's `_TGZ_ORIGPATH` and `_TGZ_NAME` CMake variables may have values specified for them to change where the file is downloaded from.

For each filter that has its `_USE_LOCALCONTENT` variable set to `ON`, CMake will search the directory specified by `TGZPATH` for the compressed file with that filter's library source code according to the default names, which can be changed by setting a filter's `_TGZ_NAME` variable. For example, to locate compressed files with zlib and libaec libraries with the names `zlib.tar.gz` and `libaec.tar.gz` in the directory `/home/user/filters`, specify the following CMake options when configuring HDF5:

  - `-DTGZPATH="/home/user/filters"`
  - `-DZLIB_TGZ_NAME="zlib.tar.gz"`
  - `-DLIBAEC_TGZ_NAME="libaec.tar.gz"`

CMake will combine `TGZPATH` with each of `ZLIB_TGZ_NAME` and `LIBAEC_TGZ_NAME` to form full paths to each file.

<b>Example:</b> Configure HDF5 with zlib-ng and libaec from downloaded tar.gz files with the default names

```bash
cmake -DHDF5_ALLOW_EXTERNAL_SUPPORT=TGZ -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_USE_ZLIB_NG=ON -DZLIB_USE_EXTERNAL=ON -DHDF5_ENABLE_SZIP_SUPPORT=ON -DSZIP_USE_EXTERNAL=ON <hdf5_source_dir>
```

<b>Example:</b> Configure HDF5 with zlib-ng and libaec from tar.gz files on the system with custom names

```bash
cmake -DHDF5_ALLOW_EXTERNAL_SUPPORT=TGZ -DTGZPATH="/path/to/files" -DHDF5_ENABLE_ZLIB_SUPPORT=ON -DHDF5_USE_ZLIB_NG=ON -DZLIB_USE_EXTERNAL=ON -DZLIB_USE_LOCALCONTENT=ON -DZLIBNG_TGZ_NAME="zlib-ng-2.3.3.tar.gz" -DHDF5_ENABLE_SZIP_SUPPORT=ON -DSZIP_USE_EXTERNAL=ON -DLIBAEC_USE_LOCALCONTENT=ON -DLIBAEC_TGZ_NAME="libaec-1.1.6.tar.gz" <hdf5_source_dir>
```

## Troubleshooting

If the HDF5 build process does not properly locate libraries for a data filter on the system, debugging of the process can be enabled when configuring HDF5 by passing the option `-DCMAKE_FIND_DEBUG_MODE=ON`. This option is intended for advanced users, but will give information about what paths on the system are or are not being searched.
