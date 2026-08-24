# 

TODO
For information on using HDF5 see the documentation, tutorials and examples
found at https://support.hdfgroup.org/documentation/hdf5/latest/.

A summary of the features included in the built HDF5 installation can be found
in the `libhdf5.settings` file in the same directory as the static and/or
shared HDF5 libraries. However CMake provides a programmable method to
determine the features of the library. The CMake installation will
provide a CMake package configuration file, located in the installation folder,
`cmake/hdf5-config.cmake`, and can be used to determine the features of the library.
The file is accessed by using the `find_package` command in your `CMakeLists.txt` file.

Set the `HDF5_ROOT` CMake variable `-DHDF5_ROOT=<install_path>`
or environment variable, `set(ENV{HDF5_ROOT} "<install_path>")`
to the installed location of HDF5. Examples:

* On Windows: `HDF5_ROOT=C:/Program Files/HDF_Group/HDF5/z.y.x/`
* On UNIX and alike: `HDF5_ROOT=<install root folder>/HDF_Group/HDF5/z.y.x/`

If you are using shared libraries, you may need to add to the path
environment variable. Set the path environment variable to the
installed location of the library files for HDF5. Examples:

* On Windows (*.dll): `PATH=%PATH%;C:/Program Files/HDF_Group/HDF5/z.y.x/bin`
* On UNIX and alike (*.so): `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<install root folder>/HDF_Group/HDF5/z.y.x/lib`

(Note there are no quote characters used on Windows and all platforms
use forward slashes)

Add the following to your `CMakeLists.txt` file:
```cmake
find_package (HDF5 CONFIG COMPONENTS C shared)
```

The components are optional and can be omitted if not needed. The components are:

  - `shared`
  - `static`
  - `C`
  - `CXX`
  - `Fortran`
  - `HL`
  - `CXX_HL`
  - `Fortran_HL`
  - `Java`
  - `Tools`

Users who want to build and run an application with HDF5 in Visual Studio
without using CMake should consult the [USING_HDF5_VS.md](./USING_HDF5_VS.md) file.
Building applications with the dynamic/shared hdf5 libraries requires that the
`H5_BUILT_AS_DYNAMIC_LIB` compile definition be used.