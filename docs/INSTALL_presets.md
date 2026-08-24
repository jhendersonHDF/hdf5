# Building HDF5 with CMake presets

[CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html) can be a convenient
way of building HDF5 using fixed sets of options for common configurations. These instructions
assume that the HDF5 source code has already been obtained and is located in a directory referred to
below as `<source_dir>`.

---

## Table of Contents

* [Using workflow presets](#using-workflow-presets)
  * [Prerequisites](#prerequisites)
* [Using configuration presets](#using-configuration-presets)

---

## Using workflow presets

Workflow presets can be used to configure, build, test and package HDF5 for common configurations
using a single command. The available workflow presets can be shown by running

```Shell
cmake --workflow --list-presets
```

from the `<source_dir>` directory.

### Prerequisites



## Using configuration presets

TODO

The available configuration presets can be shown by running

```Shell
cmake --list-presets=configure
```

from the `<source_dir>` directory.

## 

Available CMake presets for building HDF5 can be shown with

```bash
cmake --workflow --list-presets
```

TODO
## 4. Quick Start Presets

### 4.1. Solution

HDF Group provides a file in the source, `CMakePresets.json`. This file is in the HDF5 library source as well as the `HDF5Examples` source of the installed binary.

* **Library Source File:** Builds HDF5 with the options for building a typical shared library with the common languages for a platform. The features include building the tools, examples, plugins, and the shared and static libraries.

* **HDF5Examples Source File:** Builds the examples with the components that were enabled by the options selected when the install HDF5 library was built. The typical library built by HDF5 and available from the HDF5 Releases page includes C, Java, and Fortran compilers along with the tools, examples, plugins, and the shared and static libraries.

### 4.2. Discussion

The `CMakePresets.json` file is located in the root directory of the HDF5 source and the `HDF5Examples` source of the installed binary.

It is from here you will execute the `cmake` command to build HDF5. The following example shows how to build HDF5 or examples with the `CMakePresets.json` file:

```bash
# Change directory to the source folder
$ cd <source-folder>

# Execute the workflow
$ cmake --workflow --preset ci-StdShar-<compiler-type> --fresh
```
> **Note:** `<compiler-type>` should be replaced with `GNUC`, `MSVC`, or `Clang`.

The above example will create a `build` folder in the source parent directory, which will contain the results of the build, including installation package files when the library is built.

### Prerequisites
* **Ninja build system** (recommended, should be downloaded if not available)

### Quick Start (3 steps)
1. Change to the HDF5 source directory:
   ```bash
   cd /path/to/hdf5-2.x.y
   ```
2. Execute a workflow preset:
   ```bash
   cmake --workflow --preset ci-StdShar-GNUC --fresh       # Linux/Mac with GCC
   cmake --workflow --preset ci-StdShar-MSVC --fresh       # Windows with MSVC
   cmake --workflow --preset ci-StdShar-Clang --fresh      # Linux/Mac with Clang
   ```
3. Find your build artifacts in:
   ```text
   ../build/ci-StdShar-<compiler>/
   ```

That's it! The workflow preset automatically:
* Configures the build
* Compiles libraries and tools
* Runs tests
* Creates installation packages

### Available Presets
View all available presets:
```bash
cmake --list-presets
```

**Common presets:**
* **Standard Builds:**
  * `ci-StdShar-GNUC`        (Standard shared libraries - GCC)
  * `ci-StdShar-MSVC`        (Standard shared libraries - MSVC)
  * `ci-StdShar-Clang`       (Standard shared libraries - Clang)
  * `ci-MinShar-GNUC`        (Minimal shared libraries - GCC)
* **Java Builds:**
  * `ci-StdShar-GNUC-Java-FFM`     (Java FFM bindings - GCC)
  * `ci-StdShar-GNUC-Java-JNI`     (Java JNI bindings - GCC)
* **Maven Deployment (JNI default - Java 8+):**
  * `ci-MinShar-GNUC-Maven-Snapshot`               (JNI snapshots for Maven)
  * `ci-MinShar-GNUC-Maven`                        (JNI release for Maven)
* **Maven Deployment (FFM optional - Java 25+):**
  * `ci-MinShar-GNUC-Maven-FFM-Snapshot`           (FFM snapshots for Maven)
  * `ci-MinShar-GNUC-Maven-FFM`                    (FFM release for Maven)

See [Section XI]((#section-xi)) for creating custom preset configurations.

### Individual Preset Commands (Advanced)
If you prefer to run preset steps individually (where `<compiler-type>` is `GNUC`, `MSVC`, or `Clang`):

```bash
cd /path/to/hdf5-source
cmake --preset ci-StdShar-<compiler-type>                 # Configure
cmake --build --preset ci-StdShar-<compiler-type>         # Build
ctest --preset ci-StdShar-<compiler-type>                 # Test
cpack --preset ci-StdShar-<compiler-type>                 # Package
```

The workflow preset (shown in Quick Start above) runs all these steps automatically.


<a id="section-xi"></a>
## XI. Creating Custom Preset Configurations

The quickest way to customize your build is to create a `CMakeUserPresets.json` file in the HDF5 source directory.

**Basic Customization Steps:**
1. Copy `CMakePresets.json` to `CMakeUserPresets.json`.
2. Edit `CMakeUserPresets.json`. Change configuration names from `ci-*` to `my-*` and modify the fields as needed.

### Example: Using Local Support Files

To change external support files to use a local directory:

```json
{
  "name": "my-base-tgz",
  "hidden": true,
  "inherits": "ci-base",
  "cacheVariables": {
    "HDF5_ALLOW_EXTERNAL_SUPPORT": {"type": "STRING", "value": "TGZ"},
    "TGZPATH": {"type": "PATH", "value": "${sourceParentDir}/temp"}
  }
},
{
  "name": "my-StdCompression",
  "hidden": true,
  "inherits": "my-base-tgz",
  "cacheVariables": {
    ...
  }
},
{
  "name": "my-StdShar",
  "hidden": true,
  "inherits": "my-StdCompression",
  "cacheVariables": {
    ...
  }
},
{
  "name": "my-StdShar-GNUC",
  "description": "My Custom GNUC Standard Config for x64 (Release)",
  "inherits": [
    "ci-x64-Release-GNUC",
    "ci-CPP",
    "ci-Fortran",
    "ci-Java",
    "my-StdShar",
    "my-StdExamples"
  ]
}
```

Then you can change or add options for your specific case.

### Example: Maven Deployment Preset

For Maven deployment with custom repository URL:

```json
{
  "name": "my-maven-custom",
  "inherits": "ci-MinShar-GNUC-Maven-Snapshot",
  "cacheVariables": {
    "MAVEN_REPOSITORY_URL": {
      "type": "STRING",
      "value": "https://your-repo.com/maven"
    },
    "HDF5_ENABLE_ROS3_VFD": {
      "type": "BOOL",
      "value": "ON"
    }
  }
}
```

Build with:
```bash
cmake --workflow --preset my-maven-custom --fresh
```

> **Note:** This example uses JNI (default). For FFM, inherit from `ci-MinShar-GNUC-Maven-FFM-Snapshot`.

### Preset File Details

`CMakePresets.json` and `CMakeUserPresets.json`:

* Live in the project's root directory.
* Both have the same format.
* `CMakePresets.json` is for project-wide build details (**don't modify**).
* `CMakeUserPresets.json` is for local build customizations (modify to your heart's content).

The HDF Group presets require CMake 3.26 and use the Ninja build system.
Ninja may need to be installed separately on some platforms.

Hidden presets (marked `"hidden": true`) are used for inheritance and
cannot be used directly. They are defined in `config/cmake-presets/hidden-presets.json`.

For more information:
    https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html

TODO

### Preset Architecture
- **Layered inheritance**: Base presets + feature-specific + platform-specific
- **Hidden presets**: Reusable components (`ci-base`, `ci-Debug`, `ci-Release`, `ci-Maven`, `ci-Maven-Snapshot`, `ci-Maven-Minimal`, `ci-Maven-Minimal-Snapshot`)
- **Platform presets**: `ci-GNUC`, `ci-Clang`, `ci-MSVC`, `ci-macos`
- **Maven presets**: Hidden base configurations for Maven deployment support
- **Minimal Maven presets**: Streamlined configurations for Java artifact generation only
- **Build type matrix**: Debug, Release (RelWithDebInfo + docs), Maven variants

### Key Preset Patterns
```bash
# Standard shared library builds
cmake --workflow --preset ci-StdShar-GNUC --fresh      # GCC
cmake --workflow --preset ci-StdShar-Clang --fresh     # Clang
cmake --workflow --preset ci-StdShar-MSVC --fresh      # MSVC

# Maven-enabled builds (Java artifacts with deployment support)
cmake --workflow --preset ci-StdShar-GNUC-Maven --fresh          # Maven release (full build)
cmake --workflow --preset ci-StdShar-GNUC-Maven-Snapshot --fresh # Maven snapshot (full build)
cmake --workflow --preset ci-MinShar-GNUC-Maven --fresh          # Maven release (minimal build)
cmake --workflow --preset ci-MinShar-GNUC-Maven-Snapshot --fresh # Maven snapshot (minimal build)

# Multi-platform Maven presets (minimal builds for Java artifacts only)
cmake --workflow --preset ci-MinShar-MSVC-Maven --fresh          # Windows Maven
cmake --workflow --preset ci-MinShar-Clang-Maven --fresh         # macOS Maven

# Naming convention: ci-[Features]-[Compiler][-Maven[-Snapshot]]
# Features: Std (standard), Min (minimal), StdShar (standard shared), MinShar (minimal shared)
# Maven: Adds Maven deployment support with platform-specific JARs
# Snapshot: Adds -SNAPSHOT suffix for development versions
# Minimal Maven presets: Skip examples, testing, tools, C++, Fortran - Java artifacts only
# Java Examples Maven Integration: Comprehensive testing of Java examples with Maven artifacts across all platforms
```

### Preset Configuration Strategy
- **Binary directory**: `${sourceParentDir}/build/${presetName}`
- **Install directory**: `${sourceParentDir}/install/${presetName}`
- **Generator**: Ninja (default for most presets)
- **External libraries**: TGZ/GIT support for zlib, szip, libaec
