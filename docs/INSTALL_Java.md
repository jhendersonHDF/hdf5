# HDF5 Java-specific installation information

---

## Table of Contents

* [1. Java implementation selection](#1-java-implementation-selection)

---

## 1. Java implementation selection

HDF5 Java bindings support two native interface implementations:
- JNI (Java Native Interface). Default, works with Java 11+, production-stable.
- FFM (Foreign Function & Memory). Optional, requires Java 25+, modern native access.

Maven Artifacts:

- `org.hdfgroup:hdf5-java-ffm` is FFM implementation
- `org.hdfgroup:hdf5-java-jni` is JNI implementation

Both implementations use the same `hdf.hdf5lib.*` package structure for seamless migration.

To build HDF5 with Maven deployment support:

    cmake -DHDF5_BUILD_JAVA:BOOL=ON -DHDF5_ENABLE_MAVEN_DEPLOY:BOOL=ON ../hdf5

To build Maven snapshot versions for development:

    cmake -DHDF5_BUILD_JAVA:BOOL=ON -DHDF5_ENABLE_MAVEN_DEPLOY:BOOL=ON -DHDF5_MAVEN_SNAPSHOT:BOOL=ON ../hdf5

> **Note:** FFM is selected for Java 25+ if `HDF5_ENABLE_JNI` is `OFF`.
> To force JNI even with Java 25+:
>
> ```
> cmake -DHDF5_BUILD_JAVA:BOOL=ON -DHDF5_ENABLE_MAVEN_DEPLOY:BOOL=ON -DHDF5_ENABLE_JNI:BOOL=ON ../hdf5
>```

Or use the Maven-enabled CMake presets (recommended):

```bash
# Minimal build for Java artifacts only (recommended for Maven deployment)
# Linux (GCC) - JNI (default, Java 8+):
cmake --workflow --preset ci-MinShar-GNUC-Maven --fresh          # JNI Release
cmake --workflow --preset ci-MinShar-GNUC-Maven-Snapshot --fresh # JNI Snapshot

# Linux (GCC) - FFM (optional, Java 25+):
cmake --workflow --preset ci-MinShar-GNUC-Maven-FFM --fresh          # FFM Release
cmake --workflow --preset ci-MinShar-GNUC-Maven-FFM-Snapshot --fresh # FFM Snapshot

# Windows (MSVC) - JNI (default):
cmake --workflow --preset ci-MinShar-MSVC-Maven --fresh          # JNI Release
cmake --workflow --preset ci-MinShar-MSVC-Maven-Snapshot --fresh # JNI Snapshot

# Windows (MSVC) - FFM (optional):
cmake --workflow --preset ci-MinShar-MSVC-Maven-FFM --fresh          # FFM Release
cmake --workflow --preset ci-MinShar-MSVC-Maven-FFM-Snapshot --fresh # FFM Snapshot

# macOS (Clang) - JNI (default):
cmake --workflow --preset ci-MinShar-Clang-Maven --fresh          # JNI Release
cmake --workflow --preset ci-MinShar-Clang-Maven-Snapshot --fresh # JNI Snapshot

# macOS (Clang) - FFM (optional):
cmake --workflow --preset ci-MinShar-Clang-Maven-FFM --fresh         # FFM Release
cmake --workflow --preset ci-MinShar-Clang-Maven-FFM-Snapshot --fresh # FFM Snapshot

# ROS3 VFD (S3 cloud storage) - Add to any preset above:
cmake --workflow --preset ci-MinShar-GNUC-Maven --fresh -DHDF5_ENABLE_ROS3_VFD=ON
```

>**Note:** Presets are platform-specific. Use `cmake --list-presets` to see available presets for your
current platform. Minimal Maven presets skip examples, testing, tools, C++, and Fortran builds to
optimize for Java artifact generation only.

#### Java Examples Maven Integration

The HDF5 Java examples are available as a separate Maven artifact: `org.hdfgroup:hdf5-java-examples`. It
contains platform-specific dependencies to ensure compatibility with the HDF5 Java library, and
complete examples with documentation for all HDF5 Java functionality. See
[HDF5Examples/JAVA/README-MAVEN.md](../HDF5Examples/JAVA/README-MAVEN.md) for usage instructions.

#### Testing Java Examples with Maven Artifacts

- Maven staging workflow validates examples against artifacts.
- Cross-platform testing ensures compatibility on all supported platforms.
- Native library error handling validates JAR structure in Maven-only environments.
- Fork-based testing allows validation on repository forks before canonical deployment.
- Dynamic repository workflows adapt to any GitHub repository automatically.
