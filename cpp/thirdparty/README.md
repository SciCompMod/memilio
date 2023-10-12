# 3rd-Party Dependencies

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. See [the MEmilio C++ README](../README.md) for the full list of dependencies.

## Bundled Dependencies

Most dependencies of this project don't need to be installed manually. Dependencies are bundled in two different ways: cloning an external repository or as an archive included in the MEmilio project. Using the `MEMILIO_USE_BUNDLED_<XYZ>` CMake options (where `<XYZ>` is the name of the dependency), installed packages can be used instead of the bundled packages, using the usual `find_package` mechanism.

### External Repositories

The repository of the dependency is cloned during CMake configuration into the `<build>/_deps/<xyz>-src` directory. The dependency is then built together with the MEmilio project. The version of the package is set in the [thirdparty CMakeLists.txt](CMakeLists.txt). To upgrade the version, simply increase the version number there.

### Archived (Boost)

We currently bundle only a minimal extract of boost library as an archive that contains only the libraries of boost that we use. Currently, these libraries are filesystem, outcome, and optional, including transitive dependencies. The archive has been created using the boost tool `bcp`, see https://www.boost.org/doc/libs/1_72_0/tools/bcp/doc/html/index.html.

To upgrade boost, follow these steps: 

1. call `bcp` to copy the required files (see the `bcp` documentation for details)
```bash
./bcp optional outcome filesystem path_to_epi_source/cpp/thirdparty/boost_<version>
```
2. compress the folder into a `.tar.gz` archive and replace the existing archive in the repository.
3. adapt the file `cpp/cmake/BuildBoost.cmake`. At least update the version number and archive name.
   
