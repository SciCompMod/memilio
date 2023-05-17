# 3rd-Party Dependencies

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. See [the MEmilio C++ README](../README.md) for the full list of dependencies.

## Bundled Dependencies

Most dependencies of this project don't need to be installed manually. Dependencies are bundled in two different ways: cloning an external repository or as an archive included in the MEmilio project. Using the `MEMILIO_USE_BUNDLED_<XYZ>` CMake options (where `<XYZ>` is the name of the dependency), installed packages can be used instead of the bundled packages, using the usual `find_package` mechanism.

### External Repositories

The repository of the dependency is cloned during CMake configuration into the `<build>/_deps/<xyz>-src` directory. The dependency is then built together with the MEmilio project. The version of the package is set in the [thirdparty CMakeLists.txt](CMakeLists.txt). To upgrade the version, simply increase the version number there.
