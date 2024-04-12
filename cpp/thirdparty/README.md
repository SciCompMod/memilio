# 3rd-Party Dependencies

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. See [the MEmilio C++ README](../README.md) for the full list of dependencies.

## Bundled Dependencies

Most dependencies of this project don't need to be installed manually. These dependencies are bundled by cloning an external repository. Using the `MEMILIO_USE_BUNDLED_<XYZ>` CMake options (where `<XYZ>` is the name of the dependency), installed packages can be used instead of the bundled packages, using the usual `find_package` mechanism.

### External Repositories

The repository of the dependency is cloned during CMake configuration into the `<build>/_deps/<xyz>-src` directory. The dependency is then built together with the MEmilio project. The version of the package is set in the [thirdparty CMakeLists.txt](CMakeLists.txt). To upgrade the version, simply increase the version number there.

Note: Cloning the boost git repository can take a while. Especially for this dependency, it may be useful to set the `MEMILIO_USE_BUNDLED_BOOST` option to `OFF` if the package is already installed. The installed boost version must be at least version 1.76.0.
 It is planned to offer the option to use a minimal extract of boost library as an archive included in the MEmilio project. With the minimal version of boost, only limited functionality of MEmilio can be used.