# 3rd-Party Dependencies

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. See [the MEmilio C++ README](../README.md) for the full list of dependencies.

## Conan Package Manager

The Conan package manager (https://conan.io) is used to aquire all dependencies automatically during CMake configuration. If this is not desired, e.g., because Conan is not available on your system or it is not connected to the internet, Conan can be disabled using the CMake option `MEMILIO_USE_CONAN=OFF`. 

The recommended way to install Conan is with `pip`:
```bash
pip install conan
```
See [here](https://docs.conan.io/en/latest/installation.html) for other ways to install Conan. Conan is run by CMake, it is not required to run Conan manually.

You can use a virtual python environment to install Conan. This may be the same environment you used to install our python packages, but this is not required. Create and activate a virtual environment in the `venv` directory with
```bash
python -m venv venv
source ./venv/bin/activate
```

Conan caches packages in your user directory (even if using a virtual environment), so packages only need to be downloaded once. However, multiple downloads are required for different compilers and configurations, e.g. one set of packages for Debug builds and one set for Release builds. The configuration is detected automatically and logged by CMake, no action is required by the user. 

For multi-configuration CMake generators (e.g. Visual Studio, XCode), packages for all configurations in the `CMAKE_CONFIGURATION_TYPES` variable are downloaded. Note that prebuilt packages are usually only available for Release and Debug configurations. Missing packages are built by Conan on demand, which can take a long time. It is recommended to set `CMAKE_CONFIGURATION_TYPES` to only the configurations you really require, e.g., `cmake .. "-DCMAKE_CONFIGURATION_TYPES=Debug;Release"`.

There is a CMake option `MEMILIO_USE_SYSTEM_<XYZ>` for each dependency. This stops Conan from downloading and using the package for the specified dependency. This is not recommended but may be required, e.g., to use a specific version of the dependency optimized for your system. 

If you disable Conan or specific packages, CMake will try to find a compatible version of each dependency on your system instead. You may have to give CMake hints where to find the dependency on your system, e.g., use `-DHDF5_DIR=<path/to/HDF5>` to point CMake to your HDF5 installation. To be precise, the path given must contain a `HDF5Config.cmake` file which should be included in your HDF5 distribution.
