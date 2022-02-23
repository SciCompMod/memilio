# 3rd-Party Dependencies

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. The following table lists all libraries used by `memilio`:

| Library | Version  | Notes |
|---------|----------|-------|
| spdlog  | 1.5.0    | https://github.com/gabime/spdlog |
| Eigen   | 3.3.9    | http://gitlab.com/libeigen/eigen |
| Boost   | 1.75.0   | https://www.boost.org/ |
| JsonCpp | 1.9.5    | https://github.com/open-source-parsers/jsoncpp |
| HDF5    | 1.12.0   | https://www.hdfgroup.org/, package libhdf5-dev on apt (Ubuntu) |
| GoogleTest | 1.10.0| https://github.com/google/googletest |

## Conan Package Manager

The Conan package manager (https://conan.io) is used to aquire all dependencies automatically during CMake configuration. 

The recommended way to install Conan is with `pip`:
```bash
pip install conan
```
See [here](https://docs.conan.io/en/latest/installation.html) for other ways to install Conan. Conan is called by CMake, usually no action is required by the user.

You can use a virtual python environment to install Conan. This may be the same environment you used to install our python packages, but this is not required. Create and activate a virtual environment in the `venv` directory with
```bash
python -m venv venv
source ./venv/bin/activate
```

The Conan repository contains pre-built packages for many configurations, so most packages don't have to be built on your system during installation. Conan caches packages in your user directory (even if using a virtual environment), so packages only need to be downloaded and/or built once. However, multiple downloads are required for different compilers and configurations, e.g. one set of packages for Debug builds and one set for Release builds. 

By default, Conan stores the packages it builds or downloads in your home directory in `~/.conan`. On Windows, there is an additional storage directory to work around path length restrictions, usually at `C:\.conan`. The locations can be configured by editing the `~/.conan/conan.conf` file or setting the `CONAN_USER_HOME` and `CONAN_USER_HOME_SHORT` environment variables. See the [documentation](https://docs.conan.io/en/latest/reference/config_files/conan.conf.html#conan-conf) to learn more about configuring Conan.

For multi-configuration CMake generators (e.g. Visual Studio, XCode), packages for all configurations in the `CMAKE_CONFIGURATION_TYPES` variable are downloaded. Note that prebuilt packages are usually only available for Release and Debug configurations. Missing packages are built by Conan on demand, which can take a long time. It is recommended to set `CMAKE_CONFIGURATION_TYPES` to only the configurations you really require, e.g., `cmake .. "-DCMAKE_CONFIGURATION_TYPES=Debug;Release"`.

In some rare cases due to limitations of Conan and the large number of system variables, a pre-built package is installed that doesn't actually work on your system. If you get errors during compilation of `memilio` that are related to third party dependencies, it may help to force a local build of that package by setting the cmake variable `MEMILIO_CONAN_BUILD`. The variable accepts a semicolon separated list of packages or keywords. Package names can contain wildcard patterns, e.g., `b*` to build all packages that start with `b`. Accepted keywords:
- `missing`: build packages if no matching pre-built package is found in the repository (default)
- `always`: build all packages
- `cascade`: build any package that has a dependency on another package that is being built
Example: `MEMILIO_CONAN_BUILD=fmt;missing;cascade` forces Conan to build the `fmt` package, any package that depends on `fmt`, and any package that is not available pre-built. 
Check the `--build` flag of [`conan install`](https://docs.conan.io/en/latest/reference/commands/consumer/install.html) for more information.
