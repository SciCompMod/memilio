# 3rd-Party Dependencies #

This directory contains CMake configuration of the dependencies of the MEmilio C++ library. See [the MEmilio C++ README] for the full list of dependencies.

*Upgrading Boost (For developers only)*

We currently bundle only a minified boost library that contains only boost filesystem.
To update boost to a new version, the boost tool `bcp` can be used to bundle only the 
required files, see:  https://www.boost.org/doc/libs/1_72_0/tools/bcp/doc/html/index.html

To do so, call

```bash
./bcp filesystem path_to_epi_source/cpp/thirdparty/boost_<version>
```

to copy the required files. Then, this folder has to be compressed into a `.tar.gz`.

Then, the file `cpp/cmake/BuildBoost.cmake` has to be adapted accordingly.
   