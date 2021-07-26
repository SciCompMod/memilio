3rd-Party Dependencies
======================

The project requires a number of other libraries.
Some of them are bundled with the project, such that
they can be built as part of the project.

If a bundled library should not be used but another 
version installed in the system, the cmake argument

    -DEPI_USE_BUNDLED_<LIBNAME>=OFF

must be used.

Epidemiology uses the following dependencies

+------------+------------+----------------+-----------+ 
| Library    | Version    |  Used in       | Bundled   |
+============+============+================+===========+
| Eigen      | 3.3.0      | epidemiology   | Yes       |
+------------+------------+----------------+-----------+
| Spdlog     | 1.5.0      | epidemiology   | Yes       |
+------------+------------+----------------+-----------+
| Boost      |            |                |           |
| outcome    | 1.75.0     | epidemiology   | Yes       |
+------------+------------+----------------+-----------+
| Boost      |            |                |           |
| optional   | 1.75.0     | epidemiology   | Yes       |
+------------+------------+----------------+-----------+
| JSonCPP    | 1.7.4      | epidemiology-io| Yes       |
+------------+------------+----------------+-----------+
| Boost      |            |                |           |
| filesystem | 1.75.0     | epidemiology-io| Yes       |
+------------+------------+----------------+-----------+
| HDF5       | 1.10       | epidemiology-io| No        |
+------------+------------+----------------+-----------+

Updating Boost
--------------

We currently bundle only a minified boost library that contains only boost filesystem.
To update boost to a new version, the boost tool `bcp` can be used to bundle only the 
required files, see:  https://www.boost.org/doc/libs/1_72_0/tools/bcp/doc/html/index.html

To do so, call

.. code:: sh

    ./bcp filesystem \
       path_to_epi_source/cpp/thirdparty/boost_<version>

to copy the required files. Then, this folder has to be compressed into a `.tar.gz`.

Then, the file `cpp/cmake/BuildBoost.cmake` has to be adapted accordingly.
   