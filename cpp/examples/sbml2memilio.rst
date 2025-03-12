SBML integration 
================

The SBML integration into memilio works via the `sbml2memilio` executable. 
Building it requires the installation of `libsbml <https://sbml.org/software/libsbml/>`_ . Then the general build 
command will produce the executable in the `build/bin`-directory. 


Usage
------
Call the executable on a given SBML file as 

```
sbml2memilio <sbml-file>
```
If `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ is not installted, it will end with an error, but produce the necessary files nevertheless.

It will produce 
 - a folder that can be copied to the `cpp/models`-directory of the memilio repository
 - an implementation-file that can be copied to the `cpp/examples`-directory of the memilio repository
 - a file called `CMakeListsAddition.txt` that contains the necessary additions to the `cpp/examples/CMakeLists.txt`-file to include the new model in the build process.
 - a file called `CMakeListsFolderNames.txt` that contains the necessary additions to the `cpp/CMakeLists.txt` to include the new model directory in the build process (starting in line 152).

 Once the directory and example file are copied to the correct locations and the necessary changes to the CMakeLists are done, 
 calling the general build function will also build an executable for the new model.
 

Changing parameters
----------------------
As some parameters (mainly the duration of simulation) are not part of the SBML file, they are set to generic values in 
the generated example file. They can also be changed there.


Limitations
-------------
Not every feature implemented in SBML is also supported by memilio.
The following features are the most important not supported features:
    - Events not triggered by time
    - Usage of multiple SBML-compartments
    - Rules (unless they are RateRules)

In general: Please check your models after the conversion to ensure that everything is working as expected.

