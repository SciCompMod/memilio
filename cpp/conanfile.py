from conans import ConanFile, CMake

class MemilioConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "boost/1.75.0", "spdlog/1.5.0", "eigen/3.3.9"
    options = {"with_tests": [True, False], "with_io": [True, False]}
    generators = "cmake_find_package_multi"
    default_options = {"with_tests": True, "with_io": True}

    def configure(self):
        if self.options.with_tests:
           self.requires("gtest/1.10.0")
        if self.options.with_io:
            self.requires("hdf5/1.12.0")
            self.requires("jsoncpp/1.9.5")

