from conans import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout


class HilbertplotCoreConan(ConanFile):
    name = "hilbertplot-core"
    version = "0.1"

    # Optional metadata
    license = "MIT"
    author = "Dario Fragas Gonzalez dariofg98@gmail.com"
    url = "https://github.com/dfg-98/hilbertplot-core"
    description = "Generate Hilbert's curve and hilbert plots"
    topics = ("data-visualization", "data-minig")

    requires = ["cmake/3.16.3", "fftw/3.3.9"]

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "src/*", "include/*"

    scm = {
        "type": "git",
        "subfolder": "hilbertplot-core",
        "url": "https://github.com/dfg-98/hilbertplot-core.git",
        "revision": "master",
    }

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["hilbertplot-core"]
