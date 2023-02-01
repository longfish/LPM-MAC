# LPM-CPP
A multi-threaded implementation of a nonlocal lattice particle method (LPM) using an iterative solution procedure. Below shows the image of slip system activation in a SENT sample [1].

![Slip](slip_system.png)

## Prerequisites
- Linux operating system (also work in Windows or other systems, but require other Integrated Development Environment like [Visual Studio](https://visualstudio.microsoft.com/))
- Editor (such as VS Code)
- CMake version 3.25+

## Building instructions

### Intel MKL environment
1. go to the Intel oneAPI website (https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html) download and install the base and hpc toolkits (install 2021.4+ version of oneAPI to have a proper support for CMake)
2. `source /opt/intel/oneapi/setvars.sh linux64 --force` (directory may change if customize the oneAPI installation)
3. install Ninja: `sudo apt-get install ninja-build`
4. install openmp: `sudo apt-get install libomp-dev`
5. note: you may need to add the MKL include path into your editor preference (like VS Code)

### Compile and run LPM-CPP
1. `source /opt/intel/oneapi/setvars.sh linux64 `
2. `git clone https://github.com/longfish/LPM-CPP.git`  # clone the project source files into your own machine
3. `mkdir build & cd build`
4. `cmake .. -G "Ninja" -DMKL_INTERFACE=ilp64 -DCMAKE_BUILD_TYPE=Release` # change to *Debug* from *Release* if needed
5. `cmake --build . -j 8`

### Run the code
`./lpmcpp`

The results will be in the build folder.

### Default numerical example
A 3D simple cubic lattice example for elastic material is provided by default in `./src` directory.

### Examples
There are some example codes in the `./examples` folder that contains additional numerical cases in [1, 2]. They define the main() functions of the project. One should change the example file name to "lpmc_project.c" and replace it with the same file in `./src` directory then run the code (possibly require some other modifications).

### References
1. Meng C, Wei H, Chen H, et al. Modeling plasticity of cubic crystals using a nonlocal lattice particle method[J]. Computer Methods in Applied Mechanics and Engineering, 2021, 385: 114069.

2. Meng C, Liu Y. Nonlocal Damage-enhanced Plasticity Model for Ductile Fracture Analysis Using a Lattice Particle Method[J]. arXiv preprint arXiv:2108.01214, 2021.
