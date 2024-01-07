# LPM-MAC

A multi-threaded implementation of a nonlocal lattice particle method (LPM) using an iterative solution procedure. The default force unit is N, displacement unit is mm, pressure unit is MPa. Below shows the image of slip system activation in a SENT sample [1].

![Slip](slip_system.png)

## Building instructions

### Eigen + OpenBLAS

#### For macOS

1. Install Homebrew from `https://brew.sh/` 
2. Install LLVM and openmp: `brew install llvm` and `brew install libomp`, add into environment variables
4. Install OpenBLAS: `brew install openblas`, export `LIBRARY_PATH` and `CPLUS_INCLUDE_PATH` environment variables for OpenBLAS
4. Download eigen from `https://eigen.tuxfamily.org/`
5. Download boost from `https://www.boost.org/`
6. Export both eigen and boost root folder into `CPLUS_INCLUDE_PATH` environment variable

### Compile and run LPM

0. `git pull https://github.com/longfish/LPM-CPP.git` under `LPM-CPP` folder if need an updated version of the code
1. `source /opt/intel/oneapi/setvars.sh linux64 `
2. `git clone https://github.com/longfish/LPM-CPP.git`  # clone the project source files into your own machine
3. `mkdir build`
4. `cd build`
5. `cmake .. -G "Ninja" -DMKL_INTERFACE=ilp64 -DCMAKE_BUILD_TYPE=Release` # change from *Release* to *Debug* for debugging (e.g., valgrind)
6. `cmake --build . -j 8`

## Run the code

`./lpmcpp`

The results will be in the build folder.

## Examples

There are some example files in the `./examples` folder that contains additional numerical cases such as those in [1, 2]. They define the `run()` functions of the project. Please change/add the example file and also include them into `src/lmpcpp.cpp` to run the code. Please note that other code pieces, such as `assembly.h`, `lpm.h`, etc. may also need to be changed.

## Loading/Geometry

To ease some common geometries that often see in numerical simulations, customized loading and geometry codes are provided. The `loading` folder contains Jupyter code, while `geometry` folder contains cpp code.    

## References

1. Meng C, Wei H, Chen H, et al. Modeling plasticity of cubic crystals using a nonlocal lattice particle method[J]. Computer Methods in Applied Mechanics and Engineering, 2021, 385: 114069.
2. Meng C, Liu Y. Damage-augmented nonlocal lattice particle method for fracture simulation of solids[J]. International Journal of Solids and Structures, 2022, 243: 111561.

## Appendix

1. Valgrind command: `valgrind --leak-check=yes --show-leak-kinds=all --log-file=valgrind.rpt ./lpmcpp`