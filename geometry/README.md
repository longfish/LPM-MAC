# Create geometry (mainly) for polycrystal simulation
For all Fortran code, run the following:

`ulimit -s unlimited`
`export KMP_STACKSIZE=5G`
`ifort -qmkl=parallel -qopenmp create_BCC_configuration.F90 -o bcc`