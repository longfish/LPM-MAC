# Create geometry (mainly) for polycrystal simulation
For all Fortran code, run the following:

To run Fortran code, 
`ifort -qmkl=parallel -qopenmp create_BCC_configuration.F90 -o bcc`

To run Atomsk (you need to first download the Atomsk package),
`atomsk --create fcc 4.046 Al aluminium.xsf`
`atomsk --polycrystal aluminium.xsf polycrystal.txt final.lmp -wrap ` # if need to wrap the particles
