# Membrane Analyses
 Some ideas for analyses of membrane simulations


### Recalculate atomic forces in specific molecules of your system in a OpenMM trajectory
- recalc_forces_openmm.py

```
usage: recalc_forces_openmm.py [-h] -crd CRD -psf PSF -dcd DCD -toppar TOPPAR -boxsize BOXSIZE BOXSIZE BOXSIZE -mol MOL -which_forces WHICH_FORCES [-temp TEMP] [-pressure PRESSURE] [-dt DT]
                               [-is_drude IS_DRUDE]

optional arguments:
  -h, --help            show this help message and exit
  -crd CRD              Input coordinate file (.crd)
  -psf PSF              Topology file in XPLOR format (.psf)
  -dcd DCD              Trajectory file in DCD format (.dcd)
  -toppar TOPPAR        Force field stream file (ex. "toppar.str").
  -boxsize BOXSIZE BOXSIZE BOXSIZE
                        Box vector lengths (a,b,c) in nanometers
  -mol MOL              Resname of your molecule of interest.
  -which_forces WHICH_FORCES
                        All forces or just Nonbonded forces
  -temp TEMP            Target temperature, in Kelvin. Default = 298.
  -pressure PRESSURE    Target pressure, in bar. Default is 1.0.
  -dt DT                Integration step (in fs). Default = 1.
  -is_drude IS_DRUDE    If it is a Drude system or not
```


### Calculate solvent density in a membrane system (useful for water permeation)
- calc_solv_dens.py

```
python calc_solv_dens.py   topology_file    trajectory_file     output_density_grid.dx
```
