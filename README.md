# Growing Critical: Self-Organized Criticality in a Developing Neural System
*Felipe Yaroslav Kalle Kossio, Sven Goedeke, Benjamin van den Akker, Borja Ibarz, and Raoul-Martin Memmesheimer*, 
[Phys. Rev. Lett. 121, 058301](https://doi.org/10.1103/PhysRevLett.121.058301), [arXiv](https://doi.org/10.48550/arXiv.1811.02861).

To replicate Figures 1, 2 and 3. 

Install GNU Scientific Library and Octave:
```
sudo apt update
sudo apt install libgsl-dev octave
```        
Compile C code:
```
make -C src
```        
Create directories for simulation output:
```
mkdir -p data/nc data/sc data/nc_ref1 data/nc_ref4 data/ab
```        
Run simulation with different parameters, for example:
```
xargs -a params/sc.txt ./simulate
```
Create Figure X:
```
cd figures
octave figureX.m
```
Figure 1 is a schematic. 
Figure 2 requires simulation results for parameters from nc.txt.
Figure 3 requires simulation results for parameters from nc.txt and sc.txt.
Figure 4(b) requires simulation results for parameters from nc.txt, nc_ref1.txt, nc_ref4.txt and ab.txt.

Versions:
- GCC 11.4.0
- GSL 2.7.1
- GNU Octave 6.4.0
