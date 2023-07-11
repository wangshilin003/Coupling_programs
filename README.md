# Coupling-programs-between-the-two-scales
These programs are used for the paper of "Multi-scale Simulation of The Effect of Microbial Growth on The Permeability of Porous Media".
Program of "dugks1.cpp" is used for calculating the velocity  and nutrient concentration distributions of the REV computational cells in step 1. 
The representative cells for microbial growth during bioclogging are determined according to experimental observations.
The velocity and nutrient concentration of determined cells are used as the inlet conditions of the pore scale calculations, which are performed by the program of "LBM-IMB-CA.cpp" in step 2.
The permeability at the pore scale is obtained and then used for calculating the equivalent porosity, which is calculated by the program of "Equivalent porosity calculation.cpp" in step 3.
The obtained equivalent porosity at different time is fitted as a variation curve over time.
The variation curve of equivalent porosity is then used in the REV scale calculation for describing the determined cells porosity variation. The second REV scale calculation coupled with the variation curve of equivalent porosity is conducted by the program of "dugks2.cpp" in step 4, and the flow field and permeability of at the REV scale is obtained. 
