# glucoseRASpodocytes
A mathematical model of glucose-stimulated local renin-angiotensin system in podocytes

[![DOI](https://zenodo.org/badge/94033856.svg)](https://zenodo.org/badge/latestdoi/94033856)

## Overview
The code for the model is provided as is, without guarantees and without support. The corresponding manuscript describing the model is currently under review for publication.

## Glucose RAS in Podocytes Model
### Authors
Minu R. Pilvankar, Michele A. Higgins, and Ashlee N. Ford Versypt, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

## Related Publication for Model Details
M. R. Pilvankar, M. A. Higgins, A. N. Ford Versypt, Mathematical Model for Glucose Dependence of the Local Renin-Angiotensin System in Podocytes, submitted 2017.

### Main files

* param_estimation_Approach1and2.m. Runs the model to estimate the parameters using Approach 1 and Appraoch 2 (described in the manuscript). The simulation results compare the change in concentrations of ANG II with increasing glucose for Approach 1 and different scenarios of Approach 2. 
   Runs the model 

### Supplementary files
 
* .mat files. 
   Needed to run the programs to pass data, parameters, and calculated values.
    
* export_fig.m
   MATLAB package to nicely export figures.
   Download: https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
   Tips for usage: https://github.com/altmany/export_fig/blob/master/README.md

(c) Ashlee N. Ford Versypt, 2017
