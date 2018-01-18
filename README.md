# R codes for Interval Designs
This repository contains the R codes to implement several existing interval designs including CCD, mTPI, BOIN, Keyboard and UMPBI designs
# Description
The primary objective of phase I oncology trials is to identify the maximum tolerated dose (MTD), whose induced dose-limiting toxicity (DLT) probability is the closest to the target toxicity rate. Interval designs have recently attracted enormous attention in phase I clinical trials due to their simplicity and desirable finite-sample performance. Theentireprocedureofthe interval design is guided by comparing the observed toxicity probability (or the number of DLTs) with a prespecified toxicity tolerance interval. We provides R codes to implement several existing interval designs including CCD, mTPI, BOIN, Keyboard and UMPBI designs
# Functions
The repository includes two files:
* get.boundary.R: The file contains ```get.boundary()``` function to generate escalation and de-escalation boundaries for several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.
```rscript
get.boundary(target, ncohort, cohortsize=3,design=1,cutoff.eli)
```
* get.oc.R: The file contains ```get.boundary()``` function to conduct simulation studies and generate operating characteristics for several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.
```rscript
get.oc(target, p.true, ncohort, cohortsize, startdose, design, cutoff.eli, ntrial)
```
