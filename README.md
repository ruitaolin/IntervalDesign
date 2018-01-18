# R codes for Interval Designs
This repository contains the R codes to implement several existing interval designs including CCD, mTPI, BOIN, Keyboard and UMPBI designs

# Description
The primary objective of phase I oncology trials is to identify the maximum tolerated dose (MTD), whose induced dose-limiting toxicity (DLT) probability is the closest to the target toxicity rate. Interval designs have recently attracted enormous attention in phase I clinical trials due to their simplicity and desirable finite-sample performance. Theentireprocedureofthe interval design is guided by comparing the observed toxicity probability (or the number of DLTs) with a prespecified toxicity tolerance interval. We provides R codes to implement several existing interval designs including CCD, mTPI, BOIN, Keyboard and UMPBI designs

# Functions
The repository includes two files:
* get.boundary.R: The file contains ```get.boundary()``` function to generate escalation and de-escalation boundaries for several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.
```rscript
get.boundary(target, ncohort, cohortsize,design,cutoff.eli)
```
* get.oc.R: The file contains ```get.boundary()``` function to conduct simulation studies and generate operating characteristics for several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.
```rscript
get.oc(target, p.true, ncohort, cohortsize, startdose, design, cutoff.eli, ntrial)
```

#Inputs
* ```target```: The target toxicity probability, e.g., ```target<-0.30```.
* ```p.true```: The true toxicity rate for each dose, e.g., ```p.true<-c(0.1,0.2,0.3,0.4,0.5,0.6)```.
* ```ncohort```: The total number of cohorts.
* ```cohortsize```: The cohort size.
* ```startdose```: The starting dose of the trial
* ```design```: The specific interval design. 1 is the CCD design, 2 is the mTPI design, 3 is the BOIN design, 4 is the Keyboard design (or mTPI2 design), 5 is the UMPBI design. Default values for design parameters of these designs are utilized.
* ```cutoff.eli```: The cutoff to eliminate the overly toxic dose for safety monitoring ```cutoff.eli<-0.95```.
* ```ntrial```: The number of simulated trials ```ntrial<-1000```.

