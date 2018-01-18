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

# Inputs
* ```target```: The target toxicity probability, e.g., ```target<-0.30```.
* ```p.true```: The true toxicity rate for each dose, e.g., ```p.true<-c(0.1,0.2,0.3,0.4,0.5,0.6)```.
* ```ncohort```: The total number of cohorts.
* ```cohortsize```: The cohort size.
* ```startdose```: The starting dose of the trial
* ```design```: The specific interval design. 1 is the CCD design, 2 is the mTPI design, 3 is the BOIN design, 4 is the Keyboard design (or mTPI2 design), 5 is the UMPBI design. Default values for design parameters of these designs are utilized.
* ```cutoff.eli```: The cutoff to eliminate the overly toxic dose for safety monitoring, e.g., ```cutoff.eli<-0.95```.
* ```ntrial```: The number of simulated trials, e.g., ```ntrial<-1000```.

# Example
We take the UMPBI design as an example, i.e., ```design<-5```. Suppose the target toxicity rate is 0.3, the number of cohorts is 12 with three patients in a cohort, and the elimination boundary is set at 0.95. 
* We generate the escalation and de-escalation boundaries for the UMPBI design.
```rscript
 get.boundary(target=0.3,ncohort=12,cohortsize=3,design=5,cutoff.eli=0.95)
```
The output is given by 
```rscript
Number of patients treated 3 6 9 12 15 18 21 24 27 30 33 36
Eliminate if # of DLT >=   3 4 5  7  8  9 10 11 12 14 15 16
Deescalate if # of DLT >=  2 3 4  5  6  7  8  9 10 11 12 13
Escalate if # of DLT <=    0 1 2  2  3  4  5  5  6  7  8  9
```
* If we need to generate the operating characteristics of the UMPBI design based on the scenario ```p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60)```, then we can use the ```get.oc()``` function. 
```rscript
p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60)
get.oc(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,design=5,cutoff.eli=0.95)
```
The output is given by 
```rscript
selection percentage at each dose level (%):
0.2  5.4  27.1  48.5  17.4  1.3  
number of patients treated at each dose level:
4.2  6.1  10.3  10.2  4.5  0.8  
number of toxicity observed at each dose level:
0.3  0.6  2.0  3.1  2.0  0.5  
average number of toxicties: 8.5 
average number of patients: 36.0 
```
# Authors and Reference
* Ruitao Lin (ruitaolin@gmail.com) and Guosheng Yin (gyin@hku.hk)
* The R codes were adapted based on the codes written by Suyu Liu and Ying Yuan.

* Guo W, Wang SJ, Yang S, et al. A Bayesian interval dose-finding design addressing Ockham’s razor: mTPI-2. Contemp Clin Trials 2017; 58: 23–33.
* Ivanova A, Flournoy N, Chung Y. Cumulative cohort design for dose finding. J Stat Plan Inference 2007; 137: 2316–2327.
* Ji Y, Liu P, Li Y, et al. A modified toxicity probability interval method for dose-finding trials. Clin Trials 2010; 7: 653–663.
* Lin R, Yin G. Uniformly most powerful Bayesian interval design for phase I dose-Finding trials.
* Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials. J R Stat Soc Ser C 2015; 64: 507–523.
* Yan F, Mandrekar SJ and Yuan Y. Keyboard: A Novel Bayesian toxicity probability interval design for phase I clinical trials. Clin Cancer Res 2017; 23: 3994–4003



