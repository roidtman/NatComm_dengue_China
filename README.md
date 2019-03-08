# Inter-annual variation in seasonal dengue epidemics driven by multiple interacting factors in Guangzhou, China
## Rachel J. Oidtman, Shengjie Lai, Zhoujie Huang, Juan Yang, Amir S. Siraj, Robert C. Reiner, Andrew J. Tatem, T. Alex Perkins, Hongjie Yu


Here, we give a brief description of code and data used to fit the model, run analyses, and generate output for the Oidtman et al. (2019). 

## Getting started

The R scripts are written assuming you have the following folder structure:

```
DENV_china
│   README.md
└─── main_text_code
└─── output
└─── data
└─── supplemental_code
```
Where all of the MCMC code, analysis code, and processing code is in the 'code' folder. All of the figures generated and resulting .RData files feed to the 'output' folder. All requisite data is in the 'data' folder. All of the R scripts are written assuming you are inside of the 'code' folder. 

### Software and packages

We used R version 3.4.1, "Single Candle". Each script loads requisite libraries. 

\newline R packages necessary for these analyses:

* fda
* VGAM
* mvtnorm
* lubridate
* coda
* mgcv
* BayesianTools
* scam
* parallel
* RColorBrewer

```
install.packages(c('fda', 'VGAM', 'mvtnorm', 'lubridate', 'coda', 'mgcv', 'BayesianTools', 'scam', 'parallel', 'RColorBrewer'))
```

## Analysis

### 1. Fit mosquito curves and estimate transmission coefficient prior

Code to fit mosquito curves and estiamte prior distributions for the transmission coefficient (beta_0) are available in 1a_mosquito_spline_mcmc.R then 1b_beta_surface.R, respectively.

### 2. Maximum likelihood estimate of mosquito curve

Code to fit maximum likelihood estimates of the mosquito curves in 2_mosquito_optim.R. 

### 3. Fit model

We ran several parallel sequential monte carlo chains to estimate model parameters on the Notre Dame Center for Resource Computings servers. Representative code is available in 3_smc_main_text.R. 

### 4. Run factorial experiments

Code to run factorial simulation experiments and produce analyses. 

### 5. Generate main text figures 

Code to generate simulations and main text figures. 

### Supplementary code and figures

All code to run supplementary analyses, models, and figures are in /supplemental_code


## Code authors

* **Rachel J. Oidtman**
* Amir S. Siraj
* T. Alex Perkins

## Acknowledgments

* Corey Chivers (https://github.com/cjbayesian/MHadaptive/blob/master/R/positiveDefinite.R)
