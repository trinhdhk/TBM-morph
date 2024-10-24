# TBM-morph: Morphometric analysis for *Tuberculous meningitis*

This repo hosts the code for the morphometric analysis for changes in sub-cortical volumetry over time by treatment arms.

## Pre-requistes

R >= 4.4.0
RStan 2.32.6 (Stan 2.32.2)
cmdstanr

## Structure

renv.lock contains project package information

R/ contains cdoe for the project:

  - 00_volumetry_extract.R: extract volumetry information from xml files extracted by NiftyReg
  - 01_clinical_extract.R; extract clinical outcome and stratum variables.
  - 02_prepared_data.R: use the extracted information to create long dataset etc.
  - model/: model code, icv_ for CSF analysis and RV for regional subcortical anlysis, \_Normal, \_T, and \_GT for Gaussian, Student, and Generalised Student families
  - loo/: model validation
  - infer/: inference code, ploting etc.
  
## Status:

Under construction

(c) Trinh Dong, 2024
