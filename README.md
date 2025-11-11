# TND-Survival

This repository contains the R scripts and SLURM submission files used to perform the simulation study described in our manuscript *Improving Assessment of Vaccine Effectiveness by Coupling Test-Negative Design Studies with Survival Models*.

All simulations were implemented in **R (version ≥ 4.1)** and executed on a **Linux-based high-performance computing (HPC) cluster** using **SLURM** for parallel job management.

---

## File Descriptions

- **000-func-definition.R**  
  Defines all utility functions used in the simulation, including data generation, model fitting, result reading, saving, and summarization.  
  This file must be sourced before running any simulation scripts.

- **011-simu-part-one.R**  
  Runs the first part of the simulation (Constant VE scenario).  
  It loads the function definitions, sets simulation parameters, executes multiple replicates, and saves intermediate outputs.

- **011-simu-part-one.slurm**  
  SLURM batch script for submitting `011-simu-part-one.R` to an HPC cluster.  
  It specifies computational resources (e.g., CPU, memory, wall time) and executes the R script using `R CMD BATCH`.

- **012-simu-part-one-read-results.R**  
  Reads and summarizes the results from Part One simulations.  
  Generates summary tables and figures presented in the manuscript.

- **021-simu-part-two.R**  
  Runs the second part of the simulation (Waning VE scenario).  
  Similar structure to Part One but with modified model assumptions and time-varying vaccine effectiveness.

- **021-simu-part-two.slurm**  
  SLURM submission script for running `021-simu-part-two.R` on the HPC cluster.

- **022-simu-part-two-read-results.R**  
  Reads and summarizes the results from Part Two simulations.  
  Generates summary tables and figures used in the manuscript.

---

## How to Run

### 1. Prepare Environment
Install required R packages before running the scripts:
```r
install.packages(c("tidyverse", "coxme", "survival", "glue"))
```

### 2. Run on a SLURM Cluster (Recommended)
Before submission, specify your HPC account information in the `.slurm` files and set the desired save directory in the `.R` scripts.

Submit jobs using:
```bash
sbatch 011-simu-part-one.slurm
sbatch 021-simu-part-two.slurm
```

### 3. Run Locally (Optional)
To test locally, reduce the number of simulations and iterations within the R scripts.  
Modify the `for` loop parameters and specify a local save directory before execution.

---

## Contact
For questions or additional information, please contact the author:  
**Shangchen Song**  
Email: [s.song@ufl.edu]
