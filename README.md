# Genetic Analysis Pipeline: BAG Construction and WGS Analysis

This repository contains the code for a three-part genetic analysis pipeline:

1.  **BAG Construction**: Building Biological Age Gap models using proteomic data.
2.  **WGS Main Analysis**: Performing Whole Genome Sequencing association analyses using the constructed BAGs.
3.  **WGS Post Analysis**: Conducting downstream analyses including Polygenic Risk Scores (PRS), enrichment, genetic correlation, heritability, and Leave-One-Variant-Out (LOVO) tests.

## 1. System Requirements

### Hardware Requirements

- A standard high-performance computing (HPC) cluster or server with sufficient RAM to handle large genomic datasets (BED/BIM/FAM files).
- Multi-core support is required for parallel processing (configured in shell scripts).

### Software Requirements

#### Operating Systems

- Tested on **Linux** (Ubuntu/CentOS).
- Requires a Bash shell environment.

#### External Dependencies

The following software must be installed separately and linked in the `./software/` directory as referenced in the scripts:

| Software     | Version       | Description                        | Download Link                                          |
| :----------- | :------------ | :--------------------------------- | :----------------------------------------------------- |
| **REGENIE**  | v3.4.1 / v4.0 | Whole genome regression            | [REGENIE](https://rgcgithub.github.io/regenie/)        |
| **PLINK**    | v2.0          | Genomic data management            | [PLINK 2.0](https://www.coggenomics.org/plink/)        |
| **VEP**      | v2.0          | Variant Effect Predictor           | [Ensembl VEP](https://github.com/Ensembl/ensembl-vep/) |
| **MPH**      | Latest        | Mask-based Phenotype Heritability  | [MPH](https://jiang18.github.io/mph/)                  |
| **PRSice-2** | v2.3.5        | Polygenic Risk Score analysis      | [PRSice-2](https://choishingwan.github.io/PRSice/)     |
| **RICE**     | Latest        | Rare variant analysis pipeline     | [RICE](https://github.com/jwilliams10/RICE/)           |
| **GREP**     | v1.0.0        | Genetic Regression Analysis        | [GREP](https://github.com/saorisakaue/GREP)            |
| **GCTA**     | Latest        | Genome-wide Complex Trait Analysis | [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) |
| **R**        | v4.2.0+       | Statistical computing              | [R Project](https://www.r-project.org/)                |

*Note: Ensure the executable paths in the shell scripts (e.g., `./software/regenie`) match your installation directory.*

#### R Package Dependencies

The R scripts require the following packages. You can install them via CRAN, Bioconductor, or Conda.

**Core Packages:**
`dplyr`, `data.table`, `readr`, `caret`, `parallel`, `glmnet`, `doParallel`, `foreach`

**Post-Analysis Packages:**
`clusterProfiler`, `org.Hs.eg.db`, `BiocManager`

**Installation Command (R):**

```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(c("dplyr", "data.table", "readr", "caret", "glmnet", "doParallel", "foreach"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
```

## 2. Installation Guide

### Instructions

1. **Clone the Repository**:
   `git clone <YOUR_REPOSITORY_URL>`
   `cd <REPOSITORY_NAME>`

2. **Install External Tools**:
   Download the software listed in the **System Requirements** section. Create a folder named `software` in the root directory and place the executables there (or update the paths in the `.sh` scripts to point to your installed binaries).

   Example structure:

   - `./software/regenie`
   - `./software/PRSice/PRSice_linux`
   - `./software/gcta`

3. **Setup R Environment**:
   Install the required R packages using the commands provided above.

### Typical Install Time

- **R Packages**: ~5-10 minutes on a normal desktop with internet access.
- **External Tools**: Varies significantly based on download speed and compilation requirements (if building from source). Pre-compiled binaries usually take <10 minutes to download and configure.

## 3. Demo

This pipeline is designed to be run in three sequential stages. Below are the instructions to run a demonstration on a single phenotype (e.g., "Liver") and a single chromosome (e.g., Chromosome 1).

### Instructions to Run

*Prerequisite: Ensure input data files (.bed, .pheno, .covar, etc.) are placed in the expected `./data/` directories as referenced in the scripts.*

**Step 1: BAG Construction**
Run the R script to construct the Biological Age Gap for specific organs.

```shell
cd BAG_construction
Rscript Proteomic_BAG_construction.R
```

**Step 2: WGS Main Analysis**
Run the step 1 (null model) and step 2 (association) scripts.

```shell
cd ../WGS_main_analysis
bash step1.sh Liver
bash step2_genebased.sh 1 Liver LoF
bash step2_single.sh 1 Liver
```

**Step 3: WGS Post Analysis**
Run downstream analyses.

```shell
cd ../WGS_post_analysis
bash PRS.sh Liver
Rscript enrichment.R
bash genetic_correlation.sh Liver 1
bash heritability.sh Liver
bash lovo.sh 1 Liver mask mask_id gene_name 0.01 domain
```

### Expected Output

- **BAG Construction**: CSV files containing predicted biological ages and residuals (BAGs) saved in `./elg_results/`.
- **WGS Main**: REGENIE output files (.log, .regenie, .vcf.gz if applicable) containing association statistics (p-values, effect sizes) in `./result/`.
- **Post Analysis**: Plots and tables for enrichment, PRS curves, and correlation matrices.

### Expected Run Time

- **Demo Run Time**: Highly **dependent on data volume** (number of samples, variants, and chromosomes processed).
  - *Small scale demo (1 chromosome, few samples)*: Minutes to hours.
  - *Full scale analysis*: Hours to days, typically requiring an HPC cluster.

## 4. Instructions for Use

To run the software on your own data:

1. **Prepare Data**: Organize your genotype data (BED/BIM/FAM), phenotype files (.txt with FID/IID/Phenotype), and covariate files according to the formats expected by PLINK/REGENIE. Place them in the `./data/` subdirectories corresponding to each module.
2. **Configure Paths**: Edit the shell scripts (`step1.sh`, `step2_*.sh`, etc.) to match your file paths and desired parameters (e.g., number of threads, memory limits).
3. **Execute Workflow**:
   - Start with `BAG_construction/Proteomic_BAG_construction.R` to generate your traits.
   - Use the generated traits as input for `WGS_main_analysis`.
   - Feed the summary statistics from the main analysis into `WGS_post_analysis` scripts.
4. **Parallelization**: For the WGS analysis, it is recommended to submit `step2_genebased.sh` and `step2_single.sh` as array jobs across chromosomes (1-22) using a job scheduler (SLURM/PBS) to speed up computation.

### Reproduction Instructions

To reproduce the quantitative results in the manuscript:

1. Ensure the exact same version of external software (specifically REGENIE v3.4.1/v4.0 and PRSice-2 v2.3.5) is used.
2. Use the provided mask definitions and set-lists located in the `./data/` directory (if included) or generate them using the specified tools (VEP/MPH).
3. Run the full pipeline sequentially as described in the **Demo** section for all 22 chromosomes and all phenotypes listed in the R scripts (`organ_base_filename` vector).

## License

This project is provided as-is for research purposes. Please cite the associated manuscript if you use this code.
