# Dependencies

## R packages

- **caret**, **glmnet**, **ROCR** -- classification and regression.
- **FNN**, **mixtools**, **spatstat**, **transport** -- statistical
  distance/distribution analysis.
- **Rsamtools** -- Bioconductor package for reading BAM files
  (`readFragmentBam.R`); install via `BiocManager::install()`, not
  CRAN.
- **ggplot2**, **ggpubr**, **gplots**, **RColorBrewer** -- plotting.
- **readxl**, **stringr**, **seewave**, **devtools** -- supporting
  utilities.

## Python

No `requirements.txt` is included. `est_rel_entro_HJW.py` and
`kld_crc_test.py` use `numpy` and `pandas`.

## Included data

`genome_bins.bed` and `sample_reference.csv` are included reference
data files used by the analysis scripts.

## Hardcoded paths

`DELFI2_divergence.R` and `DELFI_divergence.R` are large scripts with
many hardcoded absolute paths to the original author's machine and to
a specific HPC cluster (e.g. `~/genomedk/...`). Active (non-commented)
instances are flagged with a `# EDIT:` comment directly above them --
update these before running a script.
