## FragmentomicsSubclinicalDisease

Utilization of cfDNA fragment size patterns ​for disease detection & classification ​based on low-coverage WGS data.

## Quick start

This repository contains R and Python scripts analyzing cell-free DNA
fragment length distributions for disease detection and
classification. See [DEPENDENCIES.md](DEPENDENCIES.md) for the
required R and Python packages.

## Repository contents

- `DELFI_divergence.R`, `DELFI2_divergence.R` -- the core
  Kullback-Leibler (KL) divergence fragment-length analysis.
- `normalizeFRL.R`, `readFragmentBam.R`, `regressionDelfi.R`,
  `testBoxplot.R` -- supporting analysis scripts.
- `kld_crc_test.py`, `est_rel_entro_HJW.py` -- a Python port of the KL
  divergence analysis, using the third-party HJW estimator (see LICENSE).
- `kerasMNISTconvnet.R`, `testMNISTconvnet.R`, `mnist_convet.py` --
  MNIST ConvNet examples.
- `genome_bins.bed`, `sample_reference.csv` -- included reference
  data.
- `paper/` -- the LaTeX source and compiled PDF of the associated
  paper.
- `results/` -- an insert-size histogram.
- **License:** see [LICENSE](LICENSE) -- research/educational use.

## About

We consider the relative entropy between cohorts’ cfDNA fragment lengths and test two hypotheses.

1. We can pinpoint particular lengths for which disease differs from healthy.

2. We can identify distinct differences for colorectal (CRC) as well as other cancer types (ovarian, pancreatic, gastric, breast, lung cancer and cholangiocarcinoma).

Preliminary KL divergence analysis of the Delfi data shows:

1. Cancer vs healthy:

- Healthy individuals and cancer patients exhibit differences for
particular fragment lengths (classification of new clinical samples and early detection of disease).
- We measure two to three peaks on the divergence histogram (identify the disease stage).

2. Cancer vs cancer:

- CRC patients and other cancers exhibit differences for particular
fragment lengths (identify the tissue of origin).
- At least 8% of the fragments belong to diverging populations (determine the degree of overlap between the regulation of different tumors).

For detailed information, see: https://www.researchgate.net/publication/382382448_Analysis_of_Genome-Wide_Cell-Free_DNA_Fragment_Length_Distributions_in_Colorectal_Cancer

