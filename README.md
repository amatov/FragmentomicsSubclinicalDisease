# Utilization of cfDNA fragment size patterns ​for disease detection & classification ​based on low-coverage WGS data 

We consider the relative entropy between cohorts’ cfDNA fragment lengths and test two hypotheses.

1. We can pinpoint particular lengths for which disease differs from healthy.

2. We can identify distinct differences for colorectal (CRC) as well as other cancer types (ovarian, pancreatic, gastric, breast, lung cancer and cholangiocarcinoma).

Preliminary Kullback-Leibler divergence ([PMC5812299](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5812299/)) analysis of the Delfi ([PMC6774252](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6774252/)) data shows:

1. Cancer vs healthy:

- Healthy individuals and cancer patients exhibit differences for
particular fragment lengths (classification of new clinical samples and early detection of disease).
- We measure two to three peaks (see [KLD_CRC_FRL.pdf](https://gitlab.com/amatov/dnafrl/-/blob/master/KLD_CRC_FRL.pdf)) on the divergence histogram (identify the disease stage).

2. Cancer vs cancer:

- CRC patients and other cancers exhibit differences for particular
fragment lengths (identify the tissue of origin).
- At least 8% of the fragments belong to diverging populations (determine the degree of overlap between the regulation of different tumors).
