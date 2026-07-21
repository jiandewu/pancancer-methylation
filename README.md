# TCGA Pan-Cancer Methylation Analysis

This repository contains historical workflows for comparing TCGA
HumanMethylation450 profiles and evaluating methylation associations with
clinical and survival outcomes.

The principal multi-cancer workflow targets BRCA, COAD, LUAD, and PRAD. The
repository also retains earlier paired-sample, TNBC-focused, and methylation
lollipop-plot analyses.

## Repository layout

```text
analysis/
├── multi_cancer/  # numbered multi-cancer DMP, DMR, survival, and plotting steps
├── legacy/        # earlier paired-sample methylation-survival scripts
├── tnbc/          # TNBC-specific cohort preparation scripts
└── lollipop/      # CpG and mutation lollipop visualizations
```

Private reports and manuscript drafts have been removed from the complete Git
history. Downloaded data, participant-level clinical tables, RData objects, and
generated figures are not included.

## Environment

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate tcga-pancancer-methylation
```

On platforms where Bioconda cannot resolve `ChAMP` with its historical
`DNAcopy` dependency, install the current Bioconductor release after activating
the environment:

```bash
R -q -e 'BiocManager::install("ChAMP", ask = FALSE, update = FALSE)'
```

## Multi-cancer workflow

Run scripts from the repository's `analysis` directory and execute the six
`multi_cancer/MethySurvival_step*.R` files in numerical order.

Configure the first step with environment variables:

```bash
cd analysis
export TCGA_CANCER=BRCA
export ANALYSIS_COHORT=BRCA
export PANCANCER_WORK_DIR=./multi_cancer
```

Supported historical cancer selections include `BRCA`, `COAD`, `LUAD`, and
`PRAD`. For the TNBC subset, use `TCGA_CANCER=BRCA` and
`ANALYSIS_COHORT=TNBC`; this requires an approved private cohort table at
`data/private/caseID_clinical.txt`.

Expected inputs include UCSC Xena HumanMethylation450 matrices, TCGA clinical
tables obtained through `TCGAbiolinks`, and curated survival tables. The
scripts create cohort-specific `Rdata/` and `figure/` directories.

The Venn-diagram and image-arrangement utilities consume outputs from several
completed cancer-specific runs.

## Reference inputs

The lollipop scripts expect external resources under `analysis/references/`
and `analysis/data/external/`, including gene ranges, the 450K annotation, and
mutation calls. These resources are intentionally excluded and should only be
added when their redistribution terms permit it.

## Data policy

Do not commit TCGA participant barcodes, clinical exports, methylation
matrices, survival tables, RData workspaces, or unpublished reports. Keep
controlled and participant-level files in approved storage.

This repository contains analysis code only. Review the applicable TCGA, GDC,
and UCSC Xena access and redistribution terms before obtaining or sharing data.
