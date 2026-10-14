# Reproducible 16S rRNA Amplicon Analysis Workflow

[![Validate demo analysis](https://github.com/FarhadRafiei/16s-amplicon-explorer/actions/workflows/ci.yml/badge.svg)](https://github.com/FarhadRafiei/16s-amplicon-explorer/actions/workflows/ci.yml)

This repository provides a compact, reproducible workflow for exploratory analysis of 16S rRNA gene amplicon data. It includes an executable demonstration pipeline, automated tests, a container specification, and a documented production workflow for QIIME 2 and DADA2.

The project is intended to demonstrate transparent data handling and reporting. It separates a small teaching implementation from the validated tools required for analysis of real sequencing data.

## Scope and limitations

The files in `data/demo/` and `data/reference/` are synthetic. They are included solely to test the workflow and illustrate its outputs. They must not be used to draw biological conclusions.

The Python implementation performs basic validation, exact-sequence dereplication, k-mer reference matching, and diversity calculations. It is not a substitute for read denoising, chimera detection, contaminant assessment, or a region-specific taxonomic classifier.

## Workflow overview

1. Validate sequence characters, sequence length, ambiguous-base content, and sample metadata.
2. Dereplicate retained reads into exact amplicon sequence variants (ASVs).
3. Generate an ASV count table by sample.
4. Assign provisional taxonomy using nearest-reference 8-mer Jaccard similarity.
5. Calculate alpha diversity (observed ASVs, Shannon, Simpson) and Bray–Curtis dissimilarity.
6. Write tabular outputs and a standalone HTML report.

## Requirements

- Python 3.10 or later
- `make` (optional, for convenience commands)
- Docker (optional, for containerized execution)

No third-party Python packages are required for the demonstration workflow.

## Quick start

```bash
git clone https://github.com/FarhadRafiei/16s-amplicon-explorer.git
cd 16s-amplicon-explorer
make test
make demo
open results/demo/report.html
```

## Inputs

The demonstrator accepts a demultiplexed FASTA file, a tab-separated metadata file, and a reference FASTA file.

### Read headers

Each sequence header must begin with a sample identifier followed by a vertical bar:

```text
>sample_01|read_001
ACGT...
```

### Metadata

The metadata table must be tab-separated and contain a `sample_id` column. Additional columns, such as treatment group or sampling site, are retained as metadata.

```text
sample_id	group
sample_01	control
sample_02	treatment
```

### Reference catalog

Reference sequence headers use the following format:

```text
>reference_id|d__Bacteria;p__...;g__...;s__...
ACGT...
```

## Running a custom analysis

```bash
python3 src/analyze_16s.py \
  --fasta path/to/demultiplexed.fasta \
  --metadata path/to/metadata.tsv \
  --reference path/to/reference.fasta \
  --output results/my_analysis
```

| Output | Description |
| --- | --- |
| `run_summary.tsv` | Number of retained and rejected reads; ASV count |
| `feature_table.tsv` | ASV abundance by sample |
| `taxonomy.tsv` | Provisional taxonomic assignments and similarity scores |
| `alpha_diversity.tsv` | Observed ASVs, Shannon and Simpson indices |
| `beta_diversity_bray_curtis.tsv` | Pairwise Bray–Curtis dissimilarities |
| `report.html` | Self-contained analysis summary |

## Containerized execution

```bash
docker build -t 16s-amplicon-analysis .
docker run --rm -v "$(pwd):/project" 16s-amplicon-analysis \
  --fasta data/demo/reads.fasta \
  --metadata data/demo/metadata.tsv \
  --reference data/reference/toy_16s_references.fasta \
  --output results/docker-demo
```

## Recommended workflow for real data

For paired-end Illumina reads, use QIIME 2 and DADA2 with parameters selected after quality inspection. A typical analysis includes:

1. Define metadata, inclusion criteria, controls, and the study hypothesis.
2. Import reads and inspect per-base quality profiles.
3. Trim and denoise with DADA2; merge pairs and remove chimeras.
4. Construct a phylogeny and calculate diversity metrics at a documented sampling depth.
5. Assign taxonomy using a versioned, primer-region-matched SILVA or GTDB classifier.
6. Examine negative controls, library sizes, and group dispersion before hypothesis testing.
7. Record software versions, reference database release, and all analysis parameters.

Illustrative commands only—truncation lengths and classifier choice must be determined from the study data:

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 240 --p-trunc-len-r 200 \
  --o-table table.qza --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
qiime feature-classifier classify-sklearn \
  --i-classifier region-matched-classifier.qza \
  --i-reads rep-seqs.qza --o-classification taxonomy.qza
```

## Reproducibility and data stewardship

- Do not commit raw participant data, identifiers, consent-restricted metadata, or credentials.
- Record data accessions, checksums, primer region, extraction and sequencing batches, controls, tool versions, and parameters.
- Interpret 16S results at the taxonomic resolution supported by the amplicon region and classifier.
- Cite the original dataset, analysis software, and reference database in derived work.

## Repository structure

```text
data/                 Synthetic demonstration inputs
docs/                 Study-planning and portfolio notes
src/                  Analysis implementation
tests/                Automated tests
results/              Generated output; excluded from version control
.github/workflows/    Continuous integration
```

## License and citation

This project is released under the MIT License. Update `CITATION.cff` with the author and repository details before citing or reusing the work.
