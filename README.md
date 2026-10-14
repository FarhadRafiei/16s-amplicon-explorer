# 16S Amplicon Explorer

A reproducible starter project for exploring bacterial communities from 16S rRNA gene amplicons. It includes a dependency-free, transparent FASTA demonstration pipeline and a production-oriented QIIME 2 workflow outline for raw paired-end reads.

> **Scientific note:** 16S is a DNA gene marker, although it encodes 16S rRNA. The bundled reads and reference catalog are synthetic teaching data. They validate the project mechanics only and must not be used for biological inference or publication.

## What it produces

- filtering summary and exact sequence variants (ASVs)
- feature/count table by sample
- transparent nearest-reference k-mer taxonomy assignments
- alpha-diversity metrics: observed ASVs, Shannon, and Simpson
- standalone HTML report

## Quick start

Python 3.10+ is sufficient; there are no package dependencies.

```bash
make demo
open results/demo/report.html
make test
```

The demo writes outputs to `results/demo/`:

| File | Contents |
| --- | --- |
| `run_summary.tsv` | filtering and ASV totals |
| `feature_table.tsv` | per-sample ASV counts |
| `taxonomy.tsv` | sequence, assignment, and confidence proxy |
| `alpha_diversity.tsv` | sample diversity metrics |
| `report.html` | readable summary |

## Input contract for the demo analyzer

The script accepts an already-demultiplexed FASTA file. Its headers must begin with `sample_id|`, and every sample must occur in the `sample_id` column of a tab-separated metadata file.

```bash
python3 src/analyze_16s.py \
  --fasta path/to/demultiplexed.fasta \
  --metadata path/to/metadata.tsv \
  --reference path/to/reference.fasta \
  --output results/my_run
```

Reference FASTA headers follow `reference_id|taxonomy`, where taxonomy uses semicolon-separated ranks. The included k-mer assignment is deliberately simple and auditable; replace it with a trained, versioned classifier for research results.

## Recommended production workflow (QIIME 2)

For real Illumina data, preserve raw FASTQ files outside Git and record all versions and parameters. A typical workflow is:

1. Prepare a QIIME 2 manifest and metadata file, validate metadata, and import reads.
2. Inspect demultiplexed quality profiles; use DADA2 to trim, denoise, merge pairs, and remove chimeras.
3. Build a rooted phylogeny (`mafft`, `mask`, `fasttree`, `midpoint-root`) and calculate core diversity metrics at a documented sampling depth.
4. Classify representative sequences with a region-matched, versioned SILVA or GTDB classifier; retain the classifier and database provenance.
5. Export feature tables, taxonomy, diversity results, and visualizations. Confirm negative controls and read-depth behavior before interpretation.

Illustrative commands (adjust artifact types and trim positions after inspecting your own data):

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv --output-path demux.qza --input-format PairedEndFastqManifestPhred33V2
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
qiime dada2 denoise-paired --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 240 --p-trunc-len-r 200 \
  --o-table table.qza --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
qiime feature-classifier classify-sklearn --i-classifier silva-region-classifier.qza \
  --i-reads rep-seqs.qza --o-classification taxonomy.qza
```

Use current QIIME 2 documentation and a reference database license appropriate to your study. Never make taxonomic or differential-abundance claims from the bundled toy catalog.

## Reproducibility and governance

- The GitHub workflow runs tests and the demo pipeline on every push.
- Do not commit participant data, API keys, raw sequencing data, or patient metadata.
- For a study, add a `data/README.md` with sample provenance, primer region, extraction kit, negative controls, sequencing platform, database release, and exact software versions.
- Add a license and cite all databases, tools, and original data before sharing results.

## Project layout

```text
data/          synthetic demo data and tiny test reference
src/           analysis implementation
tests/         automated checks
results/       generated analysis outputs (ignored by Git)
.github/       continuous integration
```
