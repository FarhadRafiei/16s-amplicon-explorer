# Portfolio case study: 16S Amplicon Explorer

## The question

**Can a transparent computational workflow quantify and compare microbial community composition across environmental sample groups?**

This project demonstrates the complete reasoning chain that I would apply to a small 16S rRNA gene amplicon study: define sample metadata, inspect input quality, derive amplicon sequence variants (ASVs), annotate them against a reference database, calculate diversity, and communicate limits of inference.

## Why this is suitable for bioinformatics graduate applications

It connects molecular biology and computation rather than treating either as a black box. The repository shows that I understand:

- 16S rRNA gene sequencing is a marker-gene assay, not shotgun metagenomics.
- primer region, extraction batch, contamination controls, and sequencing quality affect what may be concluded.
- ASVs are reproducible sequence features; they are not automatically species or strains.
- alpha diversity describes a sample; beta diversity compares samples.
- taxonomy is database- and classifier-dependent and must carry provenance.
- a repeatable result needs parameters, tests, version control, and an explicit data policy.

## What I built

The lightweight Python analyzer is intentionally inspectable. It validates FASTA and sample metadata, applies basic filters, dereplicates exact ASVs, builds a count table, assigns nearest-reference taxonomy by k-mer similarity, and calculates Shannon, Simpson, and Bray–Curtis values. It produces TSV files and a standalone HTML report.

For an actual study, the README provides a QIIME 2 / DADA2 production pathway. The toy analyzer does not claim to replace denoising, chimera removal, contaminant analysis, or a trained classifier; making that distinction is part of responsible analysis.

## How to discuss it in an interview

> “I built a small but reproducible 16S microbiome project to bridge my cellular and molecular biology background with bioinformatics. I implemented the transparent core analysis myself, then documented the validated QIIME 2 workflow I would use for real data. The most important design choice was separating a demonstrable teaching dataset from biological claims: the repository reports the limits of its taxonomy and includes the controls, provenance, and parameter decisions a real study needs.”

## Next research-ready extension

Use a public, appropriately licensed paired-end 16S dataset with a defined hypothesis. Record the accession, sample selection criteria, primer region, reference database release, and all QIIME 2 artifact versions. Add negative controls, compare groups with PERMANOVA only after checking dispersion, and publish a concise methods/results notebook. Do not reuse the synthetic demo data for that work.
