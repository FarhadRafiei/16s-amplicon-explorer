#!/usr/bin/env python3
"""Small, dependency-free 16S amplicon demo analysis.

This is designed for transparent demonstrations and teaching. For production
analysis of FASTQ reads, use the documented QIIME 2 workflow in README.md.
"""
import argparse
import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path

VALID = set("ACGTN")


def parse_fasta(path):
    records, header, chunks = [], None, []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def parse_metadata(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "sample_id" not in rows[0]:
        raise ValueError("Metadata must be a tab-separated file with a sample_id column")
    return {row["sample_id"]: row for row in rows}


def kmers(sequence, k=8):
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 0.0


def read_reference(path):
    references = []
    for header, sequence in parse_fasta(path):
        fields = header.split("|", 1)
        references.append((fields[0], fields[1] if len(fields) > 1 else "Unassigned", kmers(sequence)))
    return references


def assign_taxonomy(sequence, references):
    sequence_kmers = kmers(sequence)
    label, taxonomy, score = max(
        ((name, tax, jaccard(sequence_kmers, ref_kmers)) for name, tax, ref_kmers in references),
        key=lambda x: x[2],
    )
    return label, taxonomy if score >= 0.25 else "Unassigned", score


def write_tsv(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Analyze demultiplexed 16S-style FASTA reads.")
    parser.add_argument("--fasta", required=True, help="FASTA headers must start with sample_id|")
    parser.add_argument("--metadata", required=True, help="TSV with sample_id column")
    parser.add_argument("--reference", required=True, help="Reference FASTA: id|taxonomy")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-length", type=int, default=100)
    parser.add_argument("--max-n-fraction", type=float, default=0.02)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    metadata, references = parse_metadata(args.metadata), read_reference(args.reference)
    kept, rejected = [], Counter()
    for header, sequence in parse_fasta(args.fasta):
        sample = header.split("|", 1)[0]
        if sample not in metadata:
            rejected["missing_metadata"] += 1
        elif set(sequence) - VALID:
            rejected["invalid_bases"] += 1
        elif len(sequence) < args.min_length:
            rejected["too_short"] += 1
        elif sequence.count("N") / len(sequence) > args.max_n_fraction:
            rejected["too_many_N"] += 1
        else:
            kept.append((sample, sequence))
    if not kept:
        raise ValueError("No reads passed filtering")

    asv_counts, sample_asvs = Counter(seq for _, seq in kept), defaultdict(Counter)
    for sample, seq in kept:
        sample_asvs[sample][seq] += 1
    ordered_asvs = [seq for seq, _ in asv_counts.most_common()]
    asv_ids = {seq: f"ASV_{i:03d}" for i, seq in enumerate(ordered_asvs, 1)}

    taxonomy_rows = []
    for seq in ordered_asvs:
        match, taxonomy, score = assign_taxonomy(seq, references)
        taxonomy_rows.append({"feature_id": asv_ids[seq], "sequence": seq, "total_count": asv_counts[seq], "reference_match": match, "taxonomy": taxonomy, "kmer_jaccard": f"{score:.3f}"})
    write_tsv(output / "taxonomy.tsv", list(taxonomy_rows[0]), taxonomy_rows)

    feature_rows = []
    for sample in sorted(metadata):
        row = {"sample_id": sample, **{asv_ids[seq]: sample_asvs[sample][seq] for seq in ordered_asvs}}
        feature_rows.append(row)
    write_tsv(output / "feature_table.tsv", ["sample_id"] + [asv_ids[s] for s in ordered_asvs], feature_rows)

    diversity = []
    for sample in sorted(metadata):
        counts = list(sample_asvs[sample].values())
        total = sum(counts)
        proportions = [x / total for x in counts] if total else []
        diversity.append({"sample_id": sample, "group": metadata[sample].get("group", ""), "reads_retained": total, "observed_asvs": len(counts), "shannon": f"{-sum(p * math.log(p) for p in proportions):.3f}" if proportions else "0.000", "simpson": f"{1 - sum(p * p for p in proportions):.3f}" if proportions else "0.000"})
    write_tsv(output / "alpha_diversity.tsv", list(diversity[0]), diversity)

    summary = {"input_reads": len(kept) + sum(rejected.values()), "retained_reads": len(kept), "rejected_reads": sum(rejected.values()), "asvs": len(ordered_asvs), **{f"rejected_{k}": v for k, v in sorted(rejected.items())}}
    write_tsv(output / "run_summary.tsv", list(summary), [summary])
    report(output / "report.html", summary, diversity, taxonomy_rows)


def report(path, summary, diversity, taxonomy):
    def table(rows):
        cols = list(rows[0]) if rows else []
        head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
        body = "".join("<tr>" + "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols) + "</tr>" for r in rows)
        return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    metrics = "".join(f"<li><strong>{html.escape(k.replace('_', ' '))}:</strong> {v}</li>" for k, v in summary.items())
    page = f"""<!doctype html><html><head><meta charset='utf-8'><title>16S analysis report</title><style>body{{font:16px system-ui;margin:2rem;max-width:1100px;color:#17212b}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{border:1px solid #ccd6dd;padding:.45rem;text-align:left}}th{{background:#eaf3f5}}code{{background:#f4f4f4;padding:.15rem}}</style></head><body><h1>16S amplicon analysis report</h1><p>Demo-oriented exact-sequence analysis. Taxonomy assignments are nearest-reference k-mer matches, not a replacement for validated classifiers.</p><h2>Run summary</h2><ul>{metrics}</ul><h2>Alpha diversity</h2>{table(diversity)}<h2>ASV taxonomy</h2>{table(taxonomy)}</body></html>"""
    path.write_text(page, encoding="utf-8")


if __name__ == "__main__":
    main()
