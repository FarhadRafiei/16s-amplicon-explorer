PYTHON ?= python3

.PHONY: demo test clean

demo:
	$(PYTHON) src/analyze_16s.py --fasta data/demo/reads.fasta --metadata data/demo/metadata.tsv --reference data/reference/toy_16s_references.fasta --output results/demo

test:
	$(PYTHON) -m unittest discover -s tests -v

clean:
	rm -rf results/demo
