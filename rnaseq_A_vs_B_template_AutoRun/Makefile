.PHONY: env download refs A B de all clean

env:
	mamba env create -f envs/conda_env.yml || true

download:
	bash scripts/01_download_airway_fastq.sh

refs:
	bash scripts/00_setup_references.sh

A:
	bash scripts/02A_star_featurecounts.sh

B:
	bash scripts/02B_salmon_quant.sh

de:
	Rscript analysis/03_deseq2_compare.R

all: download refs A B de

clean:
	rm -rf results/approachA_star results/approachA_featurecounts results/approachB_salmon results/deseq2
