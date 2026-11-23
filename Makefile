# Convenience targets. Run from the repository root. Assumes the conda environment is active.
.PHONY: help env counts deseq2 app test test-shell test-r simulate deploy clean

help:
	@echo "make env         create the conda environment from environment.yml"
	@echo "make deseq2      run step 5 (DESeq2 + GSEA + run QC) on existing counts"
	@echo "make app         launch the Shiny dashboard"
	@echo "make test        run all automated tests"
	@echo "make simulate    t-test vs DESeq2 across replicate numbers (simulation)"
	@echo "make deploy      deploy the dashboard to shinyapps.io"
	@echo "make clean       remove generated results (keeps your input data)"

env:
	conda env create -f environment.yml

counts:
	bash scripts/04_count_featurecounts.sh

deseq2:
	bash scripts/05_run_deseq2.sh

app:
	Rscript -e 'shiny::runApp("app")'

test: test-shell test-r

test-shell:
	bash tests/test_shell_steps.sh

test-r:
	Rscript -e 'testthat::test_dir("tests/testthat")'

deploy:
	Rscript deploy.R

clean:
	rm -rf results app/data/demo_deseq2_data.RData

simulate:
	Rscript analysis/simulate_design_power.R
