# Troubleshooting

Problems that came up while building this project, with the cause and the fix.

| Symptom | Cause | Fix |
|---|---|---|
| `bash: scripts/05_run_deseq2.sh: No such file or directory` | Running from inside `scripts/`. | Run from the repo root, or use the file name alone. The scripts locate their own folder. |
| `there is no package called 'shiny'` (or DESeq2) | R packages not installed in the active environment. | `conda env create -f environment.yml`, or install the R packages listed there. |
| conda fails with `mkdir: /Volumes/BINF: Permission denied` | An environment was being created inside a path with a space; some Bioconductor data packages run an install script that does not quote paths. | Create the environment in the normal conda location (not inside the project folder). All project scripts quote paths. |
| `Cannot find app_data.RData` when launching the app | Step 5 has not run. | `bash scripts/05_run_deseq2.sh`, then launch from the repo root. |
| `unable to find an inherited method for function 'select'` | DESeq2 (through AnnotationDbi) masked dplyr's `select`. | Use `dplyr::select` explicitly (already done in the app). |
| Step 4 reports a low assigned percentage | Wrong `STRANDEDNESS`, or GTF and genome build mismatch (for example `chr1` vs `1`), or lots of intronic or genomic DNA reads. | Try `STRANDEDNESS=2` then `1`; confirm the GTF matches the index; read the `Unassigned_*` rows. Delete `gene_counts.txt` or run with `FORCE=1` to recount. |
| `Could not match count columns to sample IDs by name` | featureCounts column names do not contain the sample IDs. | Make `samples.tsv` IDs match the BAM file names, or keep the same order in both. |
| `Contrast levels not found` | `condition` values in `samples.tsv` do not match the contrast in `scripts/deseq2.R`. | Edit the `contrast` setting at the top of the script. |
| Warning: `A group has <2 replicates` | Expected for this dataset (1 late vs 2 early). | Treat all results as exploratory. |
| No package for Apple Silicon (osx-arm64) | Some bioconda tools lack ARM builds. | `CONDA_SUBDIR=osx-64 conda env create -f environment.yml`. |
| The laptop slept during a long run | No keep-awake. | `caffeinate -i bash scripts/run_all.sh`. Finished samples are skipped on rerun. |
| Dashboard shows "No pathway results" | Pathway analysis was skipped (for example `msigdbr` missing) or found no sets. | Check the message in the step 5 output; install `msigdbr` and `fgsea`, then rerun step 5. |
