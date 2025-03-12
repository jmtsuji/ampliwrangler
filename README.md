# ampliwrangler
Simple command-line utilities for enhancing QIIME2-based amplicon analyses.
Replacement for the [`qiime2-helpers`](https://github.com/jmtsuji/qiime2-helpers) repo

## Installation
`ampliwrangler` can be installed on Linux or MacOS as follows:

Method 1 (make a separate conda env) - requires miniforge/miniconda
```bash
git clone https://github.com/jmtsuji/ampliwrangler.git

conda env create -n ampliwrangler --file=ampliwrangler/environment.yml
conda activate ampliwrangler

cd ampliwrangler
pip install --editable .

# See available commands
ampliwrangler -h
```

Method 2 (just install via pip) - requires python 3. Can be done inside a conda env if desired.
```bash
git clone https://github.com/jmtsuji/ampliwrangler.git
cd ampliwrangler
pip install --editable .

# See available commands
ampliwrangler -h
```

## Included modules
See full usage instructions for each module in the Appendix at the end of this README.

### `ampliwrangler tabulate`
Adds taxonomy and representative sequence information onto a Feature Table.

Example (test data in repo):
```bash
# Assuming you are in the Github repo directory
input_dir="testing/generate_combined_feature_table/inputs"

ampliwrangler tabulate \
  -f "${input_dir}/feature_table.tsv" \
  -s "${input_dir}/representative_seqs.fasta" \
  -t "${input_dir}/taxonomy.tsv" \
  -o "test_table.tsv" \
  --parse_taxonomy
```

### `ampliwrangler count`
Summarizes total counts per sample for a QZA `FeatureCounts[Frequency]` archive. 
This script is nice if trying to decide how to rarefy your data, for example.

See full usage instructions at the end of the README.

## Testing
Run an automated end-to-end test of `ampliwrangler tabulate` via:
```bash
testing/test-tabulate.sh \
  testing/tabulate
```

## Citation
If you benefit from using qiime2helpers, please cite the repo in a way similar to the following:
> Tsuji, JM. ampliwrangler: simple command-line utilities for enhancing QIIME2-based amplicon analyses. 2025. https://github.com/jmtsuji/ampliwrangler

## Appendix: full usage instructions for modules
Copied from the command line help messages

`ampliwrangler tabulate`
```
TODO
```

`ampliwrangler count`
```
TODO
```
