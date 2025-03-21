# ampliwrangler
[![DOI](https://zenodo.org/badge/947071183.svg)](https://doi.org/10.5281/zenodo.15018024)

Simple command-line utilities for enhancing QIIME2-based amplicon analyses.

Replacement for the [`qiime2-helpers`](https://github.com/jmtsuji/qiime2-helpers) repo.

## About
When working with QIIME2 for amplicon analysis, I often run into steps in the pipeline that I wished could be
enhanced or streamlined. Ampliwrangler exists to fill in some of those gaps to make amplicon analysis smoother. Modules
of ampliwrangler can directly interface with QIIME2 and support QZA format files along with exported TSV/BIOM/FASTA file
types. The hope is that ampliwrangler can be used either during a QIIME2 analysis or post-analysis to improve your
amplicon analysis experience!

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
`ampliwrangler` consists of several loosely related modules that help speed up or enhance QIIME2-based amplicon
analyses:
- tabulate: feature table manipulations, including adding taxonomy and sequence info as "metadata" to the table for ease
  of viewing downstream
- count: get the total read counts for all samples in the feature table (or just the min count) - this is helpful for
  quickly seeing what read count to rarefy to.
- split: split a manifest file by a metadata column value such as sequencing run ID. It's recommended to run DADA2 (a
  sequence denoising tool) on one sequencing run at a time, so it's helpful to be able to split the manifest file this
  way for pre-DADA2 steps.

Basic usage instructions are included below. See full usage instructions for each module in Appendix 2 at the end of
this README.

### ampliwrangler tabulate
Manipulate feature tables, including adding taxonomy and sequence information as "metadata".

Example (test data in repo):
```bash
# Assuming you are in the Github repo directory
input_dir="testing/tabulate/inputs"

ampliwrangler tabulate \
  -f "${input_dir}/feature-table.tsv" \
  -s "${input_dir}/dna-sequences.fasta" \
  -t "${input_dir}/taxonomy.tsv" \
  -o "feature-table-with-metadata.tsv" \
  --parse_taxonomy

# Direct use of QZA files is also supported!
ampliwrangler tabulate \
  -f "${input_dir}/feature-table.qza" \
  -s "${input_dir}/dna-sequences.qza" \
  -t "${input_dir}/taxonomy.qza" \
  -o "feature-table-with-metadata.tsv" \
  --parse_taxonomy
```
The tabuilate module is quite feature rich. Other functions, including normalization of feature table data, renaming
ASVs, sorting by taxonomy, and so on are detailed in the help statement in Appendix 2 below.

### ampliwrangler count
Summarize total counts per sample for a feature table.
This script is nice if trying to decide how to rarefy your data, for example.

Example (test data in repo):
```bash
# Assuming you are in the Github repo directory
input_dir="testing/count/inputs"

ampliwrangler count \
  -i "${input_dir}/feature-table.qza" \
  -c "sample-counts.tsv" \
  -m "min-sample-count.txt"
  
# Or just write the minimum count value (for rarefaction) to stdout
min_sample_count=$(ampliwrangler count -i "${input_dir}/feature-table.qza" -m -)
echo ${min_sample_count} # to see the min sample count
```

### ampliwrangler split
Split a manifest file into multiple sub-files by a metadata column value such as sequencing run ID.

Example (test data in repo):
```bash
# Assuming you are in the Github repo directory
input_dir="testing/split/inputs"

ampliwrangler split \
  -i "${input_dir}/manifest.tsv" \
  -o "manifest-split-dir" \
  -r "run-id"
```

## Citation
If you benefit from using qiime2helpers, please cite the repo in a way similar to the following:
> Tsuji, JM. ampliwrangler: simple command-line utilities for enhancing QIIME2-based amplicon analyses. 2025. 
> Available at https://github.com/jmtsuji/ampliwrangler (doi:https://doi.org/10.5281/zenodo.15018024).

## Appendix 1: testing
Run automated end-to-end tests as follows (in the command line, assuming you are in the git repo dir):

```bash
# ampliwrangler tabulate
testing/test-tabulate.sh testing/tabulate

# ampliwrangler count
testing/test-count.sh testing/count

## ampliwrangler split
testing/test-split.sh testing/split
```

## Appendix 2: full usage instructions for modules
Copied from the command line help messages

ampliwrangler
```commandline
usage: ampliwrangler [-h] [-V] {tabulate,count,split} ...

Ampliwrangler: simple command-line utilities for enhancing QIIME2-based amplicon analyses. Copyright Jackson M. Tsuji, 2025. Version: 0.1.3

positional arguments:
  {tabulate,count,split}
                        Available modules:
    tabulate            Creates a TSV-format QIIME2 feature table with overlaid taxonomy and sequence information.
    count               Creates a tabular summary of the total read counts of all samples in a feature table.
    split               Simple script to split a QIIME2 manifest file into multiple files based on the run ID column.

options:
  -h, --help            show this help message and exit
  -V, --version         Print tool version and exit
```

ampliwrangler tabulate
```commandline
usage: ampliwrangler tabulate [-h] [-lf PATH] [-O] [-v] -f TABLE [-o TABLE] [-s SEQ] [-t TAX] [-n {percent,proportion}] [-S] [-R] [-P] [-u] [-C] [-H NUM] [-i NAME] [-I NAME]

options:
  -h, --help            show this help message and exit

Basic config settings:
  -lf PATH, --logfile PATH
                        Log filepath (default: None)
  -O, --overwrite       Overwrite existing files/directories. By setting this flag, you risk erasing old data.
  -v, --verbose         Enable verbose logging

Input/output file options:
  -f TABLE, --feature_table TABLE
                        The path to the input feature table file. TSV, BIOM, or QZA file types are supported.
  -o TABLE, --output_feature_table TABLE
                        The path to the output TSV feature table. Will write to STDOUT (-) if nothing is provided.
  -s SEQ, --sequences SEQ
                        The path to the input ASV/OTU sequence file. Sequences will be added as the "Sequences" column. You can optionally omit this flag and not have sequences added to the table. FastA and QZA file types are supported.
  -t TAX, --taxonomy TAX
                        The path to the input taxonomy file. Taxonomy will be added as the "Taxonomy" column. You can optionally omit this flag and not have taxonomy added to the table. TSV (qiime2 format) or QZA file types are supported.

Optional table manipulation options:
  -n {percent,proportion}, --normalize {percent,proportion}
                        Optionally normalize each sample by the desired normalization method (e.g., to convert to percent relative abundance).
  -S, --sort_features   Optionally sort Feature IDs roughly based on overall abundance.
  -R, --rename_features
                        Optionally rename the Feature IDs sequentially, roughly based on overall abundance. Automatically sets --sort_features.
  -P, --parse_taxonomy  Optionally parse Silva taxonomy into 7 ranks with columns "Domain", "Phylum", etc.
  -u, --fill_unresolved_taxonomy
                        Optionally add Unresolved_[taxon] labels for blank taxonomy ranks. Requires --parse_taxonomy.
  -C, --sort_columns    Sort feature table columns by sample ID (alphabetically)
  -H NUM, --header_row NUM
                        The row to use as header for the input feature table (e.g., 0 = top row). [Default: "auto" = auto detect best row]
  -i NAME, --original_feature_id_colname NAME
                        The name of the Feature ID column of the input feature table. [Default: "auto" = auto detect]
  -I NAME, --final_feature_id_colname NAME
                        The name of the first column of the output feature table. [Default: "Feature ID"]
```

ampliwrangler count
```commandline
usage: ampliwrangler count [-h] [-lf PATH] [-O] [-v] -i TABLE [-c TABLE] [-m TXT]

options:
  -h, --help            show this help message and exit

Basic config settings:
  -lf PATH, --logfile PATH
                        Log filepath (default: None)
  -O, --overwrite       Overwrite existing files/directories. By setting this flag, you risk erasing old data.
  -v, --verbose         Enable verbose logging

Input/output file options:
  -i TABLE, --input_filepath TABLE
                        The path to the input FeatureTable file, either QZA, BIOM, or TSV.
  -c TABLE, --output_counts_filepath TABLE
                        Optional path to write output counts per sample in TSV format. Provide "-" to write to stdout.
  -m TXT, --output_min_count_filepath TXT
                        Optional path to write a single-line text file with the lowest count value in the dataset. Provide "-" to write to stdout.
```

ampliwrangler split
```commandline
usage: ampliwrangler split [-h] [-lf PATH] [-O] [-v] -i TSV -o DIR [-r STR]

options:
  -h, --help            show this help message and exit

Basic config settings:
  -lf PATH, --logfile PATH
                        Log filepath (default: None)
  -O, --overwrite       Overwrite existing files/directories. By setting this flag, you risk erasing old data.
  -v, --verbose         Enable verbose logging

Input/output file options:
  -i TSV, --input_filepath TSV
                        The path to the input manifest file. Must match QIIME standards AND have a column with a unique ID for each sequencing run as specified in --run_id_column
  -o DIR, --output_dir DIR
                        The directory where output files (named "manifest_[RUN_ID].tsv") will be saved.

Other params:
  -r STR, --run_id_column STR
                        Name of the column in the input manifest file that contains the unique IDs for each run. [default: run_ID]
```
