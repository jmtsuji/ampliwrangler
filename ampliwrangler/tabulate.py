#!/usr/bin/env python
"""
ampliwrangler tabulate
Description: Creates a TSV-format QIIME2 feature table with overlaid taxonomy and sequence information
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import sys
import os
import logging
import argparse
import re

import pandas as pd
from ampliwrangler.utils import check_output_file
from ampliwrangler.load import load_feature_table, load_taxonomy, load_sequence_table

# GLOBAL VARIABLES
# TODO - move these somewhere central
TAXONOMY_RANKS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

logger = logging.getLogger(__name__)


def main(args):
    """
    Runs the workflow based on the provided command line arguments.
    """

    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
        logging.getLogger('ampliwrangler.tabulate').setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logging.getLogger('ampliwrangler.tabulate').setLevel(logging.INFO)

    # Check that taxonomy_filepath is set if parse_taxonomy is True
    if (args.parse_taxonomy is True) and (args.taxonomy is False):
        error = RuntimeError('Although --parse_taxonomy was called, no --taxonomy_filepath was supplied.')
        logger.error(error)
        raise error

    check_output_file(args.output_feature_table, overwrite=args.overwrite)

    # Startup messages
    logger.info(f'Running {os.path.basename(sys.argv[0])} tabulate')
    logger.debug('### SETTINGS ###')
    logger.debug(f'Feature table filepath: {args.feature_table}')
    logger.debug(f'Representative sequences filepath: {args.sequences}')
    logger.debug(f'Taxonomy filepath: {args.taxonomy}')
    logger.debug(f'Output table filepath: {args.output_feature_table}')
    logger.debug(f'Normalization method: {args.normalize}')
    logger.debug(f'Sort Feature IDs roughly by relative abundance?: {args.sort_features}')
    logger.debug(f'Rename Feature IDs sequentially?: {args.rename_features}')
    logger.debug(f'Parse Silva taxonomy into 7 ranks?: {args.parse_taxonomy}')
    logger.debug(f'Fill in blank taxonomy ranks with Unresolved_[taxon]?: {args.fill_unresolved_taxonomy}')
    logger.debug(f'Sort table columns alphabetically?: {args.sort_columns}')
    logger.debug(f'Header row detection for feature table: {args.header_row}')
    logger.debug(f'Feature ID column name detection for input feature table: {args.original_feature_id_colname}')
    logger.debug(f'Feature ID column name for the final output table: {args.final_feature_id_colname}')
    logger.debug(f'Overwrite existing files: {args.overwrite}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('################')

    feature_table = generate_combined_feature_table(feature_table_filepath=args.feature_table,
                                                    sequence_filepath=args.sequences,
                                                    taxonomy_filepath=args.taxonomy,
                                                    normalization_method=args.normalize,
                                                    sort_features=args.sort_features,
                                                    rename_features=args.rename_features,
                                                    parse_taxonomy=args.parse_taxonomy,
                                                    fill_unresolved_taxonomy=args.fill_unresolved_taxonomy,
                                                    sort_columns=args.sort_columns,
                                                    header_row=args.header_row,
                                                    original_feature_id_colname=args.original_feature_id_colname,
                                                    final_feature_id_colname=args.final_feature_id_colname)
    # Write output
    if args.output_feature_table == '-':
        # Write to STDOUT
        logger.info("Writing merged table to STDOUT")
        feature_table.to_csv(sys.stdout, sep='\t', index=False)
    else:
        logger.info(f'Writing merged table to {args.output_feature_table}')
        feature_table.to_csv(args.output_feature_table, sep='\t', index=False)

    logger.info(f'{os.path.basename(sys.argv[0])} tabulate: done.')


def mask_metadata_columns(feature_table: pd.DataFrame, metadata_columns: list = 'default') -> pd.DataFrame:
    """
    Set any detected metadata (non-count / Feature ID) columns as indices

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :param metadata_columns: list of column names to consider as metadata. If 'default' (string), then will
                             search for ['Feature ID', 'Taxonomy', 'Sequence'] + TAXONOMY_RANKS
    """

    default_metadata_cols = ['Feature ID', 'Taxonomy', 'Sequence'] + TAXONOMY_RANKS
    if metadata_columns == 'default':
        metadata_columns = default_metadata_cols

    # A copy of the feature table is made just in case pandas changes the main table object when these edits are made
    feature_table_masked = feature_table.copy(deep=True)

    mask_columns = []
    for metadata_column in metadata_columns:
        if metadata_column in feature_table_masked.columns:
            mask_columns.append(metadata_column)

    if len(mask_columns) > 0:
        logger.debug(f'Masking metadata columns: {mask_columns}')
        feature_table_masked = feature_table_masked.set_index(mask_columns)

    return feature_table_masked


def sort_feature_table_columns(feature_table: pd.DataFrame, sort_sample_columns: bool = False,
                               metadata_columns: list = 'default') -> pd.DataFrame:
    """
    Sort order of columns of feature table

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :param sort_sample_columns: boolean of whether to sort sample columns alphabetically or leave in their current order
    :param metadata_columns: list of column names to consider as metadata. If 'default' (string), then will
                             search for ['Taxonomy'] + TAXONOMY_RANKS + ['Sequence']
    """
    default_metadata_cols = ['Taxonomy'] + TAXONOMY_RANKS + ['Sequence']
    if metadata_columns == 'default':
        metadata_columns = default_metadata_cols

    available_metadata_columns = []
    sample_columns = []
    for feature_table_column in feature_table.columns:
        if feature_table_column in metadata_columns:
            available_metadata_columns.append(feature_table_column)
        else:
            # TODO - don't hard-code Feature ID but find a way to deal with this more programmatically
            if feature_table_column != 'Feature ID':
                sample_columns.append(feature_table_column)

    # Make the above into pandas Series and sort
    metadata_columns = pd.Series(metadata_columns)
    available_metadata_columns_sorted = metadata_columns[metadata_columns.isin(available_metadata_columns)]
    sample_columns = pd.Series(sample_columns)
    if sort_sample_columns:
        logger.debug('Sorting sample column names alphabetically')
        sample_columns = sample_columns.sort_values()

    column_order = ['Feature ID'] + sample_columns.to_list() + available_metadata_columns_sorted.to_list()
    logger.debug(f'Sorting columns in final order: {column_order}')
    feature_table = feature_table[column_order]

    return feature_table


def normalize_feature_table(feature_table: pd.DataFrame, normalization_method='percent') -> pd.DataFrame:
    """
    Normalize a feature table (e.g., of counts) per sample.

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :param normalization_method: method for normalization. Can choose 'percent' or 'proportion' currently.
    :return: QIIME2 FeatureTable[Frequency] artifact (pandas DataFrame) with normalized feature abundances
    """
    normalization_methods = ['percent', 'proportion']
    if normalization_method not in normalization_methods:
        error = ValueError(f'Input normalization method "{normalization_method}" is not one of the available '
                           f'normalization methods: "{",".join(normalization_methods)}".')
        logger.error(error)
        raise error

    logger.debug(f'Normalizing feature table using method: {normalization_method}')
    feature_table_normalized = mask_metadata_columns(feature_table)
    total_per_sample = feature_table_normalized.sum(axis=0)

    # Any samples with zero total count cannot be normalized properly, so they must be removed
    zero_count_sample_names = []
    for sample_name, read_total in total_per_sample.items():
        if read_total == 0:
            zero_count_sample_names.append(sample_name)

    if len(zero_count_sample_names) > 0:
        logger.info(f'Removing {len(zero_count_sample_names)} samples during normalization that had zero counts: '
                    f'{", ".join(zero_count_sample_names)}')
        feature_table_normalized = feature_table_normalized.drop(columns=zero_count_sample_names)
        total_per_sample = total_per_sample.drop(labels=zero_count_sample_names)

    # Perform the normalization
    if normalization_method == 'proportion':
        feature_table_normalized = feature_table_normalized.div(total_per_sample)
    elif normalization_method == 'percent':
        feature_table_normalized = feature_table_normalized.div(total_per_sample).multiply(100)
    else:
        # This in theory should never be called, because the error should be caught above.
        raise ValueError()

    feature_table_normalized = feature_table_normalized.reset_index()

    return feature_table_normalized


def add_taxonomy_to_feature_table(feature_table: pd.DataFrame, taxonomy_filepath: str) -> pd.DataFrame:
    """
    Adds taxonomy values as the Taxonomy column to a QIIME2 feature table

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :param taxonomy_filepath: Path to the taxonomy.tsv file output by the QIIME2 classifier
    :return: QIIME2 FeatureTable[Frequency] artifact with taxonomy in the Taxonomy column
    """
    # Check if Taxonomy column already exists
    if 'Taxonomy' in feature_table.columns:
        error = ValueError('"Taxonomy" column already exists in provided feature table. Cannot add taxonomy.')
        logger.error(error)
        raise error

    # Load taxonomy file
    logger.debug('Loading taxonomy file and adding taxonomy info')
    taxonomy_table = load_taxonomy(taxonomy_filepath)\
        .rename(columns={'Taxon': 'Taxonomy'})
    taxonomy_table = taxonomy_table[['Feature ID', 'Taxonomy']]

    # Merge
    feature_table = pd.merge(feature_table, taxonomy_table, how='left', on='Feature ID', validate='1:1')

    return feature_table


def add_sequences_to_feature_table(feature_table: pd.DataFrame, seq_filepath: str) -> pd.DataFrame:
    """
    Adds ASV/OTU sequences as the Sequence column to a QIIME2 feature table

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :param seq_filepath: Path to the dna-sequences.fasta file output by the QIIME2 denoising/clustering step
    :return: QIIME2 FeatureTable[Frequency] artifact with representative sequences in the ReprSequences column
    """
    # Check if ReprSequence column already exists
    if 'Sequence' in feature_table.columns:
        error = ValueError('"Sequence" column already exists in provided feature table. Cannot add representative '
                           'sequences.')
        logger.error(error)
        raise error

    logger.debug('Loading and adding representative sequences')
    seq_table = load_sequence_table(seq_filepath)

    # Merge
    feature_table = pd.merge(feature_table, seq_table, how='left', on='Feature ID', validate='1:1')

    return feature_table


def sort_feature_table(feature_table: pd.DataFrame) -> pd.DataFrame:
    """
    Roughly sorts a QIIME2 feature table by abundances of Features. The sort is performed based on the maximum
    count (or percent abundance) per feature.

    :param feature_table: QIIME2 FeatureTable[Frequency] artifact loaded as a pandas DataFrame
    :return: sorted QIIME2 FeatureTable[Frequency] artifact
    """
    # Choose a column name to store the sort info; name should not duplicate another column name in the table
    sort_column_id = 'sort'
    sort_column_id_checked = False
    while sort_column_id_checked is False:
        if sort_column_id in feature_table.columns:
            sort_column_id = sort_column_id + '-tmp'
        else:
            sort_column_id_checked = True

    logger.debug('Sorting feature table by maximum count of each feature')
    feature_table = mask_metadata_columns(feature_table)
    # Get max value per feature
    feature_table[sort_column_id] = feature_table.max(axis=1)
    feature_table = feature_table.reset_index()

    feature_table = feature_table.sort_values(by=[sort_column_id, 'Feature ID'], ascending=[False, True])\
        .drop(columns=sort_column_id)

    return feature_table


def parse_silva_taxonomy_entry(taxonomy_entry: str, resolve: bool = True) -> list:
    """
    Parse a single Silva taxonomy entry

    :param taxonomy_entry: Silva taxonomy string (see example below)
    :param resolve: Boolean of whether to resolve blank taxonomy entries with Unresolved_[last_classified_rank]
    :return: list of 7 rank entries for the taxonomy
    """
    # Example:
    # d__Bacteria; p__Margulisbacteria; c__Margulisbacteria; o__Margulisbacteria; f__Margulisbacteria;
    # g__Margulisbacteria; s__microbial_mat

    taxonomy_split = str(taxonomy_entry).split(sep='; ')

    if len(taxonomy_split) > 7:
        error = ValueError(f'Taxonomy entry is {len(taxonomy_split)} long, not 7 as expected. Full entry: '
                           f'"{taxonomy_entry}".')
        logger.error(error)
        raise error
    elif len(taxonomy_split) < 7:
        # Sometimes Silva entries are short; rest is unresolved
        for entry_index in range(len(taxonomy_split) + 1, 8):
            taxonomy_split.append('')

    # Remove taxonomy prefixes
    # TODO - confirm they are in the right order (d,p,c,o,f,g,s)
    taxonomy_split = [re.sub("[dpcofgs]__", "", level) for level in taxonomy_split]

    # Fill in empty parts, if they exist
    if '' in taxonomy_split and resolve is True:
        taxonomy_split = resolve_taxonomy_entry(taxonomy_split)

    return taxonomy_split


def resolve_taxonomy_entry(taxonomy_split: list) -> list:
    """
    Fills in blank taxonomy entries in a list of taxonomy entries

    :param taxonomy_split: list of each taxonomy rank entry for a feature, with rank prefix removed. The list must be
                           sorted so that higher taxonomy entries (e.g., Domain) are before lower entries (e.g, Genus).
    :return: taxonomy_split but with empty entries filled in with Unresolved_[last_classified_rank]
    """
    # Get the positions of the blank entries
    empty_taxa = []
    for taxonomy_level in taxonomy_split:
        if taxonomy_level == '':
            empty_taxa.append(True)
        else:
            empty_taxa.append(False)

    # Get the numeric index of the first empty taxon
    # See https://stackoverflow.com/a/9542768, accessed Sept. 18, 2019
    first_empty_taxon = empty_taxa.index(True)

    if False in empty_taxa[first_empty_taxon:]:
        error = ValueError(f'There seems to be an empty entry in the middle of your taxonomy levels. Cannot resolve. '
                           f'Full taxonomy: "{",".join(taxonomy_split)}".')
        logger.error(error)
        raise error

    filler_entry = f'Unresolved_{taxonomy_split[(first_empty_taxon - 1)]}'
    for taxonomy_level_index in range(first_empty_taxon, 7):
        taxonomy_split[taxonomy_level_index] = filler_entry

    return taxonomy_split


def generate_combined_feature_table(feature_table_filepath: str, sequence_filepath: str, taxonomy_filepath: str,
                                    normalization_method: str = None, sort_features: bool = False,
                                    rename_features: bool = False, parse_taxonomy: bool = False,
                                    fill_unresolved_taxonomy: bool = False, sort_columns: bool = False,
                                    header_row='auto', original_feature_id_colname: str = 'auto',
                                    final_feature_id_colname: str = 'Feature ID') -> pd.DataFrame:
    """
    Loads and parses a feature table, along with optional sequence and taxonomy info, to generate a combined feature
    table. Optionally parses taxonomy, sorts features, and renames features.

    :param feature_table_filepath: path to the TSV-format feature-table.tsv output by QIIME2; a FeatureTable[Frequency]
                                   artifact exported as TSV.
    :param sequence_filepath: path to the dna-sequences.fasta output by QIIME2.
    :param taxonomy_filepath: path to the taxonomy.tsv output by QIIME2.
    :param normalization_method: set to 'percent' or 'proportion' to normalize the table to percent or proportion
                                 relative abundance, per sample. If None, no normalization is performed.
    :param sort_features: whether to sort features roughly based on abundances.
    :param rename_features: whether to rename features roughly based on abundance rank. Sort will also be performed.
    :param parse_taxonomy: whether to parse taxonomy into 7-rank taxonomy columns.
    :param fill_unresolved_taxonomy: whether to fill in blank taxonomy ranks with Unresolved_[taxon]. Requires
                                     parse_taxonomy to be True.
    :param sort_columns: whether to sort columns alphabetically by sample ID
    :param header_row: for TSV format files, manually specify the header row here, as per the header argument in
                        pd.read_csv (e.g., by integer, where 0 is the top row), or specify 'auto' to try to auto-detect
                        the correct header row.
    :param original_feature_id_colname: for TSV format files, optionally specify the name of the feature ID column in
                                        the provided table. Must be the first column in the table. Auto-detection is
                                        attempted if 'auto' is provided.
    :param final_feature_id_colname: name to give the column where Feature IDs are stored in the final output table.
    :return: a combined feature table as a pandas DataFrame.
    """
    if rename_features is True:
        logger.debug('Changing sort_features to True because rename_features is True.')
        sort_features = True

    if fill_unresolved_taxonomy is True:
        logger.debug('Setting parse_taxonomy to True because fill_unresolved_taxonomy is True.')
        parse_taxonomy = True

    # Load the feature table
    logger.debug('Loading feature table')
    feature_table = load_feature_table(feature_table_filepath, header_row=header_row,
                                       original_feature_id_colname=original_feature_id_colname,
                                       final_feature_id_colname='Feature ID')

    # Normalize
    if normalization_method is not None:
        feature_table = normalize_feature_table(feature_table, normalization_method)

    # Add taxonomy
    if taxonomy_filepath is not None:
        feature_table = add_taxonomy_to_feature_table(feature_table, taxonomy_filepath)

    # Parse taxonomy
    if parse_taxonomy is True:
        logger.debug('Parsing taxonomy into 7 ranks')
        taxonomy_entries_parsed = map(lambda entry: parse_silva_taxonomy_entry(entry, resolve=fill_unresolved_taxonomy),
                                      feature_table['Taxonomy'].tolist())
        taxonomy_table_parsed = pd.DataFrame(taxonomy_entries_parsed, columns=TAXONOMY_RANKS)
        # Bind to main table in place of 'Taxonomy'
        feature_table = pd.concat([feature_table, taxonomy_table_parsed], axis=1, sort=False)
        feature_table = feature_table.drop(columns='Taxonomy')

    # Add representative sequences
    if sequence_filepath is not None:
        feature_table = add_sequences_to_feature_table(feature_table, sequence_filepath)

    # Sort Feature IDs
    if sort_features is True:
        feature_table = sort_feature_table(feature_table)

    # Rename Feature IDs
    if rename_features is True:
        logger.debug('Renaming feature IDs sequentially')
        num_rows = feature_table.shape[0]
        feature_table['Feature ID'] = range(num_rows)

    # Sort the feature table columns sensibly
    feature_table = sort_feature_table_columns(feature_table, sort_sample_columns=sort_columns)

    # Change first column to that desired by user
    if final_feature_id_colname != 'Feature ID':
        logger.debug('Changing "Feature ID" colname to "' + final_feature_id_colname + '"')
        feature_table = feature_table.rename(columns={'Feature ID': final_feature_id_colname})

    return feature_table


def subparse_cli(subparsers, parent_parser: argparse.ArgumentParser = None):
    """
    Parses the CLI arguments and adds them as a subparser to an existing parser.

    :param subparsers: A special subparser action object created from an existing parser by the add_subparsers() method.
                       For example, parser = argparse.ArgumentParser(); subparsers = parser.add_subparsers().
    :param parent_parser: An optional ArgParse object with additional arguments (e.g., shared across all modules) to
                          add to this CLI parser. This can be a unique parser and does not need to be the parent of the
                          subparsers object. If None, then no parent will be added for the subparser.
    :return: An ArgumentParser object created by subparsers.add_parser()
    """

    description = 'Creates a TSV-format QIIME2 feature table with overlaid taxonomy and sequence information.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('tabulate', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(tabulate=True)

    file_settings = subparser.add_argument_group('Input/output file options')
    table_settings = subparser.add_argument_group('Optional table manipulation options')

    file_settings.add_argument('-f', '--feature_table', metavar='TABLE', required=True,
                               help='The path to the input feature table file. TSV, BIOM, or QZA file types are '
                                    'supported.')
    file_settings.add_argument('-o', '--output_feature_table', metavar='TABLE', required=False, default='-',
                               help='The path to the output TSV feature table. Will write to STDOUT (-) if nothing '
                                    'is provided.')
    file_settings.add_argument('-s', '--sequences', metavar='SEQ', required=False, default=None,
                               help='The path to the input ASV/OTU sequence file. Sequences will be added as the '
                                    '"Sequences" column. You can optionally omit this flag and not have sequences '
                                    'added to the table. FastA and QZA file types are supported.')
    file_settings.add_argument('-t', '--taxonomy', metavar='TAX', required=False, default=None,
                               help='The path to the input taxonomy file. Taxonomy will be added as the "Taxonomy" '
                                    'column. You can optionally omit this flag and not have taxonomy added to the '
                                    'table. TSV (qiime2 format) or QZA file types are supported.')

    table_settings.add_argument('-n', '--normalize', required=False, choices=['percent', 'proportion'],
                                help='Optionally normalize each sample by the desired normalization method (e.g., to '
                                     'convert to percent relative abundance).')
    table_settings.add_argument('-S', '--sort_features', required=False, action='store_true',
                                help='Optionally sort Feature IDs roughly based on overall abundance.')
    table_settings.add_argument('-R', '--rename_features', required=False, action='store_true',
                                help='Optionally rename the Feature IDs sequentially, roughly based on overall '
                                     'abundance. Automatically sets --sort_features.')
    table_settings.add_argument('-P', '--parse_taxonomy', required=False, action='store_true',
                                help='Optionally parse Silva taxonomy into 7 ranks with columns "Domain", "Phylum", '
                                     'etc.')
    table_settings.add_argument('-u', '--fill_unresolved_taxonomy', required=False, action='store_true',
                                help='Optionally add Unresolved_[taxon] labels for blank taxonomy ranks. Requires '
                                     '--parse_taxonomy.')
    table_settings.add_argument('-C', '--sort_columns', required=False, action='store_true',
                                help='Sort feature table columns by sample ID (alphabetically)')
    table_settings.add_argument('-H', '--header_row', metavar='NUM', required=False, default='auto',
                                help='The row to use as header for the input feature table (e.g., 0 = top row). '
                                     '[Default: "auto" = auto detect best row]')
    table_settings.add_argument('-i', '--original_feature_id_colname', metavar='NAME', required=False, default='auto',
                                help='The name of the Feature ID column of the input feature table. [Default: '
                                     '"auto" = auto detect]')
    table_settings.add_argument('-I', '--final_feature_id_colname', metavar='NAME', required=False,
                                default='Feature ID',
                                help='The name of the first column of the output feature table. [Default: '
                                     '"Feature ID"]')

    return subparser
