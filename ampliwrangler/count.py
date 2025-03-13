#!/usr/bin/env python
"""
ampliwrangler count
Description: Creates a tabular summary of the total read counts of all samples in a qiime2 feature table
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import sys
import os
import shutil
import logging
import argparse
import uuid
import biom

import zipfile as zf
import pandas as pd

logger = logging.getLogger(__name__)


def main(args):
    """
    Runs the workflow based on the provided command line arguments.
    """
    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
        logging.getLogger('ampliwrangler.count').setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logging.getLogger('ampliwrangler.count').setLevel(logging.INFO)

    # Startup messages
    logger.info(f'Running {os.path.basename(sys.argv[0])}')
    logger.info(f'Input filepath: {args.input_filepath}')
    logger.info(f'Output filepath: {args.output_filepath}')
    logger.info(f'Min count filepath: {args.min_count_filepath}')
    logger.info(f'Temp directory: {args.tmp_dir}')
    logger.info(f'Verbose logging: {args.verbose}')

    process_sample_counts(args.input_filepath, args.output_filepath, args.min_count_filepath, args.tmp_dir)


def unpack_biom_from_qza(input_filepath, tmp_dir, qza_biom_path='data/feature-table.biom') -> tuple:
    """
    Unzip a BIOM file from a QIIME2 QZA archive.

    :global BIOM_PATH_PARTIAL: partial path to the biom file within the QZA archive
    :param input_filepath: Path to the QIIME2 QZA archive, FeatureTable[Frequency] format
    :param tmp_dir: The base directory to extract the ZIP file to (will create a random subfolder)
    :param qza_biom_path: expected path for the BIOM file inside the unpacked QZA file
    :return: tuple of the path to the unzipped BIOM file and the path to the random temp subfolder
    """
    logger.info('Unpacking QZA file')

    with zf.ZipFile(input_filepath, 'r') as qza_data:
        # Dump list of contents and determine path to the BIOM file
        qza_data_contents = qza_data.namelist()

        # Find the path to the biom file
        # See https://stackoverflow.com/a/12845341 (accessed Sept. 12, 2019)
        biom_filepath = list(filter(lambda x: qza_biom_path in x, qza_data_contents))
        # TODO - check this is a list of length 1
        logger.debug(f'Found biom file in QZA file at {biom_filepath[0]}')

        # Extract the biom file to a randomly generated subfolder in tmp_dir
        tmp_subdir = os.path.join(tmp_dir, uuid.uuid4().hex)
        extraction_path = qza_data.extract(biom_filepath[0], path=tmp_subdir)
        logger.debug(f'Extracted temp BIOM file to {extraction_path}')

    biom_info = (extraction_path, tmp_subdir)

    return biom_info


def load_biom_df_from_qza(input_filepath: str, tmp_dir: str) -> pd.DataFrame:
    """
    Load a BIOM file from a QIIME2 QZA archive as a Pandas dataframe.

    :param input_filepath: Path to the QIIME2 QZA archive, FeatureTable[Frequency] format
    :param tmp_dir: The base directory to extract the ZIP file to (will create a random subfolder)
    :return: pandas DataFrame of the BIOM table's count data
    """
    # Unzip the BIOM table from within the QZA (ZIP) file
    biom_path, tmp_subdir = unpack_biom_from_qza(input_filepath, tmp_dir)

    # Load biom file and convert to pandas dataframe
    logger.debug("Loading BIOM file")
    biom_data = biom.load_table(biom_path)
    biom_table = biom_data.to_dataframe()

    # Delete tmp dir
    logger.debug(f'Removing temp dir {tmp_subdir}')
    shutil.rmtree(tmp_subdir)

    return biom_table


def process_sample_counts(input_filepath: str, output_filepath: str, min_count_filepath: str = None,
                          tmp_dir: str = '.'):
    """
    Load feature table, sum columns, and write output

    :param input_filepath: path to the input QZA file
    :param output_filepath: path to the output TSV file with sample count data (or - = stdout)
    :param min_count_filepath: path to the output TXT file where just the min count info is stored
    :param tmp_dir: temporary directory for loading QZA data (a random subfolder will be temporarily created)
    :return: writes output files to output_filepath and optionally min_count_filepath
    """
    # Load the BIOM table from within the QZA (ZIP) file
    biom_table = load_biom_df_from_qza(input_filepath, tmp_dir)

    # Sum columns
    logger.info('Summarizing BIOM file')
    sample_counts = biom_table.sum().to_frame().reset_index()
    sample_counts.columns = ['sample-id', 'count']
    # Sort by count
    sample_counts = sample_counts.sort_values(by=['count'], ascending=False)

    # Write output
    if output_filepath == '-':
        # Write to STDOUT
        logger.info('Writing counts summary file to STDOUT')
        sample_counts.to_csv(sys.stdout, sep='\t', index=False)
    else:
        logger.info(f'Writing counts summary file to {output_filepath}')
        sample_counts.to_csv(output_filepath, sep='\t', index=False)

    # Get min count and write, if desired
    if min_count_filepath:
        logger.info(f'Writing min count to {min_count_filepath}')
        min_count = int(min(sample_counts['count']))
        with open(min_count_filepath, 'w') as min_count_file:
            min_count_file.write(f'{str(min_count)}\n')

    logger.info(f'{os.path.basename(sys.argv[0])}: done.')


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

    description = 'Creates a tabular summary of the total read counts of all samples in a qiime2 feature table.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('count', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(count=True)

    file_settings = subparser.add_argument_group('Input/output file options')
    workflow_settings = subparser.add_argument_group('Workflow options')

    file_settings.add_argument('-i', '--input_filepath', metavar='QZA', required=True,
                               help='The path to the input QZA FeatureTable file.')
    file_settings.add_argument('-o', '--output_filepath', metavar='TSV', required=False, default='-',
                               help='The path to the output TSV file. Will write to STDOUT (-) if nothing is provided.')
    file_settings.add_argument('-m', '--min_count_filepath', metavar='TXT', required=False, default=None,
                               help='Optional path to write a single-line text file with the lowest count value in the '
                                    'dataset.')

    workflow_settings.add_argument('-T', '--tmp_dir', metavar='DIR', required=False, default='.',
                                   help='Optional path to the temporary directory used for unpacking the QZA. A random '
                                        'subdirectory will be temporarily created here [default: .].')

    return subparser
