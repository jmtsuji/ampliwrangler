#!/usr/bin/env python
"""
ampliwrangler split
Description: Simple script to split a QIIME2 manifest file into multiple files based on the run_ID column
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import sys
import os
import logging
import argparse

import pandas as pd
from ampliwrangler.utils import check_output_file, set_up_output_directory

logger = logging.getLogger(__name__)


def main(args):
    """
    Runs the workflow based on the provided command line arguments.
    """
    # Startup checks
    if args.verbose is True:
        logger.setLevel(logging.DEBUG)
        logging.getLogger('ampliwrangler.split').setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logging.getLogger('ampliwrangler.split').setLevel(logging.INFO)

    # Startup messages
    logger.info(f'Running {os.path.basename(sys.argv[0])}')
    logger.debug(f'Input filepath: {args.input_filepath}')
    logger.debug(f'Output directory: {args.output_dir}')
    logger.debug(f'Run ID column name: {args.run_id_column}')
    logger.debug(f'Overwrite existing files: {args.overwrite}')

    set_up_output_directory(args.output_dir, args.overwrite)
    split_manifest_by_run_id(args.input_filepath, args.output_dir, args.run_id_column, args.overwrite)


def split_manifest_by_run_id(input_filepath: str, output_dir: str, run_id_column: str = 'run_ID',
                             overwrite: bool = False):
    """
    Splits an input manifest file by run based on values in the run ID column.
    Writes outputs as individual manifest files.

    :param input_filepath: path to the input manifest file
    :param output_dir: path to the directory where manifest files, split by run ID, will be written
    :param run_id_column: name of the column in the input manifest file with run IDs
    :param overwrite: whether or not to overwrite any existing files
    :return: writes output files to output_dir
    """
    # Load manifest file
    # TODO - check for second row specifying types and add as a second header row if it exists
    logger.info('Loading manifest file')
    manifest_table = pd.read_csv(input_filepath, sep='\t', header=0)

    # Does the run ID column exist?
    if run_id_column not in manifest_table.columns.values.tolist():
        error = FileNotFoundError(f'Did not find the run ID column "{run_id_column}" in the provided manifest file '
                                  f'at {input_filepath}.')
        logger.error(error)
        raise error

    # Get unique run IDs
    unique_run_ids = set(manifest_table[run_id_column])

    # Split table and write
    for run_id in unique_run_ids:
        output_filename = f'manifest_{run_id}.tsv'
        output_filepath = os.path.join(output_dir, output_filename)

        check_output_file(output_filepath, overwrite=overwrite)
        logger.info(f'Writing run "{run_id}" to file "{output_filename}".')

        single_run_table = manifest_table[manifest_table[run_id_column] == run_id]
        pd.DataFrame.to_csv(single_run_table, output_filepath, sep='\t', index=False)

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

    description = 'Simple script to split a QIIME2 manifest file into multiple files based on the run ID column.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('split', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(split=True)

    file_settings = subparser.add_argument_group('Input/output file options')
    run_settings = subparser.add_argument_group('Other params')

    file_settings.add_argument('-i', '--input_filepath', metavar='MANIFEST', required=True,
                               help='The path to the input manifest file. Must match QIIME standards AND have a column '
                                    'with a unique ID for each sequencing run as specified in --run_id_column')
    file_settings.add_argument('-o', '--output_dir', metavar='DIR', required=True,
                               help='The directory where output files (named "manifest_[RUN_ID].tsv") will be saved.')

    run_settings.add_argument('-r', '--run_id_column', metavar='STR', required=False, default='run_ID',
                               help='Name of the column in the input manifest file that contains the unique IDs for '
                                    'each run. [default: run_ID]')

    return subparser
