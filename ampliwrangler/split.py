#!/usr/bin/env python
"""
ampliwrangler split
Description: Simple script to split a QIIME2 manifest file into multiple files based on the run_ID column
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import sys
import os
import time
import logging
import argparse

import pandas as pd

logger = logging.getLogger(__name__)


def main(args):
    # Set user variables
    input_filepath = args.input_filepath
    output_dir = args.output_dir

    # Startup messages
    logger.info("Running " + os.path.basename(sys.argv[0]))
    logger.info("Version: " + SCRIPT_VERSION)
    logger.info("Input filepath: " + input_filepath)
    logger.info("Output directory: " + output_dir)

    # Load manifest file
    # TODO - check for second row specifying types and add as a second header row if it exists
    logger.info("Loading manifest file")
    manifest_table = pd.read_csv(input_filepath, sep='\t', header=0)

    # Does run_ID exist?
    if 'run_ID' not in manifest_table.columns.values.tolist():
        logger.error("Did not find the 'run_ID' column in the provided manifest file. Cannot divide the manifest file "
                     "by run. Exiting...")
        sys.exit(1)

    # Get unique run IDs
    unique_run_ids = set(manifest_table['run_ID'])

    # Split table and write
    for run_ID in unique_run_ids:
        output_filename = "manifest_" + run_ID + ".tsv"
        output_filepath = os.path.join(output_dir, output_filename)
        # TODO - check if output_filepath already exists and do not write output unless --force is specified
        logger.info("Writing run '" + run_ID + "' to file '" + output_filename + "'")

        single_run_table = manifest_table[manifest_table['run_ID'] == run_ID]
        pd.DataFrame.to_csv(single_run_table, output_filepath, sep='\t', index=False)

    logger.info(os.path.basename(sys.argv[0]) + ": done.")


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

    description = 'Simple script to split a QIIME2 manifest file into multiple files based on the run_ID column.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('count', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(count=True)

    file_settings = subparser.add_argument_group('Input/output file options')
    workflow_settings = subparser.add_argument_group('Workflow options')

    file_settings.add_argument('-i', '--input_filepath', metavar='input', required=True,
                               help='The path to the input manifest file. Must match QIIME standards AND have a column '
                                    'named "run_ID" with a unique ID for each Illumina run')
    file_settings.add_argument('-o', '--output_dir', metavar='output', required=True,
                               help='The directory where output files (named "manifest_[RUN_ID].tsv") will be saved. '
                                    'Will OVERWRITE existing files.')

    workflow_settings.add_argument('-v', '--verbose', required=False, action='store_true',
                                   help='Enable for verbose logging.')

    return subparser
