#!/usr/bin/env python
"""
ampliwrangler count
Description: Creates a tabular summary of the total read counts of all samples in a qiime2 feature table
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import sys
import os
import logging
import argparse

from ampliwrangler.params import VERSION
from ampliwrangler.utils import check_output_file
from ampliwrangler.load import load_feature_table

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

    if args.output_counts_filepath:
        check_output_file(args.output_counts_filepath, overwrite=args.overwrite)
    if args.output_min_count_filepath:
        check_output_file(args.output_min_count_filepath, overwrite=args.overwrite)

    if (args.output_counts_filepath == "-") & (args.output_min_count_filepath == "-"):
        error = RuntimeError('Cannot supply "-" (stdout) to both --output_counts_filepath and '
                             '--output_min_count_filepath')
        logger.error(error)
        raise error

    # Startup messages
    logger.info(f'Running {os.path.basename(sys.argv[0])} count, version {VERSION}')
    logger.debug(f'Input filepath: {args.input_filepath}')
    logger.debug(f'Counts-per-sample filepath: {args.output_counts_filepath}')
    logger.debug(f'Min count filepath: {args.output_min_count_filepath}')
    logger.debug(f'Overwrite existing files: {args.overwrite}')
    logger.debug(f'Verbose logging: {args.verbose}')

    process_sample_counts(args.input_filepath, args.output_counts_filepath, args.output_min_count_filepath)

    logger.info(f'{os.path.basename(sys.argv[0])} count: done.')


def process_sample_counts(feature_table_filepath: str, counts_filepath: str = None, min_count_filepath: str = None):
    """
    Load feature table, sum columns, and write output

    :param feature_table_filepath: path to the input FeatureTable (QZA, BIOM, or TSV format)
    :param counts_filepath: path to an output TSV file with sample count data (or - = stdout)
    :param min_count_filepath: path to an output TXT file where just the min count info is stored (or - = stdout)
    :return: writes output files to output_filepath and optionally min_count_filepath
    """
    feature_data = load_feature_table(feature_table_filepath)\
        .set_index('Feature ID')

    # Sum columns
    logger.info('Summarizing FeatureTable counts')
    sample_counts = feature_data.sum().to_frame().reset_index()
    sample_counts.columns = ['sample-id', 'count']
    # Sort by count
    sample_counts = sample_counts.sort_values(by=['count'], ascending=False)

    if counts_filepath:
        # Write counts-per-sample
        if counts_filepath == '-':
            # Write to STDOUT
            logger.info('Writing counts-per-sample to STDOUT')
            sample_counts.to_csv(sys.stdout, sep='\t', index=False)
        else:
            logger.info(f'Writing counts-per-sample as file: {counts_filepath}')
            sample_counts.to_csv(counts_filepath, sep='\t', index=False)

    if min_count_filepath:
        # Write minimum read count per sample
        min_count = int(min(sample_counts['count']))

        if min_count_filepath == '-':
            # Write to STDOUT
            logger.info('Writing minimum count per sample to STDOUT')
            sys.stdout.write(f'{min_count}\n')
        else:
            logger.info(f'Writing minimum count per sample as file: {min_count_filepath}')
            with open(min_count_filepath, 'w') as output_handle:
                output_handle.write(f'{min_count}\n')


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

    description = 'Creates a tabular summary of the total read counts of all samples in a feature table.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('count', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(count=True)

    file_settings = subparser.add_argument_group('Input/output file options')

    file_settings.add_argument('-i', '--input_filepath', metavar='TABLE', required=True,
                               help='The path to the input FeatureTable file, either QZA, BIOM, or TSV.')
    file_settings.add_argument('-c', '--output_counts_filepath', metavar='TABLE', required=False, default=None,
                               help='Optional path to write output counts per sample in TSV format. Provide "-" to '
                                    'write to stdout.')
    file_settings.add_argument('-m', '--output_min_count_filepath', metavar='TXT', required=False, default=None,
                               help='Optional path to write a single-line text file with the lowest count value in the '
                                    'dataset. Provide "-" to write to stdout.')

    return subparser
