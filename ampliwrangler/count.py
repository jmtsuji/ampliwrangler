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

from ampliwrangler.utils import check_output_file, load_feature_table

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

    check_output_file(args.output_filepath, overwrite=args.overwrite)
    if args.min_count_filepath:
        check_output_file(args.min_count_filepath, overwrite=args.overwrite)

    # Startup messages
    logger.info(f'Running {os.path.basename(sys.argv[0])} count')
    logger.debug(f'Input filepath: {args.input_filepath}')
    logger.debug(f'Output filepath: {args.output_filepath}')
    logger.debug(f'Min count filepath: {args.min_count_filepath}')
    logger.debug(f'Temp directory: {args.tmp_dir}')
    logger.debug(f'Overwrite existing files: {args.overwrite}')
    logger.debug(f'Verbose logging: {args.verbose}')

    process_sample_counts(args.input_filepath, args.output_filepath, args.min_count_filepath, args.tmp_dir)

    logger.info(f'{os.path.basename(sys.argv[0])} count: done.')


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
    feature_data = load_feature_table(input_filepath, tmp_dir=tmp_dir)

    # Sum columns
    logger.info('Summarizing FeatureTable counts')
    sample_counts = feature_data.sum().to_frame().reset_index()
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
    workflow_settings = subparser.add_argument_group('Workflow options')

    file_settings.add_argument('-i', '--input_filepath', metavar='TABLE', required=True,
                               help='The path to the input FeatureTable file, either QZA, BIOM, or TSV.')
    file_settings.add_argument('-o', '--output_filepath', metavar='TABLE', required=False, default='-',
                               help='The path to the output TSV file. Will write to STDOUT (-) if nothing is provided.')
    file_settings.add_argument('-m', '--min_count_filepath', metavar='TXT', required=False, default=None,
                               help='Optional path to write a single-line text file with the lowest count value in the '
                                    'dataset.')

    workflow_settings.add_argument('-T', '--tmp_dir', metavar='DIR', required=False, default='.',
                                   help='Optional path to the temporary directory used for unpacking the QZA. A random '
                                        'subdirectory will be temporarily created here [default: .].')

    return subparser
