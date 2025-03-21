#!/usr/bin/env python
"""
ampliwrangler
Description: Simple command-line utilities for enhancing QIIME2-based amplicon analyses.
Copyright: Jackson M. Tsuji, 2025
"""

import os
import sys
import logging
import argparse

from ampliwrangler.params import VERSION
from ampliwrangler.utils import check_output_file, set_up_output_directory
from ampliwrangler.tabulate import main as tabulate_main
from ampliwrangler.tabulate import subparse_cli as subparse_tabulate_cli
from ampliwrangler.count import main as count_main
from ampliwrangler.count import subparse_cli as subparse_count_cli
from ampliwrangler.split import main as split_main
from ampliwrangler.split import subparse_cli as subparse_split_cli


def main():
    """
    Collects input arguments and selects a command to perform.
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of ampliwrangler.py. Run parser inside main() so it can be called 
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    # Initialize the root logger with a stream handler
    logger = logging.getLogger()
    formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(filename)s: %(funcName)s: %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if (hasattr(args, 'version')) and (args.version is True):
        print(f'Version: {VERSION}')
        sys.exit(0)

    if (hasattr(args, 'verbose')) and (args.verbose is True):
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # File logger setup
    if (hasattr(args, 'logfile')) and (args.logfile is not None):
        check_output_file(output_filepath=args.logfile, overwrite=args.overwrite)

        # Note: The logging level set here is the minimum level that the logger can write to. The actual level is
        #       defined for the root logger above in the args.verbose conditional.
        custom_path_file_handler = logging.FileHandler(filename=args.logfile, mode='w')
        custom_path_file_handler.setFormatter(formatter)
        custom_path_file_handler.setLevel(logging.DEBUG)
        logger.addHandler(custom_path_file_handler)

    if hasattr(args, 'output_dir'):
        set_up_output_directory(output_directory_filepath=args.output_dir, overwrite=args.overwrite)

        # Start log file in the output dir
        # Note: The logging level set here is the minimum level that the logger can write to. The actual level is
        #       defined for the root logger above in the args.verbose conditional.
        output_dir_file_handler = logging.FileHandler(filename=os.path.join(args.output_dir, 'log.txt'), mode='w')
        output_dir_file_handler.setFormatter(formatter)
        output_dir_file_handler.setLevel(logging.DEBUG)
        logger.addHandler(output_dir_file_handler)

    # Select the sub-command to run.
    if hasattr(args, 'tabulate'):
        tabulate_main(args)
    elif hasattr(args, 'count'):
        count_main(args)
    elif hasattr(args, 'split'):
        split_main(args)
    else:
        parser.print_help()


def parse_cli():
    """
    Parses the CLI arguments.
    :return: An argparse parser object.
    """
    cli_title = (f'Ampliwrangler: simple command-line utilities for enhancing QIIME2-based amplicon analyses.\n'
                 f'Copyright Jackson M. Tsuji, 2025.\n'
                 f'Version: {VERSION}')
    parser = argparse.ArgumentParser(description=cli_title)
    parser.add_argument('-V', '--version', required=False, action='store_true', help='Print tool version and exit')

    # Common arguments across all modules
    parent_parser = argparse.ArgumentParser(add_help=False)
    basic_config = parent_parser.add_argument_group('Basic config settings')
    basic_config.add_argument('-lf', '--logfile', metavar='PATH', required=False, default=None, type=str,
                              help='Log filepath (default: None)')
    basic_config.add_argument('-O', '--overwrite', required=False, action='store_true',
                              help='Overwrite existing files/directories. By setting this flag, you risk erasing old '
                                   'data.')
    basic_config.add_argument('-v', '--verbose', required=False, action='store_true',
                              help='Enable verbose logging')

    # Declare sub-commands
    subparsers = parser.add_subparsers(help='Available modules:')
    subparse_tabulate_cli(subparsers, parent_parser)
    subparse_count_cli(subparsers, parent_parser)
    subparse_split_cli(subparsers, parent_parser)

    return parser


if __name__ == '__main__':
    main()
