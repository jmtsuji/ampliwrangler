#!/usr/bin/env python
"""
utils.py
Description: utility functions within ampliwrangler
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import logging
import os
import sys

logger = logging.getLogger(__name__)


def check_output_file(output_filepath: str, overwrite: bool = False):
    """
    Checks if OK to create an output file. Raises an error if the output file already exists (unless overwrite=True).

    :param output_filepath: path to the desired output file
    :param overwrite: if True, the keep going with a warning if the output file already exists
    """

    output_file_exists = os.path.isfile(output_filepath)

    if output_file_exists is True:
        if overwrite is False:
            logger.error(f'Output file already exists: "{output_filepath}". Will not continue. Set the '
                         f'--overwrite flag at your own risk if you want to overwrite existing files.')
            sys.exit(1)
        elif overwrite is True:
            logger.warning(f'Output file already exists: "{output_filepath}". File will be overwritten.')
        else:
            error = ValueError(f'overwrite must be True or False, but you provided "{overwrite}"')
            logger.error(error)
            raise error


def set_up_output_directory(output_directory_filepath: str, overwrite: bool = False):
    """
    Creates an output directory. Raises an error if a directory already exists (unless overwrite=True).

    :param output_directory_filepath: path to the desired output directory
    :param overwrite: if True, then keep going with a warning if the output directory already exists
    """

    output_dir_exists = os.path.isdir(output_directory_filepath)

    if output_dir_exists is True:
        if overwrite is False:
            error = RuntimeError(f'Output directory already exists: "{output_directory_filepath}". Will not continue. '
                                 f'Set the --overwrite flag at your own risk if you want to use an existing directory.')
            logger.error(error)
            raise error
        elif overwrite is True:
            logger.warning(f'Output directory already exists: "{output_directory_filepath}". Files may be overwritten.')
        else:
            raise ValueError(f'overwrite must be True or False, but you provided "{overwrite}"')

    os.makedirs(output_directory_filepath, exist_ok=True)
