#!/usr/bin/env python
# utils.py
# Utility functions within spokewrench
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024
import logging
import os
import sys
import uuid

import zipfile as zf
from Bio import SeqIO

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
            logger.error(f'Output directory already exists: "{output_directory_filepath}". Will not continue. Set the '
                         f'--overwrite flag at your own risk if you want to use an existing directory.')
            sys.exit(1)
        elif overwrite is True:
            logger.warning(f'Output directory already exists: "{output_directory_filepath}". Files may be overwritten.')
        else:
            raise ValueError(f'overwrite must be True or False, but you provided "{overwrite}"')

    os.makedirs(output_directory_filepath, exist_ok=True)


def load_fasta_sequences(fasta_filepath: str):
    """
    Loads an input FastA file as a generator.

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :return: generator of a SeqRecord object for the loaded sequences
    """

    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            yield record


def unpack_from_qza(qza_filepath: str, target_filename: str, qza_data_dir: str = 'data', tmp_dir: str = '.') -> tuple:
    """
    Unpack a data file from a QIIME2 QZA archive.

    :param qza_filepath: Path to the QIIME2 QZA archive
    :param target_filename: name of the target file inside the unpacked QZA
    :param qza_data_dir: expected path inside the QZA where the target file is stored
    :param tmp_dir: The base directory to extract the QZA file to (will create a random subfolder inside)
    :return: tuple: path to the unzipped file and the path to the random temp subfolder
    """
    logger.debug(f'Unpacking QZA file {qza_filepath}')

    with zf.ZipFile(qza_filepath, 'r') as qza_data:
        # Dump list of contents and determine path to the BIOM file
        qza_data_contents = qza_data.namelist()

        # Find the path to the target file
        # TODO - consider getting the exact path more carefully by grabbing the UUID from metadata.yaml
        target_path_partial = os.path.join(qza_data_dir, target_filename)
        matching_filepaths = []
        for filepath in qza_data_contents:
            if target_path_partial in filepath:
                matching_filepaths.append(filepath)

        if len(matching_filepaths) == 1:
            target_filepath = matching_filepaths[0]
            logger.debug(f'Found target file in QZA at {target_filepath}')
        else:
            error = ValueError(f'Found more than one file matching the target: {matching_filepaths}')
            logger.error(error)
            raise error

        # Extract the target file to a randomly generated subfolder in tmp_dir
        tmp_subdir = os.path.join(tmp_dir, f'.{uuid.uuid4().hex}')
        set_up_output_directory(tmp_subdir, overwrite=False)
        extraction_path = qza_data.extract(target_filepath, path=tmp_subdir)
        logger.debug(f'Extracted target file to {extraction_path}')

    extracted_info = (extraction_path, tmp_subdir)

    return extracted_info
