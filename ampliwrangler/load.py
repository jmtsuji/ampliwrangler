#!/usr/bin/env python
"""
load.py
Description: functions related to data loading within ampliwrangler
Copyright: Jackson M. Tsuji, 2025
"""

# Imports
import logging
import os
import shutil
import uuid
import biom

import pandas as pd
import zipfile as zf
from Bio import SeqIO

from ampliwrangler.utils import set_up_output_directory

logger = logging.getLogger(__name__)


def load_fasta_sequences(fasta_filepath: str):
    """
    Loads an input FastA file as a generator.

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :return: generator of a SeqRecord object for the loaded sequences
    """

    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            yield record


def unpack_from_qza(qza_filepath: str, target_filename: str, qza_data_dir: str = 'data') -> bytes:
    """
    Unpack a data file from a QIIME2 QZA archive.

    :param qza_filepath: Path to the QIIME2 QZA archive
    :param target_filename: name of the target file inside the unpacked QZA
    :param qza_data_dir: expected path inside the QZA where the target file is stored
    :return: contents of the unzipped file as bytes
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

        with qza_data.open(target_filepath) as target_handle:
            target_contents = target_handle.read()

    return target_contents


def load_tsv_feature_table(tsv_filepath: str, header_row='auto') -> pd.DataFrame:
    """
    Loads a TSV-format FeatureTable as a pandas DataFrame. Roughly tries to handle screening out a header line.
    :params tsv_filepath: path to the TSV-format FeatureTable
    :params header_row: manually specify the header row as per the header argument in pd.read_csv (e.g., by integer,
                        where 0 is the top row), or specify 'auto' to try to auto-detect the correct header row.
    :return: FeatureTable data as a pandas DataFrame
    """
    if header_row == 'auto':
        # Roughly check if the FeatureTable has a first header line or not
        with open(tsv_filepath, 'r') as input_handle:
            first_line = input_handle.readline()

        if first_line == '# Constructed from biom file\n':
            feature_data = pd.read_csv(tsv_filepath, sep='\t', header=1)
        else:
            feature_data = pd.read_csv(tsv_filepath, sep='\t', header=0)

    else:
        feature_data = pd.read_csv(tsv_filepath, sep='\t', header=header_row)

    return feature_data


def load_biom_feature_table(biom_filepath: str) -> pd.DataFrame:
    """
    Load a biom-format feature table as a pandas DataFrame
    :params biom_filepath: path to the biom-format file to load
    :return: pandas DataFrame of the biom file data
    """
    biom_table = biom.load_table(biom_filepath).to_dataframe()

    return biom_table


def load_qza_feature_table(qza_filepath: str, tmp_dir: str = '.') -> pd.DataFrame:
    """
    Load a BIOM file from a QIIME2 QZA archive as a Pandas dataframe.

    :param qza_filepath: Path to the QIIME2 QZA archive, FeatureTable[Frequency] format
    :param tmp_dir: The base directory to extract the ZIP file to (will create a random subfolder)
    :return: pandas DataFrame of the BIOM table's count data
    """
    # Load the biom file contents in binary, then write to a temp file that biom.load_table can open
    biom_contents = unpack_from_qza(qza_filepath, target_filename='feature-table.biom')
    tmp_subdir = os.path.join(tmp_dir, f'.{uuid.uuid4().hex}')
    set_up_output_directory(tmp_subdir, overwrite=False)
    biom_path = os.path.join(tmp_subdir, 'feature-table.biom')
    logger.debug(f'Writing temp biom file to {biom_path}')
    with open(biom_path, 'wb') as biom_handle:
        biom_handle.write(biom_contents)

    # Load biom file
    biom_table = load_biom_feature_table(biom_path)

    # Delete tmp dir
    logger.debug(f'Removing temp dir {tmp_subdir}')
    shutil.rmtree(tmp_subdir)

    return biom_table


def load_feature_table(feature_table_filepath, header_row='auto', tmp_dir: str = '.') -> pd.DataFrame:
    """
    Generalized function to load a feature table regardless of input format
    :params feature_table_filepath: path to the feature table to load (in TSV, BIOM, or QZA format)
    :params header_row: for TSV format files, manually specify the header row here, as per the header argument in
                        pd.read_csv (e.g., by integer, where 0 is the top row), or specify 'auto' to try to auto-detect
                        the correct header row.
    :params tmp_dir: for QZA format files, optionally provide the base directory to extract the ZIP file to (will create
                     a random subfolder here)
    :return: pandas DataFrame of the loaded feature table
    """
    try:
        # First try: assume TSV format
        logger.debug(f'Trying to load FeatureTable as TSV: {feature_table_filepath}')
        feature_data = load_tsv_feature_table(feature_table_filepath, header_row=header_row)
        logger.debug(f'TSV loading succeeded')
    except UnicodeDecodeError:
        try:
            # Second try: BIOM format
            logger.debug(f'TSV loading failed')
            logger.debug(f'Trying to load FeatureTable as BIOM: {feature_table_filepath}')
            feature_data = load_biom_feature_table(feature_table_filepath)
            logger.debug(f'BIOM loading succeeded')
        except UnicodeDecodeError:
            # Third try: QZA format
            logger.debug(f'BIOM loading failed')
            logger.debug(f'Trying to load FeatureTable as QZA: {feature_table_filepath}')
            feature_data = load_qza_feature_table(feature_table_filepath, tmp_dir=tmp_dir)
            logger.debug(f'QZA loading succeeded')

    return feature_data
