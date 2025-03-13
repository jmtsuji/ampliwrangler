#!/usr/bin/env python
"""
load.py
Description: functions related to data loading within ampliwrangler
Copyright: Jackson M. Tsuji, 2025
"""
# Imports
import logging
import os
import io
import shutil
import uuid
import biom

import pandas as pd
import zipfile as zf
from Bio import SeqIO

from ampliwrangler.utils import set_up_output_directory

logger = logging.getLogger(__name__)


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


def validate_feature_table(feature_data: pd.DataFrame, original_feature_id_colname: str = 'auto',
                           final_feature_id_colname: str = 'Feature ID') -> pd.DataFrame:
    """
    Validate and standardize a QIIME2 feature table based on the feature ID column.

    :param feature_data: FeatureTable loaded as a pandas DataFrame
    :param original_feature_id_colname: name of the feature ID column in the provided table. Must be the first column in
                                        the table. Auto-detection is attempted if 'auto' is provided.
    :param final_feature_id_colname: updated name of the feature ID column, for standardization. Column name is kept
                                     as-is if None is provided.
    :return: FeatureTable data as a pandas DataFrame
    """
    acceptable_feature_id_column_names = ['Feature ID', '#OTU ID']

    first_column_name = feature_data.columns[0]
    if original_feature_id_colname == 'auto':
        if first_column_name not in acceptable_feature_id_column_names:
            error = ValueError(f'First column of FeatureTable is "{first_column_name}", not one of '
                               f'{acceptable_feature_id_column_names}')
            logger.error(error)
            raise error
    else:
        if first_column_name != original_feature_id_colname:
            error = ValueError(f'Expected first column of FeatureTable to be "{original_feature_id_colname}", but got'
                               f'"{first_column_name}".')
            logger.error(error)
            raise error

    if (final_feature_id_colname is not None) & (first_column_name != final_feature_id_colname):
        if final_feature_id_colname not in feature_data.columns:
            logger.debug(f'Standardizing FeatureTable ID colname to {final_feature_id_colname} for analysis')
            feature_data = feature_data.rename(columns={first_column_name: final_feature_id_colname})
        else:
            error = ValueError(f'Final feature ID column "{final_feature_id_colname}" already exists in the provided '
                               f'feature table and is not the first column.')
            logger.error(error)
            raise error

    return feature_data


def load_tsv_feature_table(tsv_filepath: str, header_row='auto', original_feature_id_colname: str = 'auto',
                           final_feature_id_colname: str = 'Feature ID') -> pd.DataFrame:
    """
    Loads a TSV-format FeatureTable as a pandas DataFrame. Roughly tries to handle screening out a header line.
    :param tsv_filepath: path to the TSV-format FeatureTable
    :param header_row: manually specify the header row as per the header argument in pd.read_csv (e.g., by integer,
                        where 0 is the top row), or specify 'auto' to try to auto-detect the correct header row.
    :param original_feature_id_colname: name of the feature ID column in the provided table. Must be the first column in
                                        the table. Auto-detection is attempted if 'auto' is provided.
    :param final_feature_id_colname: updated name of the feature ID column, for standardization. Column name is kept
                                     as-is if None is provided.
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

    feature_data = validate_feature_table(feature_data, original_feature_id_colname=original_feature_id_colname,
                                          final_feature_id_colname=final_feature_id_colname)

    return feature_data


def load_biom_feature_table(biom_filepath: str, feature_id_colname: str = 'Feature ID') -> pd.DataFrame:
    """
    Load a biom-format feature table as a pandas DataFrame
    :param biom_filepath: path to the biom-format file to load
    :param feature_id_colname: name to give to the feature ID column after loading.
    :return: FeatureTable data as a pandas DataFrame
    """
    feature_data = biom.load_table(biom_filepath)\
        .to_dataframe()\
        .reset_index(names=feature_id_colname)

    feature_data = validate_feature_table(feature_data, original_feature_id_colname=feature_id_colname,
                                          final_feature_id_colname=None)

    return feature_data


def load_qza_feature_table(qza_filepath: str, feature_id_colname: str = 'Feature ID',
                           tmp_dir: str = '.') -> pd.DataFrame:
    """
    Load a BIOM file from a QIIME2 QZA archive as a Pandas dataframe.

    :param qza_filepath: Path to the QIIME2 QZA archive, FeatureTable[Frequency] format
    :param feature_id_colname: name to give to the feature ID column after loading.
    :param tmp_dir: The base directory to extract the ZIP file to (will create a random subfolder)
    :return: FeatureTable data as a pandas DataFrame
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
    feature_data = load_biom_feature_table(biom_path, feature_id_colname=feature_id_colname)

    # Delete tmp dir
    logger.debug(f'Removing temp dir {tmp_subdir}')
    shutil.rmtree(tmp_subdir)

    return feature_data


def load_feature_table(feature_table_filepath, header_row='auto', original_feature_id_colname: str = 'auto',
                       final_feature_id_colname: str = 'Feature ID', tmp_dir: str = '.') -> pd.DataFrame:
    """
    Generalized function to load a feature table regardless of input format

    :param feature_table_filepath: path to the feature table to load (in TSV, BIOM, or QZA format)
    :param header_row: for TSV format files, manually specify the header row here, as per the header argument in
                        pd.read_csv (e.g., by integer, where 0 is the top row), or specify 'auto' to try to auto-detect
                        the correct header row.
    :param original_feature_id_colname: for TSV format files, optionally specify the name of the feature ID column in
                                        the provided table. Must be the first column in the table. Auto-detection is
                                        attempted if 'auto' is provided.
    :param final_feature_id_colname: final name to give to the feature ID column, for standardization. Applies to all
                                     input file types. For TSV format files, if None is provided, then the feature ID
                                     column name is kept as-is.
    :param tmp_dir: for QZA format files, optionally provide the base directory to extract the ZIP file to (will create
                     a random subfolder here)
    :return: FeatureTable data as a pandas DataFrame
    """
    try:
        # First try: assume TSV format
        logger.debug(f'Trying to load FeatureTable as TSV: {feature_table_filepath}')
        feature_data = load_tsv_feature_table(feature_table_filepath, header_row=header_row,
                                              original_feature_id_colname=original_feature_id_colname,
                                              final_feature_id_colname=final_feature_id_colname)
        logger.debug(f'TSV loading succeeded')
    except UnicodeDecodeError:
        try:
            # Second try: BIOM format
            logger.debug(f'TSV loading failed')
            logger.debug(f'Trying to load FeatureTable as BIOM: {feature_table_filepath}')
            feature_data = load_biom_feature_table(feature_table_filepath, feature_id_colname=final_feature_id_colname)
            logger.debug(f'BIOM loading succeeded')
        except UnicodeDecodeError:
            # Third try: QZA format
            logger.debug(f'BIOM loading failed')
            logger.debug(f'Trying to load FeatureTable as QZA: {feature_table_filepath}')
            feature_data = load_qza_feature_table(feature_table_filepath, feature_id_colname=final_feature_id_colname,
                                                  tmp_dir=tmp_dir,)
            logger.debug(f'QZA loading succeeded')

    return feature_data


def load_taxonomy(taxonomy_filepath: str) -> pd.DataFrame:
    """
    Load a TSV or QZA-format taxonomy file as a Pandas dataframe.

    :param taxonomy_filepath: Path to the taxonomy file, in TSV or QZA format
    :return: pandas DataFrame of the BIOM table's count data
    """
    try:
        # Attempt 1: TSV format
        logger.debug(f'Trying to load taxonomy as TSV: {taxonomy_filepath}')
        taxonomy_data = pd.read_csv(taxonomy_filepath, sep='\t')
        logger.debug(f'TSV loading succeeded')
    except UnicodeDecodeError:
        # Attempt 2: QZA format
        logger.debug(f'TSV loading failed')
        logger.debug(f'Trying to load taxonomy as QZA: {taxonomy_filepath}')
        taxonomy_contents = unpack_from_qza(taxonomy_filepath, target_filename='taxonomy.tsv')
        taxonomy_data = pd.read_csv(io.BytesIO(taxonomy_contents), sep='\t')
        logger.debug(f'QZA loading succeeded')

    return taxonomy_data


def load_fasta_sequence_table(fasta_filepath: str) -> pd.DataFrame:
    """
    Loads a FastA file into a pandas DataFrame with two columns: Feature ID and Sequence

    :param fasta_filepath: path to the FastA file
    :return: DNA sequence info as a pandas DataFrame
    """
    fasta_ids = []
    fasta_seqs = []
    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            fasta_ids.append(record.id)
            fasta_seqs.append(str(record.seq))

    sequence_table = pd.DataFrame({'Feature ID': fasta_ids, 'Sequence': fasta_seqs})

    return sequence_table


def load_qza_sequence_table(qza_filepath: str) -> pd.DataFrame:
    """
    Load a QZA DNA sequences file as a pandas DataFrame "sequence table" with two columns: Feature ID and Sequence

    :param qza_filepath: path to the QZA DNA sequences file
    :return: DNA sequence info as a pandas DataFrame
    """
    # Load the FastA file contents in binary, then use this to load the sequence table
    fasta_contents = unpack_from_qza(qza_filepath, target_filename='dna-sequences.fasta')
    parsable_fasta_contents = io.StringIO(fasta_contents.decode())

    # Load sequence table
    fasta_ids = []
    fasta_seqs = []
    for record in SeqIO.parse(parsable_fasta_contents, 'fasta'):
        fasta_ids.append(record.id)
        fasta_seqs.append(str(record.seq))

    sequence_table = pd.DataFrame({'Feature ID': fasta_ids, 'Sequence': fasta_seqs})

    return sequence_table


def load_sequence_table(sequence_data_filepath: str) -> pd.DataFrame:
    """
    Load a FastA or QZA-format DNA sequence file into a Pandas dataframe.

    :param sequence_data_filepath: Path to the DNA sequence file, in FastA or QZA format
    :return: DNA sequence info as a pandas DataFrame
    """
    try:
        # Attempt 1: FastA format
        logger.debug(f'Trying to load DNA sequences as FastA: {sequence_data_filepath}')
        sequence_table = load_fasta_sequence_table(sequence_data_filepath)
        logger.debug(f'FastA loading succeeded')
    except UnicodeDecodeError:
        # Attempt 2: QZA format
        logger.debug(f'FastA loading failed')
        logger.debug(f'Trying to load DNA sequences as QZA: {sequence_data_filepath}')
        sequence_table = load_qza_sequence_table(sequence_data_filepath)
        logger.debug(f'QZA loading succeeded')

    return sequence_table
