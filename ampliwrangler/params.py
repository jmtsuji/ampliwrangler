#!/usr/bin/env python
"""
ampliwrangler params
Description: global variables and params for ampliwrangler
Copyright: Jackson M. Tsuji, 2025
"""
from importlib.metadata import version, PackageNotFoundError

TAXONOMY_RANKS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
DEFAULT_METADATA_COLUMNS = ['Feature ID', 'Taxonomy'] + TAXONOMY_RANKS + ['Sequence']
ACCEPTABLE_FEATURE_ID_COLUMN_NAMES = ['Feature ID', '#OTU ID']

try:
    VERSION = version('ampliwrangler')
except PackageNotFoundError:
    VERSION = 'unknown'
