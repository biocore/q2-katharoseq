import importlib
from qiime2.plugin import (Plugin, Citations, Str,
                           MetadataColumn, Categorical, Choices)
from q2_types.feature_table import (FeatureTable, Frequency)
from q2_types.feature_data import FeatureData, Taxonomy
from ._methods import read_count_threshold
import q2_katharoseq

citations = Citations.load('citations.bib', package='q2_katharoseq')

plugin = Plugin(
    name='katharoseq',
    version=q2_katharoseq.__version__,
    website='https://github.com/biocore/q2-katharoseq',
    package='q2_katharoseq',
    description=('This QIIME 2 plugin implements the KatharoSeq '
                 'protocol for analyzing low biomass samles.'),
    short_description='Plugin for KatharoSeq.',
)

plugin.visualizers.register_function(
    function=read_count_threshold,
    inputs={
        'table': FeatureTable[Frequency],
        'taxonomy': FeatureData[Taxonomy]
    },
    parameters={
        'control': Str % Choices(['zymobiomics', 'classic', 'atcc']),
        'positive_control_value': Str,
        'positive_control_column': MetadataColumn[Categorical],

    },
    input_descriptions={
        'table': (
            'A FeatureTable collapsed to the genus level (level 6) that '
            'contains the control samples.'
        ),
        'taxonomy': (
            ''
        ),
    },
    parameter_descriptions={
        'control': (
            'The type of positive control used.'
        ),
        'positive_control_value': (
            'The value in the control column that demarks which samples are '
            'the positive controls.'
        ),
        'positive_control_column': (
            'The column in the sample metadata that describes which samples '
            'are and are not controls.'
        ),
    },
    name='Methods for the application of the KatharoSeq protocol',
    description='KatharoSeq is high-throughput protocol combining laboratory '
                'and bioinformatic methods that can differentiate a true '
                'positive signal in samples with as few as 50 to 500 cells.',
    citations=[citations['minich2018']]
)
