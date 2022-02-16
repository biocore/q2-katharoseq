import importlib
from qiime2.plugin import (Plugin, Citations, Str, Int,
                           MetadataColumn, Categorical, Numeric, Choices)
from q2_types.feature_table import (FeatureTable, Frequency)
from q2_types.feature_data import FeatureData
from . import read_count_threshold

citations = Citations.load('citations.bib', package='q2_katharoseq')

plugin = Plugin(
    name='katharoseq',
    version='0.0.1',
    website='https://github.com/dpear/q2-katharoseq',
    package='q2_katharoseq',
    description=('This QIIME 2 plugin implements the katharoseq '
                 'protocol for analyzing low biomass samles.'),
    short_description='Plugin for katharoseq.',
)

plugin.visualizers.register_function(
    function=read_count_threshold, # NAME OF THE FUNCTION FOR USERS
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'control': Str % Choices(['zymobiomics', 'classic', 'atcc']),
        'threshold': Int,
        'positive_control_value': Str,
        'positive_control_column': MetadataColumn[Categorical],
        'cell_count_column': MetadataColumn[Numeric]
    },
    input_descriptions={
        'table': (
            'Feature table of frequencies.'
        )
    },
    parameter_descriptions={
        'control': (
            'Environment used to calculate reads aligning to mock '
            'community. Available options are '
            '(atcc, zymobiomics, classic).'
        ),
        'threshold': (
            'Threshold to use in calculating minimum frequency. '
            'Must be int in [0,100].'
        ),
        'positive_control_value': (
            'Value that indicates a positive control.'
        ),
        'positive_control_column': (
            'Column that contains positive controls.'
        ),
        'cell_count_column': (
            'Column that contains cell count information. '
            'Used to normalize '
        ),
    },
    name='katharoseq',
    description='Katharoseq protocol for analyzing low biomass samples.',
    citations=[]
)
