import importlib
from qiime2.plugin import (Plugin, Citations, Str,
                           MetadataColumn, Categorical, Choices)
from q2_types.feature_table import (FeatureTable, Frequency)
from q2_types.feature_data import FeatureData, Taxonomy
from ._methods import read_count_threshold

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
            ''
        ),
        'taxonomy': (
            ''
        ),
    },
    parameter_descriptions={
        'control': (
            ''
        ),
        'positive_control_value': (
            ''
        ),
        'positive_control_column': (
            ''
        ),
    },
    name='katharoseq',
    description='',
    citations=[]
)
