import importlib
from qiime2.plugin import (Plugin, Citations, Str, Int,
                           MetadataColumn, Categorical, Numeric, Choices)
from q2_types.feature_table import (FeatureTable, Frequency)
from . import (read_count_threshold, estimating_biomass, biomass_plot,
               control_type)
import q2_katharoseq
from q2_katharoseq._type import EstimatedBiomass
from q2_katharoseq._format import EstimatedBiomassFmt, EstimatedBiomassDirFmt


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


plugin.register_formats(EstimatedBiomassFmt, EstimatedBiomassDirFmt)
plugin.register_semantic_types(EstimatedBiomass)
plugin.register_semantic_type_to_format(EstimatedBiomass,
                                        artifact_format=EstimatedBiomassDirFmt)


plugin.visualizers.register_function(
    function=read_count_threshold,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'control': Str % Choices(control_type.keys()),
        'threshold': Int,
        'positive_control_value': Str,
        'positive_control_column': MetadataColumn[Categorical],
        'cell_count_column': MetadataColumn[Numeric]
    },
    input_descriptions={
        'table': (
            'A FeatureTable collapsed to the genus level (level 6) that '
            'contains the control samples.'
        ),
    },
    parameter_descriptions={
        'control': (
            'The type of positive control used.'
        ),
        'threshold': (
            'Threshold to use in calculating minimum frequency. '
            'Must be int in [0,100].'
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
    citations=[]
)

plugin.methods.register_function(
    function=estimating_biomass,
    inputs={},
    parameters={'total_reads': MetadataColumn[Numeric],
                'control_cell_extraction': MetadataColumn[Numeric],
                'positive_control_column': MetadataColumn[Categorical],
                'positive_control_value': Str,
                'extraction_mass_g': MetadataColumn[Numeric],
                'min_total_reads': Int,
                'pcr_template_vol': Int,
                'dna_extract_vol': Int},
    outputs=[('estimated_biomass', EstimatedBiomass)],
    input_descriptions={},
    parameter_descriptions={
        'total_reads': 'The total sum of the reads or ASVs for each sample.',
        'control_cell_extraction': (
            'The estimated number of cells or genomes used as input to your '
            'library prep. One may typically estimate this by determining the '
            'total number of cells from a stock solution used to make '
            'standard titrations. Each titration will have an estimated '
            'number of microbial cells put into the extraction. The final '
            'estimate will depend on the elution volume and the final volume '
            'used into the library prep (e.g. 16S PCR).'),
        'positive_control_column': (
            'The column in the sample metadata that describes which samples '
            'are and are not controls.'),
        'positive_control_value': (
            'The value in the control column that demarks which samples are '
            'the positive controls.'),
        'extraction_mass_g': (
            'The column in the sample metadata that describes the sample '
            '(e.g. stool, tissue, soil, etc) mass (in grams - typically '
            'converted from mg)'),
        'min_total_reads': 'The minimum threshold to apply.',
        'pcr_template_vol': (
            'The volume of DNA used as template in the '
            'library prep (PCR reaction)'),
        'dna_extract_vol': (
            'The final elution volume used during DNA extraction')},
    output_descriptions={
        'estimated_biomass': (
            'A dataframe containing the details on estimated biomass')
        },
    name='Estimate the biomass of samples using KatharoSeq controls.',
    description='Estimate the biomass of samples using KatharoSeq controls.',
    citations=[]
)

plugin.visualizers.register_function(
    function=biomass_plot,
    inputs={},
    parameters={'total_reads': MetadataColumn[Numeric],
                'control_cell_extraction': MetadataColumn[Numeric],
                'positive_control_column': MetadataColumn[Categorical],
                'positive_control_value': Str,
                'min_total_reads': Int
                },
    input_descriptions={},
    parameter_descriptions={
        'total_reads': 'The total sum of the reads or ASVs for each sample.',
        'control_cell_extraction': (
            'The estimated number of cells or genomes used as input to your '
            'library prep. One may typically estimate this by determining the '
            'total number of cells from a stock solution used to make '
            'standard titrations. Each titration will have an estimated '
            'number of microbial cells put into the extraction. The final '
            'estimate will depend on the elution volume and the final volume '
            'used into the library prep (e.g. 16S PCR).'),
        'positive_control_column': (
            'The column in the sample metadata that describes which samples '
            'are and are not controls.'),
        'positive_control_value': (
            'The value in the control column that demarks which samples are '
            'the positive controls.'),
        'min_total_reads': 'The minimum threshold to apply.',
    },
    name='Plot the results of estimating_biomass.',
    description='Plot the results of estimating_biomass.',
    citations=[]
)

importlib.import_module('q2_katharoseq._transformer')
