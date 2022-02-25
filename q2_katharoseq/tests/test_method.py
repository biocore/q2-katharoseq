from unittest import TestCase, main

import tempfile
import os
import pandas as pd
import qiime2
from qiime2 import CategoricalMetadataColumn
from qiime2 import NumericMetadataColumn

from q2_katharoseq import read_count_threshold, estimating_biomass
from q2_katharoseq._methods import allosteric_sigmoid
from q2_katharoseq._methods import get_threshold

from os.path import dirname, abspath, join
from inspect import currentframe, getfile


class KatharoSeqTestCase(TestCase):

    def setUp(self):
        self.output_dir = '.'
        self.positive_control_value = 'a'
        ind = pd.Index(['s1', 's2', 's3', 's4'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        self.positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [1, 2, 3, 4],
            index=ind,
            name='cell_count_column')
        self.cell_count_column = NumericMetadataColumn(cell_count_column)

        columns = [
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
            'f__Bacillaceae;g__Bacillus',
            'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
            'o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus',
            'f3',
            'f4']
        self.table = pd.DataFrame(
            [[0, 1, 2, 3],
             [0, 1, 2, 3],
             [5, 4, 3, 2],
             [7, 2, 3, 4]],
            index=ind,
            columns=columns)

        self.control = 'classic'
        self.threshold = 50

    def test_outputs_index(self):
        with tempfile.TemporaryDirectory() as output_dir:
            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                self.table,
                self.control)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

    def test_invalid_threshold(self):
        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                'Threshold must be between 0 and 100.'):

            threshold = -1
            read_count_threshold(
                output_dir,
                threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                self.table,
                self.control)

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                'Threshold must be between 0 and 100.'):

            threshold = 101
            read_count_threshold(
                output_dir,
                threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                self.table,
                self.control)

    def test_no_positive_controls_in_col(self):
        ind = pd.Index(['s1', 's2', 's3', 's4'],
                       name='sampleid')
        positive_control_column = pd.Series(
                    ['not_a', 'b', 'not_a', 'b'],  # change 'a'
                    index=ind,
                    name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                'No positive controls found '
                'in positive control column.'):

            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                positive_control_column,
                self.cell_count_column,
                self.table,
                self.control)

    def test_no_positive_controls_in_table(self):
        ind = pd.Index(
                ['s5', 's6', 's7', 's8'],
                name='sampleid')
        table = pd.DataFrame(
            [[0, 1, 2, 3],
             [0, 1, 2, 3],
             [5, 4, 3, 2],
             [7, 2, 3, 4]],
            index=ind,  # change index
            columns=['f1', 'f2', 'f3', 'f4'])

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                KeyError,
                'No positive controls found '
                'in table.'):

            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                table,
                self.control)

    def test_sigmoid(self):
        x = 1
        h = 2
        k_prime = 3
        a = allosteric_sigmoid(x, h, k_prime)
        self.assertTrue(a == .25)

    def test_threshold(self):
        r1 = 2
        r2 = 3
        thresh = 50
        min_freq = get_threshold(r1, r2, thresh)
        self.assertTrue(min_freq == 37)

    def test_estimating_biomass(self):
        fp = join(dirname(abspath(getfile(currentframe()))), 'support_files')

        data = pd.read_csv(
            f'{fp}/input_estimating_biomass.tsv', sep='\t', dtype={
                'sample_name': str, 'total_reads': float,
                'control_cell_into_extraction': float,
                'extraction_mass_g': float,
                'positive_control': str})

        data = qiime2.Metadata.load(f'{fp}/input_estimating_biomass.tsv')

        obs = estimating_biomass(
            total_reads=data.get_column('total_reads'),
            control_cell_extraction=data.get_column('control_cell_into_extraction'),  # noqa
            min_total_reads=1150,
            positive_control_value='True',
            positive_control_column=data.get_column('positive_control'),
            pcr_template_vol=5,
            dna_extract_vol=60,
            extraction_mass_g=data.get_column('extraction_mass_g')
        )
        exp = pd.read_csv(
            f'{fp}/output_estimating_biomass.tsv', sep='\t', index_col=0)
        pd.testing.assert_frame_equal(obs, exp)


    def test_biomass_plot(self):
        fp = join(dirname(abspath(getfile(currentframe()))), 'support_files')

        data = pd.read_csv(
            f'{fp}/input_estimating_biomass.tsv', sep='\t', dtype={
                'sample_name': str, 'total_reads': float,
                'control_cell_into_extraction': float,
                'extraction_mass_g': float,
                'positive_control': str})

        data = qiime2.Metadata.load(f'{fp}/input_estimating_biomass.tsv')

        biomass_plot(
            total_reads=data.get_column('total_reads'),
            control_cell_extraction=data.get_column('control_cell_into_extraction'),  # noqa
            min_total_reads=1150,
            positive_control_value='True',
            positive_control_column=data.get_column('positive_control')
        )
        self.assertTrue(os.path.exists(index_fp))


if __name__ == '__main__':
    main()
