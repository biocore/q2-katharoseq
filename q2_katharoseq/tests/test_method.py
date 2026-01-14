from unittest import TestCase, main

import tempfile
import os
import sys
from io import StringIO
import pandas as pd
import qiime2
from qiime2 import CategoricalMetadataColumn
from qiime2 import NumericMetadataColumn

from q2_katharoseq import (read_count_threshold,
                           estimating_biomass,
                           biomass_plot)
from q2_katharoseq._methods import allosteric_sigmoid
from q2_katharoseq._methods import get_threshold

from os.path import dirname, abspath, join
from inspect import currentframe, getfile


class KatharoSeqTestCase(TestCase):

    def setUp(self):
        self.output_dir = '.'
        self.positive_control_value = 'a'
        # need at least 3 positive controls with unique cell counts
        ind = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        self.positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        # s1, s3, s5 are positive controls with values 100, 1000, 10000
        cell_count_column = pd.Series(
            [100, 2, 1000, 4, 10000, 6],
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
            [[1, 1, 2, 3],
             [2, 1, 2, 3],
             [10, 4, 3, 2],
             [20, 2, 3, 4],
             [100, 5, 6, 7],
             [200, 8, 9, 10]],
            index=ind,
            columns=columns)

        self.control = 'classic'
        self.threshold = 50

        folder = '../../example'
        self.fp = join(dirname(abspath(getfile(currentframe()))), folder)

    def test_specify_asv_as_control(self):
        with tempfile.TemporaryDirectory() as output_dir:
            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                self.table,
                'asv',
                'f4')

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

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
                "No positive controls found in positive control column. "
                "Searched for 'a' but found these values:"):

            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                positive_control_column,
                self.cell_count_column,
                self.table,
                self.control)

    def test_no_positive_controls_in_table(self):
        # use sample IDs that don't overlap with metadata at all
        ind = pd.Index(
                ['x1', 'x2', 'x3', 'x4'],
                name='sampleid')
        table = pd.DataFrame(
            [[0, 1, 2, 3],
             [0, 1, 2, 3],
             [5, 4, 3, 2],
             [7, 2, 3, 4]],
            index=ind,  # completely different from metadata indices
            columns=['f1', 'f2', 'f3', 'f4'])

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                KeyError,
                "No positive controls found in feature table"):

            read_count_threshold(
                output_dir,
                self.threshold,
                self.positive_control_value,
                self.positive_control_column,
                self.cell_count_column,
                table,
                self.control)

    def test_partial_positive_controls_warning(self):
        """Test that a warning is raised when some controls are missing."""
        # create metadata with 4 positive controls
        ind = pd.Index(['s1', 's2', 's3', 's4'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'a', 'a', 'a'],  # all are positive controls
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [100, 1000, 10000, 100000],
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        # create table with only 3 of the 4 samples (s4 missing)
        columns = [
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
            'f__Bacillaceae;g__Bacillus',
            'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
            'o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus',
            'f3',
            'f4']
        table = pd.DataFrame(
            [[5, 1, 2, 3],
             [10, 1, 2, 3],
             [15, 4, 3, 2]],
            index=pd.Index(['s1', 's2', 's3'], name='sampleid'),
            columns=columns)

        with tempfile.TemporaryDirectory() as output_dir:
            # capture stderr to check for warning message
            captured_stderr = StringIO()
            old_stderr = sys.stderr
            sys.stderr = captured_stderr

            try:
                read_count_threshold(
                    output_dir,
                    self.threshold,
                    'a',
                    positive_control_column,
                    cell_count_column,
                    table,
                    self.control)
            finally:
                sys.stderr = old_stderr

            # check that warning was printed to stderr
            stderr_output = captured_stderr.getvalue()
            self.assertIn("Only 3 of 4 positive controls found", stderr_output)

    def test_insufficient_dilution_series(self):
        """Test error when fewer than 3 unique cell count values."""
        ind = pd.Index(['s1', 's2', 's3', 's4'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        # only 2 unique values for the 'a' controls
        cell_count_column = pd.Series(
            [100, 2, 100, 4],  # s1 and s3 have same value
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        columns = [
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
            'f__Bacillaceae;g__Bacillus',
            'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
            'o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus',
            'f3',
            'f4']
        table = pd.DataFrame(
            [[0, 1, 2, 3],
             [0, 1, 2, 3],
             [5, 4, 3, 2],
             [7, 2, 3, 4]],
            index=ind,
            columns=columns)

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                "Insufficient dilution series"):

            read_count_threshold(
                output_dir,
                self.threshold,
                'a',
                positive_control_column,
                cell_count_column,
                table,
                self.control)

    def test_zero_total_reads(self):
        """Test error when positive control samples have zero total reads."""
        ind = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [100, 2, 1000, 4, 10000, 6],
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        columns = [
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
            'f__Bacillaceae;g__Bacillus',
            'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
            'o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus',
            'f3',
            'f4']
        # s1 has zero total reads (all zeros)
        table = pd.DataFrame(
            [[0, 0, 0, 0],
             [2, 1, 2, 3],
             [10, 4, 3, 2],
             [20, 2, 3, 4],
             [100, 5, 6, 7],
             [200, 8, 9, 10]],
            index=ind,
            columns=columns)

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                "positive control sample.*with zero total reads"):

            read_count_threshold(
                output_dir,
                self.threshold,
                'a',
                positive_control_column,
                cell_count_column,
                table,
                self.control)

    def test_control_taxa_not_found(self):
        """Test error when control taxa are not in the feature table."""
        ind = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [100, 2, 1000, 4, 10000, 6],
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        # table has no classic control taxa columns
        table = pd.DataFrame(
            [[1, 1, 2, 3],
             [2, 1, 2, 3],
             [10, 4, 3, 2],
             [20, 2, 3, 4],
             [100, 5, 6, 7],
             [200, 8, 9, 10]],
            index=ind,
            columns=['f1', 'f2', 'f3', 'f4'])

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                "None of the classic control taxa were found"):

            read_count_threshold(
                output_dir,
                self.threshold,
                'a',
                positive_control_column,
                cell_count_column,
                table,
                'classic')

    def test_asv_zero_reads_in_controls(self):
        """Test error when ASV has zero reads in all positive controls."""
        ind = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [100, 2, 1000, 4, 10000, 6],
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        # asv column has 0 in all positive control rows (s1, s3, s5)
        table = pd.DataFrame(
            [[1, 1, 2, 0],   # s1: a (control) - asv=0
             [2, 1, 2, 10],  # s2: b - asv=10
             [10, 4, 3, 0],  # s3: a (control) - asv=0
             [20, 2, 3, 15], # s4: b - asv=15
             [100, 5, 6, 0], # s5: a (control) - asv=0
             [200, 8, 9, 20]],  # s6: b - asv=20
            index=ind,
            columns=['f1', 'f2', 'f3', 'target_asv'])

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                "specified ASV has zero reads in all.*positive control"):

            read_count_threshold(
                output_dir,
                self.threshold,
                'a',
                positive_control_column,
                cell_count_column,
                table,
                'asv',
                'target_asv')

    def test_no_variation_in_correct_assign(self):
        """Test error when all positive controls have identical correct_assign."""
        ind = pd.Index(['s1', 's2', 's3', 's4', 's5', 's6'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b', 'a', 'b'],
            index=ind,
            name='positive_control_column')
        positive_control_column = CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series(
            [100, 2, 1000, 4, 10000, 6],
            index=ind,
            name='cell_count_column')
        cell_count_column = NumericMetadataColumn(cell_count_column)

        # all positive controls (s1, s3, s5) have identical ratios
        # asv reads: 10, and total reads: 20 each -> correct_assign = 0.5
        table = pd.DataFrame(
            [[10, 10],  # s1: 10/20 = 0.5
             [5, 15],
             [20, 20],  # s3: 20/40 = 0.5
             [10, 30],
             [50, 50],  # s5: 50/100 = 0.5
             [25, 75]],
            index=ind,
            columns=['target_asv', 'other'])

        with tempfile.TemporaryDirectory() as output_dir, \
            self.assertRaisesRegex(
                ValueError,
                "identical correct assignment ratio"):

            read_count_threshold(
                output_dir,
                self.threshold,
                'a',
                positive_control_column,
                cell_count_column,
                table,
                'asv',
                'target_asv')

    def test_sigmoid(self):
        x = 1.0
        h = 2.0
        k_prime = 3.0
        a = allosteric_sigmoid(x, h, k_prime)
        self.assertTrue(a == .25)

    def test_threshold(self):
        r1 = [3.5, 2.3, 1.3, 3.4]
        r2 = [1.1, 2.2, 1.7, 2.3]
        thresh = 50.0
        min_freq = get_threshold(r1, r2, thresh)
        self.assertTrue(min_freq == 1)

    def test_estimating_biomass(self):

        data = qiime2.Metadata.load(f'{self.fp}/fmp_metadata.tsv')

        table = f'{self.fp}/fmp_collapsed_table.qza'
        table = qiime2.Artifact.load(table).view(pd.DataFrame)

        obs = estimating_biomass(
            # total_reads=data.get_column('total_reads'),
            table=table,
            control_cell_extraction=data.get_column('control_cell_into_extraction'),  # noqa
            min_total_reads=1150,
            positive_control_value='control',
            positive_control_column=data.get_column('control_rct'),
            pcr_template_vol=5,
            dna_extract_vol=60,
            extraction_mass_g=data.get_column('extraction_mass_g')
        )

        exp = pd.read_csv(f'{self.fp}/est_biomass_output.csv', index_col=0)
        pd.testing.assert_frame_equal(obs, exp)

    def test_biomass_plot(self):

        data = qiime2.Metadata.load(f'{self.fp}/fmp_metadata.tsv')
        table = qiime2.Artifact.load(f'{self.fp}/fmp_collapsed_table.qza')
        table = table.view(pd.DataFrame)

        with tempfile.TemporaryDirectory() as output_dir:
            biomass_plot(
                output_dir,
                table=table,
                control_cell_extraction=data.get_column(
                    'control_cell_into_extraction'),  # noqa
                min_total_reads=1150,
                positive_control_value='control',
                positive_control_column=data.get_column('control_rct')
            )

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))


if __name__ == '__main__':
    main()
