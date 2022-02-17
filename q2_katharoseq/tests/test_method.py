from unittest import TestCase, main

import tempfile
import os
import pandas as pd
from qiime2 import CategoricalMetadataColumn
from qiime2 import NumericMetadataColumn


from q2_katharoseq import read_count_threshold
from q2_katharoseq import get_threshold
from q2_katharoseq import allosteric_sigmoid


class KatharoSeqTestCase(TestCase):

    def setUp(self):

        self.output_dir = '.'
        self.positive_control_value = 'a'
        ind = pd.Index(['s1', 's2', 's3', 's4'],
                       name='sampleid')
        positive_control_column = pd.Series(
            ['a', 'b', 'a', 'b'],
            index=ind)
        self.positive_control_column = qiime2.CategoricalMetadataColumn(
            positive_control_column)

        cell_count_column = pd.Series([1, 2, 3, 4])
        self.cell_count_column = qiime2.NumericMetadataColumn(
            cell_count_column)

        self.table = pd.DataFrame(
            [[0, 1, 2, 3],
             [0, 1, 2, 3],
             [5, 4, 3, 2],
             [7, 2, 3, 4]],
            index=ind,
            columns=['f1', 'f2', 'f3', 'f4'])

        self.control = 'classic'
        self.threshold = 50

    def test_katharo(self):

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

        positive_control_column = pd.Series(
                    ['not_a', 'b', 'not_a', 'b'],  # change 'a'
                    index=['s1', 's2', 's3', 's4'])

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
                ValueError,
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


if __name__ == '__main__':
    main()
