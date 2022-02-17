from unittest import TestCase, main
from q2_katharoseq import read_count_threshold, estimating_biomass
import pandas as pd
import qiime2 as q2
from os.path import dirname, abspath, join
from inspect import currentframe, getfile


class KatharoSeqTestCase(TestCase):

    def test_katharo(self):

        output_dir = '.'
        positive_control_value = 'a'
        positive_control_column = pd.Series(['a', 'b', 'a', 'b'],
            index=['s1', 's2', 's3', 's4'])
        positive_control_column = CategoricalMetadataColumn(positive_control_column)

        cell_count_column = pd.Series([1,2,3,4])
        cell_count_column = NumericMetadataColumn(cell_count_column)

        table = pd.DataFrame([[0,1,2,3],[0,1,2,3],[5,4,3,2],[7,2,3,4]],
            index=['s1','s2','s3','s4'], columns=['f1','f2','f3','f4'])

        control = 'classic'
        threshold = 50

        read_count_threshold(output_dir,
                        threshold,
                        positive_control_value,
                        positive_control_column,
                        cell_count_column,
                        table,
                        control)


    def test_estimating_biomass(self):
        fp = join(dirname(abspath(getfile(currentframe()))), 'support_files')
        data = pd.read_csv(
            f'{fp}/input_estimating_biomass.tsv', sep='\t', dtype={
                'sample_name': str, 'total_reads': float,
                'control_cell_into_extraction': float,
                'extraction_mass_g': float,
                'positive_control': str})
        data.set_index('sample_name', inplace=True)

        obs = estimating_biomass(
            total_reads=q2.NumericMetadataColumn(data['total_reads']),
            control_cell_extraction=q2.NumericMetadataColumn(
                data['control_cell_into_extraction']),
            min_total_reads=1500,
            positive_control_value='True',
            positive_control_column=q2.CategoricalMetadataColumn(
                data['positive_control']),
            pcr_template_vol=5,
            dna_extract_vol=60,
            extraction_mass_g=q2.NumericMetadataColumn(
                data['extraction_mass_g'])
        )
        exp = pd.read_csv(
            f'{fp}/output_estimating_biomass.tsv', sep='\t', index_col=0)
        pd.testing.assert_frame_equal(obs, exp)


if __name__=='__main__':
    main()
