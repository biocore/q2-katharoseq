from unittest import TestCase, main
from q2_katharoseq import read_count_threshold

class KatharoSeqTestCase(TestCase):

    # def setUp(self):

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

        read_count_threshold(output_dir,
                        positive_control_value,
                        positive_control_column,
                        cell_count_column,
                        table,
                        control)



if __name__=='__main__':
    main()