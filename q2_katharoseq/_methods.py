import pandas as pd
import qiime2

def read_count_threshold(output_dir:str,
                        positive_control_value: str,
                        table: pd.DataFrame,
                        taxonomy: pd.DataFrame,
                        positive_control_column: qiime2.CategoricalMetadataColumn,
                        control: str):
    pass