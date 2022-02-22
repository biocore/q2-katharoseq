import pandas as pd

from .plugin_setup import plugin
from ._format import EstimatedBiomassFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> EstimatedBiomassFmt:
    ff = EstimatedBiomassFmt()
    data.to_csv(str(ff))
    return ff


@plugin.register_transformer
def _2(ff: EstimatedBiomassFmt) -> pd.DataFrame:
    return pd.read_csv(str(ff), index_col='sample-id')
