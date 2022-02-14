import pandas as pd
import qiime2
import math
from sklearn.linear_model import LinearRegression

def read_count_threshold(output_dir:str,
                        positive_control_value: str,
                        table: pd.DataFrame,
                        taxonomy: pd.DataFrame,
                        positive_control_column: qiime2.CategoricalMetadataColumn,
                        control: str):
    pass


def estimating_biomass(
    total_reads: qiime2.NumericMetadataColumn,
    min_total_reads: Int,
    positive_control_value: str,
    positive_control_column: qiime2.CategoricalMetadataColumn,
    ) -> pd.DataFrame:

    total_reads = total_reads.to_dataframe()
    filtered = total_reads[
        (total_reads.total_reads > 1150) &
        (total_reads.total_reads.notnull())].copy()
    filtered['log_total_reads'] = filtered.total_reads.apply(math.log10)

    control = filtered[
        filtered[positive_control_column] == positive_control_value].copy()
    control['log_control_cell_into_extraction'] = \
        control.control_cell_into_extraction.apply(math.log10)

    lm = LinearRegression()
    lm.fit(
        control.log_total_reads.values.reshape(-1, 1),
        control.log_control_cell_into_extraction)

    # estimates number of cells in the total PCR reaction volume (e.g in this
    # case: cells per 5ul) save into new column called <estimated_biomass_per_pcrrxn>
    filtered['estimated_biomass_per_pcrrxn'] = \
        10**filtered.log_total_reads*lm.coef_[0]+lm.intercept_

    # go back to original file and calculate the log cells per extraction save
    # as new column <estimated_biomass_per_pcrrxn> enter the PCR reaction
    # volume
    pcr_template_vol = 5
    # enter the DNA extraction volume
    dna_extract_vol = 60
    # normalize this volume up to the DNA extraction volume (DNA extraction
    # volume / PCR volume) (cells per 60 ul) save into new column called
    # <estimated_biomass_per_dnarxn>
    filtered['estimated_biomass_per_dnarxn'] = \
        filtered.estimated_biomass_per_pcrrxn*(
            dna_extract_vol/pcr_template_vol)

    # normalize this volume now based on the amount of actual physical tissue
    # or primary material placed into the DNA extraction tube:
    # <extraction_mass_g> save as <estimated_cells_per_g>
    filtered['estimated_cells_per_g'] = \
        filtered['estimated_biomass_per_dnarxn']/filtered['extraction_mass_g']

    # log transform biomass estimate and save as new column
    # <log_estimated_cells_per_g>
    filtered['log_estimated_cells_per_g'] = \
        filtered.estimated_cells_per_g.apply(math.log10)

    return filtered
