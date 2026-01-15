import pandas as pd
import qiime2
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import q2templates
from importlib.resources import files
import math
import sys
from sklearn.linear_model import LinearRegression

control_type = {
    'atcc': [
        'd__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;'
        'f__Clostridiaceae;g__Clostridium',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Enterobacterales;f__Enterobacteriaceae;g__',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Enterobacterales;f__Enterobacteriaceae;'
        'g__Escherichia-Shigella',
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;'
        'f__Staphylococcaceae;g__Staphylococcus'],
    'zymobiomics': [
        'd__Bacteria;p__Firmicutes;c__Bacilli;'
        'o__Lactobacillales;f__Listeriaceae;g__Listeria',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas',
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
        'f__Bacillaceae;g__Bacillus',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Enterobacterales;f__Enterobacteriaceae;__',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Enterobacterales;f__Enterobacteriaceae;'
        'g__Escherichia-Shigella',
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
        'f__Lactobacillaceae;g__Lactobacillus',
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
        'f__Enterococcaceae;g__Enterococcus',
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales'
        ';f__Staphylococcaceae;g__Staphylococcus'],
    'classic': [
        'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
        'f__Bacillaceae;g__Bacillus',
        'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
        'o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus'],
    'single': [
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Burkholderiales;f__Comamonadaceae;g__Variovorax'],
    'asv': ''
}


# Define the allosteric sigmoid equation
def allosteric_sigmoid(x, h, k_prime):
    y = x ** h / (k_prime + x ** h)
    return y


def get_threshold(r1, r2, thresh):
    # assign variables and solve for X (number of reads to pass filter)
    popt, pcov = curve_fit(allosteric_sigmoid, r1, r2, method='dogbox')
    h = popt[0]  # first value printed above graph
    k = popt[1]  # second value printed above graph
    y = thresh   # desired thresh (50%, 70%, 90%)
    min_log_reads = np.power((k/(1/y-1)), (1/h))
    min_freq = np.power(10, min_log_reads).astype(int)
    return min_freq


def fit_lm(table,
           min_total_reads,
           positive_control_column,
           positive_control_value,
           control_cell_extraction):

    total_reads = table.sum(axis=1)
    filtered = pd.DataFrame(total_reads[total_reads > min_total_reads])
    filtered = filtered.rename(columns={0: 'total_reads'})
    filtered['log_total_reads'] = filtered.total_reads.apply(math.log10)

    positive_control_column = positive_control_column.to_series().loc[
        filtered.index]
    positive_controls = positive_control_column[
        positive_control_column == positive_control_value]
    positive_controls = filtered.loc[positive_controls.index]

    positive_controls['control_cell_extraction'] = \
        control_cell_extraction.to_series().loc[positive_controls.index]
    positive_controls['log_control_cell_extraction'] = \
        positive_controls.control_cell_extraction.apply(math.log10)

    lm = LinearRegression()
    lm.fit(
        positive_controls.log_total_reads.values.reshape(-1, 1),
        positive_controls.log_control_cell_extraction)

    return lm, filtered, positive_controls


def read_count_threshold(
        output_dir: str,
        threshold: int,
        positive_control_value: str,
        positive_control_column: qiime2.CategoricalMetadataColumn,
        cell_count_column: qiime2.NumericMetadataColumn,
        table: pd.DataFrame,
        control: str,
        asv: str = None) -> None:
    if control == 'asv':
        if asv is None:
            raise ValueError("Control type set to asv but no asv provided")
        if asv not in table.columns:
            raise ValueError("asv not found in the feature table")

    # conversions
    positive_control_column = positive_control_column.to_series()
    cell_count_column = cell_count_column.to_series()

    # filter columns
    positive_controls = positive_control_column[
        positive_control_column == positive_control_value]

    if not positive_controls.shape[0]:
        unique_values = positive_control_column.unique()[:10]
        raise ValueError(
            f"No positive controls found in positive control column. "
            f"Searched for '{positive_control_value}' but found these values: "
            f"{list(unique_values)}"
        )
    positive_controls = pd.Series(positive_controls)

    # check shapes - validate overlap between metadata and feature table
    n_controls_metadata = len(positive_controls)
    inds = positive_controls.index.intersection(table.index)
    n_controls_in_table = len(inds)

    if n_controls_in_table == 0:
        missing_samples = list(positive_controls.index[:5])
        table_samples = list(table.index[:5])
        raise KeyError(
            f"No positive controls found in feature table. "
            f"Found {n_controls_metadata} controls in metadata but none match "
            f"the feature table sample IDs.\n"
            f"Example control IDs from metadata: {missing_samples}\n"
            f"Example sample IDs from table: {table_samples}\n"
            f"Check that sample IDs match between your metadata and "
            f"feature table."
        )

    if n_controls_in_table < n_controls_metadata:
        missing = positive_controls.index.difference(table.index)
        missing_cell_counts = cell_count_column.loc[missing]
        print(
            f"Warning: Only {n_controls_in_table} of {n_controls_metadata} "
            f"positive controls found in feature table. Missing "
            f"{len(missing)} samples with cell counts: "
            f"{sorted(missing_cell_counts.unique().tolist())}. "
            f"Proceeding with available controls.",
            file=sys.stderr
        )

    table_positive = table.loc[inds]

    # get cell counts only for samples that are in the table
    cell_counts = cell_count_column.loc[inds]

    # validate we have enough points for curve fitting
    unique_cell_counts = cell_counts.unique()
    if len(unique_cell_counts) < 3:
        raise ValueError(
            f"Insufficient dilution series: only {len(unique_cell_counts)} "
            f"unique cell count values found ({sorted(unique_cell_counts)}). "
            f"At least 3 different concentrations are required for "
            f"curve fitting."
        )

    if threshold > 100 or threshold < 0:
        raise ValueError('Threshold must be between 0 and 100.')
    df = table_positive

    # visual check
    max_cell_counts = cell_counts.idxmax()

    max_input = df.loc[max_cell_counts]
    max_inputT = max_input.T
    max_inputT = max_inputT.sort_values(ascending=False).head(10)

    # calculate the total number of reads per sample
    df['asv_reads'] = df.sum(axis=1)

    # validate no zero-read samples
    zero_read_samples = df[df['asv_reads'] == 0].index.tolist()
    if zero_read_samples:
        raise ValueError(
            f"Found {len(zero_read_samples)} positive control sample(s) "
            f"with zero total reads: {zero_read_samples[:5]}. Cannot "
            f"compute correct assignment ratio. These samples may have "
            f"been filtered out upstream."
        )

    # number reads aligning to mock community input
    if control == 'asv':
        df['control_reads'] = df[asv]
    else:
        # validate control taxa exist in the table
        control_taxa = control_type[control]
        missing_taxa = [t for t in control_taxa if t not in df.columns]
        if len(missing_taxa) == len(control_taxa):
            raise ValueError(
                f"None of the {control} control taxa were found in the "
                f"feature table. Expected taxa like: "
                f"{control_taxa[0][:50]}..."
            )
        # use only taxa that exist
        present_taxa = [t for t in control_taxa if t in df.columns]
        df['control_reads'] = df[present_taxa].sum(axis=1)

    # validate control has reads in at least some positive controls
    total_control_reads = df['control_reads'].sum()
    if total_control_reads == 0:
        if control == 'asv':
            raise ValueError(
                f"The specified ASV has zero reads in all {len(df)} "
                f"positive control samples. Cannot build KatharoSeq "
                f"curve. Verify the ASV sequence is correct and present "
                f"in your positive control dilution series."
            )
        else:
            raise ValueError(
                f"The {control} control taxa have zero total reads across "
                f"all {len(df)} positive control samples. Cannot build "
                f"KatharoSeq curve. Verify the control type matches your "
                f"experimental setup."
            )

    # percent correctly assigned
    df['correct_assign'] = df['control_reads'] / df['asv_reads']

    # validate there's variation in correct_assign for curve fitting
    unique_correct_assign = df['correct_assign'].nunique()
    if unique_correct_assign < 2:
        raise ValueError(
            f"All positive control samples have identical correct "
            f"assignment ratio ({df['correct_assign'].iloc[0]:.4f}). "
            f"Cannot fit curve without variation in the data. Check "
            f"that your dilution series spans a range of concentrations."
        )

    # define katharo
    katharo = df[['correct_assign', 'control_reads', 'asv_reads']].copy()
    katharo['log_asv_reads'] = np.log10(katharo['asv_reads'].values)

    # fit curve to data
    try:
        popt, pcov = curve_fit(allosteric_sigmoid,
                               katharo['log_asv_reads'],
                               katharo['correct_assign'],
                               method='dogbox')
    except RuntimeError as e:
        raise RuntimeError(
            f"Curve fitting failed to converge. This typically indicates "
            f"the data does not follow the expected sigmoid pattern. "
            f"Details: {e}\n"
            f"Data summary - log_asv_reads range: "
            f"[{katharo['log_asv_reads'].min():.2f}, "
            f"{katharo['log_asv_reads'].max():.2f}], correct_assign range: "
            f"[{katharo['correct_assign'].min():.4f}, "
            f"{katharo['correct_assign'].max():.4f}]"
        ) from e

    # plot
    x = np.linspace(0, 5, 50)
    y = allosteric_sigmoid(x, *popt)
    plt.plot(katharo['log_asv_reads'],
             katharo['correct_assign'],
             'o', label='data')
    plt.plot(x, y, label='fit')
    plt.ylim(0, 1.05)
    plt.legend(loc='best')
    plt.savefig(os.path.join(output_dir, 'fit.svg'))
    plt.close()

    # find threshold
    min_freq = get_threshold(katharo['log_asv_reads'],
                             katharo['correct_assign'],
                             threshold/100)

    # visualizer
    max_input_html = q2templates.df_to_html(max_inputT.to_frame())
    context = {'minimum_frequency': min_freq,
               'threshold': threshold,
               'table': max_input_html}
    TEMPLATES = files('q2_katharoseq') / 'read_count_threshold_assets'
    index = TEMPLATES / 'index.html'
    q2templates.render(str(index), output_dir, context=context)


def estimating_biomass(
        table: pd.DataFrame,
        control_cell_extraction: qiime2.NumericMetadataColumn,
        min_total_reads: int,
        positive_control_value: str,
        positive_control_column: qiime2.CategoricalMetadataColumn,
        pcr_template_vol: int,
        dna_extract_vol: int,
        extraction_mass_g: qiime2.NumericMetadataColumn) -> pd.DataFrame:

    lm, filtered, positive_controls = fit_lm(table,
                                             min_total_reads,
                                             positive_control_column,
                                             positive_control_value,
                                             control_cell_extraction)

    filtered['estimated_biomass_per_pcrrxn'] = \
        10**((filtered.log_total_reads*lm.coef_[0])+lm.intercept_)
    filtered['estimated_biomass_per_dnarxn'] = \
        filtered.estimated_biomass_per_pcrrxn*(
            dna_extract_vol/pcr_template_vol)

    filtered['extraction_mass_g'] = extraction_mass_g.to_series().loc[
        filtered.index]
    filtered['estimated_cells_per_g'] = \
        filtered['estimated_biomass_per_dnarxn']/filtered['extraction_mass_g']
    filtered['log_estimated_cells_per_g'] = \
        filtered.estimated_cells_per_g.apply(math.log10)

    filtered.index.rename('sample_name', inplace=True)

    return filtered


def biomass_plot(
        output_dir: str,
        table: pd.DataFrame,
        control_cell_extraction: qiime2.NumericMetadataColumn,
        min_total_reads: int,
        positive_control_value: str,
        positive_control_column: qiime2.CategoricalMetadataColumn) -> None:

    lm, filtered, positive_controls = fit_lm(table,
                                             min_total_reads,
                                             positive_control_column,
                                             positive_control_value,
                                             control_cell_extraction)

    # make plot
    y = positive_controls['log_control_cell_extraction']
    x = positive_controls['log_total_reads']
    intercept = lm.intercept_
    slope = lm.coef_[0]

    plt.clf()
    plt.scatter(x, y, color='black')
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')
    plt.xlabel('Log reads')
    plt.ylabel('Log cells')
    plt.savefig(os.path.join(output_dir, 'fit.svg'))
    plt.close()

    # visualizer
    TEMPLATES = files('q2_katharoseq') / 'estimating_biomass_assets'
    index = TEMPLATES / 'index.html'
    q2templates.render(str(index), output_dir)
