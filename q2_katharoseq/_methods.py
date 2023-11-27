import pandas as pd
import qiime2
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import q2templates
import pkg_resources
import math
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
    'gg2-2022.10': [
        'd__Bacteria;p__Firmicutes_D;c__Bacilli;o__Bacillales_B_306089;f__Bacillaceae_H_294103;g__Bacillus_P_294101',
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales_616009;f__Xanthomonadaceae_616009;g__Lysobacter_A_615995',
        'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus'
        ],
    'single': [
        'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;'
        'o__Burkholderiales;f__Comamonadaceae;g__Variovorax']
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
        control: str) -> None:

    # CONVERSIONS
    positive_control_column = positive_control_column.to_series()
    cell_count_column = cell_count_column.to_series()

    # # FILTER COLUMNS
    positive_controls = positive_control_column[
        positive_control_column == positive_control_value]

    if not positive_controls.shape[0]:
        raise ValueError('No positive controls found in ' +
                         'positive control column.')
    positive_controls = pd.Series(positive_controls)

    cell_counts = cell_count_column.loc[positive_controls.index]

    # # CHECK SHAPES
    inds = positive_controls.index.intersection(table.index)
    print(inds)
    if len(inds) == 0:
        raise KeyError('No positive controls found in table.')
    table_positive = table.loc[inds]

    if threshold > 100 or threshold < 0:
        raise ValueError('Threshold must be between 0 and 100.')
    df = table_positive

    # VISUAL CHECK: TOP 7 TAXA MAKE UP MOST OF THE
    # READS IN HIGHEST INPUT SAMPLE
    max_cell_counts = cell_counts.idxmax()

    if max_cell_counts not in df.index.values:
        raise KeyError('No positive controls found in table.')

    max_input = df.loc[max_cell_counts]
    max_inputT = max_input.T
    max_inputT = max_inputT.sort_values(ascending=False).head(10)

    # Calculate the total number of reads per sample
    df['asv_reads'] = df.sum(axis=1)

    # NUMBER READS ALIGNING TO MOCK COMMUNITY INPUT
    df['control_reads'] = df[control_type[control]].sum(axis=1)

    # PERCENT CORRECTLY ASSIGNED
    df['correct_assign'] = df['control_reads'] / df['asv_reads']

    # DEFINE KATHARO
    katharo = df[['correct_assign', 'control_reads', 'asv_reads']].copy()
    katharo['log_asv_reads'] = np.log10(katharo['asv_reads'].values)

    # FIT CURVE TO DATA
    popt, pcov = curve_fit(allosteric_sigmoid,
                           katharo['log_asv_reads'],
                           katharo['correct_assign'],
                           method='dogbox')

    # PLOT
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

    # FIND THRESHOLD
    min_freq = get_threshold(katharo['log_asv_reads'],
                             katharo['correct_assign'],
                             threshold/100)

    # VISUALIZER
    max_input_html = q2templates.df_to_html(max_inputT.to_frame())
    context = {'minimum_frequency': min_freq,
               'threshold': threshold,
               'table': max_input_html}
    TEMPLATES = pkg_resources.resource_filename(
        'q2_katharoseq', 'read_count_threshold_assets')
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context=context)


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

    # MAKE PLOT
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

    # VISUALIZER: git push origin visualizer
    TEMPLATES = pkg_resources.resource_filename(
        'q2_katharoseq', 'estimating_biomass_assets')
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir)
