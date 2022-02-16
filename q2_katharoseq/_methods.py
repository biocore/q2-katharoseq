import pandas as pd
import qiime2
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import q2templates
import pkg_resources

# define the allosteric sigmoid equation
def allosteric_sigmoid(x, h, k_prime):
    y = x ** h / (k_prime + x ** h)
    return y

def _threshold(r1,r2,thresh):
        # 50% threshold (recommended)
        # assign variables and solve for X (number of reads to pass filter)
        popt, pcov = curve_fit(allosteric_sigmoid, r1, r2, method='dogbox')
        h = popt[0]  # first value printed above graph
        k = popt[1]   # second value printed above graph
        y = thresh #0.5 ## what you want to solve for

        min_log_reads = np.power((k/(1/y-1)),(1/h))
        min_freq = np.power(10, min_log_reads).astype(int)
        return min_freq

def read_count_threshold(output_dir:str,
                        threshold: int,
                        positive_control_value: str,
                        positive_control_column: qiime2.CategoricalMetadataColumn,
                        cell_count_column: qiime2.NumericMetadataColumn,
                        table: pd.DataFrame,
                        control: str)  -> None:

    positive_control_column = positive_control_column.to_series()
    cell_count_column       = cell_count_column.to_series()

    positive_controls = positive_control_column[positive_control_column==positive_control_value]
    cell_counts       = cell_count_column.loc[positive_controls.index]

    if not positive_controls.shape[0]:
        raise ValueError('No positive controls found in positive control column.')

    table_positive = table.loc[set(positive_controls.index)]

    if not table_positive.shape[0]:
        raise ValueError('No positive control samples found in table.')

    if threshold > 100 or threshold < 0:
        raise ValueError('Threshold must be between 0 and 100.')

    df = table_positive

    # quick visual check that top 7 taxa make up most of the reads in highest input sample ("d1")
    max_cell_counts = cell_counts.idxmax()
    max_input = df.loc[max_cell_counts]
    max_inputT = max_input.T
    max_inputT = max_inputT.sort_values(ascending = False).head(10)

    control_type = {'atcc':['d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium',
                           'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;__',
                           'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella',
                           'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus'],
                    'zymobiomics':['d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Listeriaceae;g__Listeria',
                           'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas',
                           'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',
                           'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;__',
                           'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella',
                           'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus',
                           'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus',
                           'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus'],
                    'classic':['d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',
                            'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus']
                    }

    # Calculate the total number of reads per sample
    df['asv_reads'] = df.sum(axis=1)
    # calculate the number of reads aligning to the mock community input genera
    df['control_reads'] = df[control_type[control]].sum(axis=1)
    # calculate the percent correctly assigned
    df['correct_assign'] = df['control_reads'] / df['asv_reads']

    # PLOTTING
    katharo = df[['correct_assign','control_reads','asv_reads']]
    katharo['log_asv_reads'] = np.log10(katharo['asv_reads'].values)

    # fit the curve to your data
    popt, pcov = curve_fit(allosteric_sigmoid, katharo['log_asv_reads'], katharo['correct_assign'], method='dogbox')

    # plot the fit
    # When plotting 
    x = np.linspace(0, 5, 50)
    y = allosteric_sigmoid(x, *popt)
    plt.plot(katharo['log_asv_reads'], katharo['correct_assign'], 'o', label='data')
    plt.plot(x,y, label='fit')
    plt.ylim(0, 1.05)
    plt.legend(loc='best')
    plt.savefig(os.path.join(output_dir, 'fit.svg'))
    plt.close()

    min_freq = _threshold(katharo['log_asv_reads'], katharo['correct_assign'], threshold/100)

    # TODO: Put into visualizer

    max_input_html = q2templates.df_to_html(max_inputT.to_frame())
    context = {'minimum_frequency': min_freq,
               'threshold': threshold,
               'table': max_input_html
            }


    TEMPLATES = pkg_resources.resource_filename('q2_katharoseq', 'read_count_threshold_assets')
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context=context)








