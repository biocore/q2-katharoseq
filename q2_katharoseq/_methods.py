import pandas as pd
import qiime2
import pkg_resources

# from qiime2.plugins import taxa
# Cancelling this no longer using collapse
# Require the table to be collapsed

def read_count_threshold(output_dir:str,
                        positive_control_value: str,
                        positive_control_column: qiime2.CategoricalMetadataColumn,
                        cell_count_column: qiime2.NumericMetadataColumn,
                        table: pd.DataFrame,
                        control: str)  -> None:

    positive_control_column = positive_control_column.to_series()
    cell_count_column       = cell_count_column.to_series()

    positive_controls = positive_control_column[positive_control_column==positive_control_value]
    cell_counts       = cell_counts.loc[positive_controls.index]

    if not positive_controls.shape[0]:
        raise ValueError('No positive controls found in positive control column.')

    table_positive = table.loc[set(postive_controls.index)]

    if not table_positive.shape[0]:
        raise ValueError('No positive control samples found in table.')

    df = table_positive

    # quick visual check that top 7 taxa make up most of the reads in highest input sample ("d1")
    max_cell_counts = cell_counts.idxmax()
    max_input = df.loc[max_cell_counts]
    max_inputT = max_input.T
    max_inputT = max_inputT.sort_values(max_inputT.columns[0], ascending = False).head(10)

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
                    'classic':['k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',
                            'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Paracoccus']
                    }

    # Calculate the total number of reads per sample
    df['asv_reads'] = df.sum(axis=1)
    # calculate the number of reads aligning to the mock community input genera
    df['control_reads'] = df[control_type[control]].sum(axis=1)
    # calculate the percent correctly assigned
    df['correct_assign'] = df['control_reads'] / df['asv_reads']

    # PLOTTING
    katharo = df[['correct_assign','control_reads','asv_reads']]
    katharo['log_asv_reads'] = np.log10(katharo['asv_reads'])

    # define the allosteric sigmoid equation
    def allosteric_sigmoid(x, h, k_prime):
        y = x ** h / (k_prime + x ** h)
        return y

    # fit the curve to your data
    popt, pcov = curve_fit(allosteric_sigmoid, katharo['log_asv_reads'], katharo['correct_assign'], method='dogbox')
    # print(popt)
    # plot fit curve
    x = np.linspace(0, 5, 50)
    y = allosteric_sigmoid(x, *popt)

    # plot the fit
    pylab.plot(katharo['log_asv_reads'], katharo['correct_assign'], 'o', label='data')
    pylab.plot(x,y, label='fit')
    pylab.ylim(0, 1.05)
    pylab.legend(loc='best')
    plt.savefig(os.path.join(output_dir, 'fit.svg'))
    plt.close()

    # 50% threshold (recommended)

    # assign variables and solve for X (number of reads to pass filter)
    h = popt[0]  # first value printed above graph
    k = popt[1]   # second value printed above graph
    y = 0.5 ## what you want to solve for

    min_log_reads = np.power((k/(1/y-1)),(1/h))
    min_freq_50 = np.power(10, min_log_reads).astype(int)
    min_freq_50

    # TODO: Put into visualizer

    result = pd.DataFrame([min_freq_50],
        index = ['50_percent_threshold'])
    table_html = q2templates.df_to_html(result.to_frame())

    # context = {
    #     'table': table_html,
    #     'sample_size': sample_size,
    #     'mismatched_ids': mismatched_ids
    # }

    context = {
        'table': table_html,
    }

    # TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')

    # index = os.path.join(
    #     TEMPLATES, 'mantel_assets', 'index.html')

    index = output_dir + '/index.html'

    q2templates.render(index, output_dir, context=context)




    table_html = q2templates.df_to_html(result.to_frame())

    # We know the distance matrices have matching ID sets at this point, so we
    # can safely generate all pairs of IDs using one of the matrices' ID sets
    # (it doesn't matter which one).
    scatter_data = []
    for id1, id2 in itertools.combinations(dm1.ids, 2):
        scatter_data.append((dm1[id1, id2], dm2[id1, id2]))

    plt.figure()
    x = 'Pairwise Distance (%s)' % label1
    y = 'Pairwise Distance (%s)' % label2
    scatter_data = pd.DataFrame(scatter_data, columns=[x, y])
    sns.regplot(x=x, y=y, data=scatter_data, fit_reg=False)
    plt.savefig(os.path.join(output_dir, 'mantel-scatter.svg'))
    plt.close()

    context = {
        'table': table_html,
        'sample_size': sample_size,
        'mismatched_ids': mismatched_ids
    }

    TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')

    index = os.path.join(
        TEMPLATES, 'mantel_assets', 'index.html')
    q2templates.render(index, output_dir, context=context)








